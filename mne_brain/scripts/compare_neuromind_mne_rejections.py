"""Compare NeuroMIND Julia and mne_brain epoch rejection decisions.

Reads existing pipeline outputs only; it does not rerun preprocessing.

Default inputs, when run from ``mne_brain/``:
  - Julia/NeuroMIND: ``../results/subjects``
  - MNE-Python:      ``results/subjects``

Outputs:
  - ``rejection_comparison_by_recording.csv``
  - ``rejection_comparison_summary.md``
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any


TASK_TO_COND = {"eyesclosed": "EC", "eyesopen": "EO"}
COND_TO_TASK = {v: k for k, v in TASK_TO_COND.items()}


def _read_json(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}


def _read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as fh:
            return list(csv.DictReader(fh))
    except FileNotFoundError:
        return []


def _decision(n_valid: int | None, n_total: int | None, min_epochs: int) -> str:
    if n_valid is None or n_total is None:
        return "missing"
    if n_total == 0 or n_valid == 0:
        return "exclude"
    if n_valid < min_epochs:
        return "manual_review"
    return "include"


def _rate(n_rejected: int | None, n_total: int | None) -> float | None:
    if n_rejected is None or not n_total:
        return None
    return n_rejected / n_total


def _top(counter: Counter[str], n: int = 5) -> list[str]:
    return [name for name, _count in counter.most_common(n)]


def _fmt_top(counter: Counter[str], n: int = 5) -> str:
    return ";".join(f"{name}:{count}" for name, count in counter.most_common(n))


def _split_channels(text: str, sep: str) -> list[str]:
    return [item.strip() for item in text.split(sep) if item.strip()]


def load_julia_recordings(root: Path, min_epochs: int) -> dict[tuple[str, str, str], dict[str, Any]]:
    out: dict[tuple[str, str, str], dict[str, Any]] = {}
    for summary_path in root.glob("sub-*/ses-*/*/artifact_rejection_summary.json"):
        subject = summary_path.parts[-4].replace("sub-", "")
        session = summary_path.parts[-3].replace("ses-", "")
        task = summary_path.parts[-2]
        condition = TASK_TO_COND.get(task, task)
        base = summary_path.parent

        summary = _read_json(summary_path)
        n_total = _as_int(summary.get("n_total"))
        n_valid = _as_int(summary.get("n_valid"))
        n_rejected = _as_int(summary.get("n_rejected"))

        channel_counts: Counter[str] = Counter()
        channel_summary = _read_csv(base / "channel_artifact_summary.csv")
        for row in channel_summary:
            channel = row.get("channel", "").strip()
            n_bad = _as_int(row.get("n_bad"))
            if channel and n_bad:
                channel_counts[channel] += n_bad

        reason_counts: Counter[str] = Counter()
        rejected_segments = _read_csv(base / "rejected_segments.csv")
        for row in rejected_segments:
            reason = row.get("rejection_reason", "").strip()
            if reason:
                reason_counts[reason] += 1
            if not channel_summary:
                for channel in _split_channels(row.get("channels_violating", ""), ";"):
                    channel_counts[channel] += 1

        key = (subject, session, condition)
        out[key] = {
            "present": True,
            "n_total": n_total,
            "n_valid": n_valid,
            "n_rejected": n_rejected,
            "rejection_rate": _rate(n_rejected, n_total),
            "decision": _decision(n_valid, n_total, min_epochs),
            "top_channels": _top(channel_counts),
            "top_channels_with_counts": _fmt_top(channel_counts),
            "reason_counts": ";".join(f"{k}:{v}" for k, v in reason_counts.most_common()),
            "dominant_reason": reason_counts.most_common(1)[0][0] if reason_counts else "",
            "n_channels_used": _as_int(summary.get("n_channels_used")),
        }
    return out


def load_mne_recordings(root: Path, min_epochs: int) -> dict[tuple[str, str, str], dict[str, Any]]:
    out: dict[tuple[str, str, str], dict[str, Any]] = {}
    for summary_path in root.glob("sub-*/ses-*/*/tables/epoch_summary_*.json"):
        subject = summary_path.parts[-5].replace("sub-", "")
        session = summary_path.parts[-4].replace("ses-", "")
        task = summary_path.parts[-3]
        condition = TASK_TO_COND.get(task, summary_path.stem.rsplit("_", 1)[-1])
        base = summary_path.parent

        summary = _read_json(summary_path)
        n_total = _as_int(summary.get("n_initial"))
        n_valid = _as_int(summary.get("n_valid"))
        n_rejected = _as_int(summary.get("n_rejected"))

        channel_counts: Counter[str] = Counter()
        for row in _read_csv(base / f"epoch_drops_{condition}.csv"):
            for channel in _split_channels(row.get("reason", ""), ","):
                channel_counts[channel] += 1

        key = (subject, session, condition)
        out[key] = {
            "present": True,
            "n_total": n_total,
            "n_valid": n_valid,
            "n_rejected": n_rejected,
            "rejection_rate": _rate(n_rejected, n_total),
            "decision": _decision(n_valid, n_total, min_epochs),
            "top_channels": _top(channel_counts),
            "top_channels_with_counts": _fmt_top(channel_counts),
            "reason_counts": "mne_ptp_channel_drop:" + str(n_rejected or 0),
            "dominant_reason": "mne_ptp_channel_drop" if n_rejected else "",
            "quality_exclusion": (base / f"quality_exclusion_{condition}.json").is_file(),
            "n_channels": _as_int(summary.get("n_channels")),
        }
    return out


def _as_int(value: Any) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _as_pct(value: float | None) -> str:
    if value is None:
        return ""
    return f"{100.0 * value:.1f}"


def compare_recordings(
    julia: dict[tuple[str, str, str], dict[str, Any]],
    mne: dict[tuple[str, str, str], dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for subject, session, condition in sorted(set(julia) | set(mne)):
        j = julia.get((subject, session, condition), {})
        m = mne.get((subject, session, condition), {})
        j_top = set(j.get("top_channels", []))
        m_top = set(m.get("top_channels", []))
        union = j_top | m_top
        overlap = j_top & m_top
        jaccard = len(overlap) / len(union) if union else None

        j_decision = j.get("decision", "missing")
        m_decision = m.get("decision", "missing")
        rate_diff = _diff(j.get("rejection_rate"), m.get("rejection_rate"))

        rows.append(
            {
                "subject": subject,
                "session": session,
                "condition": condition,
                "julia_present": bool(j),
                "mne_present": bool(m),
                "julia_decision": j_decision,
                "mne_decision": m_decision,
                "decision_match": j_decision == m_decision,
                "julia_n_total": j.get("n_total", ""),
                "julia_n_valid": j.get("n_valid", ""),
                "julia_n_rejected": j.get("n_rejected", ""),
                "julia_rejection_pct": _as_pct(j.get("rejection_rate")),
                "mne_n_total": m.get("n_total", ""),
                "mne_n_valid": m.get("n_valid", ""),
                "mne_n_rejected": m.get("n_rejected", ""),
                "mne_rejection_pct": _as_pct(m.get("rejection_rate")),
                "rejection_pct_diff_abs": "" if rate_diff is None else f"{100.0 * abs(rate_diff):.1f}",
                "julia_dominant_reason": j.get("dominant_reason", ""),
                "mne_dominant_reason": m.get("dominant_reason", ""),
                "julia_top_channels": j.get("top_channels_with_counts", ""),
                "mne_top_channels": m.get("top_channels_with_counts", ""),
                "top_channel_overlap": ";".join(sorted(overlap)),
                "top_channel_jaccard": "" if jaccard is None else f"{jaccard:.3f}",
                "motive_similarity": _motive_similarity(j_decision, m_decision, overlap, jaccard),
            }
        )
    return rows


def _diff(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    return a - b


def _motive_similarity(
    julia_decision: str,
    mne_decision: str,
    top_overlap: set[str],
    jaccard: float | None,
) -> str:
    if julia_decision == "missing" or mne_decision == "missing":
        return "missing"
    if julia_decision != mne_decision:
        return "decision_mismatch"
    if julia_decision == "include":
        return "both_include"
    if jaccard is not None and jaccard >= 0.6:
        return "similar_top_channels"
    if len(top_overlap) >= 2:
        return "partial_channel_overlap"
    if len(top_overlap) == 1:
        return "weak_channel_overlap"
    return "different_top_channels"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path: Path, rows: list[dict[str, Any]], julia_root: Path, mne_root: Path) -> None:
    total = len(rows)
    decision_matches = sum(row["decision_match"] for row in rows)
    excluded_j = [r for r in rows if r["julia_decision"] == "exclude"]
    excluded_m = [r for r in rows if r["mne_decision"] == "exclude"]
    both_excluded = [r for r in rows if r["julia_decision"] == "exclude" and r["mne_decision"] == "exclude"]
    only_j = [r for r in rows if r["julia_decision"] == "exclude" and r["mne_decision"] != "exclude"]
    only_m = [r for r in rows if r["mne_decision"] == "exclude" and r["julia_decision"] != "exclude"]

    similarity_counts = Counter(str(r["motive_similarity"]) for r in rows)
    decision_pairs = Counter((str(r["julia_decision"]), str(r["mne_decision"])) for r in rows)

    lines = [
        "# NeuroMIND vs MNE Rejection Comparison",
        "",
        f"- Julia root: `{julia_root}`",
        f"- MNE root: `{mne_root}`",
        f"- Recordings compared: {total}",
        f"- Decision matches: {decision_matches}/{total} ({100.0 * decision_matches / max(total, 1):.1f}%)",
        f"- Julia exclusions: {len(excluded_j)}",
        f"- MNE exclusions: {len(excluded_m)}",
        f"- Excluded by both: {len(both_excluded)}",
        f"- Excluded only by Julia: {len(only_j)}",
        f"- Excluded only by MNE: {len(only_m)}",
        "",
        "## Decision Pairs",
        "",
        "| Julia | MNE | n |",
        "|---|---|---:|",
    ]
    for (j_decision, m_decision), count in sorted(decision_pairs.items()):
        lines.append(f"| {j_decision} | {m_decision} | {count} |")

    lines += [
        "",
        "## Motive Similarity",
        "",
        "| Class | n |",
        "|---|---:|",
    ]
    for label, count in similarity_counts.most_common():
        lines.append(f"| {label} | {count} |")

    lines += [
        "",
        "## Excluded By Both",
        "",
        "| Recording | Julia top channels | MNE top channels | Similarity |",
        "|---|---|---|---|",
    ]
    for row in both_excluded:
        rec = f"{row['subject']} {row['session']} {row['condition']}"
        lines.append(
            f"| {rec} | {row['julia_top_channels']} | "
            f"{row['mne_top_channels']} | {row['motive_similarity']} |"
        )

    if only_j:
        lines += ["", "## Excluded Only By Julia", "", "| Recording | Julia valid/total | MNE valid/total | Julia top channels |", "|---|---:|---:|---|"]
        for row in only_j:
            rec = f"{row['subject']} {row['session']} {row['condition']}"
            lines.append(
                f"| {rec} | {row['julia_n_valid']}/{row['julia_n_total']} | "
                f"{row['mne_n_valid']}/{row['mne_n_total']} | {row['julia_top_channels']} |"
            )

    if only_m:
        lines += ["", "## Excluded Only By MNE", "", "| Recording | Julia valid/total | MNE valid/total | MNE top channels |", "|---|---:|---:|---|"]
        for row in only_m:
            rec = f"{row['subject']} {row['session']} {row['condition']}"
            lines.append(
                f"| {rec} | {row['julia_n_valid']}/{row['julia_n_total']} | "
                f"{row['mne_n_valid']}/{row['mne_n_total']} | {row['mne_top_channels']} |"
            )

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--julia-root", type=Path, default=Path("../results/subjects"))
    parser.add_argument("--mne-root", type=Path, default=Path("results/subjects"))
    parser.add_argument("--out-dir", type=Path, default=Path("/private/tmp/neuromind_mne_rejection_comparison"))
    parser.add_argument("--min-epochs", type=int, default=20)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    julia_root = args.julia_root.resolve()
    mne_root = args.mne_root.resolve()
    out_dir = args.out_dir.resolve()

    julia = load_julia_recordings(julia_root, args.min_epochs)
    mne = load_mne_recordings(mne_root, args.min_epochs)
    rows = compare_recordings(julia, mne)

    csv_path = out_dir / "rejection_comparison_by_recording.csv"
    md_path = out_dir / "rejection_comparison_summary.md"
    write_csv(csv_path, rows)
    write_summary(md_path, rows, julia_root, mne_root)

    print(f"Compared {len(rows)} recordings")
    print(f"CSV: {csv_path}")
    print(f"Summary: {md_path}")


if __name__ == "__main__":
    main()
