"""Phase 5 — Epoching, baseline correction and AR for a single recording.

Default pilot: sub-M05 / ses-T2 / EC (eyesclosed).

Usage
-----
python scripts/run_phase5_m05.py
python scripts/run_phase5_m05.py --subject M07 --session T1 --condition EC
python scripts/run_phase5_m05.py --condition EO

The script reads whatever ICA-cleaned (or filtered) .npz cache Phase 4
produced and emits:

    results/subjects/sub-{id}/ses-{sess}/{task}/
    ├── cache/epochs_{cond}-epo.fif       ← native MNE, input for Phase 6
    └── tables/
        ├── epoch_summary_{cond}.json
        └── epoch_drops_{cond}.csv        ← only if rejections occurred
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from mne_brain import load_config
from mne_brain.processing import (
    apply_baseline,
    load_cleaned_raw,
    make_epochs,
    reject_artifacts,
    save_epoch_results,
)


def _resolve_dirs(cfg, subject: str, session: str, condition: str) -> tuple[Path, Path]:
    task = "eyesclosed" if condition == "EC" else "eyesopen"
    base = (
        Path(cfg.root)
        / str(cfg.paths.get("results", "results"))
        / "subjects"
        / f"sub-{subject}"
        / f"ses-{session}"
        / task
    )
    return base, base / "cache"


def _load_ica_suggestions(tables_dir: Path, condition: str) -> list[int]:
    """Read suggested rejected components from Phase 4 summary JSON."""
    summary_path = tables_dir / f"ica_summary_{condition}.json"
    if not summary_path.is_file():
        return []
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    return [int(c) for c in data.get("suggested_rejected_components", [])]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject", default="M05")
    parser.add_argument("--session", default="T2")
    parser.add_argument("--condition", default="EC")
    args = parser.parse_args()

    cfg = load_config()
    base_dir, cache_dir = _resolve_dirs(cfg, args.subject, args.session, args.condition)
    tables_dir = base_dir / "tables"

    # ── 1. Load ICA-cleaned signal as native mne.io.RawArray ──────────────
    rejected = _load_ica_suggestions(tables_dir, args.condition)
    raw = load_cleaned_raw(cache_dir, args.condition, rejected_components=rejected)
    print(f"Raw: {len(raw.ch_names)} ch × {raw.n_times} samples @ {raw.info['sfreq']} Hz")

    # ── 2. Segment → baseline → AR  (all native MNE) ─────────────────────
    epochs = make_epochs(raw, cfg)
    n_initial = len(epochs)
    print(f"Epochs created: {n_initial}  ({epochs.tmax - epochs.tmin + 1/epochs.info['sfreq']:.1f} s each)")

    epochs = apply_baseline(epochs, cfg)
    print("Baseline applied: mean over entire epoch")

    epochs = reject_artifacts(epochs, cfg)

    # ── 3. Save outputs ──────────────────────────────────────────────────
    summary = save_epoch_results(epochs, n_initial, base_dir, args.condition)

    print(f"\n── Phase 5 complete ────────────────────────────────────")
    print(f"  Valid epochs   : {summary['n_valid']} / {summary['n_initial']}")
    print(f"  Rejected       : {summary['n_rejected']}  ({summary['rejection_rate']:.1%})")
    print(f"  Cache          : {cache_dir / f'epochs_{args.condition}-epo.fif'}")
    print(f"  Summary        : {tables_dir / f'epoch_summary_{args.condition}.json'}")


if __name__ == "__main__":
    main()
