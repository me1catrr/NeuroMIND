"""Full mne_brain pipeline (Phases 3–7) — single subject or batch mode.

Phases
------
3  QC + filtering (BrainVision → filtered_{cond}.npz)
4  ICA (filtered_{cond}.npz → ica_{cond}-ica.fif + cleaned_{cond}.npz)
5  Epoching + baseline + AR (cleaned_{cond}.npz → epochs_{cond}-epo.fif)
6  PSD spectral analysis (epochs → band_power, spectrum PNG)
7  wPLI connectivity (epochs → wPLI matrices + heatmap PNGs)

Outputs (per condition, in results/subjects/sub-{id}/ses-{sess}/{task}/)
-------
cache/   filtered_*.npz · cleaned_*.npz · *-ica.fif · *-epo.fif · *.npz
tables/  qc_*.csv · ica_*.csv · ica_summary_*.json · psd_*.csv · band_power_*.csv ·
         wpli_*.csv · connectivity_summary_*.json
figures/ psd_spectrum_*.png · wpli_heatmap_*_*.png · ica/ica_topomaps_*.png

Usage — single subject
----------------------
python scripts/run_full_pipeline.py --subject M05 --session T2
python scripts/run_full_pipeline.py --subject M07 --session T1 --condition EC
python scripts/run_full_pipeline.py --subject M20 --session T2 --from-phase 5

Usage — batch (multiple subjects in sequence)
---------------------------------------------
python scripts/run_full_pipeline.py --all --skip-done
python scripts/run_full_pipeline.py --subjects M07 M20 MC10
python scripts/run_full_pipeline.py --all --groups HC --sessions T1 --conditions EC
python scripts/run_full_pipeline.py --all --from-phase 7 --skip-done   # only wPLI
python scripts/run_full_pipeline.py --all --limit 5 --dry-run

Batch behaviour
---------------
- One row per (subject, session, condition) appended to
  `results/batch_pipeline_log.csv` (or --log-path)
- A subject failure does NOT stop the batch — the next recording continues
- `--skip-done` checks for `tables/connectivity_summary_{COND}.json` (Phase 7
  artifact). If present and `--from-phase <= 7`, the recording is skipped.
"""

from __future__ import annotations

import argparse
import csv
import gc
import json
import re
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import mne
import numpy as np
import pandas as pd

from mne_brain import load_config
from mne_brain.bids import load_eeg_bids
from mne_brain.common.types import RecordingMeta
from mne_brain.connectivity import compute_wpli, connectivity_matrices, save_connectivity_results
from mne_brain.ica import (
    apply_ica_rejection,
    compute_ica_features,
    evaluate_ica_components,
    load_filtered_recording,
    mne_artifact_suggestions,
    run_ica,
    suggest_rejected_components,
)
from mne_brain.preprocessing import describe_filter_chain, filter_recording, qc_report
from mne_brain.processing import (
    apply_baseline,
    load_cleaned_raw,
    make_epochs,
    reject_artifacts,
    save_epoch_exclusion,
    save_epoch_results,
)
from mne_brain.spectral import compute_psd, save_spectral_results


MNE_BRAIN_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BATCH_LOG = MNE_BRAIN_ROOT / "results" / "batch_pipeline_log.csv"

CONDITION_BY_TASK = {"eyesclosed": "EC", "eyesopen": "EO"}
METADATA_PATTERN = re.compile(
    r"^sub-(?P<subject>[^_]+)_ses-(?P<session>[^_]+)_task-(?P<task>[^_]+)_run-(?P<run>\d+)_eeg_metadata\.json$"
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────


def _resolve_dirs(cfg, subject: str, session: str, condition: str) -> tuple[Path, Path, Path]:
    task = "eyesclosed" if condition == "EC" else "eyesopen"
    base = (
        Path(cfg.root)
        / str(cfg.paths.get("results", "results"))
        / "subjects"
        / f"sub-{subject}"
        / f"ses-{session}"
        / task
    )
    return base, base / "cache", base / "tables"


def _section(title: str) -> None:
    print(f"\n{'─' * 55}")
    print(f"  {title}")
    print("─" * 55)


def _elapsed(t0: float) -> str:
    return f"{time.time() - t0:.1f}s"


def _exclusion_path(cfg, subject: str, session: str, condition: str) -> Path:
    _base_dir, _cache_dir, tables_dir = _resolve_dirs(cfg, subject, session, condition)
    return tables_dir / f"quality_exclusion_{condition}.json"


def is_excluded(cfg, subject: str, session: str, condition: str) -> bool:
    """True if Phase 5 already excluded this recording because no epochs survived."""
    return _exclusion_path(cfg, subject, session, condition).is_file()


def _exclude_recording(
    cfg,
    subject: str,
    session: str,
    condition: str,
    summary: dict,
    reason: str,
) -> None:
    base_dir, _cache_dir, _tables_dir = _resolve_dirs(cfg, subject, session, condition)
    save_epoch_exclusion(base_dir, condition, summary, reason)
    print(f"  [EXCLUDE] {reason}")


def _save_matrix_csv(path: Path, matrix: np.ndarray, index: list[str] | None = None) -> None:
    df = pd.DataFrame(matrix)
    if index is not None and len(index) == matrix.shape[0]:
        df.insert(0, "channel", index)
    df.to_csv(path, index=False)


def _save_ica_topomaps(ica, figures_dir: Path, condition: str) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    saved: list[str] = []
    n = int(ica.n_components_)
    for start in range(0, n, 10):
        picks = list(range(start, min(start + 10, n)))
        try:
            figs = ica.plot_components(picks=picks, show=False)
            if not isinstance(figs, list):
                figs = [figs]
            for idx, fig in enumerate(figs):
                name = f"ica_topomaps_{condition}_{start + idx:02d}.png"
                fig.savefig(figures_dir / name, dpi=150, bbox_inches="tight")
                plt.close(fig)
                saved.append(str(figures_dir / name))
        except Exception:
            continue
    return saved


# ─────────────────────────────────────────────────────────────────────────────
# Phase 3 — QC + Filtering
# ─────────────────────────────────────────────────────────────────────────────


def run_phase3(cfg, subject: str, session: str, condition: str) -> None:
    _section(f"Phase 3 — QC + Filtering  [{condition}]")
    t0 = time.time()
    base_dir, cache_dir, tables_dir = _resolve_dirs(cfg, subject, session, condition)
    cache_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    rec = load_eeg_bids(cfg, subject, session, condition, source="brainvision")
    print(f"  Loaded: {rec.data.shape[0]} ch × {rec.data.shape[1]} samples @ {rec.meta.fs} Hz")

    qc = qc_report(rec, cfg, output_dir=tables_dir)
    bad_ch = qc.loc[qc["is_bad"], "channel"].to_list()
    print(f"  QC bad channels: {bad_ch or 'none'}")

    filtered = filter_recording(rec, cfg)
    np.savez_compressed(
        cache_dir / f"filtered_{condition}.npz",
        data=filtered.data,
        times=filtered.times,
        fs=filtered.meta.fs,
        channel_names=np.asarray(filtered.meta.channel_names, dtype=object),
    )
    (tables_dir / f"filter_chain_{condition}.json").write_text(
        json.dumps(describe_filter_chain(cfg), indent=2), encoding="utf-8"
    )
    print(f"  Filtered cache saved  [{_elapsed(t0)}]")


# ─────────────────────────────────────────────────────────────────────────────
# Phase 4 — ICA
# ─────────────────────────────────────────────────────────────────────────────


def run_phase4(cfg, subject: str, session: str, condition: str) -> list[int]:
    """Returns suggested rejected component indices."""
    _section(f"Phase 4 — ICA  [{condition}]")
    t0 = time.time()
    base_dir, cache_dir, tables_dir = _resolve_dirs(cfg, subject, session, condition)
    figures_dir = base_dir / "figures" / "ica"
    figures_dir.mkdir(parents=True, exist_ok=True)

    filtered_path = cache_dir / f"filtered_{condition}.npz"
    meta = RecordingMeta(
        subject_id=subject, session_id=session, condition=condition, run=1,
        fs=float(cfg.recording.get("fs", 500.0)),
        n_channels=0, channel_names=[], channel_positions=None,
        bids_path=str(filtered_path),
    )
    rec = load_filtered_recording(filtered_path, meta)
    fitted = run_ica(rec, cfg)

    features = compute_ica_features(
        fitted.component_maps, fitted.sources,
        rec.meta.fs, rec.meta.channel_names,
    )
    evaluated = evaluate_ica_components(features, artifact_thresh=1.5)
    custom_suggestions = suggest_rejected_components(evaluated)
    mne_sugg = mne_artifact_suggestions(fitted.ica, fitted.raw)

    evaluated.to_csv(tables_dir / f"ica_component_features_{condition}.csv", index=False)
    _save_matrix_csv(tables_dir / f"ica_component_maps_{condition}.csv",
                     fitted.component_maps, rec.meta.channel_names)
    _save_matrix_csv(tables_dir / f"ica_mixing_matrix_{condition}.csv",
                     fitted.result.mixing_matrix, rec.meta.channel_names)
    _save_matrix_csv(tables_dir / f"ica_unmixing_matrix_{condition}.csv",
                     fitted.result.unmixing_matrix)
    np.savez_compressed(
        cache_dir / f"ica_activations_{condition}.npz",
        activations=fitted.sources, times=rec.times, fs=rec.meta.fs,
    )
    fitted.ica.save(cache_dir / f"ica_{condition}-ica.fif", overwrite=True)

    cleaned = apply_ica_rejection(rec, fitted, custom_suggestions)
    np.savez_compressed(
        cache_dir / f"cleaned_{condition}.npz",
        data=cleaned.data, times=cleaned.times, fs=cleaned.meta.fs,
        channel_names=np.asarray(cleaned.meta.channel_names, dtype=object),
        rejected_components=np.asarray(custom_suggestions, dtype=int),
    )

    preview_path = tables_dir / f"signal_preview_{condition}.json"
    if preview_path.exists():
        preview = json.loads(preview_path.read_text(encoding="utf-8"))
        n_prev = len(preview.get("times", []))
        n_avail = cleaned.data.shape[1]
        preview["post_ica_uv"] = cleaned.data[:, :min(n_prev, n_avail)].round(6).tolist()
        preview_path.write_text(json.dumps(preview), encoding="utf-8")

    topomap_files = _save_ica_topomaps(fitted.ica, figures_dir, condition)

    summary = {
        "subject": subject, "session": session, "condition": condition,
        "method": str(cfg.ica.get("method", "fastica")),
        "library": "MNE-Python",
        "n_channels": rec.meta.n_channels,
        "n_components": int(fitted.sources.shape[0]),
        "random_seed": int(cfg.ica.get("random_seed", 42)),
        "suggested_rejected_components": custom_suggestions,
        "mne_suggestions": mne_sugg,
        "topomap_files": topomap_files,
    }
    (tables_dir / f"ica_summary_{condition}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"  ICA: {fitted.sources.shape[0]} components, {len(custom_suggestions)} suggested for rejection")
    print(f"  Suggested: {custom_suggestions}  [{_elapsed(t0)}]")
    return custom_suggestions


# ─────────────────────────────────────────────────────────────────────────────
# Phase 5 — Epoching + Baseline + AR
# ─────────────────────────────────────────────────────────────────────────────


def run_phase5(cfg, subject: str, session: str, condition: str) -> bool:
    _section(f"Phase 5 — Epoching + AR  [{condition}]")
    t0 = time.time()
    base_dir, cache_dir, tables_dir = _resolve_dirs(cfg, subject, session, condition)

    ica_summary_path = tables_dir / f"ica_summary_{condition}.json"
    rejected = []
    if ica_summary_path.is_file():
        data = json.loads(ica_summary_path.read_text(encoding="utf-8"))
        rejected = [int(c) for c in data.get("suggested_rejected_components", [])]

    raw = load_cleaned_raw(cache_dir, condition, rejected_components=rejected)
    print(f"  Raw: {len(raw.ch_names)} ch × {raw.n_times} samples @ {raw.info['sfreq']} Hz")

    epochs = make_epochs(raw, cfg)
    n_initial = len(epochs)
    print(f"  Epochs: {n_initial} × {len(epochs.times)} samples")

    epochs = apply_baseline(epochs, cfg)
    epochs = reject_artifacts(epochs, cfg)
    n_valid = len(epochs)
    print(f"  After AR: {n_valid}/{n_initial} kept  ({n_initial - n_valid} rejected)  [{_elapsed(t0)}]")

    summary = save_epoch_results(epochs, n_initial, base_dir, condition)
    if n_valid == 0:
        reason = (
            "0 valid epochs after artifact rejection; downstream PSD/wPLI skipped. "
            f"Inspect tables/epoch_drops_{condition}.csv and relax AR or review bad channels."
        )
        _exclude_recording(cfg, subject, session, condition, summary, reason)
        return False
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Phase 6 — PSD Spectral Analysis
# ─────────────────────────────────────────────────────────────────────────────


def run_phase6(cfg, subject: str, session: str, condition: str) -> bool:
    _section(f"Phase 6 — Spectral PSD  [{condition}]")
    t0 = time.time()
    base_dir, cache_dir, tables_dir = _resolve_dirs(cfg, subject, session, condition)

    fif_path = cache_dir / f"epochs_{condition}-epo.fif"
    epochs = mne.read_epochs(str(fif_path), preload=True, verbose=False)
    print(f"  Epochs: {len(epochs)} × {len(epochs.times)} samples, {len(epochs.ch_names)} ch")
    if len(epochs) == 0:
        summary_path = tables_dir / f"epoch_summary_{condition}.json"
        summary = {}
        if summary_path.is_file():
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        reason = (
            "0 valid epochs in saved FIF; downstream PSD/wPLI skipped. "
            f"Inspect tables/epoch_drops_{condition}.csv and rerun from Phase 5 after QC changes."
        )
        _exclude_recording(cfg, subject, session, condition, summary, reason)
        return False

    spectrum = compute_psd(epochs, cfg)
    save_spectral_results(spectrum, cfg, base_dir, condition)
    print(f"  PSD computed, figures + tables saved  [{_elapsed(t0)}]")
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Phase 7 — wPLI Connectivity
# ─────────────────────────────────────────────────────────────────────────────


def run_phase7(cfg, subject: str, session: str, condition: str) -> dict | None:
    _section(f"Phase 7 — wPLI Connectivity  [{condition}]")
    t0 = time.time()
    base_dir, cache_dir, tables_dir = _resolve_dirs(cfg, subject, session, condition)

    fif_path = cache_dir / f"epochs_{condition}-epo.fif"
    epochs = mne.read_epochs(str(fif_path), preload=True, verbose=False)
    print(f"  Epochs: {len(epochs)} × {len(epochs.ch_names)} ch @ {epochs.info['sfreq']} Hz")
    if len(epochs) == 0:
        summary_path = tables_dir / f"epoch_summary_{condition}.json"
        summary = {}
        if summary_path.is_file():
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        reason = (
            "0 valid epochs in saved FIF; wPLI skipped. "
            f"Inspect tables/epoch_drops_{condition}.csv and rerun from Phase 5 after QC changes."
        )
        _exclude_recording(cfg, subject, session, condition, summary, reason)
        return None

    con, band_names = compute_wpli(epochs, cfg)
    matrices = connectivity_matrices(con, band_names)
    summary = save_connectivity_results(con, band_names, cfg, base_dir, condition)

    print(f"  wPLI computed for {len(band_names)} bands  [{_elapsed(t0)}]")
    print(f"\n  Mean wPLI per band:")
    for band, stats in summary["band_summary"].items():
        lo, hi = cfg.bands[band]
        print(f"    {band:<12} [{lo:4.1f}–{hi:5.1f} Hz]  "
              f"mean={stats['mean_wpli']:.4f}  max={stats['max_wpli']:.4f}")
    return summary


# ─────────────────────────────────────────────────────────────────────────────
# Per-recording orchestration
# ─────────────────────────────────────────────────────────────────────────────


def run_one_recording(
    cfg, subject: str, session: str, condition: str, *, from_phase: int,
) -> str:
    """Run phases [from_phase..7] for one (subject, session, condition)."""
    base_dir, _cache_dir, _tables_dir = _resolve_dirs(cfg, subject, session, condition)
    base_dir.mkdir(parents=True, exist_ok=True)
    if from_phase <= 3:
        run_phase3(cfg, subject, session, condition)
    if from_phase <= 4:
        run_phase4(cfg, subject, session, condition)
    if from_phase <= 5:
        if not run_phase5(cfg, subject, session, condition):
            return "excluded"
    if from_phase <= 6:
        if not run_phase6(cfg, subject, session, condition):
            return "excluded"
    if from_phase <= 7:
        if run_phase7(cfg, subject, session, condition) is None:
            return "excluded"
    return "ok"


def is_done(cfg, subject: str, session: str, condition: str) -> bool:
    """True if Phase 7 already produced its summary for this recording."""
    _base_dir, _cache_dir, tables_dir = _resolve_dirs(cfg, subject, session, condition)
    return (tables_dir / f"connectivity_summary_{condition}.json").is_file()


# ─────────────────────────────────────────────────────────────────────────────
# Discovery & filtering (batch mode)
# ─────────────────────────────────────────────────────────────────────────────


def discover_recordings(cfg) -> pd.DataFrame:
    """Scan `data/BIDS/raw/` and return one row per metadata JSON."""
    bids_root = Path(cfg.root) / "data" / "BIDS"
    raw_dir = bids_root / "raw"
    if not raw_dir.is_dir():
        raise FileNotFoundError(f"BIDS raw dir not found: {raw_dir}")

    rows: list[dict[str, object]] = []
    for path in sorted(raw_dir.glob("sub-*_eeg_metadata.json")):
        m = METADATA_PATTERN.match(path.name)
        if not m:
            continue
        task = m.group("task")
        condition = CONDITION_BY_TASK.get(task)
        if condition is None:
            continue
        try:
            with path.open("r", encoding="utf-8") as fh:
                meta = json.load(fh)
        except Exception:
            meta = {}
        rows.append(
            {
                "subject": m.group("subject"),
                "session": m.group("session"),
                "task": task,
                "condition": condition,
                "run": int(m.group("run")),
                "group": meta.get("group", ""),
                "json_path": str(path),
            }
        )
    return pd.DataFrame(rows)


def filter_recordings(
    df: pd.DataFrame, *,
    groups: list[str] | None,
    sessions: list[str] | None,
    conditions: list[str] | None,
    subjects: list[str] | None,
    start_from: str | None,
    limit: int | None,
) -> pd.DataFrame:
    out = df.copy()
    if groups:
        out = out[out["group"].isin(groups)]
    if sessions:
        out = out[out["session"].isin(sessions)]
    if conditions:
        out = out[out["condition"].isin(conditions)]
    if subjects:
        out = out[out["subject"].isin(subjects)]
    if start_from:
        out = out[out["subject"] >= start_from]
    out = out.sort_values(["subject", "session", "condition"]).reset_index(drop=True)

    if limit is not None and limit > 0:
        keep_keys: set[tuple[str, str]] = set()
        kept_rows: list[int] = []
        for idx, row in out.iterrows():
            key = (row["subject"], row["session"])
            if key not in keep_keys:
                if len(keep_keys) >= limit:
                    continue
                keep_keys.add(key)
            kept_rows.append(idx)
        out = out.loc[kept_rows].reset_index(drop=True)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Logging (batch mode)
# ─────────────────────────────────────────────────────────────────────────────

LOG_FIELDS = [
    "timestamp",
    "subject",
    "session",
    "task",
    "condition",
    "group",
    "from_phase",
    "status",
    "elapsed_s",
    "error_type",
    "error_message",
]


def append_log_row(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.is_file()
    with path.open("a", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=LOG_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow({k: row.get(k, "") for k in LOG_FIELDS})


def _print_batch_summary(counts: dict[str, int], t_total: float,
                         *, log_path: Path | None = None) -> None:
    total_elapsed = time.time() - t_total
    print(f"\n{'═' * 60}")
    print(f"  Batch summary  [{total_elapsed/60:.1f} min total]")
    print(f"{'═' * 60}")
    print(f"  ok      : {counts.get('ok', 0)}")
    print(f"  excluded: {counts.get('excluded', 0)}")
    print(f"  skipped : {counts.get('skip', 0)}")
    print(f"  errors  : {counts.get('error', 0)}")
    if log_path is not None:
        print(f"\n  Log written to: {log_path}")


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry points
# ─────────────────────────────────────────────────────────────────────────────


def run_single(args) -> None:
    """Process a single (subject, session) for one or both conditions."""
    cfg = load_config()
    conditions = [args.condition] if args.condition else ["EC", "EO"]
    t_total = time.time()

    print(f"\n{'═' * 55}")
    print(f"  mne_brain Full Pipeline — sub-{args.subject} / ses-{args.session}")
    print(f"  Conditions: {conditions}")
    print(f"  Starting from Phase {args.from_phase}")
    print(f"{'═' * 55}")

    for cond in conditions:
        base_dir, _, _ = _resolve_dirs(cfg, args.subject, args.session, cond)
        base_dir.mkdir(parents=True, exist_ok=True)
        try:
            status = run_one_recording(
                cfg, args.subject, args.session, cond, from_phase=args.from_phase
            )
            if status == "excluded":
                print(f"\n  Condition {cond} excluded after epoch QC")
        except Exception as exc:
            print(f"\n  ✗ ERROR in condition {cond}: {exc}")
            traceback.print_exc()
            continue

    results_base = (
        Path(cfg.root) / str(cfg.paths.get("results", "results"))
        / "subjects" / f"sub-{args.subject}" / f"ses-{args.session}"
    )
    print(f"\n{'═' * 55}")
    print(f"  Pipeline complete  [{time.time() - t_total:.1f}s total]")
    print(f"{'═' * 55}")
    for cond in conditions:
        task = "eyesclosed" if cond == "EC" else "eyesopen"
        out_dir = results_base / task
        print(f"\n  [{cond}] → {out_dir}")
        if (out_dir / "figures").exists():
            figures = list((out_dir / "figures").glob("*.png"))
            print(f"    figures/   : {len(figures)} PNG(s)")
        if (out_dir / "tables").exists():
            tables = list((out_dir / "tables").glob("*.csv")) + list((out_dir / "tables").glob("*.json"))
            print(f"    tables/    : {len(tables)} file(s)")
        if (out_dir / "cache").exists():
            caches = list((out_dir / "cache").glob("*"))
            print(f"    cache/     : {len(caches)} file(s)")


def run_batch(args) -> None:
    """Process many recordings in sequence, with per-row CSV log."""
    mne.set_log_level("ERROR")
    cfg = load_config()
    log_path = Path(args.log_path)

    print(f"\n{'═' * 60}")
    print(f"  mne_brain — Batch pipeline")
    print(f"{'═' * 60}")
    print(f"  Filters:")
    print(f"    groups       = {args.groups or 'all'}")
    print(f"    sessions     = {args.sessions or 'all'}")
    print(f"    conditions   = {args.conditions or ([args.condition] if args.condition else 'all')}")
    print(f"    subjects     = {args.subjects or 'all'}")
    print(f"    start_from   = {args.start_from or '(none)'}")
    print(f"    limit pairs  = {args.limit or '(none)'}")
    print(f"  Pipeline:")
    print(f"    from-phase   = {args.from_phase}")
    print(f"    skip-done    = {args.skip_done}")
    print(f"    dry-run      = {args.dry_run}")
    print(f"  Log file: {log_path}")

    all_recs = discover_recordings(cfg)
    if all_recs.empty:
        print("\n  No metadata JSON files found in data/BIDS/raw. Run "
              "`python scripts/build_bids_full.py` first.")
        return

    conditions_filter = args.conditions
    if not conditions_filter and args.condition:
        conditions_filter = [args.condition]

    selected = filter_recordings(
        all_recs,
        groups=args.groups,
        sessions=args.sessions,
        conditions=conditions_filter,
        subjects=args.subjects,
        start_from=args.start_from,
        limit=args.limit,
    )

    pairs = (
        selected[["subject", "session"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    print(f"\n  Discovered : {len(all_recs)} recording(s)")
    print(f"  Selected   : {len(selected)} recording(s) across "
          f"{len(pairs)} (subject, session) pair(s)")

    if args.dry_run:
        print("\n  Dry-run preview:")
        for _, row in selected.iterrows():
            print(f"    sub-{row['subject']} ses-{row['session']} {row['condition']} "
                  f"[{row['group']}]")
        return

    counts = {"ok": 0, "excluded": 0, "skip": 0, "error": 0}
    t_total = time.time()

    for pair_idx, (_, pair) in enumerate(pairs.iterrows(), start=1):
        subject = str(pair["subject"])
        session = str(pair["session"])
        rec_for_pair = selected[
            (selected["subject"] == subject) & (selected["session"] == session)
        ].sort_values("condition")
        n_pair = len(rec_for_pair)

        print(f"\n{'═' * 60}")
        print(f"  [{pair_idx}/{len(pairs)}] sub-{subject} ses-{session} "
              f"— {n_pair} condition(s): "
              f"{', '.join(rec_for_pair['condition'].tolist())}")
        print(f"{'═' * 60}")

        for _, rec in rec_for_pair.iterrows():
            condition = str(rec["condition"])
            group = str(rec["group"])
            t0 = time.time()
            row: dict[str, object] = {
                "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "subject": subject,
                "session": session,
                "task": rec["task"],
                "condition": condition,
                "group": group,
                "from_phase": args.from_phase,
            }

            if (
                args.skip_done
                and args.from_phase <= 7
                and (
                    is_done(cfg, subject, session, condition)
                    or is_excluded(cfg, subject, session, condition)
                )
            ):
                row.update(status="skip", elapsed_s=f"{time.time() - t0:.2f}")
                append_log_row(log_path, row)
                counts["skip"] += 1
                reason = "connectivity_summary present"
                if is_excluded(cfg, subject, session, condition):
                    reason = "quality_exclusion present"
                print(f"  [SKIP] {condition} — {reason}")
                continue

            try:
                status = run_one_recording(
                    cfg, subject, session, condition, from_phase=args.from_phase
                )
                elapsed = time.time() - t0
                if status == "excluded":
                    row.update(status="excluded", elapsed_s=f"{elapsed:.2f}")
                    counts["excluded"] += 1
                    print(f"  [EXC] {condition} excluded after {elapsed:.1f}s")
                else:
                    row.update(status="ok", elapsed_s=f"{elapsed:.2f}")
                    counts["ok"] += 1
                    print(f"  [OK ] {condition} done in {elapsed:.1f}s")
            except KeyboardInterrupt:
                row.update(
                    status="error",
                    elapsed_s=f"{time.time() - t0:.2f}",
                    error_type="KeyboardInterrupt",
                    error_message="user interrupted",
                )
                append_log_row(log_path, row)
                print("\n  Interrupted by user. Stopping batch.")
                _print_batch_summary(counts, t_total)
                sys.exit(130)
            except Exception as exc:
                elapsed = time.time() - t0
                tb_last = traceback.format_exc().strip().splitlines()[-1]
                row.update(
                    status="error",
                    elapsed_s=f"{elapsed:.2f}",
                    error_type=type(exc).__name__,
                    error_message=tb_last,
                )
                counts["error"] += 1
                print(f"  [ERR] {condition} failed after {elapsed:.1f}s: {tb_last}")
                traceback.print_exc()

            append_log_row(log_path, row)
            gc.collect()

    _print_batch_summary(counts, t_total, log_path=log_path)


# ─────────────────────────────────────────────────────────────────────────────
# Argument parser
# ─────────────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Single-subject options
    single = parser.add_argument_group("single subject (default mode)")
    single.add_argument("--subject", help="BIDS subject id (e.g. M05, MC10)")
    single.add_argument("--session", help="Session id (T1 or T2)")
    single.add_argument(
        "--condition", choices=["EC", "EO"], default=None,
        help="EC, EO, or omit for both",
    )

    # Batch options
    batch = parser.add_argument_group("batch mode")
    batch.add_argument(
        "--all", action="store_true",
        help="Process every recording in data/BIDS/raw (with optional filters)",
    )
    batch.add_argument(
        "--subjects", nargs="*",
        help="Restrict to specific bids ids (e.g. M07 MC10). "
             "Implies batch mode when more than one is given.",
    )
    batch.add_argument(
        "--groups", nargs="*", choices=["MS", "HC"],
        help="Restrict to one or more groups",
    )
    batch.add_argument(
        "--sessions", nargs="*", choices=["T1", "T2"],
        help="Restrict to one or more sessions",
    )
    batch.add_argument(
        "--conditions", nargs="*", choices=["EC", "EO"],
        help="Restrict to one or more conditions",
    )
    batch.add_argument(
        "--start-from", default=None,
        help="Only process subjects whose bids id is >= START_FROM "
             "(lexicographic)",
    )
    batch.add_argument(
        "--limit", type=int, default=None,
        help="Limit to the first N (subject, session) pairs after filtering",
    )
    batch.add_argument(
        "--skip-done", action="store_true",
        help="Skip recordings that already have a Phase 7 connectivity_summary",
    )
    batch.add_argument(
        "--dry-run", action="store_true",
        help="List what would be processed without running anything",
    )
    batch.add_argument(
        "--log-path", default=str(DEFAULT_BATCH_LOG),
        help="Path to append per-recording log rows (batch mode)",
    )

    # Common
    common = parser.add_argument_group("common")
    common.add_argument(
        "--from-phase", type=int, default=3, choices=[3, 4, 5, 6, 7],
        help="Start each recording from this phase",
    )

    return parser


def is_batch_mode(args) -> bool:
    """Decide whether the user wants batch behaviour."""
    if args.all:
        return True
    if args.subjects and len(args.subjects) > 1:
        return True
    if args.groups or args.sessions or args.conditions:
        return True
    if args.start_from is not None or args.limit is not None:
        return True
    if args.dry_run or args.skip_done:
        return True
    return False


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if is_batch_mode(args):
        run_batch(args)
        return

    if not args.subject or not args.session:
        if args.subjects and len(args.subjects) == 1 and args.session:
            args.subject = args.subjects[0]
        else:
            parser.error(
                "Single-subject mode needs --subject and --session.  "
                "For all subjects use --all, for a subset use --subjects ... or "
                "filters such as --groups/--sessions/--conditions."
            )

    run_single(args)


if __name__ == "__main__":
    main()
