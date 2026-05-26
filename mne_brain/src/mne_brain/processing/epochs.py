"""Phase 5 — Epoching, baseline correction and artifact rejection.

Uses native MNE objects throughout:
  mne.make_fixed_length_epochs  → segmentation
  epochs.apply_baseline         → baseline correction
  epochs.drop_bad               → amplitude-based AR

Mirrors NeuroMIND/src/segmentation/Epochs.jl (profile "eeg_julia"):
  - 1 s fixed-length epochs, no overlap
  - mean baseline over the whole epoch
  - ±70 µV amplitude threshold

Threshold note
--------------
MNE's drop_bad() uses **peak-to-peak** (max − min) in Volts.
NeuroMIND checks absolute amplitude (|sample| > 70 µV).
For a mean-baseline-corrected signal centred at 0 the equivalence is:

    NeuroMIND ±70 µV  ≡  MNE p2p 140 µV = 140e-6 V

This file uses the p2p conversion so results are as close as
possible to NeuroMIND's rejection rate.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import mne
import numpy as np
import pandas as pd

from mne_brain.common.types import PipelineConfig


# ── Segment ────────────────────────────────────────────────────────────────────

def make_epochs(raw: mne.io.BaseRaw, cfg: PipelineConfig) -> mne.Epochs:
    """Segment a continuous Raw into fixed-length epochs.

    Parameters
    ----------
    raw:
        MNE Raw object (typically the ICA-cleaned signal in Volts).
    cfg:
        Pipeline configuration.  Profile "eeg_julia" forces 1 s / no overlap.

    Returns
    -------
    mne.Epochs  (preloaded, no baseline applied yet)
    """
    seg = cfg.segmentation
    profile = str(seg.get("profile", "default"))

    if profile == "eeg_julia":
        duration = 1.0
        overlap = 0.0
    else:
        duration = float(seg.get("epoch_length_s", 1.0))
        overlap = float(seg.get("epoch_overlap", 0.0))

    epochs = mne.make_fixed_length_epochs(
        raw,
        duration=duration,
        overlap=overlap,
        preload=True,
        reject_by_annotation=True,
        verbose=False,
    )
    return epochs


# ── Baseline ───────────────────────────────────────────────────────────────────

def apply_baseline(epochs: mne.Epochs, cfg: PipelineConfig) -> mne.Epochs:
    """Apply mean-baseline correction using epochs.apply_baseline().

    ``(None, None)`` instructs MNE to subtract the mean of the *entire epoch*,
    matching NeuroMIND's ``method = "mean"``.
    """
    bl_cfg = cfg.baseline
    if not bl_cfg.get("apply", True):
        return epochs
    epochs.apply_baseline((None, None), verbose=False)
    return epochs


# ── Artifact rejection ─────────────────────────────────────────────────────────

def reject_artifacts(epochs: mne.Epochs, cfg: PipelineConfig) -> mne.Epochs:
    """Drop epochs that exceed the amplitude threshold via epochs.drop_bad().

    Converts NeuroMIND's ±amplitude threshold to MNE's peak-to-peak criterion:

        ptp_V = 2 × amp_uV × 1e-6

    Profile "eeg_julia"  →  ±70 µV  →  reject = {'eeg': 140e-6}
    Profile "default"    →  ±100 µV →  reject = {'eeg': 200e-6}
    """
    ar = cfg.artifact_rejection
    if not ar.get("enabled", True):
        return epochs

    profile = str(ar.get("profile", "default"))

    if profile == "eeg_julia":
        amp_uv = float(ar.get("max_amplitude_uv", 70.0))
    else:
        amp_uv = float(ar.get("amplitude_threshold_uv", 100.0))

    ptp_v = 2.0 * amp_uv * 1e-6          # ±amp_uv → peak-to-peak
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All epochs were dropped!", category=RuntimeWarning)
        epochs.drop_bad(reject={"eeg": ptp_v}, verbose=False)
    return epochs


# ── I/O helpers ────────────────────────────────────────────────────────────────

def load_cleaned_raw(
    cache_dir: Path,
    condition: str,
    *,
    ica_fif_path: Path | None = None,
    rejected_components: list[int] | None = None,
) -> mne.io.RawArray:
    """Return a native MNE RawArray of the ICA-cleaned signal.

    Priority:
    1. ``cleaned_{condition}.npz``  — already cleaned by run_phase4
    2. ``filtered_{condition}.npz`` + apply saved ``ica_{condition}-ica.fif``
    3. ``filtered_{condition}.npz`` alone  (ICA not yet applied)
    """
    cleaned_path = cache_dir / f"cleaned_{condition}.npz"
    filtered_path = cache_dir / f"filtered_{condition}.npz"

    if cleaned_path.is_file():
        return _npz_to_raw(cleaned_path)

    raw = _npz_to_raw(filtered_path)

    fif = ica_fif_path or (cache_dir / f"ica_{condition}-ica.fif")
    if fif.is_file() and rejected_components is not None:
        ica = mne.preprocessing.read_ica(str(fif), verbose=False)
        ica.exclude = list(rejected_components)
        ica.apply(raw, verbose=False)
        print(f"  ICA applied: excluded components {rejected_components}")

    return raw


def _npz_to_raw(npz_path: Path) -> mne.io.RawArray:
    """Load a filtered/cleaned .npz cache as a native mne.io.RawArray (Volts)."""
    with np.load(npz_path, allow_pickle=True) as npz:
        data_uv = np.asarray(npz["data"], dtype=np.float64)
        fs = float(npz["fs"])
        ch_names = [str(c) for c in npz["channel_names"].tolist()]

    info = mne.create_info(
        ch_names=ch_names,
        sfreq=fs,
        ch_types=["eeg"] * len(ch_names),
    )
    raw = mne.io.RawArray(data_uv * 1e-6, info, verbose=False)   # µV → V
    try:
        raw.set_montage("standard_1020", on_missing="ignore", verbose=False)
    except Exception:
        pass
    return raw


# ── QC summary & persistence ───────────────────────────────────────────────────

def epoch_summary(epochs: mne.Epochs, n_initial: int) -> dict:
    """Build a QC summary dict from a processed Epochs object."""
    n_valid = len(epochs)
    n_rejected = n_initial - n_valid

    drop_rows = [
        {"epoch_idx": i, "reasons": list(log)}
        for i, log in enumerate(epochs.drop_log)
        if log
    ]

    return {
        "n_initial": n_initial,
        "n_valid": n_valid,
        "n_rejected": n_rejected,
        "rejection_rate": round(n_rejected / max(n_initial, 1), 4),
        "epoch_duration_s": round(
            float(epochs.tmax - epochs.tmin + 1.0 / epochs.info["sfreq"]), 6
        ),
        "sfreq_hz": float(epochs.info["sfreq"]),
        "n_channels": len(epochs.ch_names),
        "drop_log_preview": drop_rows[:20],
    }


def save_epoch_results(
    epochs: mne.Epochs,
    n_initial: int,
    out_dir: Path,
    condition: str,
    *,
    save_fif: bool = True,
) -> dict:
    """Persist epoch QC outputs and optionally the native MNE .fif cache.

    Outputs
    -------
    tables/epoch_summary_{condition}.json
    tables/epoch_drops_{condition}.csv   (only if rejections occurred)
    cache/epochs_{condition}-epo.fif     (native MNE, used by Phase 6)
    """
    out_dir = Path(out_dir)
    tables_dir = out_dir / "tables"
    cache_dir = out_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    summary = epoch_summary(epochs, n_initial)

    (tables_dir / f"epoch_summary_{condition}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    drop_rows = [
        {"epoch_idx": i, "reason": ", ".join(log)}
        for i, log in enumerate(epochs.drop_log)
        if log
    ]
    if drop_rows:
        pd.DataFrame(drop_rows).to_csv(
            tables_dir / f"epoch_drops_{condition}.csv", index=False
        )

    if save_fif and len(epochs) > 0:
        fif_path = cache_dir / f"epochs_{condition}-epo.fif"
        epochs.save(str(fif_path), overwrite=True, verbose=False)
        print(f"  Epochs saved → {fif_path.name}")
    elif save_fif:
        print("  Epochs FIF not saved because no epochs survived AR")

    return summary


def save_epoch_exclusion(
    out_dir: Path,
    condition: str,
    summary: dict,
    reason: str,
) -> dict:
    """Persist a machine-readable marker for recordings excluded after epoch QC."""
    out_dir = Path(out_dir)
    tables_dir = out_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    exclusion = {
        "condition": condition,
        "reason": reason,
        "n_initial": int(summary.get("n_initial", 0)),
        "n_valid": int(summary.get("n_valid", 0)),
        "n_rejected": int(summary.get("n_rejected", 0)),
        "rejection_rate": float(summary.get("rejection_rate", 0.0)),
        "drop_log_preview": summary.get("drop_log_preview", []),
    }
    (tables_dir / f"quality_exclusion_{condition}.json").write_text(
        json.dumps(exclusion, indent=2), encoding="utf-8"
    )
    return exclusion
