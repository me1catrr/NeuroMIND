"""Initial EEG quality-control utilities mirrored from NeuroMIND."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from mne_brain.common.types import EEGRecording, PipelineConfig


def compute_channel_stats(rec: EEGRecording) -> pd.DataFrame:
    """Compute per-channel QC metrics in microvolts.

    Mirrors `NeuroMIND/src/qc/QualityControl.jl::compute_channel_stats`.
    """
    data = np.asarray(rec.data, dtype=np.float64)
    mean_uv = data.mean(axis=1)
    rms_uv = np.sqrt(np.mean(data**2, axis=1))
    std_uv = data.std(axis=1, ddof=1)
    range_uv = data.max(axis=1) - data.min(axis=1)
    rms_zscore = (rms_uv - rms_uv.mean()) / (rms_uv.std(ddof=1) + np.finfo(float).eps)

    return pd.DataFrame(
        {
            "channel": rec.meta.channel_names,
            "mean_uv": np.round(mean_uv, 3),
            "rms_uv": np.round(rms_uv, 3),
            "std_uv": np.round(std_uv, 3),
            "range_uv": np.round(range_uv, 3),
            "rms_zscore": np.round(rms_zscore, 3),
        }
    )


def flag_bad_channels(rec: EEGRecording, *, z_threshold: float = 3.0) -> list[str]:
    """Flag channels with absolute RMS z-score above threshold."""
    stats = compute_channel_stats(rec)
    bad = stats.loc[stats["rms_zscore"].abs() > z_threshold, "channel"]
    return [str(ch) for ch in bad.to_list()]


def qc_report(
    rec: EEGRecording,
    cfg: PipelineConfig | None = None,
    *,
    output_dir: str | Path | None = None,
    z_threshold: float = 3.0,
) -> pd.DataFrame:
    """Generate a QC table and optionally save it as CSV."""
    stats = compute_channel_stats(rec)
    bad_channels = set(flag_bad_channels(rec, z_threshold=z_threshold))
    stats["is_bad"] = stats["channel"].isin(bad_channels)

    if output_dir is None and cfg is not None:
        output_dir = (
            Path(cfg.root)
            / str(cfg.paths.get("results", "results"))
            / "subjects"
            / f"sub-{rec.meta.subject_id}"
            / f"ses-{rec.meta.session_id}"
            / _task_from_condition(rec.meta.condition)
            / "tables"
        )

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        stats.to_csv(out / f"qc_channels_{rec.meta.condition}.csv", index=False)

    return stats


def _task_from_condition(condition: str) -> str:
    if condition == "EC":
        return "eyesclosed"
    if condition == "EO":
        return "eyesopen"
    return condition

