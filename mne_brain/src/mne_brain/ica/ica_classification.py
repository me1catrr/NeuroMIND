"""ICA component feature extraction and artifact suggestions."""

from __future__ import annotations

import numpy as np
import pandas as pd


def compute_ica_features(
    component_maps: np.ndarray,
    sources: np.ndarray,
    fs: float,
    channel_names: list[str],
) -> pd.DataFrame:
    """Compute NeuroMIND-style ICA features for each component."""
    n_channels, n_components = component_maps.shape
    frontal_idx = [
        idx
        for idx, ch in enumerate(channel_names)
        if ch.upper().startswith("FP") or ch.upper().startswith("F")
    ]
    temporal_idx = [
        idx
        for idx, ch in enumerate(channel_names)
        if ch.upper().startswith("T") or ch.upper().startswith("TP")
    ]
    if not frontal_idx:
        frontal_idx = list(range(min(4, n_channels)))
    if not temporal_idx:
        temporal_idx = list(range(max(0, n_channels - 4), n_channels))

    all_idx = np.arange(n_channels)
    nonfrontal_idx = np.setdiff1d(all_idx, frontal_idx)
    nontemporal_idx = np.setdiff1d(all_idx, temporal_idx)

    rows = []
    for comp in range(n_components):
        spatial = component_maps[:, comp]
        source = sources[comp]

        frontal_ratio = _safe_mean_abs(spatial[frontal_idx]) / (
            _safe_mean_abs(spatial[nonfrontal_idx]) + 1e-12
        )
        temporal_ratio = _safe_mean_abs(spatial[temporal_idx]) / (
            _safe_mean_abs(spatial[nontemporal_idx]) + 1e-12
        )
        blink_ratio = _bandpower(source, fs, 0.5, 4.0) / (_bandpower(source, fs, 4.0, 40.0) + 1e-12)
        emg_ratio = _bandpower(source, fs, 30.0, 80.0) / (_bandpower(source, fs, 1.0, 30.0) + 1e-12)
        line_ratio = _bandpower(source, fs, 48.0, 52.0) / (_bandpower(source, fs, 0.5, 100.0) + 1e-12)
        kurt = _kurtosis(source)
        std = source.std(ddof=1)
        extreme_frac = 0.0 if std == 0 else float(np.mean(np.abs(source - source.mean()) > 5 * std))

        rows.append(
            {
                "component": comp,
                "frontal_ratio": round(float(frontal_ratio), 4),
                "temporal_ratio": round(float(temporal_ratio), 4),
                "blink_ratio": round(float(blink_ratio), 4),
                "emg_ratio": round(float(emg_ratio), 4),
                "line_ratio": round(float(line_ratio), 6),
                "kurtosis": round(float(kurt), 4),
                "extreme_frac": round(float(extreme_frac), 6),
            }
        )
    return pd.DataFrame(rows)


def evaluate_ica_components(features: pd.DataFrame, *, artifact_thresh: float = 1.5) -> pd.DataFrame:
    """Add artifact scores and labels to an ICA feature table."""
    result = features.copy()
    fz = _zscore(result["frontal_ratio"].to_numpy())
    tz = _zscore(result["temporal_ratio"].to_numpy())
    blz = _zscore(result["blink_ratio"].to_numpy())
    emz = _zscore(result["emg_ratio"].to_numpy())
    lnz = _zscore(result["line_ratio"].to_numpy())
    kz = _zscore(result["kurtosis"].to_numpy())
    exz = _zscore(result["extreme_frac"].to_numpy())

    ocular = 0.4 * fz + 0.4 * blz + 0.2 * kz
    muscle = 0.5 * emz + 0.3 * tz + 0.2 * kz
    line = lnz
    jump = 0.5 * kz + 0.5 * exz

    scores = np.vstack([ocular, muscle, line, jump]).T
    labels = []
    global_score = []
    label_names = ["eye/blink", "muscle", "line_noise", "jump"]
    for row in scores:
        max_idx = int(np.argmax(row))
        score = float(row[max_idx])
        global_score.append(score)
        labels.append(label_names[max_idx] if score > artifact_thresh else "brain")

    result["ocular_score"] = np.round(ocular, 3)
    result["muscle_score"] = np.round(muscle, 3)
    result["line_score"] = np.round(line, 3)
    result["jump_score"] = np.round(jump, 3)
    result["artifact_score"] = np.round(global_score, 3)
    result["artifact_type"] = labels
    return result


def suggest_rejected_components(evaluated: pd.DataFrame) -> list[int]:
    """Return non-brain component indices suggested for rejection."""
    bad = evaluated.loc[evaluated["artifact_type"] != "brain", "component"]
    return [int(idx) for idx in bad.to_list()]


def mne_artifact_suggestions(ica, raw) -> dict[str, list[int]]:
    """Try MNE's built-in EOG and muscle detectors when possible."""
    suggestions: dict[str, list[int]] = {"eog": [], "muscle": []}
    try:
        if "Fp2" in raw.ch_names:
            eog_idx, _ = ica.find_bads_eog(raw, ch_name="Fp2", verbose="ERROR")
            suggestions["eog"] = [int(idx) for idx in eog_idx]
    except Exception:
        suggestions["eog"] = []

    try:
        muscle_idx, _ = ica.find_bads_muscle(raw, verbose="ERROR")
        suggestions["muscle"] = [int(idx) for idx in muscle_idx]
    except Exception:
        suggestions["muscle"] = []
    return suggestions


def _safe_mean_abs(values: np.ndarray) -> float:
    return float(np.mean(np.abs(values))) if values.size else 0.0


def _zscore(values: np.ndarray) -> np.ndarray:
    std = values.std(ddof=1)
    if std == 0 or np.isnan(std):
        return np.zeros_like(values, dtype=np.float64)
    return (values - values.mean()) / std


def _kurtosis(values: np.ndarray) -> float:
    centered = values - values.mean()
    variance = np.mean(centered**2)
    if variance == 0:
        return 0.0
    return float(np.mean(centered**4) / (variance**2))


def _bandpower(values: np.ndarray, fs: float, f_lo: float, f_hi: float) -> float:
    x = values - values.mean()
    freqs = np.fft.rfftfreq(x.size, d=1.0 / fs)
    psd = np.abs(np.fft.rfft(x)) ** 2 / (fs * x.size)
    idx = (freqs >= f_lo) & (freqs <= f_hi)
    return float(psd[idx].sum())

