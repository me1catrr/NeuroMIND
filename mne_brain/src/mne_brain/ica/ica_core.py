"""MNE-based ICA fitting utilities for the mne_brain pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mne_brain.common.types import EEGRecording, ICAResult, PipelineConfig, RecordingMeta


@dataclass(frozen=True)
class MNEICAResult:
    """Container for an MNE ICA object and its exported numeric data."""

    ica: object
    raw: object
    result: ICAResult
    component_maps: np.ndarray
    sources: np.ndarray
    explained_variance: np.ndarray


def load_filtered_recording(npz_path: str | Path, meta: RecordingMeta) -> EEGRecording:
    """Load a filtered recording saved by `scripts/run_phase3_m05.py`."""
    path = Path(npz_path).expanduser().resolve()
    with np.load(path, allow_pickle=True) as npz:
        data = np.asarray(npz["data"], dtype=np.float64)
        times = np.asarray(npz["times"], dtype=np.float64)
        channel_names = [str(ch) for ch in npz["channel_names"].tolist()]
        fs = float(npz["fs"])

    loaded_meta = RecordingMeta(
        subject_id=meta.subject_id,
        session_id=meta.session_id,
        condition=meta.condition,
        run=meta.run,
        fs=fs,
        n_channels=len(channel_names),
        channel_names=channel_names,
        channel_positions=meta.channel_positions,
        bids_path=str(path),
    )
    return EEGRecording(loaded_meta, data, times)


def recording_to_raw(rec: EEGRecording):
    """Convert an `EEGRecording` in microvolts to an MNE `RawArray` in volts."""
    try:
        import mne
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "MNE is required for ICA. Install with `pip install -e .`."
        ) from exc

    info = mne.create_info(
        ch_names=rec.meta.channel_names,
        sfreq=rec.meta.fs,
        ch_types=["eeg"] * len(rec.meta.channel_names),
    )
    raw = mne.io.RawArray(rec.data * 1e-6, info, verbose="ERROR")
    try:
        raw.set_montage("standard_1020", on_missing="ignore", verbose="ERROR")
    except Exception:
        pass
    return raw


def run_ica(rec: EEGRecording, cfg: PipelineConfig) -> MNEICAResult:
    """Fit MNE FastICA to a continuous filtered recording."""
    try:
        from mne.preprocessing import ICA
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "MNE is required for ICA. Install with `pip install -e .`."
        ) from exc

    raw = recording_to_raw(rec)
    ica_cfg = cfg.ica
    n_components = min(int(ica_cfg.get("n_components", 30)), rec.meta.n_channels)
    max_iter = int(ica_cfg.get("max_iter", 500))
    random_state = int(ica_cfg.get("random_seed", ica_cfg.get("seed", 42)))
    tolerance = float(ica_cfg.get("tolerance", ica_cfg.get("tol", 1e-5)))
    method = str(ica_cfg.get("method", "fastica"))

    ica = ICA(
        n_components=n_components,
        method=method,
        random_state=random_state,
        max_iter=max_iter,
        fit_params={"tol": tolerance},
    )
    ica.fit(raw, verbose="ERROR")

    sources = np.asarray(ica.get_sources(raw).get_data(), dtype=np.float64)
    component_maps = np.asarray(ica.get_components(), dtype=np.float64)
    unmixing = np.asarray(ica.unmixing_matrix_, dtype=np.float64)
    mixing = np.asarray(ica.mixing_matrix_, dtype=np.float64)
    explained = _explained_variance(ica, n_components)

    result = ICAResult(
        meta=rec.meta,
        mixing_matrix=component_maps,
        unmixing_matrix=unmixing,
        activations=sources,
        rejected_components=[],
        variance_explained=explained,
    )
    return MNEICAResult(
        ica=ica,
        raw=raw,
        result=result,
        component_maps=component_maps,
        sources=sources,
        explained_variance=explained,
    )


def apply_ica_rejection(rec: EEGRecording, fitted: MNEICAResult, rejected: list[int]) -> EEGRecording:
    """Apply selected ICA component rejection and return cleaned data in microvolts."""
    raw_clean = fitted.raw.copy()
    fitted.ica.exclude = list(rejected)
    fitted.ica.apply(raw_clean, verbose="ERROR")
    data_uv = np.asarray(raw_clean.get_data() * 1e6, dtype=np.float64)
    return EEGRecording(meta=rec.meta, data=data_uv, times=rec.times.copy())


def _explained_variance(ica, n_components: int) -> np.ndarray:
    values = getattr(ica, "pca_explained_variance_", None)
    if values is None:
        return np.full(n_components, np.nan, dtype=np.float64)
    values = np.asarray(values[:n_components], dtype=np.float64)
    total = np.nansum(values)
    if total <= 0:
        return np.full(n_components, np.nan, dtype=np.float64)
    return values / total
