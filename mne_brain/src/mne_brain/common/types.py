"""Core scientific entities mirrored from NeuroMIND/src/types.jl."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]


@dataclass(frozen=True)
class PipelineConfig:
    project: dict[str, Any]
    study: dict[str, Any]
    paths: dict[str, Any]
    recording: dict[str, Any]
    filtering: dict[str, Any]
    segmentation: dict[str, Any]
    baseline: dict[str, Any]
    artifact_rejection: dict[str, Any]
    ica: dict[str, Any]
    spectral: dict[str, Any]
    bands: dict[str, tuple[float, float]]
    connectivity: dict[str, Any]
    surrogates: dict[str, Any]
    graph: dict[str, Any]
    clinical: dict[str, Any]
    longitudinal: dict[str, Any]
    statistics: dict[str, Any]
    export_cfg: dict[str, Any]
    root: Path


@dataclass(frozen=True)
class ClinicalData:
    EDSS: float | None = None
    disease_duration_y: float | None = None
    medication: str | None = None
    fatigue_score: float | None = None
    cognition_score: float | None = None
    lesion_load: float | None = None


@dataclass(frozen=True)
class RecordingMeta:
    subject_id: str
    session_id: str
    condition: str
    run: int
    fs: float
    n_channels: int
    channel_names: list[str]
    channel_positions: dict[str, tuple[float, float]] | None
    bids_path: str


@dataclass(frozen=True)
class EEGRecording:
    meta: RecordingMeta
    data: FloatArray
    times: FloatArray

    @property
    def n_channels(self) -> int:
        return int(self.data.shape[0])

    @property
    def n_samples(self) -> int:
        return int(self.data.shape[1])

    @property
    def duration(self) -> float:
        return self.n_samples / self.meta.fs


@dataclass(frozen=True)
class EpochSet:
    meta: RecordingMeta
    data: FloatArray
    epoch_length_s: float
    n_valid: int
    rejected_idx: list[int]
    rejection_reasons: list[str] = field(default_factory=list)

    @property
    def n_epochs(self) -> int:
        return int(self.data.shape[2])

    @property
    def n_samples_epoch(self) -> int:
        return int(self.data.shape[1])

    @property
    def rejection_rate(self) -> float:
        total = self.n_valid + len(self.rejected_idx)
        return len(self.rejected_idx) / total if total else 0.0


@dataclass(frozen=True)
class ICAResult:
    meta: RecordingMeta
    mixing_matrix: FloatArray
    unmixing_matrix: FloatArray
    activations: FloatArray
    rejected_components: list[int]
    variance_explained: FloatArray


@dataclass(frozen=True)
class SpectralResult:
    meta: RecordingMeta
    psd: FloatArray
    freqs: FloatArray
    band_power: dict[str, FloatArray]
    n_epochs_used: int
    params: dict[str, Any]


@dataclass(frozen=True)
class ConnectivityMatrix:
    meta: RecordingMeta
    method: str
    matrices: dict[str, FloatArray]
    channel_names: list[str]
    space: str
    n_epochs_used: int
    params: dict[str, Any]

    @property
    def n_channels(self) -> int:
        return len(self.channel_names)


@dataclass(frozen=True)
class SurrogateResult:
    connectivity: ConnectivityMatrix
    band: str
    observed: FloatArray
    null_distribution: FloatArray
    p_values: FloatArray
    sig_mask: BoolArray
    fdr_threshold: float
    n_surrogates: int


@dataclass(frozen=True)
class GraphMetrics:
    band: str
    threshold: float
    density: float
    strength: FloatArray
    clustering: FloatArray
    path_length: float
    efficiency: float
    modularity: float
    channel_names: list[str]


@dataclass(frozen=True)
class StatResult:
    test_name: str
    statistic: float
    p_value: float
    p_adjusted: float
    effect_size: float
    significant: bool
    group_a_mean: float
    group_b_mean: float
    n_a: int
    n_b: int


@dataclass
class Session:
    id: str
    visit_number: int
    recordings: dict[str, EEGRecording] = field(default_factory=dict)
    epochs: dict[str, EpochSet] = field(default_factory=dict)
    ica: dict[str, ICAResult] = field(default_factory=dict)
    spectra: dict[str, SpectralResult] = field(default_factory=dict)
    connectivity: dict[str, ConnectivityMatrix] = field(default_factory=dict)
    surrogates: dict[str, SurrogateResult] = field(default_factory=dict)
    graph_metrics: dict[str, GraphMetrics] = field(default_factory=dict)


@dataclass
class Subject:
    id: str
    group: str
    age: int | None = None
    sex: str | None = None
    clinical: ClinicalData = field(default_factory=ClinicalData)
    sessions: dict[str, Session] = field(default_factory=dict)


@dataclass(frozen=True)
class GroupAnalysis:
    session_id: str
    condition: str
    band: str
    group_ms: str
    group_ctrl: str
    mean_connectivity_ms: FloatArray
    mean_connectivity_ctrl: FloatArray
    stat_results: list[list[StatResult]]
    channel_names: list[str]
    n_ms: int
    n_ctrl: int


@dataclass(frozen=True)
class LongitudinalAnalysis:
    subject_id: str
    condition: str
    band: str
    visits: list[str]
    connectivity_over_time: list[FloatArray]
    graph_metrics_over_time: list[GraphMetrics]
    channel_names: list[str]

