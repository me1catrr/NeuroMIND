"""Configuration loader for the mirrored NeuroMIND pipeline config."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .types import PipelineConfig


CONFIG_SECTIONS = (
    "project",
    "study",
    "paths",
    "recording",
    "filtering",
    "segmentation",
    "baseline",
    "artifact_rejection",
    "ica",
    "spectral",
    "bands",
    "connectivity",
    "surrogates",
    "graph",
    "clinical",
    "longitudinal",
    "statistics",
    "export",
)


def default_config_path() -> Path:
    return Path(__file__).resolve().parents[3] / "config" / "pipeline_config.yaml"


def load_config(path: str | Path | None = None) -> PipelineConfig:
    """Load `pipeline_config.yaml` into a typed `PipelineConfig` dataclass."""
    config_path = Path(path) if path is not None else default_config_path()
    config_path = config_path.expanduser().resolve()

    try:
        import yaml
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "PyYAML is required to load pipeline_config.yaml. Install with "
            "`pip install -e .` or `pip install pyyaml`."
        ) from exc

    with config_path.open("r", encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle) or {}

    missing = [section for section in CONFIG_SECTIONS if section not in raw]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"Missing required config section(s): {joined}")

    bands = _normalize_bands(raw["bands"])

    return PipelineConfig(
        project=dict(raw["project"]),
        study=dict(raw["study"]),
        paths=dict(raw["paths"]),
        recording=dict(raw["recording"]),
        filtering=dict(raw["filtering"]),
        segmentation=dict(raw["segmentation"]),
        baseline=dict(raw["baseline"]),
        artifact_rejection=dict(raw["artifact_rejection"]),
        ica=dict(raw["ica"]),
        spectral=dict(raw["spectral"]),
        bands=bands,
        connectivity=dict(raw["connectivity"]),
        surrogates=dict(raw["surrogates"]),
        graph=dict(raw["graph"]),
        clinical=dict(raw["clinical"]),
        longitudinal=dict(raw["longitudinal"]),
        statistics=dict(raw["statistics"]),
        export_cfg=dict(raw["export"]),
        root=config_path.parents[1],
    )


def _normalize_bands(raw_bands: dict[str, Any]) -> dict[str, tuple[float, float]]:
    bands: dict[str, tuple[float, float]] = {}
    for name, bounds in raw_bands.items():
        if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
            raise ValueError(f"Band {name!r} must contain exactly two numeric bounds")
        bands[str(name)] = (float(bounds[0]), float(bounds[1]))
    return bands
