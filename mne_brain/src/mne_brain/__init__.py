"""MNE-Python replication pipeline for NeuroMIND/BRAIN EEG analyses."""

from .common.config import load_config
from .common.types import PipelineConfig

__all__ = ["PipelineConfig", "load_config"]
