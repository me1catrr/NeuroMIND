"""Shared configuration and data types for mne_brain."""

from .config import load_config
from .types import PipelineConfig

__all__ = ["PipelineConfig", "load_config"]

