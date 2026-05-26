"""ICA fitting and artifact classification utilities."""

from .ica_classification import (
    compute_ica_features,
    evaluate_ica_components,
    mne_artifact_suggestions,
    suggest_rejected_components,
)
from .ica_core import MNEICAResult, apply_ica_rejection, load_filtered_recording, recording_to_raw, run_ica

__all__ = [
    "MNEICAResult",
    "apply_ica_rejection",
    "compute_ica_features",
    "evaluate_ica_components",
    "load_filtered_recording",
    "mne_artifact_suggestions",
    "recording_to_raw",
    "run_ica",
    "suggest_rejected_components",
]
