"""Epoching, baseline correction, and artifact rejection — native MNE."""

from mne_brain.processing.epochs import (
    apply_baseline,
    epoch_summary,
    load_cleaned_raw,
    make_epochs,
    reject_artifacts,
    save_epoch_exclusion,
    save_epoch_results,
)

__all__ = [
    "make_epochs",
    "apply_baseline",
    "reject_artifacts",
    "load_cleaned_raw",
    "epoch_summary",
    "save_epoch_exclusion",
    "save_epoch_results",
]
