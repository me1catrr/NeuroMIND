"""Unit tests for mne_brain.processing.epochs — native MNE Phase 5."""

from __future__ import annotations

import numpy as np
import pytest

mne = pytest.importorskip("mne")

from mne_brain.common.config import load_config
from mne_brain.processing import (
    apply_baseline,
    epoch_summary,
    make_epochs,
    reject_artifacts,
)


# ── Fixtures ───────────────────────────────────────────────────────────────────

def _make_raw(
    n_ch: int = 5,
    duration: float = 10.0,
    sfreq: float = 500.0,
    amplitude_uv: float = 10.0,
    seed: int = 42,
) -> "mne.io.RawArray":
    """Synthetic RawArray in Volts, centred at 0, well below ±70 µV."""
    rng = np.random.default_rng(seed)
    data_v = rng.standard_normal((n_ch, int(duration * sfreq))) * (amplitude_uv * 1e-6)
    info = mne.create_info(
        ch_names=[f"EEG{i:03d}" for i in range(n_ch)],
        sfreq=sfreq,
        ch_types=["eeg"] * n_ch,
    )
    return mne.io.RawArray(data_v, info, verbose=False)


def _inject_artifact(raw: "mne.io.RawArray", epoch_idx: int, amplitude_v: float = 200e-6) -> None:
    """Inject a bipolar spike inside the given epoch.

    A constant offset would be zeroed by baseline correction (p2p → 0),
    so we use a +/- spike pair which survives mean subtraction and yields
    p2p ≈ 2 × amplitude_v >> the 140 µV MNE threshold.
    """
    sfreq = raw.info["sfreq"]
    epoch_samp = int(sfreq)                         # 1 s = 500 samples
    start = epoch_idx * epoch_samp
    mid = start + epoch_samp // 2
    raw._data[:, mid]     =  amplitude_v            # positive spike
    raw._data[:, mid + 1] = -amplitude_v            # negative spike (next sample)


# ── make_epochs ────────────────────────────────────────────────────────────────

def test_make_epochs_count_and_duration():
    """10 s / 1 s epochs → exactly 10 epochs of 500 samples each."""
    cfg = load_config()
    raw = _make_raw(duration=10.0, sfreq=500.0)
    epochs = make_epochs(raw, cfg)

    assert len(epochs) == 10
    assert epochs.get_data().shape == (10, 5, 500)   # (n_epochs, n_ch, n_times)


def test_make_epochs_no_overlap():
    """Consecutive epochs must not share samples."""
    cfg = load_config()
    raw = _make_raw(duration=5.0, sfreq=500.0)
    epochs = make_epochs(raw, cfg)

    data = epochs.get_data()          # (epochs, ch, times)
    # Last sample of epoch k ≠ first sample of epoch k+1 for a random signal
    for k in range(len(epochs) - 1):
        assert not np.allclose(data[k, :, -1], data[k + 1, :, 0])


# ── apply_baseline ─────────────────────────────────────────────────────────────

def test_baseline_zeroes_epoch_mean():
    """Mean over time axis must be ~0 for every epoch/channel after baseline."""
    cfg = load_config()
    raw = _make_raw(duration=5.0, amplitude_uv=30.0)
    epochs = make_epochs(raw, cfg)
    epochs_bl = apply_baseline(epochs.copy(), cfg)

    data = epochs_bl.get_data()           # (n_ep, n_ch, n_times)
    epoch_means = data.mean(axis=-1)      # (n_ep, n_ch)
    assert np.abs(epoch_means).max() < 1e-14


def test_baseline_idempotent_on_zero_mean():
    """Applying baseline to already zero-mean data should be a no-op."""
    cfg = load_config()
    raw = _make_raw(duration=3.0)
    epochs = make_epochs(raw, cfg)
    ep1 = apply_baseline(epochs.copy(), cfg)
    ep2 = apply_baseline(ep1.copy(), cfg)
    np.testing.assert_array_almost_equal(ep1.get_data(), ep2.get_data())


# ── reject_artifacts ──────────────────────────────────────────────────────────

def test_reject_drops_contaminated_epoch():
    """A 200 µV artefact epoch must be rejected; clean epochs must survive."""
    cfg = load_config()
    raw = _make_raw(n_ch=5, duration=10.0, amplitude_uv=5.0)   # safe signal
    _inject_artifact(raw, epoch_idx=3, amplitude_v=200e-6)      # epoch 3 = bad

    epochs = make_epochs(raw, cfg)
    n_before = len(epochs)
    epochs_clean = reject_artifacts(epochs.copy(), cfg)

    assert len(epochs_clean) < n_before, "Contaminated epoch was not rejected"
    assert len(epochs_clean) >= n_before - 1, "Too many epochs were rejected"


def test_clean_signal_survives_ar():
    """A signal well below ±70 µV must produce zero rejections."""
    cfg = load_config()
    raw = _make_raw(duration=10.0, amplitude_uv=5.0)   # ~5 µV RMS, safe

    epochs = make_epochs(raw, cfg)
    n_before = len(epochs)
    epochs_clean = reject_artifacts(epochs.copy(), cfg)

    assert len(epochs_clean) == n_before


# ── epoch_summary ─────────────────────────────────────────────────────────────

def test_summary_fields_no_rejection():
    """Summary for a clean recording must report 0 rejected, rate=0."""
    cfg = load_config()
    raw = _make_raw(duration=5.0)
    epochs = make_epochs(raw, cfg)
    n_initial = len(epochs)
    epochs = apply_baseline(epochs, cfg)

    summary = epoch_summary(epochs, n_initial)

    assert summary["n_initial"] == n_initial
    assert summary["n_valid"] == n_initial
    assert summary["n_rejected"] == 0
    assert summary["rejection_rate"] == 0.0
    assert summary["epoch_duration_s"] == pytest.approx(1.0, abs=1e-3)
    assert summary["n_channels"] == 5


def test_summary_counts_rejections():
    """Summary must correctly count rejected epochs."""
    cfg = load_config()
    raw = _make_raw(n_ch=3, duration=10.0, amplitude_uv=5.0)
    _inject_artifact(raw, epoch_idx=1, amplitude_v=200e-6)
    _inject_artifact(raw, epoch_idx=7, amplitude_v=200e-6)

    epochs = make_epochs(raw, cfg)
    n_initial = len(epochs)
    epochs = apply_baseline(epochs, cfg)
    epochs = reject_artifacts(epochs, cfg)

    summary = epoch_summary(epochs, n_initial)

    assert summary["n_rejected"] >= 2
    assert summary["rejection_rate"] == pytest.approx(
        summary["n_rejected"] / n_initial, abs=1e-6
    )
