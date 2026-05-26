"""Unit tests for mne_brain.spectral.psd — native MNE Phase 6."""

from __future__ import annotations

import numpy as np
import pytest

mne = pytest.importorskip("mne")

from mne_brain.common.config import load_config
from mne_brain.spectral import compute_band_power, compute_psd, mean_spectrum


# ── Fixtures ───────────────────────────────────────────────────────────────────

def _make_epochs(
    n_ch: int = 4,
    n_epochs: int = 10,
    sfreq: float = 500.0,
    seed: int = 0,
) -> "mne.Epochs":
    """Synthetic preloaded Epochs in Volts (1 s each, mean-baselined)."""
    rng = np.random.default_rng(seed)
    n_samp = int(sfreq)
    data_v = rng.standard_normal((n_epochs, n_ch, n_samp)) * 20e-6
    info = mne.create_info(
        ch_names=[f"EEG{i:03d}" for i in range(n_ch)],
        sfreq=sfreq,
        ch_types=["eeg"] * n_ch,
    )
    epochs = mne.EpochsArray(data_v, info, verbose=False)
    epochs.apply_baseline((None, None), verbose=False)
    return epochs


def _make_alpha_epochs(sfreq: float = 500.0, n_epochs: int = 20) -> "mne.Epochs":
    """Epochs dominated by a 10 Hz sine (alpha band)."""
    n_samp = int(sfreq)
    t = np.linspace(0, 1.0, n_samp, endpoint=False)
    sine = np.sin(2 * np.pi * 10 * t) * 30e-6              # 30 µV @ 10 Hz
    ch_names = ["F3", "F4"]
    n_ch = len(ch_names)
    # shape must be (n_epochs, n_channels, n_samples)
    data_v = np.tile(sine[np.newaxis, np.newaxis, :], (n_epochs, n_ch, 1))
    info = mne.create_info(ch_names, sfreq, "eeg")
    epochs = mne.EpochsArray(data_v, info, verbose=False)
    epochs.apply_baseline((None, None), verbose=False)
    return epochs


# ── compute_psd ────────────────────────────────────────────────────────────────

def test_psd_output_shape():
    """EpochsSpectrum must have shape (n_epochs, n_ch, n_freqs)."""
    cfg = load_config()
    epochs = _make_epochs(n_ch=4, n_epochs=10)
    sp = compute_psd(epochs, cfg)

    n_ep, n_ch, n_freqs = sp.get_data().shape
    assert n_ep == 10
    assert n_ch == 4
    assert n_freqs > 0


def test_psd_rejects_empty_epochs():
    """PSD should fail early with a clear message when Phase 5 kept no epochs."""
    cfg = load_config()
    epochs = _make_epochs(n_epochs=1)
    epochs.drop([0], reason="TEST")

    with pytest.raises(ValueError, match="no valid epochs"):
        compute_psd(epochs, cfg)


def test_psd_frequency_resolution():
    """df must equal sfreq / n_fft = 500 / 512 ≈ 0.9766 Hz."""
    cfg = load_config()
    epochs = _make_epochs()
    sp = compute_psd(epochs, cfg)

    df = float(sp.freqs[1] - sp.freqs[0])
    assert df == pytest.approx(500.0 / 512.0, rel=1e-3)


def test_psd_nonnegative():
    """All PSD values must be ≥ 0 (power is non-negative)."""
    cfg = load_config()
    epochs = _make_epochs()
    sp = compute_psd(epochs, cfg)

    assert (sp.get_data() >= 0).all()


def test_psd_alpha_peak():
    """For a 10 Hz signal the alpha band must dominate over delta."""
    cfg = load_config()
    epochs = _make_alpha_epochs()
    sp = compute_psd(epochs, cfg)
    mean_sp = mean_spectrum(sp)

    alpha_power = mean_sp.get_data(fmin=7.8, fmax=11.7).mean()
    delta_power = mean_sp.get_data(fmin=0.5, fmax=4.0).mean()
    assert alpha_power > delta_power * 10, "Alpha should dominate for a 10 Hz signal"


# ── mean_spectrum ──────────────────────────────────────────────────────────────

def test_mean_spectrum_shape():
    """mean_spectrum must collapse the epoch axis → (n_ch, n_freqs)."""
    cfg = load_config()
    epochs = _make_epochs(n_ch=4, n_epochs=8)
    sp = compute_psd(epochs, cfg)
    mean_sp = mean_spectrum(sp)

    assert mean_sp.get_data().ndim == 2
    assert mean_sp.get_data().shape[0] == 4


def test_mean_spectrum_equals_epoch_mean():
    """mean_spectrum should equal the element-wise mean of per-epoch PSDs."""
    cfg = load_config()
    epochs = _make_epochs(n_epochs=6)
    sp = compute_psd(epochs, cfg)

    per_epoch = sp.get_data()                  # (n_ep, n_ch, n_freqs)
    expected = per_epoch.mean(axis=0)          # (n_ch, n_freqs)
    got = mean_spectrum(sp).get_data()

    np.testing.assert_allclose(got, expected, rtol=1e-6)


# ── compute_band_power ─────────────────────────────────────────────────────────

def test_band_power_keys_match_config():
    """compute_band_power must return one entry per configured band."""
    cfg = load_config()
    epochs = _make_epochs()
    sp = compute_psd(epochs, cfg)
    mean_sp = mean_spectrum(sp)

    bp = compute_band_power(mean_sp, cfg.bands)
    assert set(bp.keys()) == set(cfg.bands.keys())


def test_band_power_shape_and_positive():
    """Each band entry must be (n_channels,) and all values > 0."""
    cfg = load_config()
    epochs = _make_epochs(n_ch=4)
    sp = compute_psd(epochs, cfg)
    mean_sp = mean_spectrum(sp)
    bp = compute_band_power(mean_sp, cfg.bands)

    for band_name, values in bp.items():
        assert values.shape == (4,), f"{band_name}: wrong shape {values.shape}"
        assert (values > 0).all(), f"{band_name}: non-positive power"
