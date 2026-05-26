"""Unit tests for mne_brain.connectivity.wpli — native mne-connectivity Phase 7."""

from __future__ import annotations

import numpy as np
import pytest

mne = pytest.importorskip("mne")
pytest.importorskip("mne_connectivity")

from mne_brain.common.config import load_config
from mne_brain.connectivity import compute_wpli, connectivity_matrices


# ── Fixtures ───────────────────────────────────────────────────────────────────

def _make_epochs(n_ch: int = 4, n_epochs: int = 20, sfreq: float = 500.0, seed: int = 0):
    rng = np.random.default_rng(seed)
    data = rng.standard_normal((n_epochs, n_ch, int(sfreq))) * 20e-6
    info = mne.create_info(
        ch_names=[f"EEG{i:02d}" for i in range(n_ch)],
        sfreq=sfreq,
        ch_types=["eeg"] * n_ch,
    )
    epochs = mne.EpochsArray(data, info, verbose=False)
    epochs.apply_baseline((None, None), verbose=False)
    return epochs


def _make_coupled_epochs(n_epochs: int = 40, sfreq: float = 500.0):
    """Two channels sharing a strong 10 Hz component → high alpha wPLI."""
    rng = np.random.default_rng(7)
    t = np.linspace(0, 1.0, int(sfreq), endpoint=False)
    shared = np.sin(2 * np.pi * 10 * t) * 50e-6    # strong 10 Hz

    # ch0: shared alpha + noise
    # ch1: shared alpha (slightly phase-shifted) + noise
    # ch2: pure noise (no coupling)
    noise_scale = 5e-6
    data = np.zeros((n_epochs, 3, int(sfreq)))
    for ep in range(n_epochs):
        noise = rng.standard_normal((3, int(sfreq))) * noise_scale
        data[ep, 0] = shared + noise[0]
        data[ep, 1] = np.sin(2 * np.pi * 10 * t + 0.3) * 50e-6 + noise[1]
        data[ep, 2] = noise[2]

    info = mne.create_info(["F3", "F4", "Pz"], sfreq, "eeg")
    epochs = mne.EpochsArray(data, info, verbose=False)
    epochs.apply_baseline((None, None), verbose=False)
    return epochs


# ── compute_wpli ───────────────────────────────────────────────────────────────

def test_wpli_returns_spectral_connectivity():
    from mne_connectivity import SpectralConnectivity

    cfg = load_config()
    epochs = _make_epochs()
    con, band_names = compute_wpli(epochs, cfg)

    assert isinstance(con, SpectralConnectivity)
    assert band_names == list(cfg.bands.keys())


def test_wpli_dense_shape():
    """Dense output must be (n_ch, n_ch, n_bands)."""
    cfg = load_config()
    n_ch = 4
    epochs = _make_epochs(n_ch=n_ch)
    con, band_names = compute_wpli(epochs, cfg)

    dense = con.get_data(output="dense")
    assert dense.shape == (n_ch, n_ch, len(band_names))


def test_wpli_band_names_match_config():
    cfg = load_config()
    epochs = _make_epochs()
    _, band_names = compute_wpli(epochs, cfg)

    assert set(band_names) == set(cfg.bands.keys())
    assert band_names == list(cfg.bands.keys())


# ── connectivity_matrices ──────────────────────────────────────────────────────

def test_matrices_symmetric():
    """Each per-band wPLI matrix must be symmetric."""
    cfg = load_config()
    epochs = _make_epochs()
    con, band_names = compute_wpli(epochs, cfg)
    matrices = connectivity_matrices(con, band_names)

    for band, mat in matrices.items():
        np.testing.assert_allclose(
            mat, mat.T, atol=1e-10,
            err_msg=f"Band {band}: matrix is not symmetric"
        )


def test_matrices_diagonal_zero():
    """Self-connectivity (diagonal) must be 0."""
    cfg = load_config()
    epochs = _make_epochs()
    con, band_names = compute_wpli(epochs, cfg)
    matrices = connectivity_matrices(con, band_names)

    for band, mat in matrices.items():
        np.testing.assert_array_equal(
            np.diag(mat), 0,
            err_msg=f"Band {band}: diagonal is not zero"
        )


def test_matrices_values_in_range():
    """wPLI values must be in [0, 1]."""
    cfg = load_config()
    epochs = _make_epochs(n_epochs=30)
    con, band_names = compute_wpli(epochs, cfg)
    matrices = connectivity_matrices(con, band_names)

    for band, mat in matrices.items():
        assert mat.min() >= -1e-9, f"{band}: values below 0"
        assert mat.max() <= 1.0 + 1e-9, f"{band}: values above 1"


def test_coupled_channels_have_higher_alpha_wpli():
    """Channels sharing a 10 Hz component must show higher alpha wPLI."""
    cfg = load_config()
    epochs = _make_coupled_epochs(n_epochs=40)
    con, band_names = compute_wpli(epochs, cfg)
    matrices = connectivity_matrices(con, band_names)

    alpha_mat = matrices["ALPHA"]
    # F3–F4 (coupled) vs F3–Pz or F4–Pz (uncoupled)
    wpli_coupled = alpha_mat[0, 1]
    wpli_uncoupled = (alpha_mat[0, 2] + alpha_mat[1, 2]) / 2
    assert wpli_coupled > wpli_uncoupled, (
        f"Expected coupled wPLI ({wpli_coupled:.4f}) > "
        f"uncoupled ({wpli_uncoupled:.4f}) in alpha band"
    )
