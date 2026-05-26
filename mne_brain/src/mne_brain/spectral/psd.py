"""Phase 6 — Power spectral density and band-power estimation.

Uses native MNE:
  epochs.compute_psd()   →  EpochsSpectrum   (per-epoch Welch PSD)
  spectrum.average()     →  Spectrum         (mean across epochs)
  spectrum.get_data()    →  band-limited NumPy arrays
  spectrum.plot()        →  native MNE figure saved as PNG

Parameters matched to NeuroMIND/config/pipeline.toml:
  n_fft      = 512   →  df = 500/512 ≈ 0.9766 Hz  (identical resolution)
  n_per_seg  = 500   →  one FFT per 1-s epoch (= one Welch segment)
  window     = 'hamming'

Difference from NeuroMIND
--------------------------
NeuroMIND applies a custom Hamming *taper* (10 % of each edge only).
MNE's ``window='hamming'`` applies a full Hamming window over the whole segment.
Both produce the same frequency resolution; absolute power values differ slightly.
This difference is documented in validation/reports/ for the comparison panel.
"""

from __future__ import annotations

import json
from pathlib import Path

import mne
import numpy as np
import pandas as pd

from mne_brain.common.types import PipelineConfig


# ── PSD computation ────────────────────────────────────────────────────────────

def compute_psd(epochs: mne.Epochs, cfg: PipelineConfig) -> mne.time_frequency.EpochsSpectrum:
    """Compute per-epoch PSD using epochs.compute_psd() (Welch method).

    Parameters
    ----------
    epochs:
        Clean, baseline-corrected, AR-filtered MNE Epochs (preloaded).
    cfg:
        Pipeline config.  Key spectral params: nfft=512, window='hamming'.

    Returns
    -------
    EpochsSpectrum  shape (n_epochs, n_channels, n_freqs)
    """
    sp_cfg = cfg.spectral
    n_fft = int(sp_cfg.get("nfft", 512))
    window = str(sp_cfg.get("window", "hamming"))

    # n_per_seg = epoch length in samples → 1 Welch segment per epoch,
    # matching NeuroMIND's single-FFT-per-epoch approach.
    sfreq = epochs.info["sfreq"]
    if len(epochs) == 0:
        raise ValueError("Cannot compute PSD: no valid epochs available after artifact rejection")
    n_per_seg = len(epochs.times)   # samples in one epoch

    spectrum = epochs.compute_psd(
        method="welch",
        fmin=0.5,
        fmax=sfreq / 2.0,          # full spectrum up to Nyquist
        remove_dc=True,
        n_fft=n_fft,
        n_per_seg=n_per_seg,
        n_overlap=0,
        window=window,
        average="mean",            # average within-epoch windows (only 1 here)
        verbose=False,
    )
    return spectrum


def mean_spectrum(spectrum: mne.time_frequency.EpochsSpectrum) -> mne.time_frequency.Spectrum:
    """Average an EpochsSpectrum across epochs → Spectrum (n_channels, n_freqs)."""
    return spectrum.average()


# ── Band power ─────────────────────────────────────────────────────────────────

def compute_band_power(
    mean_sp: mne.time_frequency.Spectrum,
    bands: dict[str, tuple[float, float]],
) -> dict[str, np.ndarray]:
    """Return mean PSD power per channel for each frequency band.

    Uses spectrum.get_data(fmin, fmax) — native MNE band slicing.

    Returns
    -------
    dict  band_name → ndarray shape (n_channels,)  [µV²]
    """
    band_power: dict[str, np.ndarray] = {}
    for band_name, (f_lo, f_hi) in bands.items():
        try:
            data = mean_sp.get_data(fmin=f_lo, fmax=f_hi)   # (n_ch, n_freqs_in_band)
            # Convert V² → µV²  (MNE stores PSD in V²/Hz; we want µV² for comparison)
            band_power[band_name] = data.mean(axis=-1) * 1e12   # V²/Hz → µV²/Hz
        except Exception:
            band_power[band_name] = np.full(mean_sp.get_data().shape[0], np.nan)
    return band_power


# ── Save outputs ───────────────────────────────────────────────────────────────

def save_spectral_results(
    spectrum: mne.time_frequency.EpochsSpectrum,
    cfg: PipelineConfig,
    out_dir: Path,
    condition: str,
    *,
    save_figure: bool = True,
) -> dict:
    """Persist PSD tables, band-power CSV, params JSON, and optional figure.

    Outputs
    -------
    tables/psd_by_channel_{cond}.csv        rows=channels, cols=freq_bins (µV²/Hz)
    tables/band_power_summary_{cond}.csv    rows=channels, cols=bands
    tables/spectral_params_{cond}.json
    figures/psd_spectrum_{cond}.png         (if save_figure=True)
    cache/spectrum_{cond}.npz               raw PSD array for Phase 10 validation
    """
    out_dir = Path(out_dir)
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    cache_dir = out_dir / "cache"
    for d in (tables_dir, figures_dir, cache_dir):
        d.mkdir(parents=True, exist_ok=True)

    mean_sp = mean_spectrum(spectrum)
    ch_names = mean_sp.ch_names
    freqs = mean_sp.freqs

    # ── PSD table  (µV²/Hz) ──────────────────────────────────────────────────
    psd_v2 = mean_sp.get_data()                  # (n_ch, n_freqs)  V²/Hz
    psd_uv2 = psd_v2 * 1e12                      # → µV²/Hz

    psd_df = pd.DataFrame(
        psd_uv2,
        index=ch_names,
        columns=[f"{f:.4f}" for f in freqs],
    )
    psd_df.index.name = "channel"
    psd_df.to_csv(tables_dir / f"psd_by_channel_{condition}.csv")

    # ── Band power table ──────────────────────────────────────────────────────
    band_power = compute_band_power(mean_sp, cfg.bands)
    bp_df = pd.DataFrame(band_power, index=ch_names)
    bp_df.index.name = "channel"
    bp_df.to_csv(tables_dir / f"band_power_summary_{condition}.csv")

    # ── Params JSON ───────────────────────────────────────────────────────────
    # n_per_seg is the time-domain segment length (= epoch length in samples),
    # not the frequency axis length; derive it from sfreq and epoch duration.
    sfreq = float(spectrum.info["sfreq"])
    n_per_seg_samples = int(round(sfreq * float(cfg.segmentation.get("epoch_length_s", 1.0))))
    params = {
        "method": "welch",
        "n_fft": int(cfg.spectral.get("nfft", 512)),
        "n_per_seg": n_per_seg_samples,
        "window": str(cfg.spectral.get("window", "hamming")),
        "df_hz": round(float(freqs[1] - freqs[0]), 6) if len(freqs) > 1 else None,
        "fmin_hz": float(freqs[0]),
        "fmax_hz": float(freqs[-1]),
        "n_freqs": len(freqs),
        "n_epochs_used": len(spectrum),
        "n_channels": len(ch_names),
        "units": "uV2_per_Hz",
        "note": (
            "Full Hamming window (MNE native). "
            "NeuroMIND uses 10 % edge Hamming taper — minor absolute power difference."
        ),
    }
    (tables_dir / f"spectral_params_{condition}.json").write_text(
        json.dumps(params, indent=2), encoding="utf-8"
    )

    # ── Cache (numpy) ─────────────────────────────────────────────────────────
    np.savez_compressed(
        cache_dir / f"spectrum_{condition}.npz",
        psd_uv2=psd_uv2,
        freqs=freqs,
        channel_names=np.asarray(ch_names, dtype=object),
        n_epochs=len(spectrum),
    )

    # ── MNE native figure ─────────────────────────────────────────────────────
    if save_figure:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig = mean_sp.plot(
                picks="eeg",
                amplitude=False,        # plot power, not amplitude
                spatial_colors=True,
                show=False,
            )
            fig.savefig(figures_dir / f"psd_spectrum_{condition}.png", dpi=150, bbox_inches="tight")
            plt.close(fig)
        except Exception as exc:
            print(f"  [warn] PSD figure not saved: {exc}")

    return params
