"""Phase 7 — wPLI functional connectivity.

Uses native mne-connectivity:
  spectral_connectivity_epochs(method='wpli', faverage=True)
  → SpectralConnectivity  get_data(output='dense')  (n_ch, n_ch, n_bands)

Difference from NeuroMIND
--------------------------
NeuroMIND:  Butterworth BP (order 8) + scipy Hilbert, across-segments aggregation
            (Vinck 2011 time-domain formulation).
MNE:        Cross-spectral density estimated with multitaper/Fourier in the
            frequency domain; wPLI computed from Im(CSD).

Both estimators are mathematically equivalent in expectation but differ in
numerical precision, windowing and tapering. Expected Pearson r > 0.90
between the two for large epoch counts (validated in Phase 10).
"""

from __future__ import annotations

import json
from pathlib import Path

import mne
import numpy as np
import pandas as pd

from mne_brain.common.types import PipelineConfig


# ── Heatmap helper ─────────────────────────────────────────────────────────────

def _plot_wpli_heatmap(
    matrices: dict[str, np.ndarray],
    ch_names: list[str],
    condition: str,
    figures_dir: Path,
) -> list[str]:
    """Save one wPLI heatmap PNG per band.  Returns list of saved paths."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib import cm
    except ImportError:
        return []

    saved: list[str] = []
    n_ch = len(ch_names)
    tick_step = max(1, n_ch // 10)

    for band, mat in matrices.items():
        fig, ax = plt.subplots(figsize=(9, 7.5))
        im = ax.imshow(mat, cmap="RdYlBu_r", vmin=0.0, vmax=mat.max() or 0.5,
                       aspect="auto", interpolation="nearest")
        fig.colorbar(im, ax=ax, label="wPLI")

        ticks = list(range(0, n_ch, tick_step))
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([ch_names[t] for t in ticks], rotation=90, fontsize=8)
        ax.set_yticklabels([ch_names[t] for t in ticks], fontsize=8)
        ax.set_title(f"wPLI  —  {band}  [{condition}]", fontsize=11, fontweight="bold")
        ax.set_xlabel("Channel", fontsize=9)
        ax.set_ylabel("Channel", fontsize=9)
        plt.tight_layout()

        out_path = figures_dir / f"wpli_heatmap_{band}_{condition}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        saved.append(str(out_path))

    # ── Combined figure: all bands in one grid ────────────────────────────────
    n_bands = len(matrices)
    ncols = min(4, n_bands)
    nrows = (n_bands + ncols - 1) // ncols
    fig_all, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.5, nrows * 4.2))
    axes = np.array(axes).flatten()

    for k, (band, mat) in enumerate(matrices.items()):
        ax = axes[k]
        vmax = float(mat.max()) or 0.5
        im = ax.imshow(mat, cmap="RdYlBu_r", vmin=0.0, vmax=vmax,
                       aspect="auto", interpolation="nearest")
        fig_all.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(f"{band}", fontsize=9, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])

    for ax in axes[n_bands:]:
        ax.set_visible(False)

    fig_all.suptitle(f"wPLI Connectivity — {condition}", fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    combined_path = figures_dir / f"wpli_all_bands_{condition}.png"
    fig_all.savefig(combined_path, dpi=150, bbox_inches="tight")
    plt.close(fig_all)
    saved.append(str(combined_path))

    return saved


# ── Connectivity computation ───────────────────────────────────────────────────

def compute_wpli(
    epochs: mne.Epochs,
    cfg: PipelineConfig,
) -> "mne_connectivity.SpectralConnectivity":  # noqa: F821
    """Compute wPLI for all configured frequency bands.

    Parameters
    ----------
    epochs:
        Clean, baseline-corrected, AR-filtered MNE Epochs.
    cfg:
        Pipeline config.  Bands are read from ``cfg.bands``.

    Returns
    -------
    SpectralConnectivity
        ``get_data(output='dense')`` → (n_ch, n_ch, n_bands)
        ``freqs``                    → band centre frequencies
        ``names``                    → channel names
    """
    import warnings

    try:
        from mne_connectivity import spectral_connectivity_epochs
    except ImportError as exc:
        raise ImportError(
            "mne-connectivity is required. Install with: pip install mne-connectivity"
        ) from exc

    bands = cfg.bands
    band_names = list(bands.keys())
    fmins = tuple(float(bands[b][0]) for b in band_names)
    fmaxs = tuple(float(bands[b][1]) for b in band_names)

    # Warn once about delta reliability (< 5 cycles @ 1 s epochs), then suppress
    # the MNE repeat. The limitation is shared with NeuroMIND (same epoch length).
    epoch_s = float(epochs.tmax - epochs.tmin + 1.0 / epochs.info["sfreq"])
    min_reliable_hz = 5.0 / epoch_s
    if fmins[0] < min_reliable_hz:
        print(
            f"  [note] DELTA fmin={fmins[0]} Hz < {min_reliable_hz:.1f} Hz "
            f"({5} cycles / {epoch_s:.1f} s epoch). "
            "Delta connectivity is a lower-bound estimate — same caveat applies to NeuroMIND."
        )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*fmin.*cycles.*", category=RuntimeWarning)
        warnings.filterwarnings("ignore", message=".*no Annotations.*", category=RuntimeWarning)
        con = spectral_connectivity_epochs(
            epochs,
            method="wpli",
            fmin=fmins,
            fmax=fmaxs,
            faverage=True,          # one value per band per connection
            sfreq=float(epochs.info["sfreq"]),
            verbose=False,
        )
    return con, band_names


def connectivity_matrices(
    con: "mne_connectivity.SpectralConnectivity",  # noqa: F821
    band_names: list[str],
) -> dict[str, np.ndarray]:
    """Return one symmetric (n_ch, n_ch) wPLI matrix per band.

    MNE returns an upper-triangular result; this function mirrors it to produce
    a full symmetric matrix (diagonal = 0), matching NeuroMIND's output.
    """
    dense = con.get_data(output="dense")      # (n_ch, n_ch, n_bands)
    matrices: dict[str, np.ndarray] = {}
    for k, name in enumerate(band_names):
        mat = dense[:, :, k].copy()
        # Mirror upper → lower so the matrix is symmetric
        mat = mat + mat.T
        np.fill_diagonal(mat, 0.0)
        matrices[name] = mat
    return matrices


# ── Save outputs ───────────────────────────────────────────────────────────────

def save_connectivity_results(
    con: "mne_connectivity.SpectralConnectivity",  # noqa: F821
    band_names: list[str],
    cfg: PipelineConfig,
    out_dir: Path,
    condition: str,
) -> dict:
    """Persist connectivity matrices, summary and heatmap figures for one recording.

    Outputs  (matches NeuroMIND result structure for Phase 10 comparison)
    -------
    tables/wpli_{band}_{cond}.csv              31×31 wPLI matrix per band
    tables/connectivity_edges_{cond}.csv       long format: ch1, ch2, band, wpli
    tables/connectivity_summary_{cond}.json
    cache/connectivity_{cond}.npz
    figures/wpli_heatmap_{band}_{cond}.png     one heatmap per band
    figures/wpli_all_bands_{cond}.png          all-bands combined grid
    """
    out_dir = Path(out_dir)
    tables_dir = out_dir / "tables"
    cache_dir = out_dir / "cache"
    figures_dir = out_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    ch_names = list(con.names)
    matrices = connectivity_matrices(con, band_names)

    # ── Per-band CSV matrices ────────────────────────────────────────────────
    for band, mat in matrices.items():
        df = pd.DataFrame(mat, index=ch_names, columns=ch_names)
        df.index.name = "channel"
        df.to_csv(tables_dir / f"wpli_{band}_{condition}.csv")

    # ── Edges long-format CSV ────────────────────────────────────────────────
    rows = []
    n_ch = len(ch_names)
    for band, mat in matrices.items():
        for i in range(n_ch):
            for j in range(i + 1, n_ch):
                rows.append({
                    "ch1": ch_names[i],
                    "ch2": ch_names[j],
                    "band": band,
                    "wpli": round(float(mat[i, j]), 6),
                })
    edges_df = pd.DataFrame(rows)
    edges_df.to_csv(tables_dir / f"connectivity_edges_{condition}.csv", index=False)

    # ── Summary JSON ─────────────────────────────────────────────────────────
    summary: dict = {
        "method": "wpli",
        "library": "mne-connectivity",
        "n_channels": n_ch,
        "n_epochs_used": int(con.n_epochs_used),
        "condition": condition,
        "band_summary": {},
        "note": (
            "Frequency-domain CSD-based wPLI (MNE). "
            "NeuroMIND uses time-domain Hilbert across-segments. "
            "Expected r > 0.90 — validated in Phase 10."
        ),
    }
    for band, mat in matrices.items():
        upper = mat[np.triu_indices(n_ch, k=1)]
        summary["band_summary"][band] = {
            "mean_wpli": round(float(upper.mean()), 6),
            "max_wpli": round(float(upper.max()), 6),
            "n_connections": int(len(upper)),
        }

    # ── Heatmap figures ───────────────────────────────────────────────────────
    figure_paths = _plot_wpli_heatmap(matrices, ch_names, condition, figures_dir)
    summary["figure_paths"] = figure_paths

    (tables_dir / f"connectivity_summary_{condition}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    # ── NumPy cache (for Phase 10 validation) ────────────────────────────────
    stack = np.stack([matrices[b] for b in band_names], axis=-1)   # (n_ch, n_ch, n_bands)
    np.savez_compressed(
        cache_dir / f"connectivity_{condition}.npz",
        wpli=stack,
        band_names=np.asarray(band_names, dtype=object),
        channel_names=np.asarray(ch_names, dtype=object),
        n_epochs=int(con.n_epochs_used),
    )

    return summary
