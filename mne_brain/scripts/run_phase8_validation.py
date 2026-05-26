"""Phase 8 — Cross-pipeline validation: NeuroMIND (Julia) vs mne_brain (Python).

Compares results for one subject/session/condition between the two pipelines and
produces a statistical report + publication-quality figures.

Usage
-----
python scripts/run_phase8_validation.py                          # M05 / T2 / EC
python scripts/run_phase8_validation.py --condition EO
python scripts/run_phase8_validation.py --subject M05 --session T2 --condition EC

Outputs  → validation/reports/
-----------
comparison_{subj}_{sess}_{cond}_wpli.png        scatter + Bland-Altman grid (per band)
comparison_{subj}_{sess}_{cond}_heatmaps.png    side-by-side heatmaps (NeuroMIND vs mne_brain)
comparison_{subj}_{sess}_{cond}_psd.png         PSD overlay + log-ratio
comparison_{subj}_{sess}_{cond}_band_power.png  grouped bar chart per band
comparison_{subj}_{sess}_{cond}_summary.json    all metrics (r, ρ, MAE, RMSE, bias)

Metric notes
------------
- wPLI: upper-triangle only (n*(n-1)/2 pairs); Fp2 removed from mne_brain when comparing
  because NeuroMIND flagged it as bad and excluded it.
- PSD: NeuroMIND stores long-format (channel, freq_hz, power_uv2); mne_brain stores wide.
  Frequencies are interpolated onto a shared grid before computing correlation.
- Band power: both store µV²/Hz; direct channel-wise comparison after aligning channels.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy import stats

from mne_brain import load_config


# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────

def _get_paths(cfg, subject: str, session: str, condition: str) -> tuple[Path, Path, Path]:
    task = "eyesclosed" if condition == "EC" else "eyesopen"
    mne_brain_root = Path(cfg.root)                  # …/NeuroMIND/mne_brain
    neuromind_root = mne_brain_root.parent            # …/NeuroMIND

    nm_base = neuromind_root / "results" / "subjects" / f"sub-{subject}" / f"ses-{session}" / task
    mb_base = mne_brain_root / "results" / "subjects" / f"sub-{subject}" / f"ses-{session}" / task
    report_dir = mne_brain_root / "validation" / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    return nm_base, mb_base, report_dir


# ─────────────────────────────────────────────────────────────────────────────
# Data loading helpers
# ─────────────────────────────────────────────────────────────────────────────

def _load_nm_wpli(nm_base: Path, band: str) -> pd.DataFrame:
    """Load NeuroMIND wPLI matrix CSV (channel × channel, index=channel)."""
    p = nm_base / f"wpli_{band}.csv"
    df = pd.read_csv(p, index_col=0)
    df.index = df.index.str.strip()
    df.columns = df.columns.str.strip()
    return df


def _load_mb_wpli(mb_base: Path, band: str, condition: str) -> pd.DataFrame:
    """Load mne_brain wPLI matrix CSV (channel × channel)."""
    p = mb_base / "tables" / f"wpli_{band}_{condition}.csv"
    df = pd.read_csv(p, index_col=0)
    df.index.name = "channel"
    return df


def _align_channels(df_nm: pd.DataFrame, df_mb: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Restrict both matrices to their common channels (same order)."""
    common = [c for c in df_nm.index if c in df_mb.index]
    return df_nm.loc[common, common], df_mb.loc[common, common]


def _upper_triangle(df: pd.DataFrame) -> np.ndarray:
    """Flatten upper triangle (k=1) to 1-D array."""
    mat = df.values.astype(float)
    idx = np.triu_indices(mat.shape[0], k=1)
    return mat[idx]


def _load_nm_psd(nm_base: Path) -> pd.DataFrame:
    """Load NeuroMIND PSD from long format → wide (channel × freq)."""
    p = nm_base / "psd_by_channel.csv"
    long = pd.read_csv(p)
    long.columns = long.columns.str.strip()
    # Pivot: index=channel, columns=freq_hz, values=power_uv2
    wide = long.pivot(index="channel", columns="freq_hz", values="power_uv2")
    wide.columns = wide.columns.astype(float)
    return wide


def _load_mb_psd(mb_base: Path, condition: str) -> pd.DataFrame:
    """Load mne_brain PSD wide format (channel × freq_str)."""
    p = mb_base / "tables" / f"psd_by_channel_{condition}.csv"
    df = pd.read_csv(p, index_col=0)
    df.index.name = "channel"
    df.columns = df.columns.astype(float)
    return df


def _load_band_power(nm_base: Path, mb_base: Path, condition: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load band power summary tables from both pipelines."""
    nm = pd.read_csv(nm_base / "band_power_summary.csv", index_col=0)
    mb = pd.read_csv(mb_base / "tables" / f"band_power_summary_{condition}.csv", index_col=0)
    nm.index = nm.index.str.strip()
    mb.index = mb.index.str.strip()
    return nm, mb


# ─────────────────────────────────────────────────────────────────────────────
# Statistics
# ─────────────────────────────────────────────────────────────────────────────

def _metrics(x: np.ndarray, y: np.ndarray) -> dict:
    """Compute Pearson r, Spearman ρ, MAE, RMSE, mean bias between two vectors."""
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 3:
        return {"pearson_r": None, "spearman_r": None, "mae": None, "rmse": None, "bias": None, "n": len(x)}
    r, p_r = stats.pearsonr(x, y)
    rho, p_rho = stats.spearmanr(x, y)
    diff = y - x
    return {
        "pearson_r": round(float(r), 4),
        "pearson_p": round(float(p_r), 6),
        "spearman_r": round(float(rho), 4),
        "spearman_p": round(float(p_rho), 6),
        "mae": round(float(np.mean(np.abs(diff))), 6),
        "rmse": round(float(np.sqrt(np.mean(diff ** 2))), 6),
        "bias": round(float(np.mean(diff)), 6),      # mean(MNE − Julia)
        "bias_pct": round(float(np.mean(diff) / (np.mean(x) + 1e-12) * 100), 2),
        "n": int(len(x)),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — wPLI scatter + Bland-Altman per band
# ─────────────────────────────────────────────────────────────────────────────

_BAND_ORDER = ["DELTA", "THETA", "ALPHA", "BETA_LOW", "BETA_MID", "BETA_HIGH", "GAMMA"]
_BAND_HZ = {
    "DELTA": "0.5–4 Hz", "THETA": "4–8 Hz", "ALPHA": "7.8–11.7 Hz",
    "BETA_LOW": "12–15 Hz", "BETA_MID": "15–18 Hz", "BETA_HIGH": "18–30 Hz", "GAMMA": "30–50 Hz",
}
_CMAP_BAND = ["#7B68EE", "#4169E1", "#228B22", "#DAA520", "#FF8C00", "#DC143C", "#9400D3"]


def _fig_wpli_scatter_ba(
    nm_base: Path, mb_base: Path, condition: str, report_dir: Path, tag: str
) -> dict:
    """Scatter + Bland-Altman grid for all 7 bands. Returns per-band metrics dict."""
    n_bands = len(_BAND_ORDER)
    fig, axes = plt.subplots(2, n_bands, figsize=(n_bands * 3.0, 6.5))
    fig.suptitle(
        f"wPLI comparison: NeuroMIND (Julia) vs mne_brain (MNE-Python)\n"
        f"sub-{tag.split('_')[0]} / {condition} — upper-triangle pairs only",
        fontsize=10, fontweight="bold", y=1.01,
    )

    all_metrics: dict[str, dict] = {}

    for k, (band, color) in enumerate(zip(_BAND_ORDER, _CMAP_BAND)):
        df_nm = _load_nm_wpli(nm_base, band)
        df_mb = _load_mb_wpli(mb_base, band, condition)
        df_nm, df_mb = _align_channels(df_nm, df_mb)
        nm_v = _upper_triangle(df_nm)
        mb_v = _upper_triangle(df_mb)
        m = _metrics(nm_v, mb_v)
        all_metrics[band] = m

        # ── Scatter (row 0) ──────────────────────────────────────────────────
        ax = axes[0, k]
        ax.scatter(nm_v, mb_v, s=4, alpha=0.35, color=color, rasterized=True)
        lim = max(nm_v.max(), mb_v.max()) * 1.05
        ax.plot([0, lim], [0, lim], "k--", lw=0.8, alpha=0.6)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        r_txt = f"r={m['pearson_r']:.3f}" if m['pearson_r'] is not None else "r=n/a"
        ax.set_title(f"{band}\n{_BAND_HZ[band]}", fontsize=7.5, fontweight="bold")
        ax.text(0.05, 0.93, r_txt, transform=ax.transAxes, fontsize=7, va="top",
                color="black", bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))
        if k == 0:
            ax.set_ylabel("mne_brain (MNE)", fontsize=8)
        ax.set_xlabel("NeuroMIND (Julia)", fontsize=7)
        ax.tick_params(labelsize=6)

        # ── Bland-Altman (row 1) ──────────────────────────────────────────────
        ax2 = axes[1, k]
        mean_v = (nm_v + mb_v) / 2.0
        diff_v = mb_v - nm_v
        md = float(np.mean(diff_v))
        sd = float(np.std(diff_v))
        ax2.scatter(mean_v, diff_v, s=4, alpha=0.35, color=color, rasterized=True)
        ax2.axhline(md, color="red", lw=1.0, linestyle="--", label=f"bias={md:.3f}")
        ax2.axhline(md + 1.96 * sd, color="gray", lw=0.8, linestyle=":")
        ax2.axhline(md - 1.96 * sd, color="gray", lw=0.8, linestyle=":")
        ax2.axhline(0, color="black", lw=0.5)
        ax2.set_xlabel("Mean wPLI", fontsize=7)
        if k == 0:
            ax2.set_ylabel("Diff (MNE − Julia)", fontsize=8)
        ax2.tick_params(labelsize=6)
        bias_txt = f"bias={md:+.3f}\n±1.96σ={1.96*sd:.3f}"
        ax2.text(0.05, 0.97, bias_txt, transform=ax2.transAxes, fontsize=6, va="top",
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7))

    plt.tight_layout()
    out = report_dir / f"comparison_{tag}_{condition}_wpli.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")
    return all_metrics


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — Side-by-side heatmaps (ALPHA band, then full grid)
# ─────────────────────────────────────────────────────────────────────────────

def _fig_heatmaps(
    nm_base: Path, mb_base: Path, condition: str, report_dir: Path, tag: str
) -> None:
    n_bands = len(_BAND_ORDER)
    fig, axes = plt.subplots(2, n_bands, figsize=(n_bands * 3.0, 6.0))
    fig.suptitle(
        f"wPLI Heatmaps — NeuroMIND (top) vs mne_brain (bottom)  [{condition}]",
        fontsize=10, fontweight="bold",
    )

    for k, band in enumerate(sorted(_BAND_ORDER, key=lambda b: _BAND_ORDER.index(b))):
        df_nm = _load_nm_wpli(nm_base, band)
        df_mb = _load_mb_wpli(mb_base, band, condition)
        df_nm_a, df_mb_a = _align_channels(df_nm, df_mb)
        vmax = max(df_nm_a.values.max(), df_mb_a.values.max())

        for row, (mat, label) in enumerate([(df_nm_a.values, "NeuroMIND"), (df_mb_a.values, "mne_brain")]):
            ax = axes[row, k]
            im = ax.imshow(mat, cmap="RdYlBu_r", vmin=0, vmax=vmax, aspect="auto")
            ax.set_xticks([]); ax.set_yticks([])
            title_pre = f"{band}" if row == 0 else ""
            if title_pre:
                ax.set_title(title_pre, fontsize=7.5, fontweight="bold")
            if k == 0:
                ax.set_ylabel(label, fontsize=7.5)

    plt.tight_layout()
    out = report_dir / f"comparison_{tag}_{condition}_heatmaps.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — PSD overlay + log-ratio
# ─────────────────────────────────────────────────────────────────────────────

def _fig_psd(
    nm_base: Path, mb_base: Path, condition: str, report_dir: Path, tag: str
) -> dict:
    nm_wide = _load_nm_psd(nm_base)
    mb_wide = _load_mb_psd(mb_base, condition)

    # Common channels
    common_ch = [c for c in nm_wide.index if c in mb_wide.index]
    nm_wide = nm_wide.loc[common_ch]
    mb_wide = mb_wide.loc[common_ch]

    # Mean PSD across channels
    nm_mean = nm_wide.mean(axis=0)   # Series, index = freq
    mb_mean = mb_wide.mean(axis=0)

    # Interpolate mne_brain onto NeuroMIND freq grid (0.5–150 Hz overlap)
    nm_freqs = nm_mean.index.values.astype(float)
    mb_freqs = mb_mean.index.values.astype(float)
    f_common = nm_freqs[(nm_freqs >= mb_freqs.min()) & (nm_freqs <= mb_freqs.max())]
    mb_interp = np.interp(f_common, mb_freqs, mb_mean.values)
    nm_sel = np.interp(f_common, nm_freqs, nm_mean.values)

    # Log ratio (avoid div-by-zero)
    log_ratio = 10 * np.log10((mb_interp + 1e-15) / (nm_sel + 1e-15))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle(
        f"PSD comparison — NeuroMIND vs mne_brain  [{condition}] — mean across {len(common_ch)} channels",
        fontsize=10, fontweight="bold",
    )

    # Top: overlay (log scale)
    ax1.semilogy(nm_mean.index, nm_mean.values, lw=1.2, color="#1f77b4", label="NeuroMIND (Julia)")
    ax1.semilogy(mb_mean.index, mb_mean.values, lw=1.2, color="#ff7f0e", alpha=0.85, label="mne_brain (MNE)")
    for band, (lo, hi) in [("δ", (0.5, 4)), ("θ", (4, 8)), ("α", (7.8, 11.7)),
                            ("β", (12, 30)), ("γ", (30, 50))]:
        ax1.axvspan(lo, hi, alpha=0.06, color="gray")
        ax1.text((lo + hi) / 2, ax1.get_ylim()[0] if ax1.get_ylim()[0] > 0 else 1e-4,
                 band, ha="center", fontsize=7, color="gray")
    ax1.set_xlim(0.5, 50)
    ax1.set_ylabel("PSD (µV²/Hz)", fontsize=9)
    ax1.legend(fontsize=8)
    ax1.grid(True, which="both", alpha=0.3)
    ax1.tick_params(labelsize=8)

    # Bottom: log-ratio
    ax2.plot(f_common, log_ratio, lw=1.0, color="#2ca02c")
    ax2.axhline(0, color="black", lw=0.8)
    ax2.axhline(3, color="red", lw=0.6, linestyle="--", alpha=0.6, label="+3 dB")
    ax2.axhline(-3, color="red", lw=0.6, linestyle="--", alpha=0.6, label="−3 dB")
    ax2.set_xlim(0.5, 50)
    ax2.set_ylabel("Log ratio (dB)\n10·log₁₀(MNE/Julia)", fontsize=8)
    ax2.set_xlabel("Frequency (Hz)", fontsize=9)
    ax2.legend(fontsize=7)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(labelsize=8)

    plt.tight_layout()
    out = report_dir / f"comparison_{tag}_{condition}_psd.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")

    # Metrics in 0.5–50 Hz range
    m = _metrics(nm_sel, mb_interp)
    return {"psd_pearson_r": m["pearson_r"], "psd_mae_uv2": m["mae"],
            "psd_bias_db_mean": round(float(log_ratio.mean()), 3)}


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — Band power bar chart
# ─────────────────────────────────────────────────────────────────────────────

def _fig_band_power(
    nm_base: Path, mb_base: Path, condition: str, report_dir: Path, tag: str
) -> dict:
    nm_bp, mb_bp = _load_band_power(nm_base, mb_base, condition)
    # Common channels and bands
    common_ch = [c for c in nm_bp.index if c in mb_bp.index]
    common_bands = [b for b in _BAND_ORDER if b in nm_bp.columns and b in mb_bp.columns]

    nm_bp = nm_bp.loc[common_ch, common_bands]
    mb_bp = mb_bp.loc[common_ch, common_bands]

    # Mean across channels
    nm_means = nm_bp.mean(axis=0)
    mb_means = mb_bp.mean(axis=0)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        f"Band Power comparison — NeuroMIND vs mne_brain  [{condition}]\n"
        f"Mean across {len(common_ch)} channels (µV²/Hz)",
        fontsize=10, fontweight="bold",
    )

    # Panel A: grouped bar chart
    ax = axes[0]
    x = np.arange(len(common_bands))
    w = 0.38
    ax.bar(x - w/2, nm_means.values, w, label="NeuroMIND (Julia)", color="#1f77b4", alpha=0.85)
    ax.bar(x + w/2, mb_means.values, w, label="mne_brain (MNE)", color="#ff7f0e", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(common_bands, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("Mean Band Power (µV²/Hz)", fontsize=9)
    ax.legend(fontsize=8)
    ax.set_title("Mean band power per pipeline", fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    # Panel B: ratio MNE / Julia per band
    ax2 = axes[1]
    ratios = mb_means.values / (nm_means.values + 1e-15)
    colors_b = [c for c in _CMAP_BAND[:len(common_bands)]]
    bars = ax2.bar(x, ratios, color=colors_b, alpha=0.85)
    ax2.axhline(1.0, color="black", lw=1.0, linestyle="--")
    ax2.set_xticks(x)
    ax2.set_xticklabels(common_bands, rotation=35, ha="right", fontsize=8)
    ax2.set_ylabel("Ratio (mne_brain / NeuroMIND)", fontsize=9)
    ax2.set_title("Power ratio per band (1.0 = perfect match)", fontsize=9)
    ax2.grid(axis="y", alpha=0.3)
    for bar, r in zip(bars, ratios):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                 f"{r:.2f}×", ha="center", va="bottom", fontsize=7)

    plt.tight_layout()
    out = report_dir / f"comparison_{tag}_{condition}_band_power.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.name}")

    bp_metrics = {
        band: {
            "nm_mean": round(float(nm_means[band]), 6),
            "mb_mean": round(float(mb_means[band]), 6),
            "ratio": round(float(mb_means[band] / (nm_means[band] + 1e-15)), 4),
        }
        for band in common_bands
    }
    return bp_metrics


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject", default="M05")
    parser.add_argument("--session", default="T2")
    parser.add_argument("--condition", default="EC")
    args = parser.parse_args()

    cfg = load_config()
    nm_base, mb_base, report_dir = _get_paths(cfg, args.subject, args.session, args.condition)

    # Verify inputs exist
    for label, path in [("NeuroMIND", nm_base), ("mne_brain", mb_base)]:
        if not path.exists():
            raise FileNotFoundError(
                f"{label} results not found: {path}\n"
                "Run the respective pipeline first."
            )

    tag = f"{args.subject}_{args.session}"
    cond = args.condition

    print(f"\n{'═'*55}")
    print(f"  Phase 8 — Validation: NeuroMIND vs mne_brain")
    print(f"  sub-{args.subject} / ses-{args.session} / {cond}")
    print(f"{'═'*55}\n")

    # ── wPLI scatter + Bland-Altman ────────────────────────────────────────
    print("  [1/4] wPLI scatter + Bland-Altman…")
    wpli_metrics = _fig_wpli_scatter_ba(nm_base, mb_base, cond, report_dir, tag)

    # ── Heatmap side-by-side ───────────────────────────────────────────────
    print("  [2/4] Heatmap comparison…")
    _fig_heatmaps(nm_base, mb_base, cond, report_dir, tag)

    # ── PSD overlay ───────────────────────────────────────────────────────
    print("  [3/4] PSD overlay…")
    psd_metrics = _fig_psd(nm_base, mb_base, cond, report_dir, tag)

    # ── Band power ────────────────────────────────────────────────────────
    print("  [4/4] Band power comparison…")
    bp_metrics = _fig_band_power(nm_base, mb_base, cond, report_dir, tag)

    # ── Summary JSON ──────────────────────────────────────────────────────
    summary = {
        "subject": args.subject,
        "session": args.session,
        "condition": cond,
        "nm_pipeline": "NeuroMIND (Julia) — Hilbert time-domain wPLI",
        "mb_pipeline": "mne_brain (MNE-Python) — CSD frequency-domain wPLI",
        "note": (
            "NeuroMIND uses 30 channels (Fp2 excluded as bad); "
            "mne_brain uses 31 channels — comparison restricted to 30 shared channels. "
            "wPLI estimators differ mathematically: Hilbert vs CSD, but r > 0.85 is expected."
        ),
        "wpli_per_band": wpli_metrics,
        "psd": psd_metrics,
        "band_power_per_band": bp_metrics,
    }
    out_json = report_dir / f"comparison_{tag}_{cond}_summary.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"  → {out_json.name}")

    # ── Console summary ───────────────────────────────────────────────────
    print(f"\n{'─'*55}")
    print("  wPLI Pearson r per band:")
    for band, m in wpli_metrics.items():
        r = m.get("pearson_r")
        bias = m.get("bias")
        r_str = f"{r:.3f}" if r is not None else "n/a"
        b_str = f"{bias:+.3f}" if bias is not None else "n/a"
        print(f"    {band:<12}  r={r_str}   bias(MNE−Julia)={b_str}")

    print(f"\n  PSD:  r={psd_metrics.get('psd_pearson_r', 'n/a')}  "
          f"mean_ratio={psd_metrics.get('psd_bias_db_mean', 'n/a')} dB")

    print(f"\n  Reports saved → {report_dir}")
    print(f"{'═'*55}\n")


if __name__ == "__main__":
    main()
