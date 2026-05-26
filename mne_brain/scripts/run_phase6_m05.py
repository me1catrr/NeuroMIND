"""Phase 6 — PSD and band-power for a single recording.

Default pilot: sub-M05 / ses-T2 / EC.
Reads the native .fif epochs produced by Phase 5.

Usage
-----
python scripts/run_phase6_m05.py
python scripts/run_phase6_m05.py --subject M07 --session T1 --condition EC

Outputs
-------
results/subjects/sub-{id}/ses-{sess}/{task}/
├── cache/spectrum_{cond}.npz
└── tables/
    ├── psd_by_channel_{cond}.csv
    ├── band_power_summary_{cond}.csv
    └── spectral_params_{cond}.json
figures/psd_spectrum_{cond}.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running directly from the project root without `pip install -e .`
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import mne

from mne_brain import load_config
from mne_brain.spectral import compute_psd, save_spectral_results


def _resolve_dirs(cfg, subject: str, session: str, condition: str) -> tuple[Path, Path]:
    task = "eyesclosed" if condition == "EC" else "eyesopen"
    base = (
        Path(cfg.root)
        / str(cfg.paths.get("results", "results"))
        / "subjects"
        / f"sub-{subject}"
        / f"ses-{session}"
        / task
    )
    return base, base / "cache"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject", default="M05")
    parser.add_argument("--session", default="T2")
    parser.add_argument("--condition", default="EC")
    args = parser.parse_args()

    cfg = load_config()
    base_dir, cache_dir = _resolve_dirs(cfg, args.subject, args.session, args.condition)

    # ── 1. Load epochs from Phase 5 native .fif ──────────────────────────────
    fif_path = cache_dir / f"epochs_{args.condition}-epo.fif"
    if not fif_path.is_file():
        raise FileNotFoundError(
            f"Epochs cache not found: {fif_path}\n"
            "Run run_phase5_m05.py first."
        )
    epochs = mne.read_epochs(str(fif_path), preload=True, verbose=False)
    print(f"Epochs loaded: {len(epochs)} × {epochs.get_data().shape[-1]} samples "
          f"@ {epochs.info['sfreq']} Hz  ({len(epochs.ch_names)} channels)")

    # ── 2. Compute PSD (native MNE Welch) ────────────────────────────────────
    spectrum = compute_psd(epochs, cfg)
    print(f"EpochsSpectrum: {spectrum.get_data().shape}  "
          f"df={spectrum.freqs[1]-spectrum.freqs[0]:.4f} Hz  "
          f"[{spectrum.freqs[0]:.2f}–{spectrum.freqs[-1]:.2f} Hz]")

    # ── 3. Save outputs ───────────────────────────────────────────────────────
    params = save_spectral_results(spectrum, cfg, base_dir, args.condition)

    # ── 4. Quick band-power summary to stdout ────────────────────────────────
    mean_sp = spectrum.average()
    print(f"\n── Phase 6 complete ────────────────────────────────────")
    print(f"  Epochs used    : {params['n_epochs_used']}")
    print(f"  Freq resolution: {params['df_hz']} Hz")
    print(f"  n_freqs        : {params['n_freqs']}  ({params['fmin_hz']}–{params['fmax_hz']} Hz)")
    print(f"\n  Band-power means across channels (µV²/Hz):")
    for band, (lo, hi) in cfg.bands.items():
        try:
            bp = mean_sp.get_data(fmin=lo, fmax=hi).mean() * 1e12
            print(f"    {band:<10} [{lo:4.1f}–{hi:5.1f} Hz]  {bp:.4e}")
        except Exception:
            print(f"    {band:<10}  (error)")
    tables = base_dir / "tables"
    print(f"\n  PSD table      : {tables / f'psd_by_channel_{args.condition}.csv'}")
    print(f"  Band power     : {tables / f'band_power_summary_{args.condition}.csv'}")
    print(f"  Figure         : {base_dir / 'figures' / f'psd_spectrum_{args.condition}.png'}")


if __name__ == "__main__":
    main()
