"""Phase 7 — wPLI functional connectivity for a single recording.

Default pilot: sub-M05 / ses-T2 / EC.
Reads the native .fif epochs from Phase 5.

Usage
-----
python scripts/run_phase7_m05.py
python scripts/run_phase7_m05.py --subject M07 --session T1 --condition EC

Outputs
-------
results/subjects/sub-{id}/ses-{sess}/{task}/
├── cache/connectivity_{cond}.npz
└── tables/
    ├── wpli_{BAND}_{cond}.csv      (one 31×31 matrix per band)
    ├── connectivity_edges_{cond}.csv
    └── connectivity_summary_{cond}.json
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import mne

from mne_brain import load_config
from mne_brain.connectivity import compute_wpli, connectivity_matrices, save_connectivity_results


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
    print(f"Epochs: {len(epochs)} × {epochs.get_data().shape[-1]} samples "
          f"({len(epochs.ch_names)} ch @ {epochs.info['sfreq']} Hz)")

    # ── 2. Compute wPLI (native mne-connectivity) ─────────────────────────────
    print("Computing wPLI across all bands…  (may take ~30 s for 100 epochs × 31 ch)")
    con, band_names = compute_wpli(epochs, cfg)
    matrices = connectivity_matrices(con, band_names)
    print(f"SpectralConnectivity: {con.get_data(output='dense').shape}  "
          f"(n_ch × n_ch × n_bands)")

    # ── 3. Save outputs ───────────────────────────────────────────────────────
    summary = save_connectivity_results(con, band_names, cfg, base_dir, args.condition)

    # ── 4. Print band summary ─────────────────────────────────────────────────
    print(f"\n── Phase 7 complete ────────────────────────────────────")
    print(f"  n_epochs used  : {summary['n_epochs_used']}")
    print(f"  n_channels     : {summary['n_channels']}")
    print(f"\n  Mean wPLI per band:")
    for band, stats in summary["band_summary"].items():
        lo, hi = cfg.bands[band]
        print(f"    {band:<10} [{lo:4.1f}–{hi:5.1f} Hz]  "
              f"mean={stats['mean_wpli']:.4f}  max={stats['max_wpli']:.4f}")

    tables = base_dir / "tables"
    print(f"\n  Matrices       : {tables}/wpli_{{BAND}}_{args.condition}.csv")
    print(f"  Edges          : {tables / f'connectivity_edges_{args.condition}.csv'}")
    print(f"  Summary        : {tables / f'connectivity_summary_{args.condition}.json'}")


if __name__ == "__main__":
    main()
