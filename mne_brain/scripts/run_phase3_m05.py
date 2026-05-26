"""Run Phase 3 QC + filtering for the M05/T2/EC pilot recording."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from mne_brain import load_config
from mne_brain.bids import load_eeg_bids
from mne_brain.preprocessing import describe_filter_chain, filter_recording, qc_report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject", default="M05")
    parser.add_argument("--session", default="T2")
    parser.add_argument("--condition", default="EC")
    parser.add_argument(
        "--source",
        default="brainvision",
        choices=["auto", "brainvision", "tsv"],
        help="Use brainvision for the strict MNE path; tsv is only a fallback/debug route.",
    )
    args = parser.parse_args()

    cfg = load_config()
    rec = load_eeg_bids(
        cfg,
        args.subject,
        args.session,
        args.condition,
        source=args.source,
    )

    task = "eyesclosed" if args.condition == "EC" else "eyesopen"
    out_dir = (
        Path(cfg.root)
        / str(cfg.paths.get("results", "results"))
        / "subjects"
        / f"sub-{args.subject}"
        / f"ses-{args.session}"
        / task
    )
    tables_dir = out_dir / "tables"
    cache_dir = out_dir / "cache"
    tables_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    qc = qc_report(rec, cfg, output_dir=tables_dir)
    filtered = filter_recording(rec, cfg)

    np.savez_compressed(
        cache_dir / f"filtered_{args.condition}.npz",
        data=filtered.data,
        times=filtered.times,
        fs=filtered.meta.fs,
        channel_names=np.asarray(filtered.meta.channel_names, dtype=object),
    )

    # Save full recording for all stages (raw + pre-ICA filtered)
    # post_ica_uv will be added by run_phase4_m05.py after ICA rejection
    preview_cols = filtered.data.shape[1]   # full recording (e.g. 50000 = 100s @ 500Hz)
    preview = {
        "channel_names": filtered.meta.channel_names,
        "times": filtered.times[:preview_cols].round(6).tolist(),
        "raw_uv":      rec.data[:, :preview_cols].round(6).tolist(),
        "pre_ica_uv":  filtered.data[:, :preview_cols].round(6).tolist(),
    }
    (tables_dir / f"signal_preview_{args.condition}.json").write_text(
        json.dumps(preview),
        encoding="utf-8",
    )
    (tables_dir / f"filter_chain_{args.condition}.json").write_text(
        json.dumps(describe_filter_chain(cfg), indent=2),
        encoding="utf-8",
    )

    print(f"Loaded: {rec.data.shape[0]} ch x {rec.data.shape[1]} samples @ {rec.meta.fs} Hz")
    print(f"QC rows: {len(qc)}; bad channels: {qc.loc[qc['is_bad'], 'channel'].to_list()}")
    print(f"Filtered cache: {cache_dir / f'filtered_{args.condition}.npz'}")


if __name__ == "__main__":
    main()

