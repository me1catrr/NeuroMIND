"""Build lightweight BIDS metadata for every BrainVision recording.

Reads the original BrainVision `.vhdr` files (described in
`NeuroMIND/data/full_data/inventory.csv`) with `mne.io.read_raw_brainvision`
and writes one JSON metadata sidecar per recording into
`mne_brain/data/BIDS/raw/`. Binary EEG data is **not** copied — the JSON only
points to the original `.vhdr` via `vhdr_path`.

It also mirrors the subject-level support files from
`NeuroMIND/data/BIDS/` (participants.tsv, groups.csv, longitudinal_pairs.csv,
dataset_description.json) into `mne_brain/data/BIDS/` so the mne_brain dataset
is self-contained.

Usage
-----
python scripts/build_bids_full.py
python scripts/build_bids_full.py --groups MS              # only MS subjects
python scripts/build_bids_full.py --limit 10 --dry-run
python scripts/build_bids_full.py --overwrite              # force rewrite

Outputs
-------
mne_brain/data/BIDS/
  raw/sub-{bids_id}_ses-{bids_session}_task-{task}_run-01_eeg_metadata.json
  participants.tsv          (copied)
  groups.csv                (copied)
  longitudinal_pairs.csv    (copied)
  dataset_description.json  (copied)
results/bids_build_log.csv  (one row per recording, status + reason)
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

import mne


# ─────────────────────────────────────────────────────────────────────────────
# Paths and constants
# ─────────────────────────────────────────────────────────────────────────────

MNE_BRAIN_ROOT = Path(__file__).resolve().parents[1]
NEUROMIND_ROOT = MNE_BRAIN_ROOT.parent

INVENTORY_CSV = NEUROMIND_ROOT / "data" / "full_data" / "inventory.csv"
RAW_DATA_ROOT = (
    NEUROMIND_ROOT
    / "data"
    / "full_data"
    / "Pacientes MINDEM_IMIBIC_27 03 25"
)
NEUROMIND_BIDS_DIR = NEUROMIND_ROOT / "data" / "BIDS"

BIDS_OUT_DIR = MNE_BRAIN_ROOT / "data" / "BIDS"
RAW_OUT_DIR = BIDS_OUT_DIR / "raw"
BUILD_LOG = MNE_BRAIN_ROOT / "results" / "bids_build_log.csv"

TASK_BY_CONDITION = {"EC": "eyesclosed", "EO": "eyesopen"}

# Hardware/software template applied to every recording. Reflects the values
# documented in NeuroMIND/data/BIDS/BrainVision_hardware_software_metadata.md
# for the 210/212 recordings whose [Comment] block is present.
HARDWARE_TEMPLATE: dict[str, Any] = {
    "Manufacturer": "Brain Products",
    "ManufacturersModelName": "actiCHamp",
    "SoftwareVersions": "BrainVision Recorder Professional V. 1.21.0303",
    "AmplifierBaseUnit": "actiCHamp Base Unit (5001), S/N 16080673",
    "AmplifierModule": "actiCHamp 32 CH Module, Module 1 (5010), S/N 16091350",
    "PowerLineFrequency": 50,
    "EEGChannelCount": 31,
    "EEGReference": "REF physical channel 2",
    "EEGPlacementScheme": "10-20",
    "EEGGround": "Fpz",
    "RecordingType": "continuous",
    "HardwareFilters": {
        "HighPassFilter": {"Type": "DC"},
        "LowPassFilter": {"Cutoff": 140, "Units": "Hz"},
        "NotchFilter": {"Type": "Off"},
    },
    "SoftwareFilters": {
        "HighPassFilter": {"Cutoff": 1.59155, "Units": "s"},
        "LowPassFilter": {"Cutoff": 70, "Units": "Hz"},
        "NotchFilter": {"Cutoff": 50, "Units": "Hz"},
    },
    "ChannelResolution_uV": 0.0488281,
}

# Recordings whose extended [Comment] block is missing — listed in the
# hardware metadata report. We still treat them as standard, but flag it.
MISSING_COMMENT_BLOCK = {
    "M5_T2_MLLERGAS_011221_OJOS CERRADOS.vhdr",
    "M12_T2_MPCP_180122_ojoscerrados.vhdr",
}

# Recordings with [Comment] block but no stored impedance values.
MISSING_IMPEDANCES = {
    "M11_T2_RSB_190122_OJOSACIERTOS.vhdr",
    "M11_T2_rsb_190122_ojoscerrados.vhdr",
    "M14_T2_RAF_190122_OJOS ABIERTOS.vhdr",
    "M14_T2_raf_190122_ojos cerrados.vhdr",
    "MC27_Imc_260323_ojoscerrados.vhdr",
}


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────


def _section(title: str) -> None:
    print(f"\n{'─' * 60}")
    print(f"  {title}")
    print("─" * 60)


def load_participants(bids_dir: Path) -> dict[str, dict[str, Any]]:
    """Map `sub-{bids_id}` → dict with group / sex / age / education_level."""
    path = bids_dir / "participants.tsv"
    if not path.is_file():
        return {}
    df = pd.read_csv(path, sep="\t")
    info: dict[str, dict[str, Any]] = {}
    for _, row in df.iterrows():
        pid = str(row["participant_id"])
        info[pid] = {
            "group": str(row.get("group", "")),
            "sex": str(row.get("sex", "")),
            "age": int(row["age"]) if pd.notna(row.get("age")) else None,
            "education_level": (
                int(row["education_level"])
                if pd.notna(row.get("education_level"))
                else None
            ),
            "has_t1": bool(row.get("has_t1", False)),
            "has_t2": bool(row.get("has_t2", False)),
        }
    return info


def load_inventory(csv_path: Path) -> pd.DataFrame:
    """Read inventory.csv and normalise dtypes."""
    df = pd.read_csv(csv_path, dtype=str).fillna("")
    df["excluded_bool"] = df["excluded"].str.lower() == "true"
    return df


def resolve_vhdr_path(raw_data_root: Path, filename: str, fallback: str) -> Path | None:
    """Resolve a `.vhdr` file under raw_data_root by filename."""
    candidate = raw_data_root / filename
    if candidate.is_file():
        return candidate.resolve()
    if fallback:
        fb = Path(fallback)
        if fb.is_file():
            return fb.resolve()
    return None


def read_vhdr_with_mne(vhdr_path: Path) -> tuple[float, int, list[str]]:
    """Read header via MNE without loading the full data."""
    raw = mne.io.read_raw_brainvision(
        str(vhdr_path), preload=False, scale=1.0, verbose="ERROR"
    )
    fs = float(raw.info["sfreq"])
    ch_names = list(raw.ch_names)
    n_channels = len(ch_names)
    return fs, n_channels, ch_names


def build_metadata_dict(
    row: pd.Series,
    vhdr_path: Path,
    fs: float,
    n_channels: int,
    ch_names: list[str],
    participant: dict[str, Any],
    *,
    task: str,
    run: int = 1,
) -> dict[str, Any]:
    """Construct the metadata sidecar dictionary for a single recording."""
    bids_id = str(row["bids_id"])
    bids_session = str(row["bids_session"])
    filename = vhdr_path.name

    meta: dict[str, Any] = {
        "data_format": "brainvision",
        "vhdr_path": str(vhdr_path),
        "fs": fs,
        "n_channels": n_channels,
        "channel_names": ch_names,
        "subject": bids_id,
        "session": bids_session,
        "task": task,
        "run": run,
        "group": str(row["group"]) or participant.get("group", ""),
        "session_sub": str(row.get("session_sub", "")) or None,
        "condition_raw": str(row.get("condition_raw", "")),
        "date_str": str(row.get("date_str", "")) or None,
        "initials": str(row.get("initials", "")) or None,
        "SamplingFrequency": fs,
        **HARDWARE_TEMPLATE,
        "Participant": {
            "sex": participant.get("sex"),
            "age": participant.get("age"),
            "education_level": participant.get("education_level"),
        },
        "Source": {
            "original_filename": filename,
            "loader": "mne.io.read_raw_brainvision",
            "loader_version": mne.__version__,
            "comment_block_present": filename not in MISSING_COMMENT_BLOCK,
            "impedance_values_available": (
                filename not in MISSING_COMMENT_BLOCK
                and filename not in MISSING_IMPEDANCES
            ),
            "inventory_note": str(row.get("note", "")) or None,
        },
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    return meta


def write_metadata(meta: dict[str, Any], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(meta, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def copy_support_files(src_dir: Path, dst_dir: Path) -> list[str]:
    """Copy participants.tsv, groups.csv, longitudinal_pairs.csv, dataset_description.json."""
    copied: list[str] = []
    candidates = [
        "participants.tsv",
        "groups.csv",
        "longitudinal_pairs.csv",
        "dataset_description.json",
        "BrainVision_hardware_software_metadata.md",
    ]
    dst_dir.mkdir(parents=True, exist_ok=True)
    for name in candidates:
        src = src_dir / name
        if src.is_file():
            dst = dst_dir / name
            shutil.copy2(src, dst)
            copied.append(name)
    return copied


def write_build_log(rows: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "timestamp",
        "subject",
        "session",
        "task",
        "condition",
        "status",
        "vhdr_filename",
        "fs",
        "n_channels",
        "out_json",
        "elapsed_s",
        "reason",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


# ─────────────────────────────────────────────────────────────────────────────
# Main builder
# ─────────────────────────────────────────────────────────────────────────────


def build_all(
    *,
    groups: list[str] | None,
    sessions: list[str] | None,
    conditions: list[str] | None,
    subjects: list[str] | None,
    limit: int | None,
    overwrite: bool,
    dry_run: bool,
) -> tuple[int, int, int]:
    """Iterate the inventory and write metadata JSON files. Returns counts."""
    if not INVENTORY_CSV.is_file():
        raise FileNotFoundError(f"inventory.csv not found at {INVENTORY_CSV}")

    _section("mne_brain — build BIDS metadata for full dataset")
    print(f"  Inventory:      {INVENTORY_CSV}")
    print(f"  Raw data root:  {RAW_DATA_ROOT}")
    print(f"  Output dir:     {RAW_OUT_DIR}")
    print(f"  Overwrite:      {overwrite}")
    print(f"  Dry run:        {dry_run}")

    inv = load_inventory(INVENTORY_CSV)
    participants = load_participants(NEUROMIND_BIDS_DIR)

    if not dry_run:
        copied = copy_support_files(NEUROMIND_BIDS_DIR, BIDS_OUT_DIR)
        print(f"\n  Copied support files: {', '.join(copied) or 'none'}")

    filt = inv[~inv["excluded_bool"]]
    filt = filt[filt["condition"].isin(["EC", "EO"])]
    if groups:
        filt = filt[filt["group"].isin(groups)]
    if sessions:
        filt = filt[filt["bids_session"].isin(sessions)]
    if conditions:
        filt = filt[filt["condition"].isin(conditions)]
    if subjects:
        filt = filt[filt["bids_id"].isin(subjects)]
    if limit:
        filt = filt.head(limit)

    print(f"\n  Recordings to process: {len(filt)}")

    log_rows: list[dict[str, Any]] = []
    n_ok = 0
    n_skip = 0
    n_err = 0

    for _, row in filt.iterrows():
        t0 = time.time()
        bids_id = str(row["bids_id"])
        bids_session = str(row["bids_session"])
        condition = str(row["condition"])
        task = TASK_BY_CONDITION[condition]
        out_name = (
            f"sub-{bids_id}_ses-{bids_session}_task-{task}_run-01_eeg_metadata.json"
        )
        out_path = RAW_OUT_DIR / out_name

        log_row: dict[str, Any] = {
            "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "subject": bids_id,
            "session": bids_session,
            "task": task,
            "condition": condition,
            "vhdr_filename": str(row["filename"]),
            "out_json": str(out_path.relative_to(MNE_BRAIN_ROOT)),
        }

        if out_path.exists() and not overwrite:
            elapsed = time.time() - t0
            log_row.update(status="skip", reason="exists", elapsed_s=f"{elapsed:.3f}")
            log_rows.append(log_row)
            n_skip += 1
            continue

        vhdr = resolve_vhdr_path(
            RAW_DATA_ROOT, str(row["filename"]), str(row.get("filepath", ""))
        )
        if vhdr is None or not vhdr.is_file():
            elapsed = time.time() - t0
            log_row.update(
                status="error",
                reason=f".vhdr not found: {row['filename']}",
                elapsed_s=f"{elapsed:.3f}",
            )
            log_rows.append(log_row)
            n_err += 1
            print(f"  [!] {bids_id} {bids_session} {condition}: .vhdr not found")
            continue

        try:
            fs, n_channels, ch_names = read_vhdr_with_mne(vhdr)
        except Exception as exc:
            elapsed = time.time() - t0
            log_row.update(
                status="error",
                reason=f"read_raw_brainvision: {exc}",
                elapsed_s=f"{elapsed:.3f}",
            )
            log_rows.append(log_row)
            n_err += 1
            print(f"  [!] {bids_id} {bids_session} {condition}: MNE read error: {exc}")
            continue

        meta = build_metadata_dict(
            row,
            vhdr,
            fs,
            n_channels,
            ch_names,
            participants.get(f"sub-{bids_id}", {}),
            task=task,
        )

        if not dry_run:
            try:
                write_metadata(meta, out_path)
            except Exception as exc:
                elapsed = time.time() - t0
                log_row.update(
                    status="error",
                    reason=f"write_metadata: {exc}",
                    elapsed_s=f"{elapsed:.3f}",
                )
                log_rows.append(log_row)
                n_err += 1
                print(f"  [!] {bids_id} {bids_session} {condition}: write error: {exc}")
                continue

        elapsed = time.time() - t0
        log_row.update(
            status="ok" if not dry_run else "dry-run",
            fs=fs,
            n_channels=n_channels,
            elapsed_s=f"{elapsed:.3f}",
        )
        log_rows.append(log_row)
        n_ok += 1
        if n_ok % 20 == 0 or n_ok == 1:
            print(
                f"  [{n_ok:3d}] sub-{bids_id} ses-{bids_session} {condition}  "
                f"({n_channels} ch @ {fs:.0f} Hz)  [{elapsed:.2f}s]"
            )

    if not dry_run:
        write_build_log(log_rows, BUILD_LOG)
        print(f"\n  Build log: {BUILD_LOG.relative_to(MNE_BRAIN_ROOT)}")

    return n_ok, n_skip, n_err


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--groups",
        nargs="*",
        choices=["MS", "HC"],
        help="Restrict to one or more groups",
    )
    parser.add_argument(
        "--sessions",
        nargs="*",
        choices=["T1", "T2"],
        help="Restrict to one or more sessions",
    )
    parser.add_argument(
        "--conditions",
        nargs="*",
        choices=["EC", "EO"],
        help="Restrict to one or more conditions",
    )
    parser.add_argument(
        "--subjects",
        nargs="*",
        help="Restrict to specific bids_ids (e.g. M07 MC10)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process the first N matching recordings",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing metadata JSON files",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Do not write files, only print what would be done",
    )
    args = parser.parse_args()

    t_total = time.time()
    n_ok, n_skip, n_err = build_all(
        groups=args.groups,
        sessions=args.sessions,
        conditions=args.conditions,
        subjects=args.subjects,
        limit=args.limit,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
    )

    print(f"\n{'═' * 60}")
    print(f"  Summary  [{time.time() - t_total:.1f}s total]")
    print(f"{'═' * 60}")
    print(f"  Written / dry-run : {n_ok}")
    print(f"  Skipped (exists)  : {n_skip}")
    print(f"  Errors            : {n_err}")
    if n_err:
        print(f"\n  See {BUILD_LOG.relative_to(MNE_BRAIN_ROOT)} for details")


if __name__ == "__main__":
    main()
