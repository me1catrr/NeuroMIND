"""Loader for the NeuroMIND BIDS-EEG dataset (`NeuroMIND/data/bids/`)."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from mne_brain.common.types import EEGRecording, PipelineConfig, RecordingMeta


TASK_BY_CONDITION = {
    "EC": "eyesclosed",
    "EO": "eyesopen",
    "ec": "eyesclosed",
    "eo": "eyesopen",
    "eyesclosed": "eyesclosed",
    "eyesopen": "eyesopen",
}

CONDITION_BY_TASK = {
    "eyesclosed": "EC",
    "eyesopen": "EO",
}


@dataclass(frozen=True)
class BrainVisionHeader:
    n_channels: int
    fs: float
    ch_names: list[str]
    resolutions: list[float]
    binary_format: str
    orientation: str
    eeg_file: str


def bids_eeg_dir(bids_root: str | Path, subject_id: str, session_id: str) -> Path:
    """`{bids_root}/sub-X/ses-Y/eeg`, where each session keeps its recordings."""
    return Path(bids_root) / f"sub-{subject_id}" / f"ses-{session_id}" / "eeg"


def bids_vhdr_path(
    bids_root: str | Path, subject_id: str, session_id: str, task: str, run: int = 1
) -> Path:
    prefix = f"sub-{subject_id}_ses-{session_id}_task-{task}_run-{run:02d}"
    return bids_eeg_dir(bids_root, subject_id, session_id) / f"{prefix}_eeg.vhdr"


def bids_electrodes_path(bids_root: str | Path, subject_id: str, session_id: str) -> Path:
    return bids_eeg_dir(bids_root, subject_id, session_id) / (
        f"sub-{subject_id}_ses-{session_id}_electrodes.tsv"
    )


def load_eeg_bids(
    cfg: PipelineConfig,
    subject_id: str,
    session_id: str,
    condition: str,
    *,
    run: int = 1,
) -> EEGRecording:
    """Load a recording of the BIDS-EEG dataset as an `EEGRecording` in microvolts.

    The BrainVision trio lives in the BIDS tree itself, so no side metadata or
    absolute paths are needed. Recordings cut from a single raw file (M04, M05
    and M06 at T1) are already cropped in the dataset.
    """
    task = normalize_task(condition)
    bids_root = bids_root_dir(cfg)
    vhdr_path = bids_vhdr_path(bids_root, subject_id, session_id, task, run)
    if not vhdr_path.is_file():
        raise FileNotFoundError(f"BIDS EEG recording not found: {vhdr_path}")

    channel_positions = load_electrode_positions(
        bids_electrodes_path(bids_root, subject_id, session_id)
    )
    return load_eeg_brainvision(
        vhdr_path,
        subject_id,
        session_id,
        task,
        run=run,
        channel_positions=channel_positions,
    )


def load_eeg_brainvision(
    vhdr_path: str | Path,
    subject_id: str,
    session_id: str,
    task: str,
    *,
    run: int = 1,
    channel_positions: dict[str, tuple[float, float]] | None = None,
) -> EEGRecording:
    """Load BrainVision data with MNE and return data in microvolts.

    MNE's BrainVision reader handles the `.vhdr/.eeg/.vmrk` trio and returns EEG
    data in volts. NeuroMIND stores data in microvolts, so this function
    multiplies by `1e6` before constructing `EEGRecording`.
    """
    vhdr = Path(vhdr_path).expanduser().resolve()
    if not vhdr.is_file():
        raise FileNotFoundError(f"BrainVision .vhdr not found: {vhdr}")

    try:
        import mne
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "MNE is required to load BrainVision files. Install with `pip install -e .` "
            "or `pip install mne`."
        ) from exc

    raw = mne.io.read_raw_brainvision(vhdr, preload=True, scale=1.0, verbose="ERROR")
    raw.set_channel_types({name: "eeg" for name in raw.ch_names}, verbose="ERROR")
    try:
        raw.set_montage("standard_1020", on_missing="ignore", verbose="ERROR")
    except Exception:
        pass

    data_uv = np.asarray(raw.get_data() * 1e6, dtype=np.float64)
    fs = float(raw.info["sfreq"])
    ch_names = list(raw.ch_names)
    condition = CONDITION_BY_TASK.get(task, task)

    meta = RecordingMeta(
        subject_id=subject_id,
        session_id=session_id,
        condition=condition,
        run=run,
        fs=fs,
        n_channels=len(ch_names),
        channel_names=ch_names,
        channel_positions=channel_positions,
        bids_path=str(vhdr),
    )
    times = np.arange(data_uv.shape[1], dtype=np.float64) / fs
    return EEGRecording(meta=meta, data=data_uv, times=times)


def read_vhdr_header(vhdr_path: str | Path) -> BrainVisionHeader:
    """Parse the BrainVision header fields used by NeuroMIND."""
    path = Path(vhdr_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"BrainVision .vhdr not found: {path}")

    n_channels = 0
    fs = 500.0
    ch_names: list[str] = []
    resolutions: list[float] = []
    binary_format = "IEEE_FLOAT_32"
    orientation = "MULTIPLEXED"
    eeg_file = ""

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue

        if line.startswith("DataFile="):
            eeg_file = line.split("=", 1)[1].strip()
        elif line.startswith("NumberOfChannels="):
            n_channels = int(line.split("=", 1)[1].strip())
        elif line.startswith("SamplingInterval="):
            interval_us = float(line.split("=", 1)[1].strip())
            fs = 1_000_000.0 / interval_us
        elif line.startswith("BinaryFormat="):
            binary_format = line.split("=", 1)[1].strip()
        elif line.startswith("DataOrientation="):
            orientation = line.split("=", 1)[1].strip()
        else:
            match = re.match(r"^Ch(\d+)=(.+)$", line)
            if match:
                parts = match.group(2).split(",")
                ch_names.append(parts[0])
                resolution = float(parts[2].strip()) if len(parts) >= 3 and parts[2].strip() else 1.0
                resolutions.append(resolution)

    if n_channels > 0:
        while len(ch_names) < n_channels:
            ch_names.append(f"CH{len(ch_names) + 1}")
        while len(resolutions) < n_channels:
            resolutions.append(1.0)
        ch_names = ch_names[:n_channels]
        resolutions = resolutions[:n_channels]

    return BrainVisionHeader(
        n_channels=n_channels,
        fs=fs,
        ch_names=ch_names,
        resolutions=resolutions,
        binary_format=binary_format,
        orientation=orientation,
        eeg_file=eeg_file,
    )


def load_electrode_positions(tsv_path: str | Path) -> dict[str, tuple[float, float]] | None:
    """Load 2D electrode positions from a BIDS `*_electrodes.tsv`.

    Rows whose optional `type` column is not EEG (e.g. legacy REF/GND rows) are skipped.
    """
    path = Path(tsv_path).expanduser().resolve()
    if not path.is_file():
        return None

    positions: dict[str, tuple[float, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("type", "EEG") != "EEG":
                continue
            name = (row.get("name") or "").strip().upper()
            if not name:
                continue
            positions[name] = (float(row["x"]), float(row["y"]))
    return positions


def normalize_task(condition: str) -> str:
    try:
        return TASK_BY_CONDITION[condition]
    except KeyError as exc:
        raise ValueError(f"Unknown condition/task: {condition!r}") from exc


def bids_root_dir(cfg: PipelineConfig) -> Path:
    """Resolve the BIDS root, accepting paths relative to mne_brain or NeuroMIND.

    The default is `../data/bids`, the NeuroMIND dataset next to this project.
    """
    configured = Path(str(cfg.paths["bids_root"]))
    candidates = (
        [configured]
        if configured.is_absolute()
        else [cfg.root / configured, cfg.root.parent / configured]
    )
    for candidate in candidates:
        if (candidate / "dataset_description.json").is_file():
            return candidate.resolve()
    return candidates[0].resolve()
