"""BIDS and BrainVision loaders compatible with NeuroMIND outputs."""

from __future__ import annotations

import csv
import json
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

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


def load_eeg_bids(
    cfg: PipelineConfig,
    subject_id: str,
    session_id: str,
    condition: str,
    *,
    run: int = 1,
    source: str = "auto",
) -> EEGRecording:
    """Load a NeuroMIND BIDS recording as an `EEGRecording` in microvolts.

    The project currently has two supported layouts:

    - lightweight BIDS metadata pointing to the original BrainVision `.vhdr`
    - TSV exports with rows as channels and columns as samples
    """
    task = normalize_task(condition)
    bids_root = bids_root_dir(cfg)
    raw_dir = bids_root / "raw"
    prefix = f"sub-{subject_id}_ses-{session_id}_task-{task}_run-{run:02d}"

    metadata_path = _first_existing(
        raw_dir / f"{prefix}_eeg_metadata.json",
        raw_dir / f"{prefix}_metadata.json",
    )
    data_tsv_path = raw_dir / f"{prefix}_eeg_data.tsv"

    metadata: dict[str, Any] = {}
    if metadata_path is not None:
        metadata = _read_json(metadata_path)

    electrodes_path = bids_root / "electrodes" / f"sub-{subject_id}_ses-{session_id}_electrodes.tsv"
    channel_positions = load_electrode_positions(electrodes_path)

    source = source.lower()
    if source not in {"auto", "brainvision", "tsv"}:
        raise ValueError("source must be one of: auto, brainvision, tsv")

    vhdr_path = resolve_vhdr_path(cfg, metadata.get("vhdr_path"))
    if source in {"auto", "brainvision"} and vhdr_path:
        try:
            return load_eeg_brainvision(
                vhdr_path,
                subject_id,
                session_id,
                task,
                run=run,
                channel_positions=channel_positions,
            )
        except ModuleNotFoundError:
            if source == "brainvision" or not data_tsv_path.is_file():
                raise
            warnings.warn(
                "MNE is not installed; falling back to NeuroMIND TSV export.",
                RuntimeWarning,
                stacklevel=2,
            )

    if source in {"auto", "tsv"} and data_tsv_path.is_file():
        return load_eeg_tsv(
            data_tsv_path,
            metadata,
            subject_id,
            session_id,
            task,
            run=run,
            channel_positions=channel_positions,
        )

    searched = [str(data_tsv_path)]
    if metadata_path is not None:
        searched.append(str(metadata_path))
    raise FileNotFoundError(f"No BIDS EEG data found for {prefix}; searched: {searched}")


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


def load_eeg_tsv(
    data_path: str | Path,
    metadata: dict[str, Any] | None,
    subject_id: str,
    session_id: str,
    task: str,
    *,
    run: int = 1,
    channel_positions: dict[str, tuple[float, float]] | None = None,
) -> EEGRecording:
    """Load NeuroMIND TSV EEG exports with rows as channels."""
    path = Path(data_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"EEG TSV not found: {path}")

    rows: list[list[str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or not header or header[0] != "Channel":
            raise ValueError(f"Invalid NeuroMIND TSV header in {path}")
        rows = [row for row in reader if row]

    ch_names = [row[0] for row in rows]
    data = np.asarray([[float(value) for value in row[1:]] for row in rows], dtype=np.float64)
    meta_raw = metadata or {}
    fs = float(meta_raw.get("fs", meta_raw.get("SamplingFrequency", 500.0)))
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
        bids_path=str(path),
    )
    times = np.arange(data.shape[1], dtype=np.float64) / fs
    return EEGRecording(meta=meta, data=data, times=times)


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


def resolve_vhdr_path(cfg: PipelineConfig, metadata_vhdr_path: str | None) -> Path | None:
    """Resolve a metadata `vhdr_path`, preferring the configured raw data root.

    NeuroMIND metadata may point to a local mirror. For mne_brain, the canonical
    BrainVision source can be configured as `paths.raw_data_root`; matching is by
    basename, and files under an `EXCLUIDOS` folder are ignored.
    """
    if not metadata_vhdr_path:
        return None

    metadata_path = Path(metadata_vhdr_path).expanduser()
    raw_root = raw_data_root_dir(cfg)
    if raw_root is not None:
        match = find_vhdr_by_name(raw_root, metadata_path.name)
        if match is not None:
            return match

    if metadata_path.is_file():
        return metadata_path.resolve()
    return metadata_path


def find_vhdr_by_name(raw_data_root: str | Path, filename: str) -> Path | None:
    """Find a `.vhdr` by filename under `raw_data_root`, skipping excluded data."""
    root = Path(raw_data_root).expanduser().resolve()
    if not root.is_dir():
        return None
    for path in root.rglob(filename):
        if any(part.upper() == "EXCLUIDOS" for part in path.parts):
            continue
        if path.is_file():
            return path.resolve()
    return None


def load_electrode_positions(tsv_path: str | Path) -> dict[str, tuple[float, float]] | None:
    """Load 2D electrode positions from a NeuroMIND electrodes TSV."""
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
    """Resolve the BIDS root, accepting paths relative to mne_brain or NeuroMIND."""
    configured = Path(str(cfg.paths["bids_root"]))
    candidates = []
    if configured.is_absolute():
        candidates.append(configured)
    else:
        candidates.extend(
            [
                cfg.root / configured,
                cfg.root.parent / configured,
                cfg.root / str(configured).replace("BIDS", "bids"),
                cfg.root.parent / str(configured).replace("BIDS", "bids"),
            ]
        )

    existing = [candidate for candidate in candidates if candidate.is_dir()]
    for candidate in existing:
        raw_dir = candidate / "raw"
        if raw_dir.is_dir() and any(raw_dir.iterdir()):
            return candidate.resolve()
    for candidate in existing:
        if (candidate / "electrodes").is_dir() and any((candidate / "electrodes").iterdir()):
            return candidate.resolve()
    if existing:
        return existing[0].resolve()

    return candidates[0].resolve()


def raw_data_root_dir(cfg: PipelineConfig) -> Path | None:
    configured = cfg.paths.get("raw_data_root")
    if not configured:
        return None

    path = Path(str(configured)).expanduser()
    candidates = [path] if path.is_absolute() else [cfg.root / path, cfg.root.parent / path]
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    return candidates[0].resolve() if candidates else None


def _read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _first_existing(*paths: Path) -> Path | None:
    for path in paths:
        if path.is_file():
            return path
    return None
