import numpy as np
import pytest

from mne_brain.bids.loader import (
    bids_eeg_dir,
    bids_root_dir,
    bids_vhdr_path,
    load_eeg_bids,
    load_electrode_positions,
    read_vhdr_header,
)
from mne_brain.common.types import PipelineConfig


def _cfg(root, bids_root="../data/bids"):
    return PipelineConfig(
        project={},
        study={},
        paths={"bids_root": bids_root},
        recording={"fs": 500.0},
        filtering={},
        segmentation={},
        baseline={},
        artifact_rejection={},
        ica={},
        spectral={},
        bands={},
        connectivity={},
        surrogates={},
        graph={},
        clinical={},
        longitudinal={},
        statistics={},
        export_cfg={},
        root=root,
    )


def _write_bids_recording(bids_root, subject, session, task, adc):
    """Write a minimal BIDS-EEG BrainVision trio (2 channels, 0.5 uV/bit, float32)."""
    eeg_dir = bids_eeg_dir(bids_root, subject, session)
    eeg_dir.mkdir(parents=True)
    base = f"sub-{subject}_ses-{session}_task-{task}_run-01_eeg"
    (bids_root / "dataset_description.json").write_text('{"Name": "t"}', encoding="utf-8")
    (eeg_dir / f"{base}.vhdr").write_text(
        "\n".join(
            [
                "Brain Vision Data Exchange Header File Version 1.0",
                "[Common Infos]",
                "Codepage=UTF-8",
                f"DataFile={base}.eeg",
                f"MarkerFile={base}.vmrk",
                "DataFormat=BINARY",
                "DataOrientation=MULTIPLEXED",
                "NumberOfChannels=2",
                "SamplingInterval=2000",
                "[Binary Infos]",
                "BinaryFormat=IEEE_FLOAT_32",
                "[Channel Infos]",
                "Ch1=Fz,,0.5,µV",
                "Ch2=Cz,,0.5,µV",
            ]
        ),
        encoding="utf-8",
    )
    (eeg_dir / f"{base}.vmrk").write_text(
        "\n".join(
            [
                "Brain Vision Data Exchange Marker File, Version 1.0",
                "[Common Infos]",
                "Codepage=UTF-8",
                f"DataFile={base}.eeg",
                "[Marker Infos]",
                "Mk1=New Segment,,1,1,0,20210804120048946581",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    np.asarray(adc, dtype="<f4").T.tofile(eeg_dir / f"{base}.eeg")  # MULTIPLEXED: sample-major
    return eeg_dir


def test_load_eeg_bids_reads_brainvision_in_microvolts(tmp_path):
    project_root = tmp_path / "NeuroMIND"
    bids_root = project_root / "data" / "bids"
    bids_root.mkdir(parents=True)
    eeg_dir = _write_bids_recording(
        bids_root, "M05", "T2", "eyesclosed", [[2.0, 4.0, 6.0], [-2.0, 0.0, 2.0]]
    )
    # Electrodes live next to the recordings; BIDS has no `type` column requirement.
    (eeg_dir / "sub-M05_ses-T2_electrodes.tsv").write_text(
        "name\tx\ty\tz\nFz\t0.0\t0.7\t0.7\nCz\t0.0\t0.0\t1.0\n", encoding="utf-8"
    )

    cfg = _cfg(project_root / "mne_brain")
    rec = load_eeg_bids(cfg, "M05", "T2", "EC")

    assert rec.meta.subject_id == "M05"
    assert rec.meta.session_id == "T2"
    assert rec.meta.condition == "EC"
    assert rec.meta.fs == 500.0
    assert rec.meta.channel_names == ["Fz", "Cz"]
    assert rec.meta.channel_positions == {"FZ": (0.0, 0.7), "CZ": (0.0, 0.0)}
    # resolution 0.5 uV/bit is applied: [2, 4, 6] ADC -> [1, 2, 3] uV
    np.testing.assert_allclose(rec.data, [[1.0, 2.0, 3.0], [-1.0, 0.0, 1.0]], atol=1e-6)
    np.testing.assert_allclose(rec.times, [0.0, 0.002, 0.004])


def test_load_eeg_bids_missing_recording_raises(tmp_path):
    project_root = tmp_path / "NeuroMIND"
    (project_root / "data" / "bids").mkdir(parents=True)
    (project_root / "data" / "bids" / "dataset_description.json").write_text("{}", encoding="utf-8")
    with pytest.raises(FileNotFoundError):
        load_eeg_bids(_cfg(project_root / "mne_brain"), "M05", "T2", "EC")


def test_bids_root_dir_resolves_sibling_dataset(tmp_path):
    project_root = tmp_path / "NeuroMIND"
    bids_root = project_root / "data" / "bids"
    bids_root.mkdir(parents=True)
    (bids_root / "dataset_description.json").write_text("{}", encoding="utf-8")

    assert bids_root_dir(_cfg(project_root / "mne_brain")) == bids_root.resolve()
    assert bids_vhdr_path(bids_root, "M05", "T2", "eyesclosed").name == (
        "sub-M05_ses-T2_task-eyesclosed_run-01_eeg.vhdr"
    )


def test_read_vhdr_header(tmp_path):
    vhdr = tmp_path / "sample.vhdr"
    vhdr.write_text(
        "\n".join(
            [
                "Brain Vision Data Exchange Header File Version 1.0",
                "[Common Infos]",
                "DataFile=sample.eeg",
                "DataOrientation=MULTIPLEXED",
                "NumberOfChannels=2",
                "SamplingInterval=2000",
                "[Binary Infos]",
                "BinaryFormat=IEEE_FLOAT_32",
                "[Channel Infos]",
                "Ch1=Fz,,0.0488281,µV",
                "Ch2=Cz,,0.5,µV",
            ]
        ),
        encoding="utf-8",
    )

    header = read_vhdr_header(vhdr)

    assert header.n_channels == 2
    assert header.fs == 500.0
    assert header.ch_names == ["Fz", "Cz"]
    assert header.resolutions == [0.0488281, 0.5]
    assert header.binary_format == "IEEE_FLOAT_32"
    assert header.orientation == "MULTIPLEXED"
    assert header.eeg_file == "sample.eeg"


def test_load_electrode_positions_missing(tmp_path):
    assert load_electrode_positions(tmp_path / "missing.tsv") is None
