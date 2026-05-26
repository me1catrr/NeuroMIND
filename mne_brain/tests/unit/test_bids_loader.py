import json

import numpy as np

from mne_brain.bids.loader import (
    find_vhdr_by_name,
    load_eeg_bids,
    load_electrode_positions,
    read_vhdr_header,
    resolve_vhdr_path,
)
from mne_brain.common.types import PipelineConfig


def _cfg(root):
    return PipelineConfig(
        project={},
        study={},
        paths={"bids_root": "data/BIDS"},
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


def test_load_eeg_bids_from_tsv(tmp_path):
    project_root = tmp_path / "NeuroMIND"
    bids_root = project_root / "data" / "BIDS"
    raw_dir = bids_root / "raw"
    electrodes_dir = bids_root / "electrodes"
    raw_dir.mkdir(parents=True)
    electrodes_dir.mkdir(parents=True)

    prefix = "sub-M05_ses-T2_task-eyesclosed_run-01"
    (raw_dir / f"{prefix}_metadata.json").write_text(
        json.dumps({"fs": 500.0, "channel_names": ["Fz", "Cz"]}),
        encoding="utf-8",
    )
    (raw_dir / f"{prefix}_eeg_data.tsv").write_text(
        "Channel\tT1\tT2\tT3\nFz\t1.0\t2.0\t3.0\nCz\t-1.0\t0.0\t1.0\n",
        encoding="utf-8",
    )
    (electrodes_dir / "sub-M05_ses-T2_electrodes.tsv").write_text(
        "name\tx\ty\tz\ttype\nFz\t0.0\t0.7\t0.7\tEEG\nCz\t0.0\t0.0\t1.0\tEEG\n",
        encoding="utf-8",
    )

    cfg = _cfg(project_root / "mne_brain")
    rec = load_eeg_bids(cfg, "M05", "T2", "EC")

    assert rec.meta.subject_id == "M05"
    assert rec.meta.session_id == "T2"
    assert rec.meta.condition == "EC"
    assert rec.meta.fs == 500.0
    assert rec.meta.channel_names == ["Fz", "Cz"]
    assert rec.meta.channel_positions == {"FZ": (0.0, 0.7), "CZ": (0.0, 0.0)}
    np.testing.assert_allclose(rec.data, [[1.0, 2.0, 3.0], [-1.0, 0.0, 1.0]])
    np.testing.assert_allclose(rec.times, [0.0, 0.002, 0.004])


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


def test_find_vhdr_by_name_skips_excluded(tmp_path):
    root = tmp_path / "raw"
    excluded = root / "EXCLUIDOS"
    root.mkdir()
    excluded.mkdir()
    (excluded / "sample.vhdr").write_text("excluded", encoding="utf-8")
    expected = root / "sample.vhdr"
    expected.write_text("active", encoding="utf-8")

    assert find_vhdr_by_name(root, "sample.vhdr") == expected.resolve()


def test_resolve_vhdr_path_prefers_raw_data_root(tmp_path):
    project_root = tmp_path / "NeuroMIND"
    raw_root = tmp_path / "raw_source"
    raw_root.mkdir()
    preferred = raw_root / "sample.vhdr"
    preferred.write_text("preferred", encoding="utf-8")
    old = project_root / "old" / "sample.vhdr"
    old.parent.mkdir(parents=True)
    old.write_text("old", encoding="utf-8")

    cfg = _cfg(project_root / "mne_brain")
    cfg.paths["raw_data_root"] = str(raw_root)

    assert resolve_vhdr_path(cfg, str(old)) == preferred.resolve()
