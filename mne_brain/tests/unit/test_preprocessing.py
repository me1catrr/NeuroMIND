import numpy as np

from mne_brain.common.types import EEGRecording, PipelineConfig, RecordingMeta
from mne_brain.preprocessing import (
    compute_channel_stats,
    describe_filter_chain,
    filter_recording,
    flag_bad_channels,
)


def _meta(fs=500.0):
    return RecordingMeta(
        subject_id="M05",
        session_id="T2",
        condition="EC",
        run=1,
        fs=fs,
        n_channels=4,
        channel_names=["Fz", "Cz", "Pz", "Oz"],
        channel_positions=None,
        bids_path="synthetic",
    )


def _cfg():
    return PipelineConfig(
        project={},
        study={},
        paths={"results": "results"},
        recording={"fs": 500.0},
        filtering={
            "profile": "eeg_julia",
            "highpass_hz": 0.5,
            "lowpass_hz": 150.0,
            "notch_hz": 50.0,
            "notch_bw_hz": 1.0,
            "bandreject_lo": 99.5,
            "bandreject_hi": 100.5,
            "filter_order": 4,
        },
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
        root=".",
    )


def test_channel_stats_and_bad_channel_flagging():
    data = np.vstack(
        [
            np.ones(1000),
            np.ones(1000) * 1.1,
            np.ones(1000) * 0.9,
            np.ones(1000) * 20.0,
        ]
    )
    rec = EEGRecording(_meta(), data, np.arange(1000) / 500.0)

    stats = compute_channel_stats(rec)

    assert stats["channel"].to_list() == ["Fz", "Cz", "Pz", "Oz"]
    assert stats.loc[0, "mean_uv"] == 1.0
    assert stats.loc[3, "rms_uv"] == 20.0
    assert flag_bad_channels(rec, z_threshold=1.0) == ["Oz"]


def test_describe_filter_chain_eeg_julia_order():
    chain = describe_filter_chain(_cfg())

    assert [step["name"] for step in chain] == [
        "Notch",
        "Bandreject",
        "High-pass",
        "Low-pass",
    ]
    assert [step["method"] for step in chain] == ["filt", "filt", "filtfilt", "filtfilt"]


def test_filter_recording_reduces_50hz_component():
    fs = 500.0
    t = np.arange(0, 4, 1 / fs)
    x = np.sin(2 * np.pi * 10 * t) + 0.8 * np.sin(2 * np.pi * 50 * t)
    meta = RecordingMeta(
        subject_id="M05",
        session_id="T2",
        condition="EC",
        run=1,
        fs=fs,
        n_channels=2,
        channel_names=["Fz", "Cz"],
        channel_positions=None,
        bids_path="synthetic",
    )
    rec = EEGRecording(meta, np.vstack([x, x]), t)

    filtered = filter_recording(rec, _cfg())

    freqs = np.fft.rfftfreq(x.size, d=1 / fs)
    before = np.abs(np.fft.rfft(rec.data[0]))
    after = np.abs(np.fft.rfft(filtered.data[0]))
    idx_50 = np.argmin(np.abs(freqs - 50.0))
    assert after[idx_50] < before[idx_50] * 0.5
    assert filtered.data.shape == rec.data.shape
