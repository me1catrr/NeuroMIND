from mne_brain.common.config import load_config
from mne_brain.common.types import PipelineConfig


def test_load_default_config():
    cfg = load_config()

    assert isinstance(cfg, PipelineConfig)
    assert cfg.recording["fs"] == 500.0
    assert cfg.filtering["profile"] == "eeg_julia"
    assert cfg.filtering["filter_order"] == 4
    assert cfg.spectral["nfft"] == 512
    assert cfg.surrogates["n_surrogates"] == 200
    assert cfg.bands["ALPHA"] == (7.8, 11.7)

