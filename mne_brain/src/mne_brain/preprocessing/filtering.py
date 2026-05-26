"""Butterworth filtering compatible with NeuroMIND's `eeg_julia` profile."""

from __future__ import annotations

from typing import Any

import numpy as np

from mne_brain.common.types import EEGRecording, PipelineConfig


def filter_recording(rec: EEGRecording, cfg: PipelineConfig) -> EEGRecording:
    """Apply the configured NeuroMIND filter chain to an EEG recording.

    For `profile = "eeg_julia"` the exact intended order is:

    1. 50 Hz notch, causal `sosfilt`
    2. 99.5-100.5 Hz bandreject, causal `sosfilt`
    3. 0.5 Hz high-pass, zero-phase `sosfiltfilt`
    4. 150 Hz low-pass, zero-phase `sosfiltfilt`
    """
    filt_cfg = cfg.filtering
    data = np.array(rec.data, dtype=np.float64, copy=True)
    fs = float(rec.meta.fs)
    profile = str(filt_cfg.get("profile", "default"))
    order = int(filt_cfg.get("filter_order", 4))
    hp = float(filt_cfg.get("highpass_hz", 0.5))
    lp = float(filt_cfg.get("lowpass_hz", 150.0))
    notch = float(filt_cfg.get("notch_hz", 50.0))
    notch_bw = float(filt_cfg.get("notch_bw_hz", 1.0))
    br_lo = float(filt_cfg.get("bandreject_lo", 0.0))
    br_hi = float(filt_cfg.get("bandreject_hi", 0.0))

    if profile == "eeg_julia":
        if notch > 0.0:
            data = _apply_sos(
                data,
                fs,
                "bandstop",
                (notch - notch_bw / 2, notch + notch_bw / 2),
                order,
                method="filt",
            )
        if br_lo > 0.0 and br_hi > br_lo:
            data = _apply_sos(data, fs, "bandstop", (br_lo, br_hi), order, method="filt")
        data = _apply_sos(data, fs, "highpass", hp, order, method="filtfilt")
        data = _apply_sos(data, fs, "lowpass", lp, order, method="filtfilt")
    else:
        data = _apply_sos(data, fs, "highpass", hp, order, method="filtfilt")
        data = _apply_sos(data, fs, "lowpass", lp, order, method="filtfilt")
        if notch > 0.0:
            data = _apply_sos(
                data,
                fs,
                "bandstop",
                (notch - notch_bw / 2, notch + notch_bw / 2),
                order,
                method="filtfilt",
            )
        if br_lo > 0.0 and br_hi > br_lo:
            data = _apply_sos(data, fs, "bandstop", (br_lo, br_hi), order, method="filtfilt")

    return EEGRecording(meta=rec.meta, data=data, times=rec.times.copy())


def describe_filter_chain(cfg: PipelineConfig) -> list[dict[str, Any]]:
    """Return the filter chain as dashboard-friendly dictionaries."""
    filt_cfg = cfg.filtering
    profile = str(filt_cfg.get("profile", "default"))
    order = int(filt_cfg.get("filter_order", 4))
    hp = float(filt_cfg.get("highpass_hz", 0.5))
    lp = float(filt_cfg.get("lowpass_hz", 150.0))
    notch = float(filt_cfg.get("notch_hz", 50.0))
    notch_bw = float(filt_cfg.get("notch_bw_hz", 1.0))
    br_lo = float(filt_cfg.get("bandreject_lo", 0.0))
    br_hi = float(filt_cfg.get("bandreject_hi", 0.0))

    chain: list[dict[str, Any]] = []
    if profile == "eeg_julia":
        if notch > 0.0:
            chain.append(
                _step(
                    len(chain) + 1,
                    "Notch",
                    f"{notch - notch_bw / 2}-{notch + notch_bw / 2} Hz",
                    order,
                    "filt",
                )
            )
        if br_lo > 0.0 and br_hi > br_lo:
            chain.append(_step(len(chain) + 1, "Bandreject", f"{br_lo}-{br_hi} Hz", order, "filt"))
        chain.append(_step(len(chain) + 1, "High-pass", f"{hp} Hz", order, "filtfilt"))
        chain.append(_step(len(chain) + 1, "Low-pass", f"{lp} Hz", order, "filtfilt"))
        return chain

    chain.append(_step(1, "High-pass", f"{hp} Hz", order, "filtfilt"))
    chain.append(_step(2, "Low-pass", f"{lp} Hz", order, "filtfilt"))
    if notch > 0.0:
        chain.append(
            _step(
                len(chain) + 1,
                "Notch",
                f"{notch - notch_bw / 2}-{notch + notch_bw / 2} Hz",
                order,
                "filtfilt",
            )
        )
    if br_lo > 0.0 and br_hi > br_lo:
        chain.append(_step(len(chain) + 1, "Bandreject", f"{br_lo}-{br_hi} Hz", order, "filtfilt"))
    return chain


def apply_highpass(rec: EEGRecording, cutoff_hz: float, *, order: int = 4) -> EEGRecording:
    data = _apply_sos(rec.data, rec.meta.fs, "highpass", cutoff_hz, order, method="filtfilt")
    return EEGRecording(rec.meta, data, rec.times.copy())


def apply_lowpass(rec: EEGRecording, cutoff_hz: float, *, order: int = 4) -> EEGRecording:
    data = _apply_sos(rec.data, rec.meta.fs, "lowpass", cutoff_hz, order, method="filtfilt")
    return EEGRecording(rec.meta, data, rec.times.copy())


def apply_notch(
    rec: EEGRecording,
    center_hz: float,
    bw_hz: float = 2.0,
    *,
    order: int = 4,
) -> EEGRecording:
    band = (center_hz - bw_hz / 2, center_hz + bw_hz / 2)
    data = _apply_sos(rec.data, rec.meta.fs, "bandstop", band, order, method="filtfilt")
    return EEGRecording(rec.meta, data, rec.times.copy())


def apply_bandreject(rec: EEGRecording, lo: float, hi: float, *, order: int = 4) -> EEGRecording:
    data = _apply_sos(rec.data, rec.meta.fs, "bandstop", (lo, hi), order, method="filtfilt")
    return EEGRecording(rec.meta, data, rec.times.copy())


def _apply_sos(
    data: np.ndarray,
    fs: float,
    kind: str,
    freq: float | tuple[float, float],
    order: int,
    *,
    method: str,
) -> np.ndarray:
    from scipy.signal import butter, sosfilt, sosfiltfilt

    btype = {
        "highpass": "highpass",
        "lowpass": "lowpass",
        "bandpass": "bandpass",
        "bandstop": "bandstop",
    }[kind]
    sos = butter(order, freq, btype=btype, fs=fs, output="sos")
    if method == "filt":
        return sosfilt(sos, data, axis=1)
    if method == "filtfilt":
        return sosfiltfilt(sos, data, axis=1)
    raise ValueError(f"Unknown filter method: {method}")


def _step(step: int, name: str, freq: str, order: int, method: str) -> dict[str, Any]:
    return {
        "step": step,
        "name": name,
        "type": "Butterworth",
        "freq": freq,
        "order": order,
        "method": method,
        "applied": True,
    }
