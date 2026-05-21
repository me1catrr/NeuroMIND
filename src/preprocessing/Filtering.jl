# NeuroMIND/src/preprocessing/Filtering.jl
# Cadena de filtrado conforme al protocolo del estudio BRAIN (BRAIN.pdf):
# HP(0.5Hz) → LP(48Hz) → Notch(50Hz) → Band-reject(100-120Hz)

"""
    filter_recording(rec, cfg) -> EEGRecording

Aplica la cadena completa de filtrado del protocolo:
  1. Paso alto  (highpass_hz)
  2. Paso bajo  (lowpass_hz)
  3. Notch      (notch_hz, notch_bw_hz)
  4. Band-reject (bandreject_lo, bandreject_hi)  — si está definido
"""
function filter_recording(rec::EEGRecording, cfg::PipelineConfig)::EEGRecording
    f    = cfg.filtering
    data = copy(rec.data)
    fs   = rec.meta.fs

    data = _filt(data, fs, :highpass; freq=f["highpass_hz"],  order=get(f,"filter_order",4))
    data = _filt(data, fs, :lowpass;  freq=f["lowpass_hz"],   order=get(f,"filter_order",4))

    if get(f, "notch_hz", 0.0) > 0.0
        data = _filt(data, fs, :notch; freq=f["notch_hz"], bw=get(f,"notch_bw_hz",2.0))
    end

    lo = get(f, "bandreject_lo", 0.0)
    hi = get(f, "bandreject_hi", 0.0)
    if lo > 0.0 && hi > lo
        data = _filt(data, fs, :bandstop; freq=(lo, hi), order=get(f,"filter_order",4))
    end

    return EEGRecording(rec.meta, data, rec.times)
end

function apply_highpass(rec::EEGRecording, cutoff_hz::Real; order::Int=4)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :highpass; freq=cutoff_hz, order), rec.times)
end

function apply_lowpass(rec::EEGRecording, cutoff_hz::Real; order::Int=4)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :lowpass; freq=cutoff_hz, order), rec.times)
end

function apply_notch(rec::EEGRecording, center_hz::Real, bw_hz::Real=2.0)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :notch; freq=center_hz, bw=bw_hz), rec.times)
end

function apply_bandpass(rec::EEGRecording, lo::Real, hi::Real; order::Int=8)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :bandpass; freq=(lo,hi), order), rec.times)
end

function apply_bandreject(rec::EEGRecording, lo::Real, hi::Real; order::Int=4)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :bandstop; freq=(lo,hi), order), rec.times)
end

# ─── Kernel ───────────────────────────────────────────────────

function _filt(data::Matrix{Float64}, fs::Real, type::Symbol;
               freq=nothing, bw::Real=2.0, order::Int=4)::Matrix{Float64}
    filt = _build(type, fs, freq, bw, order)
    out  = similar(data)
    @inbounds for ch in 1:size(data,1)
        out[ch, :] = filtfilt(filt, data[ch, :])
    end
    return out
end

function _build(type::Symbol, fs::Real, freq, bw::Real, order::Int)
    nyq = fs / 2.0
    if type === :highpass
        return digitalfilter(Highpass(Float64(freq)/nyq), Butterworth(order))
    elseif type === :lowpass
        return digitalfilter(Lowpass(Float64(freq)/nyq), Butterworth(order))
    elseif type === :bandpass
        return digitalfilter(Bandpass(Float64(freq[1])/nyq, Float64(freq[2])/nyq), Butterworth(order))
    elseif type === :notch || type === :bandstop
        lo = isa(freq, Tuple) ? Float64(freq[1])/nyq : (Float64(freq) - bw/2)/nyq
        hi = isa(freq, Tuple) ? Float64(freq[2])/nyq : (Float64(freq) + bw/2)/nyq
        return digitalfilter(Bandstop(lo, hi), Butterworth(order))
    else
        error("Tipo de filtro desconocido: $type")
    end
end
