# NeuroMIND/src/preprocessing/Filtering.jl
#
# Cadena de filtrado con soporte para dos perfiles:
#
#   profile = "eeg_julia"  (reproducibilidad estricta con EEG_Julia)
#     Orden:  Notch(filt) → Bandreject(filt) → Highpass(filtfilt) → Lowpass(filtfilt)
#     Refs:   ../EEG_Julia/src/Preprocessing/filtering.jl
#             ../EEG_Julia/config/default_config.jl
#
#   profile = "default"  (legacy NeuroMIND / protocolo BrainVision)
#     Orden:  Highpass(filtfilt) → Lowpass(filtfilt) → Notch(filtfilt) → Bandreject(filtfilt)
#
# El perfil activo se lee de cfg.filtering["profile"].
# Si la clave no existe se asume "default".

"""
    filter_recording(rec, cfg) -> EEGRecording

Aplica la cadena completa de filtrado según `cfg.filtering["profile"]`.

Perfil `"eeg_julia"` (reproducibilidad estricta):
  1. Notch      50 Hz, bw 1 Hz, order 4, método: `filt` (causal)
  2. Bandreject 100 Hz, bw 1 Hz → (99.5–100.5 Hz), order 4, método: `filt` (causal)
  3. High-pass   0.5 Hz, order 4, método: `filtfilt` (zero-phase)
  4. Low-pass  150.0 Hz, order 4, método: `filtfilt` (zero-phase)

Perfil `"default"` (legacy):
  1. High-pass  (filtfilt)
  2. Low-pass   (filtfilt)
  3. Notch      (filtfilt, si notch_hz > 0)
  4. Bandreject (filtfilt, si bandreject_lo > 0 && hi > lo)
"""
function filter_recording(rec::EEGRecording, cfg::PipelineConfig)::EEGRecording
    f       = cfg.filtering
    data    = copy(rec.data)
    fs      = rec.meta.fs
    profile = get(f, "profile", "default")
    ord     = Int(get(f, "filter_order", 4))
    hp      = Float64(get(f, "highpass_hz",    0.5))
    lp      = Float64(get(f, "lowpass_hz",    48.0))
    nz      = Float64(get(f, "notch_hz",      50.0))
    nbw     = Float64(get(f, "notch_bw_hz",    2.0))
    lo      = Float64(get(f, "bandreject_lo",  0.0))
    hi      = Float64(get(f, "bandreject_hi",  0.0))

    if profile == "eeg_julia"
        # ── Reproducibilidad EEG_Julia ─────────────────────────
        # Paso 1 — Notch causal (filt, order 4, bw 1 Hz)
        nz > 0.0 && (data = _filt(data, fs, :notch;    freq=nz, bw=nbw,       order=ord, method=:filt))
        # Paso 2 — Bandreject causal (filt, order 4, 99.5–100.5 Hz)
        lo > 0.0 && hi > lo && (data = _filt(data, fs, :bandstop; freq=(lo,hi), order=ord, method=:filt))
        # Paso 3 — Highpass zero-phase (filtfilt, order 4)
        data = _filt(data, fs, :highpass; freq=hp, order=ord, method=:filtfilt)
        # Paso 4 — Lowpass zero-phase (filtfilt, order 4)
        data = _filt(data, fs, :lowpass;  freq=lp, order=ord, method=:filtfilt)
    else
        # ── Protocolo BrainVision / legacy ─────────────────────
        data = _filt(data, fs, :highpass; freq=hp,  order=ord)
        data = _filt(data, fs, :lowpass;  freq=lp,  order=ord)
        nz > 0.0 && (data = _filt(data, fs, :notch; freq=nz, bw=nbw))
        lo > 0.0 && hi > lo &&
            (data = _filt(data, fs, :bandstop; freq=(lo, hi), order=ord))
    end

    return EEGRecording(rec.meta, data, rec.times)
end

"""
    describe_filter_chain(cfg) -> Vector{NamedTuple}

Devuelve la cadena de filtros que aplicará `filter_recording` según `cfg.filtering["profile"]`.
Cada elemento: `(step, name, freq, order, method)`.

Útil para logging, dashboard y validación.
"""
function describe_filter_chain(cfg::PipelineConfig)
    f       = cfg.filtering
    profile = get(f, "profile", "default")
    ord     = Int(get(f, "filter_order", 4))
    hp      = Float64(get(f, "highpass_hz",    0.5))
    lp      = Float64(get(f, "lowpass_hz",    48.0))
    nz      = Float64(get(f, "notch_hz",      50.0))
    nbw     = Float64(get(f, "notch_bw_hz",    2.0))
    lo      = Float64(get(f, "bandreject_lo",  0.0))
    hi      = Float64(get(f, "bandreject_hi",  0.0))

    if profile == "eeg_julia"
        chain = NamedTuple[]
        step  = 1
        nz > 0.0 && (push!(chain, (step=step, name="Notch",
            freq="$(nz - nbw/2)–$(nz + nbw/2) Hz",
            order=ord, method="filt")); step += 1)
        lo > 0.0 && hi > lo && (push!(chain, (step=step, name="Bandreject",
            freq="$(lo)–$(hi) Hz",
            order=ord, method="filt")); step += 1)
        push!(chain, (step=step, name="High-pass",
            freq="$(hp) Hz",   order=ord, method="filtfilt")); step += 1
        push!(chain, (step=step, name="Low-pass",
            freq="$(lp) Hz",   order=ord, method="filtfilt"))
        return chain
    else
        chain = NamedTuple[
            (step=1, name="High-pass", freq="$(hp) Hz", order=ord, method="filtfilt"),
            (step=2, name="Low-pass",  freq="$(lp) Hz", order=ord, method="filtfilt"),
        ]
        s = 3
        nz > 0.0 && (push!(chain, (step=s, name="Notch",
            freq="$(nz - nbw/2)–$(nz + nbw/2) Hz",
            order=ord, method="filtfilt")); s += 1)
        lo > 0.0 && hi > lo && push!(chain, (step=s, name="Bandreject",
            freq="$(lo)–$(hi) Hz",
            order=ord, method="filtfilt"))
        return chain
    end
end

# ─── Filtros individuales (API pública) ───────────────────────

function apply_highpass(rec::EEGRecording, cutoff_hz::Real; order::Int=4)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :highpass; freq=cutoff_hz, order=order), rec.times)
end

function apply_lowpass(rec::EEGRecording, cutoff_hz::Real; order::Int=4)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :lowpass; freq=cutoff_hz, order=order), rec.times)
end

function apply_notch(rec::EEGRecording, center_hz::Real, bw_hz::Real=2.0)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :notch; freq=center_hz, bw=bw_hz), rec.times)
end

function apply_bandpass(rec::EEGRecording, lo::Real, hi::Real; order::Int=8)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :bandpass; freq=(lo, hi), order=order), rec.times)
end

function apply_bandreject(rec::EEGRecording, lo::Real, hi::Real; order::Int=4)::EEGRecording
    EEGRecording(rec.meta, _filt(rec.data, rec.meta.fs, :bandstop; freq=(lo, hi), order=order), rec.times)
end

# ─── Kernel interno ───────────────────────────────────────────

"""
    _filt(data, fs, type; freq, bw, order, method) -> Matrix{Float64}

Kernel de filtrado sobre todas las filas (canales) de una matriz.
`method = :filtfilt` (zero-phase, por defecto) | `:filt` (causal, compatible EEG_Julia).
"""
function _filt(data::Matrix{Float64}, fs::Real, type::Symbol;
               freq=nothing, bw::Real=2.0, order::Int=4,
               method::Symbol=:filtfilt)::Matrix{Float64}
    flt = _build(type, fs, freq, bw, order)
    out = similar(data)
    @inbounds for ch in 1:size(data, 1)
        out[ch, :] = method === :filt ? filt(flt, data[ch, :]) : filtfilt(flt, data[ch, :])
    end
    return out
end

function _build(type::Symbol, fs::Real, freq, bw::Real, order::Int)
    nyq = fs / 2.0
    if type === :highpass
        return digitalfilter(Highpass(Float64(freq) / nyq), Butterworth(order))
    elseif type === :lowpass
        return digitalfilter(Lowpass(Float64(freq) / nyq), Butterworth(order))
    elseif type === :bandpass
        return digitalfilter(Bandpass(Float64(freq[1]) / nyq, Float64(freq[2]) / nyq), Butterworth(order))
    elseif type === :notch || type === :bandstop
        lo = isa(freq, Tuple) ? Float64(freq[1]) / nyq : (Float64(freq) - bw / 2) / nyq
        hi = isa(freq, Tuple) ? Float64(freq[2]) / nyq : (Float64(freq) + bw / 2) / nyq
        return digitalfilter(Bandstop(lo, hi), Butterworth(order))
    else
        error("Tipo de filtro desconocido: $type")
    end
end
