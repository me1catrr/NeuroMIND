# NeuroMIND/src/segmentation/Epochs.jl
# Segmentación de señal EEG continua en epochs + baseline + artifact rejection.

"""
    segment_recording(rec::EEGRecording, cfg::PipelineConfig) -> EpochSet

Segmenta la señal en epochs de longitud fija (sin solapamiento por defecto).

- `profile = "eeg_julia"` → forza `epoch_length_s = 1.0`, `overlap = 0.0`
  (ignora `epoch_length_s` / `epoch_overlap` del config; replica EEG_Julia exacto)
- `profile = "default"`  → usa los valores configurados
"""
function segment_recording(rec::EEGRecording, cfg::PipelineConfig)::EpochSet
    seg     = cfg.segmentation
    profile = String(get(seg, "profile", "default"))

    if profile == "eeg_julia"
        epoch_s = 1.0
        overlap = 0.0
    else
        epoch_s = Float64(get(seg, "epoch_length_s", 1.0))
        overlap = Float64(get(seg, "epoch_overlap",  0.0))
    end

    epoch_samp = round(Int, epoch_s * rec.meta.fs)
    step_samp  = round(Int, epoch_samp * (1.0 - overlap))

    n_ch, n_total = size(rec.data)
    starts  = collect(1:step_samp:(n_total - epoch_samp + 1))
    n_ep    = length(starts)

    data = Array{Float64,3}(undef, n_ch, epoch_samp, n_ep)
    for (k, s) in enumerate(starts)
        data[:, :, k] = rec.data[:, s:(s + epoch_samp - 1)]
    end

    return EpochSet(rec.meta, data, epoch_s, n_ep, Int[])
end

"""
    apply_baseline(epochs::EpochSet, cfg::PipelineConfig) -> EpochSet

Corrección de baseline por canal y epoch.

| `method`              | Comportamiento |
|-----------------------|----------------|
| `"mean"`              | Resta la media de todo el epoch (NeuroMIND default) |
| `"median"`            | Resta la mediana de todo el epoch |
| `"first_window_mean"` | Resta la media de las muestras `[1 .. round(baseline_end_s × fs)]` (EEG_Julia: 0.00–0.10 s → muestras 1–50 a 500 Hz) |

Parámetro adicional para `"first_window_mean"`:
- `cfg.baseline["baseline_end_s"]` (default 0.10 s)
"""
function apply_baseline(epochs::EpochSet, cfg::PipelineConfig)::EpochSet
    get(cfg.baseline, "apply", true) || return epochs

    method = String(get(cfg.baseline, "method", "mean"))
    data   = copy(epochs.data)
    fs     = epochs.meta.fs
    n_ch, n_samp, n_ep = size(data)

    if method == "first_window_mean"
        bl_end_s   = Float64(get(cfg.baseline, "baseline_end_s", 0.10))
        bl_end_idx = min(max(round(Int, bl_end_s * fs), 1), n_samp)
        @inbounds for ep in 1:n_ep, ch in 1:n_ch
            bl = mean(@view data[ch, 1:bl_end_idx, ep])
            data[ch, :, ep] .-= bl
        end
    else
        @inbounds for ep in 1:n_ep, ch in 1:n_ch
            seg = @view data[ch, :, ep]
            bl  = method == "median" ? median(seg) : mean(seg)
            data[ch, :, ep] .-= bl
        end
    end

    return EpochSet(epochs.meta, data, epochs.epoch_length_s, epochs.n_valid, epochs.rejected_idx)
end

"""
    reject_artifacts(epochs::EpochSet, cfg::PipelineConfig) -> EpochSet

Rechaza epochs que superan umbrales de amplitud (y opcionalmente gradiente).
Devuelve un `EpochSet` con solo los epochs válidos.

| `profile`      | Criterio |
|----------------|----------|
| `"eeg_julia"`  | ±70 µV, primeros `n_channels_used=30` canales, sin gradiente |
| `"default"`    | `amplitude_threshold_uv=100` µV + gradiente `gradient_threshold_uv=50` µV/muestra, todos los canales |
"""
function reject_artifacts(epochs::EpochSet, cfg::PipelineConfig)::EpochSet
    ar = cfg.artifact_rejection
    get(ar, "enabled", true) || return epochs

    profile = String(get(ar, "profile", "default"))

    local min_uv::Float64, max_uv::Float64, n_ch_used::Int, use_grad::Bool, grad_thresh::Float64
    if profile == "eeg_julia"
        min_uv      = Float64(get(ar, "min_amplitude_uv", -70.0))
        max_uv      = Float64(get(ar, "max_amplitude_uv",  70.0))
        n_ch_used   = Int(get(ar, "n_channels_used", 30))
        use_grad    = false
        grad_thresh = 0.0
    else
        amp_thr     = Float64(get(ar, "amplitude_threshold_uv", 100.0))
        min_uv      = -amp_thr
        max_uv      =  amp_thr
        n_ch_used   = size(epochs.data, 1)   # todos los canales
        use_grad    = Bool(get(ar, "use_gradient", true))
        grad_thresh = Float64(get(ar, "gradient_threshold_uv", 50.0))
    end

    n_ch, n_samp, n_ep = size(epochs.data)
    ch_max   = min(n_ch_used, n_ch)
    rejected = Int[]
    kept     = Int[]

    @inbounds for ep in 1:n_ep
        bad = false
        for ch in 1:ch_max
            seg = @view epochs.data[ch, :, ep]
            if minimum(seg) < min_uv || maximum(seg) > max_uv
                bad = true; break
            end
            if use_grad && n_samp > 1 && maximum(abs.(diff(seg))) > grad_thresh
                bad = true; break
            end
        end
        bad ? push!(rejected, ep) : push!(kept, ep)
    end

    n_valid = length(kept)
    min_req = get(cfg.segmentation, "min_epochs", 20)
    n_valid < min_req &&
        @warn "$(epochs.meta.subject_id)/$(epochs.meta.session_id): solo $n_valid epochs válidos (mínimo $min_req)"

    clean = epochs.data[:, :, kept]
    return EpochSet(epochs.meta, clean, epochs.epoch_length_s, n_valid, rejected)
end

"""
    compute_epoch_quality_report(epochs, cfg) -> DataFrame

Calcula métricas de calidad para TODOS los epochs (pre-AR).
Columnas: epoch, start_s, end_s, duration_s, quality, status,
          rejection_reason, max_amp_uv, max_grad_uv, worst_channel, p2p_uv.

Respeta el perfil AR activo (`cfg.artifact_rejection["profile"]`):
- `"eeg_julia"` → umbral ±70 µV, primeros 30 canales, sin gradiente
- `"default"`   → umbral 100 µV, todos los canales, gradiente incluido

quality ∈ [0, 1] = clamp(1 − 0.5 × max(amp_ratio, grad_ratio), 0, 1)
worst_channel  = canal con mayor amplitud absoluta en el epoch.
p2p_uv         = pico a pico global (max − min) a través de todos los canales.
"""
function compute_epoch_quality_report(epochs::EpochSet, cfg::PipelineConfig)::DataFrame
    ar      = cfg.artifact_rejection
    profile = String(get(ar, "profile", "default"))
    epoch_s = epochs.epoch_length_s
    overlap = Float64(get(cfg.segmentation, "epoch_overlap", 0.0))
    step_s  = epoch_s * (1.0 - overlap)
    n_ch, n_samp, n_ep = size(epochs.data)
    ch_names = epochs.meta.channel_names

    local amp_thresh::Float64, grad_thresh::Float64, n_ch_used::Int,
          use_grad::Bool, min_uv::Float64, max_uv::Float64
    if profile == "eeg_julia"
        min_uv      = Float64(get(ar, "min_amplitude_uv", -70.0))
        max_uv      = Float64(get(ar, "max_amplitude_uv",  70.0))
        amp_thresh  = max_uv   # para ratio de calidad (umbral simétrico)
        n_ch_used   = Int(get(ar, "n_channels_used", 30))
        use_grad    = false
        grad_thresh = 1.0      # dummy, no usado
    else
        amp_thresh  = Float64(get(ar, "amplitude_threshold_uv", 100.0))
        grad_thresh = Float64(get(ar, "gradient_threshold_uv",  50.0))
        min_uv      = -amp_thresh
        max_uv      =  amp_thresh
        n_ch_used   = n_ch
        use_grad    = Bool(get(ar, "use_gradient", true))
    end

    ch_max = min(n_ch_used, n_ch)

    epoch_v    = Int[];    start_v   = Float64[]; end_v    = Float64[]
    quality_v  = Float64[]; status_v  = String[];  reason_v = String[]
    amp_v      = Float64[]; grad_v    = Float64[]
    worst_ch_v = String[];  p2p_v     = Float64[]

    @inbounds for ep in 1:n_ep
        ma = 0.0; mg = 0.0; reason = ""; bad = false
        worst_ch_idx = 1
        ep_min = Inf; ep_max = -Inf

        for ch in 1:ch_max
            seg = @view epochs.data[ch, :, ep]
            seg_min = minimum(seg); seg_max = maximum(seg)
            v = max(abs(seg_min), abs(seg_max))
            if v > ma; ma = v; worst_ch_idx = ch; end
            if (seg_min < min_uv || seg_max > max_uv) && !bad
                reason = "amplitude"; bad = true
            end
            seg_min < ep_min && (ep_min = seg_min)
            seg_max > ep_max && (ep_max = seg_max)
            if use_grad && n_samp > 1
                dv = maximum(abs.(diff(seg)))
                mg = dv > mg ? dv : mg
                if dv > grad_thresh && !bad
                    reason = "gradient"; bad = true
                end
            end
        end

        grad_ratio = use_grad ? mg / grad_thresh : 0.0
        q   = clamp(1.0 - 0.5 * max(ma / amp_thresh, grad_ratio), 0.0, 1.0)
        t0  = round((ep - 1) * step_s, digits=3)
        p2p = isinf(ep_min) ? 0.0 : round(ep_max - ep_min, digits=2)
        worst_name = (worst_ch_idx <= length(ch_names)) ?
                     ch_names[worst_ch_idx] : "CH$(worst_ch_idx)"

        push!(epoch_v,    ep)
        push!(start_v,    t0)
        push!(end_v,      round(t0 + epoch_s, digits=3))
        push!(quality_v,  round(q, digits=4))
        push!(status_v,   bad ? "rejected" : "valid")
        push!(reason_v,   reason)
        push!(amp_v,      round(ma, digits=2))
        push!(grad_v,     round(mg, digits=2))
        push!(worst_ch_v, worst_name)
        push!(p2p_v,      p2p)
    end

    return DataFrame(
        epoch            = epoch_v,
        start_s          = start_v,
        end_s            = end_v,
        duration_s       = fill(epoch_s, n_ep),
        quality          = quality_v,
        status           = status_v,
        rejection_reason = reason_v,
        max_amp_uv       = amp_v,
        max_grad_uv      = grad_v,
        worst_channel    = worst_ch_v,
        p2p_uv           = p2p_v,
    )
end

"""
    compute_channel_coverage(epochs, cfg) -> DataFrame

Proporción de épocas (%) en que cada canal pasa los umbrales de AR.
Respeta el perfil AR activo (eeg_julia / default).
"""
function compute_channel_coverage(epochs::EpochSet, cfg::PipelineConfig)::DataFrame
    ar      = cfg.artifact_rejection
    profile = String(get(ar, "profile", "default"))
    n_ch, n_samp, n_ep = size(epochs.data)
    ch_names = epochs.meta.channel_names

    local min_uv::Float64, max_uv::Float64, use_grad::Bool, grad_thresh::Float64
    if profile == "eeg_julia"
        min_uv      = Float64(get(ar, "min_amplitude_uv", -70.0))
        max_uv      = Float64(get(ar, "max_amplitude_uv",  70.0))
        use_grad    = false
        grad_thresh = 0.0
    else
        amp_thr     = Float64(get(ar, "amplitude_threshold_uv", 100.0))
        min_uv      = -amp_thr
        max_uv      =  amp_thr
        use_grad    = Bool(get(ar, "use_gradient", true))
        grad_thresh = Float64(get(ar, "gradient_threshold_uv", 50.0))
    end

    pct = Float64[]
    @inbounds for ch in 1:n_ch
        ok = 0
        for ep in 1:n_ep
            seg = @view epochs.data[ch, :, ep]
            seg_min = minimum(seg); seg_max = maximum(seg)
            pass = seg_min >= min_uv && seg_max <= max_uv
            if pass && use_grad && n_samp > 1
                pass = maximum(abs.(diff(seg))) ≤ grad_thresh
            end
            pass && (ok += 1)
        end
        push!(pct, round(100.0 * ok / n_ep, digits=1))
    end
    return DataFrame(channel = ch_names, coverage_pct = pct)
end
