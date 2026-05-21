# NeuroMIND/src/segmentation/Epochs.jl
# Segmentación de señal EEG continua en epochs + baseline + artifact rejection.

"""
    segment_recording(rec::EEGRecording, cfg::PipelineConfig) -> EpochSet

Segmenta la señal en epochs de longitud fija (sin solapamiento por defecto).
"""
function segment_recording(rec::EEGRecording, cfg::PipelineConfig)::EpochSet
    seg = cfg.segmentation
    epoch_s = get(seg, "epoch_length_s", 1.0)
    overlap  = get(seg, "epoch_overlap", 0.0)

    epoch_samp = round(Int, epoch_s * rec.meta.fs)
    step_samp  = round(Int, epoch_samp * (1.0 - overlap))

    n_ch, n_total = size(rec.data)
    starts = collect(1:step_samp:(n_total - epoch_samp + 1))
    n_epochs = length(starts)

    data = Array{Float64,3}(undef, n_ch, epoch_samp, n_epochs)
    for (k, s) in enumerate(starts)
        data[:, :, k] = rec.data[:, s:(s + epoch_samp - 1)]
    end

    return EpochSet(rec.meta, data, epoch_s, n_epochs, Int[])
end

"""
    apply_baseline(epochs::EpochSet, cfg::PipelineConfig) -> EpochSet

Corrección de baseline: resta la media de cada epoch por canal.
"""
function apply_baseline(epochs::EpochSet, cfg::PipelineConfig)::EpochSet
    get(cfg.baseline, "apply", true) || return epochs

    method = get(cfg.baseline, "method", "mean")
    data   = copy(epochs.data)

    n_ch, n_samp, n_ep = size(data)
    @inbounds for ep in 1:n_ep, ch in 1:n_ch
        seg = @view data[ch, :, ep]
        bl  = method == "median" ? median(seg) : mean(seg)
        data[ch, :, ep] .-= bl
    end

    return EpochSet(epochs.meta, data, epochs.epoch_length_s, epochs.n_valid, epochs.rejected_idx)
end

"""
    reject_artifacts(epochs::EpochSet, cfg::PipelineConfig) -> EpochSet

Rechaza epochs que superan el umbral de amplitud o gradiente.
Devuelve un `EpochSet` con solo los epochs válidos.
"""
function reject_artifacts(epochs::EpochSet, cfg::PipelineConfig)::EpochSet
    ar = cfg.artifact_rejection
    get(ar, "enabled", true) || return epochs

    amp_thresh  = get(ar, "amplitude_threshold_uv", 100.0)
    grad_thresh = get(ar, "gradient_threshold_uv", 50.0)

    n_ch, n_samp, n_ep = size(epochs.data)
    rejected = Int[]
    kept     = Int[]

    @inbounds for ep in 1:n_ep
        bad = false
        for ch in 1:n_ch
            seg = @view epochs.data[ch, :, ep]
            if maximum(abs.(seg)) > amp_thresh
                bad = true; break
            end
            if maximum(abs.(diff(seg))) > grad_thresh
                bad = true; break
            end
        end
        bad ? push!(rejected, ep) : push!(kept, ep)
    end

    n_valid  = length(kept)
    min_req  = get(cfg.segmentation, "min_epochs", 20)
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

quality ∈ [0, 1] = clamp(1 − 0.5 × max(amp/thresh, grad/thresh), 0, 1)
worst_channel  = canal con mayor amplitud máxima absoluta en el epoch.
p2p_uv         = pico a pico global (max − min) a través de todos los canales.
"""
function compute_epoch_quality_report(epochs::EpochSet, cfg::PipelineConfig)::DataFrame
    ar          = cfg.artifact_rejection
    amp_thresh  = Float64(get(ar, "amplitude_threshold_uv", 100.0))
    grad_thresh = Float64(get(ar, "gradient_threshold_uv",  50.0))
    epoch_s     = epochs.epoch_length_s
    overlap     = Float64(get(cfg.segmentation, "epoch_overlap", 0.0))
    step_s      = epoch_s * (1.0 - overlap)
    n_ch, n_samp, n_ep = size(epochs.data)
    ch_names = epochs.meta.channel_names

    epoch_v    = Int[];    start_v   = Float64[]; end_v    = Float64[]
    quality_v  = Float64[]; status_v  = String[];  reason_v = String[]
    amp_v      = Float64[]; grad_v    = Float64[]
    worst_ch_v = String[];  p2p_v     = Float64[]

    @inbounds for ep in 1:n_ep
        ma = 0.0; mg = 0.0; reason = ""; bad = false
        worst_ch_idx = 1
        ep_min = Inf; ep_max = -Inf
        for ch in 1:n_ch
            seg = @view epochs.data[ch, :, ep]
            v   = maximum(abs.(seg))
            if v > ma; ma = v; worst_ch_idx = ch; end
            if v > amp_thresh && !bad
                reason = "amplitude"; bad = true
            end
            ch_min = minimum(seg); ch_max = maximum(seg)
            ch_min < ep_min && (ep_min = ch_min)
            ch_max > ep_max && (ep_max = ch_max)
            if n_samp > 1
                dv = maximum(abs.(diff(seg)))
                mg = dv > mg ? dv : mg
                if dv > grad_thresh && !bad
                    reason = "gradient"; bad = true
                end
            end
        end
        q  = clamp(1.0 - 0.5 * max(ma / amp_thresh, mg / (grad_thresh > 0 ? grad_thresh : 1.0)), 0.0, 1.0)
        t0 = round((ep - 1) * step_s, digits=3)
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
"""
function compute_channel_coverage(epochs::EpochSet, cfg::PipelineConfig)::DataFrame
    ar          = cfg.artifact_rejection
    amp_thresh  = Float64(get(ar, "amplitude_threshold_uv", 100.0))
    grad_thresh = Float64(get(ar, "gradient_threshold_uv",  50.0))
    n_ch, n_samp, n_ep = size(epochs.data)
    ch_names = epochs.meta.channel_names

    pct = Float64[]
    @inbounds for ch in 1:n_ch
        ok = 0
        for ep in 1:n_ep
            seg = @view epochs.data[ch, :, ep]
            if maximum(abs.(seg)) ≤ amp_thresh &&
               (n_samp ≤ 1 || maximum(abs.(diff(seg))) ≤ grad_thresh)
                ok += 1
            end
        end
        push!(pct, round(100.0 * ok / n_ep, digits=1))
    end
    return DataFrame(channel = ch_names, coverage_pct = pct)
end
