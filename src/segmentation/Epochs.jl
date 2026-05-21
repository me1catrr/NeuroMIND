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
