# NeuroMIND/src/qc/QualityControl.jl
# Control de calidad de señales EEG: estadísticas por canal, detección de canales malos.

"""
    compute_channel_stats(rec::EEGRecording) -> DataFrame

Calcula estadísticas básicas por canal: media, RMS, varianza, rango, z-score de RMS.
"""
function compute_channel_stats(rec::EEGRecording)::DataFrame
    n_ch = rec.meta.n_channels
    ch   = rec.meta.channel_names

    μ     = [mean(rec.data[c, :]) for c in 1:n_ch]
    rms   = [sqrt(mean(rec.data[c, :].^2)) for c in 1:n_ch]
    σ     = [std(rec.data[c, :]) for c in 1:n_ch]
    rng   = [maximum(rec.data[c, :]) - minimum(rec.data[c, :]) for c in 1:n_ch]

    rms_z  = (rms .- mean(rms)) ./ (std(rms) + eps())

    return DataFrame(
        channel   = ch,
        mean_uv   = round.(μ,   digits=3),
        rms_uv    = round.(rms, digits=3),
        std_uv    = round.(σ,   digits=3),
        range_uv  = round.(rng, digits=3),
        rms_zscore = round.(rms_z, digits=3)
    )
end

"""
    flag_bad_channels(rec::EEGRecording; z_threshold=3.0) -> Vector{String}

Identifica canales sospechosos por z-score de RMS fuera del umbral.
"""
function flag_bad_channels(rec::EEGRecording; z_threshold::Real=3.0)::Vector{String}
    stats = compute_channel_stats(rec)
    bad   = stats[abs.(stats.rms_zscore) .> z_threshold, :channel]
    return convert(Vector{String}, bad)
end

"""
    qc_report(rec::EEGRecording, cfg::PipelineConfig) -> DataFrame

Genera un informe completo de QC y lo guarda como CSV en results/.
"""
function qc_report(rec::EEGRecording, cfg::PipelineConfig)::DataFrame
    stats     = compute_channel_stats(rec)
    bad_ch    = flag_bad_channels(rec)

    stats[!, :is_bad] = [ch in bad_ch for ch in stats.channel]

    out_dir = ensure_dirs(cfg, rec.meta.subject_id, rec.meta.session_id)
    path    = joinpath(out_dir, "tables", "qc_channels_$(rec.meta.condition).csv")
    CSV.write(path, stats)

    return stats
end
