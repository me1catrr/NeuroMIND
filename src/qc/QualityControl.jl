# NeuroMIND/src/qc/QualityControl.jl
# Control de calidad de señales EEG — Fase [2/8], sobre señal raw (pre-filtrado).
#
# Estadísticas temporales por canal (P1), métricas espectrales HFNoise/SNR (P2)
# y resumen de correlación intercanal. Definiciones alineadas con
# Report_Pre/chapters/03_raw.tex y feat/qc-extended-channel-stats.

# ─── P1: estadísticas básicas por canal (dominio temporal) ────────

"""
    compute_channel_stats(rec::EEGRecording) -> DataFrame

Estadísticas por canal (señal raw, dominio temporal):

- `mean_uv`, `rms_uv`, `std_uv`, `range_uv`, `rms_zscore` (legacy)
- `min_uv`, `max_uv`
- `skewness` — asimetría γ₁ = (1/N) Σ ((x−μ)/σ)³  (divisor N; NaN si σ≈0)
- `kurtosis` — curtosis **bruta** γ₂ = (1/N) Σ ((x−μ)/σ)⁴  (no en exceso; NaN si σ≈0)

Coincide con Report_Pre §Asimetría / §Curtosis.
"""
function compute_channel_stats(rec::EEGRecording)::DataFrame
    n_ch = rec.meta.n_channels
    ch   = rec.meta.channel_names

    μ     = [mean(rec.data[c, :]) for c in 1:n_ch]
    rms   = [sqrt(mean(rec.data[c, :].^2)) for c in 1:n_ch]
    σ     = [std(rec.data[c, :]) for c in 1:n_ch]
    x_min = [minimum(rec.data[c, :]) for c in 1:n_ch]
    x_max = [maximum(rec.data[c, :]) for c in 1:n_ch]
    rng   = x_max .- x_min

    rms_z = (rms .- mean(rms)) ./ (std(rms) + eps())

    skew_v = Vector{Float64}(undef, n_ch)
    kurt_v = Vector{Float64}(undef, n_ch)
    for c in 1:n_ch
        if σ[c] < 1e-9
            @warn "compute_channel_stats: canal plano (σ≈0); asimetría/curtosis indefinidas" channel=ch[c]
            skew_v[c] = NaN
            kurt_v[c] = NaN
        else
            sig = rec.data[c, :]
            # StatsBase.skewness usa divisor N (población), como el informe.
            skew_v[c] = skewness(sig)
            # StatsBase.kurtosis es EN EXCESO (γ₂−3); el informe usa γ₂ bruta.
            kurt_v[c] = kurtosis(sig) + 3.0
        end
    end

    return DataFrame(
        channel    = ch,
        mean_uv    = round.(μ,      digits=3),
        rms_uv     = round.(rms,    digits=3),
        std_uv     = round.(σ,      digits=3),
        range_uv   = round.(rng,    digits=3),
        rms_zscore = round.(rms_z,  digits=3),
        min_uv     = round.(x_min,  digits=3),
        max_uv     = round.(x_max,  digits=3),
        skewness   = round.(skew_v, digits=3),
        kurtosis   = round.(kurt_v, digits=3),
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

# ─── P2: métricas espectrales (HFNoise, SNR) y correlación intercanal ──

"""
    welch_psd_raw(sig, fs, nfft) -> (freqs, psd)

PSD Welch (Hann, segmentos sin solape), unilateral, en unidades²/Hz.
Usada por HFNoise y SNR. Devuelve vectores vacíos si `sig` es más corta que `nfft`.
"""
function welch_psd_raw(sig::AbstractVector{<:Real}, fs::Float64, nfft::Int)
    nfft > length(sig) && return Float64[], Float64[]
    win     = Float64[0.5 - 0.5 * cos(2π * (i - 1) / (nfft - 1)) for i in 1:nfft]
    win_pow = sum(win .^ 2)
    nhalf   = div(nfft, 2) + 1
    ps      = zeros(nhalf)
    n_seg   = div(length(sig), nfft)
    n_seg == 0 && return Float64[], Float64[]
    for i in 1:n_seg
        seg = Float64.(sig[(i - 1) * nfft + 1 : i * nfft]) .* win
        S   = abs.(rfft(seg)) .^ 2
        ps .+= S
    end
    ps ./= (n_seg * win_pow * fs)
    ps[2:end-1] .*= 2
    freqs = Float64[(k - 1) * fs / nfft for k in 1:nhalf]
    return freqs, ps
end

"""
    compute_channel_spectral_qc(rec, cfg; ...) -> DataFrame

Métricas espectrales de QC por canal (señal **raw**, un Welch por canal):

- `hfnoise` — ``P_{25–45 Hz} / P_{1–45 Hz}`` (Report_Pre §HFNoise)
- `snr_db`  — ``10·log10(P_signal / P_noise)``; bandas por nombre en `cfg.bands`
  (default ALPHA vs GAMMA)
"""
function compute_channel_spectral_qc(
    rec::EEGRecording,
    cfg::PipelineConfig;
    welch_nfft::Int = 1024,
    hfnoise_band::Tuple{Float64,Float64} = (25.0, 45.0),
    hfnoise_ref_band::Tuple{Float64,Float64} = (1.0, 45.0),
    snr_signal_band::String = "ALPHA",
    snr_noise_band::String = "GAMMA",
)::DataFrame
    n_ch = rec.meta.n_channels
    ch   = rec.meta.channel_names
    fs   = Float64(rec.meta.fs)

    sig_band   = get(cfg.bands, snr_signal_band, (7.8, 11.7))
    noise_band = get(cfg.bands, snr_noise_band,  (30.0, 50.0))

    hfnoise_v = Vector{Float64}(undef, n_ch)
    snr_v     = Vector{Float64}(undef, n_ch)

    for c in 1:n_ch
        freqs, psd = welch_psd_raw(rec.data[c, :], fs, welch_nfft)
        if isempty(freqs)
            @warn "compute_channel_spectral_qc: señal demasiado corta para nfft=$welch_nfft" channel=ch[c]
            hfnoise_v[c] = NaN
            snr_v[c]     = NaN
            continue
        end

        bandpower(f1, f2) = sum(psd[i] for (i, f) in enumerate(freqs) if f1 <= f <= f2; init=0.0)

        p_ref = bandpower(hfnoise_ref_band...)
        if p_ref <= 0
            @warn "compute_channel_spectral_qc: potencia de referencia HFNoise nula" channel=ch[c]
            hfnoise_v[c] = NaN
        else
            hfnoise_v[c] = bandpower(hfnoise_band...) / p_ref
        end

        p_sig   = bandpower(sig_band...)
        p_noise = bandpower(noise_band...)
        if p_sig <= 0 || p_noise <= 0
            @warn "compute_channel_spectral_qc: SNR indefinido (potencia nula)" channel=ch[c]
            snr_v[c] = NaN
        else
            snr_v[c] = 10 * log10(p_sig / p_noise)
        end
    end

    return DataFrame(
        channel = ch,
        hfnoise = round.(hfnoise_v, digits=4),
        snr_db  = round.(snr_v, digits=2),
    )
end

"""
    compute_correlation_summary(rec::EEGRecording) -> DataFrame

Resumen por canal de la matriz de correlación de Pearson ρᵢⱼ (señal raw):
`mean_abs_corr`, `max_corr`, `min_corr` frente al resto del montaje
(diagonal excluida). Report_Pre §Correlación intercanal.
"""
function compute_correlation_summary(rec::EEGRecording)::DataFrame
    n_ch = rec.meta.n_channels
    ch   = rec.meta.channel_names

    # cor() trata columnas como variables → transponer (canales × muestras)
    C = cor(permutedims(rec.data))

    mean_abs = Vector{Float64}(undef, n_ch)
    max_c    = Vector{Float64}(undef, n_ch)
    min_c    = Vector{Float64}(undef, n_ch)

    for c in 1:n_ch
        others = [C[c, j] for j in 1:n_ch if j != c]
        valid  = filter(!isnan, others)
        if length(valid) < length(others)
            @warn "compute_correlation_summary: correlación indefinida (posible canal plano)" channel=ch[c]
        end
        if isempty(valid)
            mean_abs[c] = NaN
            max_c[c]    = NaN
            min_c[c]    = NaN
        else
            mean_abs[c] = mean(abs.(valid))
            max_c[c]    = maximum(valid)
            min_c[c]    = minimum(valid)
        end
    end

    return DataFrame(
        channel       = ch,
        mean_abs_corr = round.(mean_abs, digits=3),
        max_corr      = round.(max_c,    digits=3),
        min_corr      = round.(min_c,    digits=3),
    )
end

# ─── Envoltorio legacy ─────────────────────────────────────────

"""
    qc_report(rec::EEGRecording, cfg::PipelineConfig) -> DataFrame

Informe QC completo (stats + is_bad) guardado como CSV en results/.
"""
function qc_report(rec::EEGRecording, cfg::PipelineConfig)::DataFrame
    stats  = compute_channel_stats(rec)
    bad_ch = flag_bad_channels(rec)
    stats[!, :is_bad] = [ch in bad_ch for ch in stats.channel]

    out_dir = ensure_dirs(cfg, rec.meta.subject_id, rec.meta.session_id)
    path    = joinpath(out_dir, "tables", "qc_channels_$(rec.meta.condition).csv")
    CSV.write(path, stats)
    return stats
end
