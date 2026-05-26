# NeuroMIND/src/ica/ICAClassification.jl
#
# Clasificación automática de componentes ICA.
# Portado de EEG_Julia/src/ICA/ICA_cleaning.jl (sin Plots, Julia puro).
#
# Features por componente:
#   frontal_ratio, temporal_ratio, blink_ratio, emg_ratio, line_ratio,
#   kurtosis, extreme_frac
# Scores:
#   ocular_score, muscle_score, line_score, jump_score → artifact_score
# Etiquetas: "brain", "eye/blink", "muscle", "line_noise", "jump", "unknown"

# ─── Auxiliares ────────────────────────────────────────────────

function _zscore(v::AbstractVector)
    m = mean(v); s = std(v)
    s == 0.0 && return zeros(length(v))
    return (v .- m) ./ s
end

function _kurtosis_simple(x::AbstractVector)
    m  = mean(x); xc = x .- m
    s2 = mean(xc .^ 2)
    s2 == 0.0 && return 0.0
    return mean(xc .^ 4) / (s2^2)
end

function _bandpower(x::AbstractVector, fs::Real, f_lo::Real, f_hi::Real)
    x  = x .- mean(x)
    N  = length(x)
    X  = rfft(x)
    freqs = (0:length(X)-1) .* (fs / N)
    psd = abs.(X) .^ 2 ./ (fs * N)
    return sum(psd[(freqs .>= f_lo) .& (freqs .<= f_hi)])
end

# ─── Features por IC ───────────────────────────────────────────

"""
    compute_ica_features(A, S, fs, ch_names) -> DataFrame

Calcula 7 features por componente ICA.
  - frontal_ratio, temporal_ratio: distribución espacial (A)
  - blink_ratio, emg_ratio, line_ratio: potencias de banda (S)
  - kurtosis, extreme_frac: estadísticos temporales (S)
"""
function compute_ica_features(
    A::AbstractMatrix{Float64},
    S::AbstractMatrix{Float64},
    fs::Real,
    ch_names::Vector{String}
)::DataFrame

    n_ch, n_ic = size(A)
    all_idx = 1:n_ch

    # Índices de canales frontales y temporales (nomenclatura 10-20)
    frontal_idx  = [i for (i, ch) in enumerate(ch_names)
                    if occursin(r"^F[PC0-9ZCT]?"i, uppercase(ch)) ||
                       startswith(uppercase(ch), "FP")]
    temporal_idx = [i for (i, ch) in enumerate(ch_names)
                    if startswith(uppercase(ch), "T") ||
                       startswith(uppercase(ch), "TP")]

    if isempty(frontal_idx);  frontal_idx  = collect(1:min(4, n_ch)); end
    if isempty(temporal_idx); temporal_idx = collect(max(1, n_ch-3):n_ch); end

    nonfrontal_idx  = setdiff(all_idx, frontal_idx)
    nontemporal_idx = setdiff(all_idx, temporal_idx)

    fr   = zeros(n_ic); tr   = zeros(n_ic)
    blnk = zeros(n_ic); emg  = zeros(n_ic)
    lnr  = zeros(n_ic); kurt = zeros(n_ic)
    exfr = zeros(n_ic)

    n_samp = size(S, 2)

    for k in 1:n_ic
        map_k = A[:, k]
        s_k   = S[k, :]

        # Spatial ratios
        fr[k] = mean(abs.(map_k[frontal_idx]))  /
                (mean(abs.(map_k[nonfrontal_idx]))  + 1e-12)
        tr[k] = mean(abs.(map_k[temporal_idx])) /
                (mean(abs.(map_k[nontemporal_idx])) + 1e-12)

        # Spectral ratios
        P04   = _bandpower(s_k, fs, 0.5,  4.0)
        P440  = _bandpower(s_k, fs, 4.0, 40.0)
        P130  = _bandpower(s_k, fs, 1.0, 30.0)
        P3080 = _bandpower(s_k, fs, 30.0, 80.0)
        Ptot  = _bandpower(s_k, fs, 0.5, 100.0)
        P4852 = _bandpower(s_k, fs, 48.0, 52.0)

        blnk[k] = P04   / (P440  + 1e-12)
        emg[k]  = P3080 / (P130  + 1e-12)
        lnr[k]  = P4852 / (Ptot  + 1e-12)

        # Statistical
        kurt[k] = _kurtosis_simple(s_k)
        μ = mean(s_k); σ = std(s_k)
        exfr[k] = σ == 0.0 ? 0.0 :
            count(t -> abs(t - μ) > 5σ, s_k) / n_samp
    end

    return DataFrame(
        component       = 1:n_ic,
        frontal_ratio   = round.(fr,   digits=4),
        temporal_ratio  = round.(tr,   digits=4),
        blink_ratio     = round.(blnk, digits=4),
        emg_ratio       = round.(emg,  digits=4),
        line_ratio      = round.(lnr,  digits=6),
        kurtosis        = round.(kurt, digits=4),
        extreme_frac    = round.(exfr, digits=6),
    )
end

# ─── Scores y clasificación ────────────────────────────────────

"""
    evaluate_ica_components(feat_df; artifact_thresh=1.5) -> DataFrame

Añade scores y etiqueta al DataFrame de features.
Etiquetas: "brain", "eye/blink", "muscle", "line_noise", "jump", "unknown"
"""
function evaluate_ica_components(
    feat_df::DataFrame;
    artifact_thresh::Real = 1.5
)::DataFrame

    fz  = _zscore(feat_df.frontal_ratio)
    tz  = _zscore(feat_df.temporal_ratio)
    blz = _zscore(feat_df.blink_ratio)
    emz = _zscore(feat_df.emg_ratio)
    lnz = _zscore(feat_df.line_ratio)
    kz  = _zscore(feat_df.kurtosis)
    exz = _zscore(feat_df.extreme_frac)

    n_ic = nrow(feat_df)

    ocular  = @. 0.4*fz + 0.4*blz + 0.2*kz
    muscle  = @. 0.5*emz + 0.3*tz + 0.2*kz
    line    = lnz
    jump    = @. 0.5*kz + 0.5*exz

    artifact_global = zeros(n_ic)
    labels          = Vector{String}(undef, n_ic)

    for k in 1:n_ic
        scores_k = (ocular[k], muscle[k], line[k], jump[k])
        mx = maximum(scores_k)
        artifact_global[k] = mx
        if mx > artifact_thresh
            labels[k] = if mx == ocular[k]; "eye/blink"
                        elseif mx == muscle[k]; "muscle"
                        elseif mx == line[k]; "line_noise"
                        else "jump"
                        end
        else
            labels[k] = "brain"
        end
    end

    result = copy(feat_df)
    result[!, :ocular_score]   = round.(ocular, digits=3)
    result[!, :muscle_score]   = round.(muscle, digits=3)
    result[!, :line_score]     = round.(line,   digits=3)
    result[!, :jump_score]     = round.(jump,   digits=3)
    result[!, :artifact_score] = round.(artifact_global, digits=3)
    result[!, :artifact_type]  = labels
    return result
end
