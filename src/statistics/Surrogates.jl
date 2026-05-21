# NeuroMIND/src/statistics/Surrogates.jl
# Inferencia estadística mediante surrogates de fase para wPLI.

"""
    surrogate_test(epochs::EpochSet, conn::ConnectivityMatrix,
                   band::String, cfg::PipelineConfig) -> SurrogateResult

Prueba de significancia mediante surrogates de fase para una banda.
Genera `n_surrogates` permutaciones de fase y construye distribución nula.
"""
function surrogate_test(
    epochs::EpochSet,
    conn::ConnectivityMatrix,
    band::String,
    cfg::PipelineConfig
)::SurrogateResult

    sur_cfg     = cfg.surrogates
    n_sur       = get(sur_cfg, "n_surrogates", 200)
    alpha       = get(sur_cfg, "alpha", 0.05)
    method      = get(sur_cfg, "method", "phase_shuffle")

    haskey(conn.matrices, band) || error("Banda '$band' no encontrada")
    W_obs = conn.matrices[band]

    (f1, f2) = cfg.bands[band]
    n_ch, n_samp, n_seg = size(epochs.data)
    fs       = epochs.meta.fs
    order    = get(cfg.connectivity, "filter_order", 8)

    # Distribución nula: wPLI de surrogates
    W_null = zeros(Float64, n_ch, n_ch, n_sur)

    for k in 1:n_sur
        epochs_sur = _phase_shuffle_epochs(epochs, f1, f2, fs, order)
        W_null[:, :, k] = _fast_wpli_band(epochs_sur, fs, f1, f2, order, n_ch, n_samp, n_seg)
    end

    # p-valor por par: proporción de surrogates ≥ observado
    p_values = zeros(Float64, n_ch, n_ch)
    @inbounds for i in 1:n_ch, j in (i+1):n_ch
        p = mean(W_null[i, j, :] .>= W_obs[i, j])
        p_values[i, j] = p
        p_values[j, i] = p
    end

    # FDR sobre el triángulo superior
    upper_idx = [(i, j) for i in 1:n_ch for j in (i+1):n_ch]
    p_vec     = [p_values[i, j] for (i, j) in upper_idx]
    fdr_thr   = fdr_correction(p_vec; alpha, method=get(sur_cfg, "fdr_method", "bh"))

    sig_mask = falses(n_ch, n_ch)
    for (k, (i, j)) in enumerate(upper_idx)
        if p_vec[k] <= fdr_thr
            sig_mask[i, j] = true
            sig_mask[j, i] = true
        end
    end

    return SurrogateResult(conn, band, W_obs, W_null, p_values, sig_mask, fdr_thr, n_sur)
end

# ─── Helpers privados ─────────────────────────────────────────

function _phase_shuffle_epochs(
    epochs::EpochSet,
    f1::Real, f2::Real,
    fs::Real, order::Int
)::Array{Float64,3}

    n_ch, n_samp, n_seg = size(epochs.data)
    out = copy(epochs.data)
    nyq = fs / 2.0
    bp  = digitalfilter(Bandpass(f1/nyq, f2/nyq), Butterworth(order))

    @inbounds for seg in 1:n_seg
        # Generar desplazamiento de fase aleatorio (mismo para todos los canales → preserva estructura)
        phase_shift = rand() * 2π
        for ch in 1:n_ch
            xf = filtfilt(bp, @view epochs.data[ch, :, seg])
            N  = length(xf)
            X  = fft(xf)
            X .*= exp.(1im .* phase_shift)
            out[ch, :, seg] = real.(ifft(X))
        end
    end
    return out
end

function _fast_wpli_band(
    data::Array{Float64,3},
    fs::Real, f1::Real, f2::Real,
    order::Int, n_ch::Int, n_samp::Int, n_seg::Int
)::Matrix{Float64}

    nyq = fs / 2.0
    bp  = digitalfilter(Bandpass(f1/nyq, f2/nyq), Butterworth(order))

    z = Vector{Matrix{ComplexF64}}(undef, n_ch)
    @inbounds for c in 1:n_ch
        Zc = Matrix{ComplexF64}(undef, n_samp, n_seg)
        for s in 1:n_seg
            xf = filtfilt(bp, @view data[c, :, s])
            Zc[:, s] = _analytic_signal(xf)
        end
        z[c] = Zc
    end

    W = zeros(Float64, n_ch, n_ch)
    @inbounds for i in 1:n_ch, j in (i+1):n_ch
        imv = imag.(z[i] .* conj.(z[j]))
        w   = abs(sum(imv)) / (sum(abs.(imv)) + eps())
        W[i, j] = w; W[j, i] = w
    end
    return W
end
