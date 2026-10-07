# NeuroMIND/src/statistics/Surrogates.jl
# Inferencia estadística mediante surrogates de desplazamiento circular para wPLI.
#
# MÉTODO: circular_shift — desplazamiento circular independiente por canal y época.
#
# Propiedades de la nulidad generada:
#   ✓ Destruye la sincronía de fase inter-canal (lo que mide wPLI)
#   ✓ Preserva el espectro de potencia de cada canal
#   ✓ Preserva la autocorrelación temporal dentro de cada canal/época
#   ✓ El filtrado bandpass ocurre UNA SOLA VEZ (igual que en el observado)
#   ✓ p-valores con corrección Monte Carlo (+1) → nunca p = 0
#   ✓ Semilla reproducible vía índice estable (no hash(String))
#
# Referencias:
#   Vinck et al. (2011) NeuroImage — wPLI
#   Theiler et al. (1992) Physica D — surrogates en series temporales

"""
    surrogate_test(epochs, conn, band, cfg; on_progress=nothing) -> SurrogateResult

Prueba de significancia wPLI mediante surrogates de desplazamiento circular.

Genera `n_surrogates` permutaciones y construye la distribución nula del wPLI.
La máscara de significancia se obtiene aplicando FDR-BH al triángulo superior.

p-valor: (#{W_null >= W_obs} + 1) / (n_sur + 1)  — nunca puede ser 0.

`on_progress` (opcional): callback `(k::Int, n_sur::Int) -> Nothing` llamado
tras cada permutación (para barras de progreso en terminal).
"""
function surrogate_test(
    epochs::EpochSet,
    conn::ConnectivityMatrix,
    band::String,
    cfg::PipelineConfig;
    on_progress = nothing,
)::SurrogateResult

    sur_cfg   = cfg.surrogates
    n_sur     = Int(get(sur_cfg, "n_surrogates", 200))
    alpha     = Float64(get(sur_cfg, "alpha", 0.05))
    seed_v    = Int(get(sur_cfg, "seed", 42))

    haskey(conn.matrices, band) || error("Banda '$band' no encontrada en ConnectivityMatrix")
    W_obs = conn.matrices[band]

    (f1, f2) = cfg.bands[band]
    n_ch, n_samp, n_seg = size(epochs.data)
    fs    = epochs.meta.fs

    # El estimador de la distribución nula DEBE ser el mismo que el observado.
    estimator = _build_estimator(cfg.connectivity)

    # Semilla reproducible: índice estable basado en posición alfabética de la banda.
    # NO usamos hash(band) porque no es estable entre versiones de Julia.
    band_names = sort(collect(keys(cfg.bands)))
    band_idx   = something(findfirst(==(band), band_names), 1)
    rng        = MersenneTwister(seed_v + band_idx * 1000)

    # ─── Distribución nula ──────────────────────────────────────
    # Circular shift independiente por canal: destruye fase inter-canal
    # sin alterar espectro ni autocorrelación de cada canal/época.
    # _compute_band_matrix usa el mismo estimador que la corrida observada
    # → distribución nula metodológicamente coherente (Hilbert/FourierCSD/Multitaper).
    W_null = zeros(Float64, n_ch, n_ch, n_sur)

    for k in 1:n_sur
        data_sur = _circular_shift_surrogate(epochs, rng)
        W_null[:, :, k] = _compute_band_matrix(
            estimator, data_sur, fs, f1, f2, n_ch, n_samp, n_seg
        )
        if on_progress !== nothing
            on_progress(k, n_sur)
        end
    end

    # ─── p-valores con corrección Monte Carlo ───────────────────
    # Fórmula: (#{W_null >= W_obs} + 1) / (n_sur + 1)
    # Garantiza p ∈ [1/(n_sur+1), 1.0], nunca p = 0.
    # Diagonal = 1.0: no hay test de auto-conectividad (wPLI[i,i] = 0 por diseño).
    p_values = ones(Float64, n_ch, n_ch)
    @inbounds for i in 1:n_ch, j in (i+1):n_ch
        cnt = count(W_null[i, j, :] .>= W_obs[i, j])
        p   = (cnt + 1) / (n_sur + 1)
        p_values[i, j] = p
        p_values[j, i] = p
    end

    # ─── FDR Benjamini-Hochberg (triángulo superior) ────────────
    # La familia de hipótesis es el triángulo superior de una banda.
    # FDR por banda es apropiado si cada banda se interpreta como
    # familia independiente. Para corrección global (todas las bandas),
    # apilar los p-vec de todas las bandas antes de llamar fdr_correction.
    upper_idx = [(i, j) for i in 1:n_ch for j in (i+1):n_ch]
    p_vec     = [p_values[i, j] for (i, j) in upper_idx]
    fdr_thr   = fdr_correction(p_vec; alpha, method=String(get(sur_cfg, "fdr_method", "bh")))

    sig_mask = falses(n_ch, n_ch)
    for (k, (i, j)) in enumerate(upper_idx)
        if p_vec[k] <= fdr_thr
            sig_mask[i, j] = true
            sig_mask[j, i] = true
        end
    end

    return SurrogateResult(conn, band, W_obs, W_null, p_values, sig_mask, fdr_thr, n_sur)
end

# ─── Validaciones ─────────────────────────────────────────────

"""
    validate_connectivity_matrix(W, band="?")

Emite advertencias `@warn` si la matriz wPLI contiene NaN/Inf,
valores fuera de [0,1], asimetría o diagonal no nula.
No lanza error — devuelve `nothing` siempre.
"""
function validate_connectivity_matrix(W::Matrix{Float64}, band::String="?")
    any(isnan, W)  && @warn "[$band] wPLI contiene NaN"
    any(isinf, W)  && @warn "[$band] wPLI contiene Inf"
    any(x -> x < -1e-6 || x > 1.0 + 1e-6, W) &&
        @warn "[$band] Valores de wPLI fuera de [0,1]: min=$(minimum(W)), max=$(maximum(W))"
    (W ≈ W') || @warn "[$band] Matriz wPLI no es simétrica (max_asym=$(maximum(abs.(W .- W'))))"
    any(!iszero, diag(W)) &&
        @warn "[$band] Diagonal de wPLI no es cero (max_diag=$(maximum(abs.(diag(W)))))"
    return nothing
end

"""
    validate_surrogate_result(sr)

Emite advertencias si los p-valores de un `SurrogateResult` están fuera del
rango teórico [(1/(n_sur+1)), 1.0]. Detecta la falta de corrección +1.
"""
function validate_surrogate_result(sr::SurrogateResult)
    p_min_theory = 1.0 / (sr.n_surrogates + 1)
    if any(x -> x < p_min_theory - 1e-10 && x > 0.0, sr.p_values)
        @warn "[$(sr.band)] p-valores menores que el mínimo teórico " *
              "$(round(p_min_theory, digits=5)) — posible falta de corrección +1"
    end
    if any(x -> x <= 0.0, sr.p_values)
        @warn "[$(sr.band)] p-valores ≤ 0 — falta corrección Monte Carlo (+1)"
    end
    if any(x -> x > 1.0 + 1e-10, sr.p_values)
        @warn "[$(sr.band)] p-valores > 1"
    end
    return nothing
end

# ─── Helpers privados ─────────────────────────────────────────

"""
    _circular_shift_surrogate(epochs, rng) -> Array{Float64,3}

Genera datos surrogate mediante desplazamiento circular INDEPENDIENTE
por canal y época. Devuelve el array de datos brutos (sin filtrar).

El cálculo de wPLI posterior usa `_compute_band_matrix(estimator, ...)` con
el mismo estimador que la corrida observada — distribución nula comparable.

El desplazamiento mínimo es `n_samp ÷ 5` muestras para evitar que
un shift muy pequeño deje prácticamente intacta la relación de fase.
"""
function _circular_shift_surrogate(
    epochs::EpochSet,
    rng::AbstractRNG
)::Array{Float64,3}
    n_ch, n_samp, n_seg = size(epochs.data)
    out       = similar(epochs.data)
    min_shift = max(1, n_samp ÷ 5)
    max_shift = max(min_shift, n_samp - min_shift)

    @inbounds for seg in 1:n_seg, ch in 1:n_ch
        shift = (min_shift == max_shift) ? min_shift : rand(rng, min_shift:max_shift)
        out[ch, :, seg] = circshift(view(epochs.data, ch, :, seg), shift)
    end
    return out
end
