# NeuroMIND/src/statistics/FDR.jl
# Corrección de tasa de falsos descubrimientos (FDR).

"""
    fdr_correction(p_values::Vector{Float64}; alpha=0.05, method="bh") -> Float64

Aplica corrección FDR Benjamini-Hochberg (default) o Bonferroni.
Devuelve el umbral de p-valor ajustado.
"""
function fdr_correction(
    p_values::Vector{Float64};
    alpha::Real = 0.05,
    method::String = "bh"
)::Float64

    method == "bonferroni" && return alpha / length(p_values)
    return _bh_threshold(p_values, Float64(alpha))
end

"""
    threshold_connectivity(conn::ConnectivityMatrix, surr::SurrogateResult) -> Matrix{Float64}

Devuelve la matriz de conectividad enmascarada: pone a 0 los edges no significativos.
"""
function threshold_connectivity(
    conn::ConnectivityMatrix,
    surr::SurrogateResult
)::Matrix{Float64}
    W = get(conn.matrices, surr.band, nothing)
    W === nothing && error("Banda '$(surr.band)' no encontrada")
    return W .* surr.sig_mask
end

# ─── BH privado ───────────────────────────────────────────────

function _bh_threshold(p::Vector{Float64}, alpha::Float64)::Float64
    m     = length(p)
    order = sortperm(p)
    ps    = p[order]

    # Buscar mayor k tal que p_{(k)} ≤ k/m * α
    k_max = 0
    for k in 1:m
        ps[k] <= k / m * alpha && (k_max = k)
    end
    k_max == 0 && return 0.0
    return ps[k_max]
end
