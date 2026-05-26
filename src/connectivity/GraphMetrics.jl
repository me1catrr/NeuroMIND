# NeuroMIND/src/connectivity/GraphMetrics.jl
# Métricas de teoría de grafos sobre matrices de conectividad wPLI.
#
# Diseño:
#   - Las métricas ponderadas (strength, clustering) usan la matriz wPLI continua.
#   - La umbralización para el grafo binario (path_length, efficiency) usa
#     la densidad configurada en [graph], NO el umbral de surrogates/FDR.
#   - Esto garantiza que las métricas son comparables entre sujetos,
#     independientemente del número de surrogates usado.

"""
    compute_graph_metrics(conn, band, cfg) -> GraphMetrics

Calcula métricas de red para una banda de frecuencia.

El umbral para el grafo binario viene de `cfg.graph["density"]` (proporción de
aristas más fuertes a mantener) o `cfg.graph["threshold"]` (valor absoluto).
Las métricas ponderadas usan la matriz wPLI sin umbralizar.
"""
function compute_graph_metrics(
    conn::ConnectivityMatrix,
    band::String,
    cfg::PipelineConfig
)::GraphMetrics

    haskey(conn.matrices, band) || error("Banda '$band' no encontrada en ConnectivityMatrix")

    g_cfg    = cfg.graph
    method   = String(get(g_cfg, "threshold_method", "proportional"))
    density  = Float64(get(g_cfg, "density", 0.1))

    W = conn.matrices[band]
    n = size(W, 1)

    # ─── Grafo binario (para path_length y efficiency) ────────
    A            = _threshold_matrix(W, method, density)
    actual_density = sum(A) / max(1, n * (n - 1))
    threshold_val  = _get_threshold(W, method, density)

    # ─── Métricas ponderadas (usan W directamente) ────────────
    # Strength: suma de pesos de todas las aristas de un nodo
    strength   = vec(sum(W, dims=2)) .- diag(W)   # excluir diagonal

    # Clustering coefficient ponderado (Barrat et al. 2004)
    clustering = _weighted_clustering(W, A, n)

    # Eficiencia y longitud de camino (sobre grafo binario)
    path_length, efficiency = _path_length_efficiency(A, n)

    # Modularidad (placeholder — requiere algoritmo Louvain)
    modularity = 0.0

    return GraphMetrics(
        band, threshold_val, actual_density,
        strength, clustering, path_length, efficiency, modularity,
        conn.channel_names
    )
end

# ─── Helpers privados ─────────────────────────────────────────

function _threshold_matrix(W::Matrix{Float64}, method::String, density::Real)::Matrix{Float64}
    n   = size(W, 1)
    A   = zeros(Float64, n, n)
    thr = _get_threshold(W, method, density)
    @inbounds for i in 1:n, j in 1:n
        i != j && W[i,j] > thr && (A[i,j] = 1.0)
    end
    return A
end

function _get_threshold(W::Matrix{Float64}, method::String, density::Real)::Float64
    if method == "proportional"
        # Mantener las `density` × 100% aristas más fuertes del triángulo superior
        vals = sort([W[i,j] for i in 1:size(W,1) for j in (i+1):size(W,1)
                     if W[i,j] > 0.0], rev=true)
        n_keep = round(Int, density * length(vals))
        return (n_keep > 0 && n_keep <= length(vals)) ? vals[n_keep] : 0.0
    else  # "absolute"
        return Float64(density)
    end
end

function _weighted_clustering(W::Matrix{Float64}, A::Matrix{Float64}, n::Int)::Vector{Float64}
    cc = zeros(Float64, n)
    @inbounds for i in 1:n
        nbrs = findall(j -> A[i,j] > 0, 1:n)
        k = length(nbrs)
        k < 2 && continue
        s_i = sum(W[i, j] for j in nbrs)
        s_i ≈ 0 && continue
        tri = 0.0
        for u in nbrs, v in nbrs
            u == v && continue
            tri += (W[i,u] * W[i,v] * A[u,v])^(1/3)
        end
        cc[i] = tri / (s_i * (k - 1))
    end
    return cc
end

function _path_length_efficiency(A::Matrix{Float64}, n::Int)
    dist = fill(Inf, n, n)
    for i in 1:n
        dist[i, i] = 0.0
        queue = [i]
        while !isempty(queue)
            u = popfirst!(queue)
            for v in 1:n
                A[u,v] > 0 && isinf(dist[i,v]) && (dist[i,v] = dist[i,u] + 1; push!(queue, v))
            end
        end
    end
    finite_d   = dist[isfinite.(dist) .& (dist .> 0)]
    path_length = isempty(finite_d) ? Inf : mean(finite_d)
    eff_vals    = [1.0/d for d in dist if d > 0 && isfinite(d)]
    efficiency  = isempty(eff_vals) ? 0.0 : mean(eff_vals)
    return path_length, efficiency
end
