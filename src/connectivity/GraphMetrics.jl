# NeuroMIND/src/connectivity/GraphMetrics.jl
# Métricas de teoría de grafos sobre matrices de conectividad.

"""
    compute_graph_metrics(conn::ConnectivityMatrix, band::String,
                          cfg::PipelineConfig) -> GraphMetrics

Umbraliza la matriz de conectividad y calcula métricas de red.
"""
function compute_graph_metrics(
    conn::ConnectivityMatrix,
    band::String,
    cfg::PipelineConfig
)::GraphMetrics

    haskey(conn.matrices, band) || error("Banda '$band' no encontrada en ConnectivityMatrix")

    g_cfg  = cfg.graph
    method = get(g_cfg, "threshold_method", "proportional")
    density = get(g_cfg, "density", 0.1)

    W = conn.matrices[band]
    n = size(W, 1)

    # Umbralización
    A = _threshold_matrix(W, method, density)
    actual_density = sum(A) / (n * (n - 1))

    # Fuerza nodal (weighted degree)
    strength = vec(sum(W .* A, dims=2))

    # Clustering coefficient (Barrat et al. 2004, weighted)
    clustering = _weighted_clustering(W, A, n)

    # Distancia media y eficiencia (sobre grafo binario)
    path_length, efficiency = _path_length_efficiency(A, n)

    # Modularidad (Louvain simplificado — placeholder)
    modularity = 0.0  # TODO: integrar LightGraphs cuando esté disponible

    threshold_val = _get_threshold(W, method, density)

    return GraphMetrics(band, threshold_val, actual_density,
                        strength, clustering, path_length, efficiency, modularity)
end

# ─── Helpers privados ─────────────────────────────────────────

function _threshold_matrix(W::Matrix{Float64}, method::String, density::Real)::Matrix{Float64}
    n    = size(W, 1)
    A    = zeros(Float64, n, n)
    thr  = _get_threshold(W, method, density)
    @inbounds for i in 1:n, j in 1:n
        i != j && W[i,j] > thr && (A[i,j] = 1.0)
    end
    return A
end

function _get_threshold(W::Matrix{Float64}, method::String, density::Real)::Float64
    if method == "proportional"
        vals = sort(W[W .> 0.0], rev=true)
        n_keep = round(Int, density * length(vals))
        return n_keep > 0 ? vals[n_keep] : 0.0
    else  # "absolute"
        return density
    end
end

function _weighted_clustering(W::Matrix{Float64}, A::Matrix{Float64}, n::Int)::Vector{Float64}
    cc = zeros(Float64, n)
    @inbounds for i in 1:n
        nbrs = findall(j -> A[i,j] > 0, 1:n)
        k = length(nbrs)
        k < 2 && continue
        s_i = sum(W[i, nbrs])
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
    # BFS para distancias — O(n²)
    dist = fill(Inf, n, n)
    for i in 1:n
        dist[i, i] = 0.0
        queue = [i]
        while !isempty(queue)
            u = popfirst!(queue)
            for v in 1:n
                A[u, v] > 0 && isinf(dist[i, v]) && (dist[i, v] = dist[i, u] + 1; push!(queue, v))
            end
        end
    end

    finite_d = dist[isfinite.(dist) .& (dist .> 0)]
    path_length = isempty(finite_d) ? Inf : mean(finite_d)
    efficiency  = mean(1.0 ./ dist[dist .> 0])
    isnan(efficiency) && (efficiency = 0.0)

    return path_length, efficiency
end
