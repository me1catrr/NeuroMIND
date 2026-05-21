# NeuroMIND/src/connectivity/CSD.jl
# Current Source Density (CSD) — módulo OPCIONAL.
# Solo se aplica si cfg.connectivity["use_csd"] == true.

"""
    apply_csd(epochs::EpochSet, cfg::PipelineConfig) -> EpochSet

Aplica la transformación CSD (spline Laplaciano) a los epochs si está habilitado.
Si `use_csd = false` en la configuración, devuelve los epochs sin cambios.
"""
function apply_csd(epochs::EpochSet, cfg::PipelineConfig)::EpochSet
    get(cfg.connectivity, "use_csd", false) || return epochs

    ch_pos = epochs.meta.channel_positions
    if ch_pos === nothing
        @warn "CSD activado pero no hay posiciones de electrodos — omitiendo CSD"
        return epochs
    end

    n_ch, n_samp, n_seg = size(epochs.data)
    data_csd = similar(epochs.data)

    # Laplaciano de superficie (spline esférico simplificado)
    G, H = _spline_matrices(ch_pos, epochs.meta.channel_names)

    @inbounds for seg in 1:n_seg
        snapshot = @view epochs.data[:, :, seg]   # (ch × samp)
        data_csd[:, :, seg] = _apply_surface_laplacian(snapshot, G, H)
    end

    return EpochSet(
        epochs.meta, data_csd, epochs.epoch_length_s,
        epochs.n_valid, epochs.rejected_idx
    )
end

# ─── Helpers privados ─────────────────────────────────────────

function _spline_matrices(
    ch_pos::Dict{String,Tuple{Float64,Float64}},
    ch_names::Vector{String}
)
    n = length(ch_names)
    G = zeros(n, n)
    H = zeros(n, n)

    xs = [get(ch_pos, uppercase(ch), (0.0, 0.0))[1] for ch in ch_names]
    ys = [get(ch_pos, uppercase(ch), (0.0, 0.0))[2] for ch in ch_names]

    for i in 1:n, j in 1:n
        i == j && continue
        dx   = xs[i] - xs[j]
        dy   = ys[i] - ys[j]
        dist = sqrt(dx^2 + dy^2 + 1e-10)
        # Función de Green del spline esférico (aproximación 2D)
        g_val = dist^2 * log(dist + 1e-10)
        G[i, j] = g_val
        H[i, j] = -2.0 * log(dist + 1e-10) - 1.0
    end
    return G, H
end

function _apply_surface_laplacian(
    data::AbstractMatrix{Float64},  # (ch × samples)
    G::Matrix{Float64},
    H::Matrix{Float64}
)::Matrix{Float64}

    n_ch, n_samp = size(data)
    # Resolver sistema G·C = data' → C: coeficientes spline
    C = G \ data'   # (ch × samples)
    # Laplaciano: H·C
    return (H * C)'  # → (ch × samples)
end
