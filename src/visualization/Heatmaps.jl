# NeuroMIND/src/visualization/Heatmaps.jl
# Heatmaps de matrices de conectividad wPLI.

using CairoMakie

"""
    plot_connectivity_heatmap(conn::ConnectivityMatrix, band::String;
                              title="", colormap=:viridis, clims=nothing) -> Figure

Genera un heatmap de la matriz wPLI para la banda especificada.
"""
function plot_connectivity_heatmap(
    conn::ConnectivityMatrix,
    band::String;
    title::String = "",
    colormap = :viridis,
    clims::Union{Nothing,Tuple{Float64,Float64}} = nothing
)::CairoMakie.Figure

    haskey(conn.matrices, band) || error("Banda '$band' no encontrada")
    W = conn.matrices[band]
    ch = conn.channel_names
    n  = length(ch)

    fig = CairoMakie.Figure(size=(800, 700))
    ax  = CairoMakie.Axis(fig[1, 1];
        title  = isempty(title) ? "wPLI — $(band)" : title,
        xlabel = "Canal",
        ylabel = "Canal",
        xticks = (1:n, ch),
        yticks = (1:n, ch),
        xticklabelrotation = π/3,
        xticklabelsize = 8,
        yticklabelsize = 8,
    )

    hm = if clims !== nothing
        CairoMakie.heatmap!(ax, W; colormap, colorrange=clims)
    else
        CairoMakie.heatmap!(ax, W; colormap)
    end
    CairoMakie.Colorbar(fig[1, 2], hm; label="wPLI")
    return fig
end

"""
    plot_group_comparison(ga::GroupAnalysis; band="", colormap=:RdBu_r) -> Figure

Heatmap de diferencia de conectividad entre grupos (A - B) con máscara de significancia.
"""
function plot_group_comparison(
    ga::GroupAnalysis;
    colormap = :RdBu_r
)::CairoMakie.Figure

    diff  = ga.mean_connectivity_ms .- ga.mean_connectivity_ctrl
    ch    = ga.channel_names
    n     = length(ch)

    lim = maximum(abs.(diff))

    fig = CairoMakie.Figure(size=(1400, 600))

    for (col, (mat, ttl)) in enumerate([
        (ga.mean_connectivity_ms,   ga.group_ms),
        (ga.mean_connectivity_ctrl, ga.group_ctrl),
        (diff, "$(ga.group_ms) − $(ga.group_ctrl)")
    ])
        ax = CairoMakie.Axis(fig[1, col];
            title  = "wPLI $(ga.band) — $ttl\n(n=$(col ≤ 2 ? (col==1 ? ga.n_ms : ga.n_ctrl) : ga.n_ms + ga.n_ctrl))",
            xticks = (1:n, ch), yticks = (1:n, ch),
            xticklabelrotation = π/3,
            xticklabelsize = 7, yticklabelsize = 7,
        )
        cr = col == 3 ? (-lim, lim) : (0.0, maximum(mat))
        cm = col == 3 ? colormap : :viridis
        hm = CairoMakie.heatmap!(ax, mat; colormap=cm, colorrange=cr)
        CairoMakie.Colorbar(fig[2, col], hm; vertical=false, label="wPLI")
    end
    return fig
end
