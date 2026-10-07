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
