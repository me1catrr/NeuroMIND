# NeuroMIND/src/visualization/GraphPlots.jl
# Visualización de métricas de red y comparación grupal de grafos.

using CairoMakie

"""
    plot_graph_metrics(gm_ms, gm_ctrl, band; save_path=nothing) -> Figure

Barplot comparando métricas de red (strength, clustering, efficiency)
entre grupo MS y controles para una banda de frecuencia.
"""
function plot_graph_metrics(
    gm_ms::Vector{GraphMetrics},
    gm_ctrl::Vector{GraphMetrics},
    band::String;
    save_path::Union{Nothing,String} = nothing
)::CairoMakie.Figure

    metrics_fn = [
        ("Strength media",     gm -> mean(gm.strength)),
        ("Clustering medio",   gm -> mean(gm.clustering)),
        ("Path Length",        gm -> gm.path_length),
        ("Efficiency global",  gm -> gm.efficiency),
    ]

    fig = CairoMakie.Figure(size=(1100, 320))
    CairoMakie.Label(fig[0, :], "Graph Metrics — Banda $(band): MS vs Control";
                     fontsize=14, font=:bold)

    for (col, (name, fn)) in enumerate(metrics_fn)
        vals_ms   = [fn(g) for g in gm_ms   if isfinite(fn(g))]
        vals_ctrl = [fn(g) for g in gm_ctrl if isfinite(fn(g))]
        isempty(vals_ms) || isempty(vals_ctrl) && continue

        ax = CairoMakie.Axis(fig[1, col]; title=name, titlesize=11,
                              xticks=([1,2], ["MS","Control"]),
                              xticklabelsize=10, yticklabelsize=9)

        # Boxplots
        CairoMakie.boxplot!(ax, fill(1, length(vals_ms)),   vals_ms;
                            color=(:tomato, 0.6), strokewidth=1)
        CairoMakie.boxplot!(ax, fill(2, length(vals_ctrl)), vals_ctrl;
                            color=(:steelblue, 0.6), strokewidth=1)

        # Puntos individuales
        CairoMakie.scatter!(ax, fill(1, length(vals_ms)),   vals_ms;
                            color=(:darkred, 0.5), markersize=5)
        CairoMakie.scatter!(ax, fill(2, length(vals_ctrl)), vals_ctrl;
                            color=(:navy, 0.5), markersize=5)
    end

    save_path !== nothing && CairoMakie.save(save_path, fig; px_per_unit=2)
    return fig
end

"""
    plot_surrogate_distribution(surr, ch_i, ch_j; save_path=nothing) -> Figure

Histograma de la distribución nula de surrogates para un par de canales,
con el valor observado de wPLI marcado.
"""
function plot_surrogate_distribution(
    surr::SurrogateResult,
    ch_i::Int,
    ch_j::Int;
    save_path::Union{Nothing,String} = nothing
)::CairoMakie.Figure

    null   = vec(surr.null_distribution[ch_i, ch_j, :])
    obs    = surr.observed[ch_i, ch_j]
    p_val  = surr.p_values[ch_i, ch_j]
    ch_a   = surr.connectivity.channel_names[ch_i]
    ch_b   = surr.connectivity.channel_names[ch_j]

    fig = CairoMakie.Figure(size=(550, 380))
    ax  = CairoMakie.Axis(fig[1, 1];
        title  = "Surrogates: $(surr.band) — $(ch_a) ↔ $(ch_b)",
        xlabel = "wPLI (distribución nula)",
        ylabel = "Frecuencia",
        titlesize = 12,
    )

    CairoMakie.hist!(ax, null; bins=30, color=(:steelblue, 0.6), strokewidth=0.5)
    CairoMakie.vlines!(ax, [obs]; color=:tomato, linewidth=2,
                       label="Observado = $(round(obs, digits=3))\np = $(round(p_val, digits=4))")
    CairoMakie.vlines!(ax, [surr.fdr_threshold]; color=:orange, linewidth=1.5,
                       linestyle=:dash, label="FDR thr = $(round(surr.fdr_threshold, digits=4))")
    CairoMakie.axislegend(ax; position=:rt, labelsize=9)

    save_path !== nothing && CairoMakie.save(save_path, fig; px_per_unit=2)
    return fig
end

"""
    save_figure(fig, path; dpi=300)

Guarda una figura CairoMakie con resolución configurada.
"""
function save_figure(fig::CairoMakie.Figure, path::String; dpi::Int=300)
    px_per_unit = dpi / 96
    CairoMakie.save(path, fig; px_per_unit)
end
