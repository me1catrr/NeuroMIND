# NeuroMIND/src/visualization/Topomaps.jl
# Mapas topográficos de potencia espectral y conectividad nodal.

using CairoMakie

"""
    plot_topomap(values::Vector{Float64}, ch_names::Vector{String},
                 ch_pos::Dict; title="", colormap=:jet, grid_res=200) -> Figure

Genera un topomap 2D por IDW (Inverse Distance Weighting) estilo BrainVision.
"""
function plot_topomap(
    values::Vector{Float64},
    ch_names::Vector{String},
    ch_pos::Dict{String,Tuple{Float64,Float64}};
    title::String = "",
    colormap = :jet,
    grid_res::Int = 200,
    clims::Union{Nothing,Tuple{Float64,Float64}} = nothing
)::CairoMakie.Figure

    xs = [get(ch_pos, uppercase(ch), (NaN, NaN))[1] for ch in ch_names]
    ys = [get(ch_pos, uppercase(ch), (NaN, NaN))[2] for ch in ch_names]

    # IDW sobre grid circular
    xg = range(-1.0, 1.0; length=grid_res)
    yg = range(-1.0, 1.0; length=grid_res)
    Z  = fill(NaN, grid_res, grid_res)

    valid = findall(i -> !isnan(xs[i]) && !isnan(ys[i]), eachindex(xs))
    xv = xs[valid]; yv = ys[valid]; vv = values[valid]

    @inbounds for j in 1:grid_res, i in 1:grid_res
        xx, yy = xg[i], yg[j]
        xx^2 + yy^2 > 1.0 && continue
        d2 = (xx .- xv).^2 .+ (yy .- yv).^2
        w  = 1.0 ./ (d2 .+ 1e-6)
        Z[j, i] = sum(w .* vv) / sum(w)
    end

    fig = CairoMakie.Figure(size=(550, 600))
    ax  = CairoMakie.Axis(fig[1, 1];
        title  = title,
        aspect = DataAspect(),
        xlabel = "", ylabel = "",
    )

    cr = clims !== nothing ? clims : (minimum(vv), maximum(vv))
    hm = CairoMakie.heatmap!(ax, xg, yg, Z; colormap, colorrange=cr)
    CairoMakie.Colorbar(fig[2, 1], hm; vertical=false, label="μV²")

    # Contorno de la cabeza
    θ = range(0, 2π; length=200)
    CairoMakie.lines!(ax, sin.(θ), cos.(θ); color=:black, linewidth=2)

    # Nariz
    θn = range(-π/6, π/6; length=30)
    CairoMakie.lines!(ax, 0.12 .* sin.(θn), 1.0 .+ 0.10 .* cos.(θn); color=:black, linewidth=2)

    # Electrodos
    inside = [i for i in valid if xs[i]^2 + ys[i]^2 <= 1.0]
    CairoMakie.scatter!(ax, xs[inside], ys[inside];
                        color=:black, markersize=5,
                        strokecolor=:white, strokewidth=1)
    return fig
end

"""
    plot_longitudinal_evolution(la::LongitudinalAnalysis; colormap=:viridis) -> Figure

Visualiza la evolución temporal de la conectividad media de un sujeto.
"""
function plot_longitudinal_evolution(la::LongitudinalAnalysis; colormap=:viridis)::CairoMakie.Figure
    n_visits = length(la.visits)
    n_visits == 0 && error("Sin datos longitudinales")

    fig = CairoMakie.Figure(size=(300 * n_visits, 350))

    global_max = maximum(maximum(m) for m in la.connectivity_over_time)

    for (k, (visit, W)) in enumerate(zip(la.visits, la.connectivity_over_time))
        ax = CairoMakie.Axis(fig[1, k];
            title  = visit,
            aspect = DataAspect(),
        )
        CairoMakie.heatmap!(ax, W; colormap, colorrange=(0.0, global_max))
    end

    CairoMakie.Label(fig[0, :],
        "$(la.subject_id) — wPLI $(la.band) ($(la.condition)) evolución temporal";
        fontsize=14)
    return fig
end
