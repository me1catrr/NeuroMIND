# NeuroMIND/src/visualization/Topomaps.jl
# Mapas topográficos de potencia espectral, ICA y conectividad nodal.

using CairoMakie

"""
    plot_topomap(values, ch_names, ch_pos; kwargs...) -> Figure

Topomap 2D por IDW (Inverse Distance Weighting), estilo publicación EEG.

# Argumentos
- `values`: valor por canal (p.ej. columna de mixing matrix ICA o potencia).
- `ch_names`: nombres de canal (BIDS / electrodes.tsv).
- `ch_pos`: `Dict{String,Tuple{Float64,Float64}}` con claves en mayúsculas.

# Keywords
- `style`: `:auto` | `:ica` | `:power` — presets de colormap / colorbar / clims.
- `colormap`, `colorbar_label`, `clims`: si se pasan, anulan el preset.
- `show_sensors`, `show_channel_labels`, `label_fontsize`.
- `contours`: trazar contorno de nivel 0 (útil en ICA).
- `idw_power`: exponente IDW (default 3).
- `grid_res`: resolución del grid (default 200).
"""
function plot_topomap(
    values::Vector{Float64},
    ch_names::Vector{String},
    ch_pos::Dict{String,Tuple{Float64,Float64}};
    title::String = "",
    colormap = nothing,
    grid_res::Int = 200,
    clims::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    colorbar_label = nothing,
    style::Symbol = :auto,
    show_sensors::Bool = true,
    show_channel_labels::Bool = true,
    label_fontsize::Real = 8.5,
    contours::Bool = false,
    idw_power::Real = 3.0,
)::CairoMakie.Figure

    xs = [get(ch_pos, uppercase(ch), (NaN, NaN))[1] for ch in ch_names]
    ys = [get(ch_pos, uppercase(ch), (NaN, NaN))[2] for ch in ch_names]

    valid = findall(i -> !isnan(xs[i]) && !isnan(ys[i]) && isfinite(values[i]), eachindex(xs))
    isempty(valid) && error("plot_topomap: ningún canal con posición y valor válidos")
    xv = Float64[xs[i] for i in valid]
    yv = Float64[ys[i] for i in valid]
    vv = Float64[values[i] for i in valid]
    names_v = String[ch_names[i] for i in valid]

    # ── Presets ───────────────────────────────────────────────
    cmap, cblabel, cr = _topomap_resolve_style(style, vv, colormap, colorbar_label, clims)

    # ── IDW sobre disco unitario ───────────────────────────────
    xg = range(-1.0, 1.0; length=grid_res)
    yg = range(-1.0, 1.0; length=grid_res)
    Z  = fill(NaN, grid_res, grid_res)
    p  = Float64(idw_power)
    eps2 = 1e-8

    @inbounds for j in 1:grid_res, i in 1:grid_res
        xx, yy = xg[i], yg[j]
        r2 = xx^2 + yy^2
        r2 > 1.0 && continue
        # Extrapolación suave al rim: proyectar query hacia anillo interior
        qx, qy = xx, yy
        if r2 > 0.95^2
            r = sqrt(r2)
            s = 0.95 / r
            qx *= s; qy *= s
        end
        d2 = (qx .- xv).^2 .+ (qy .- yv).^2
        w  = 1.0 ./ (d2 .+ eps2).^(p / 2)
        Z[j, i] = sum(w .* vv) / sum(w)
    end

    # ── Figura estilo publicación ─────────────────────────────
    fig = CairoMakie.Figure(size=(480, 520), backgroundcolor=:white)
    ax  = CairoMakie.Axis(fig[1, 1];
        title  = title,
        aspect = DataAspect(),
        backgroundcolor = :white,
    )
    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidespines!(ax)
    CairoMakie.xlims!(ax, -1.28, 1.28)
    CairoMakie.ylims!(ax, -1.22, 1.32)

    hm = CairoMakie.heatmap!(ax, xg, yg, Z; colormap=cmap, colorrange=cr)

    # Contorno cero (ICA)
    if contours
        try
            CairoMakie.contour!(ax, xg, yg, Z; levels=[0.0], color=(:black, 0.55),
                                linewidth=1.0, linestyle=:dash)
        catch
            # contour puede fallar si Z es constante
        end
    end

    # Contorno de la cabeza
    θ = range(0, 2π; length=256)
    CairoMakie.lines!(ax, sin.(θ), cos.(θ); color=:black, linewidth=2.2)

    # Nariz (superior)
    θn = range(-π/6, π/6; length=40)
    CairoMakie.lines!(ax, 0.12 .* sin.(θn), 1.0 .+ 0.10 .* cos.(θn);
                      color=:black, linewidth=2.0)

    # Orejas (laterales)
    θe = range(-π/2, π/2; length=48)
    CairoMakie.lines!(ax, -1.0 .- 0.08 .* cos.(θe), 0.12 .* sin.(θe);
                      color=:black, linewidth=1.6)
    CairoMakie.lines!(ax,  1.0 .+ 0.08 .* cos.(θe), 0.12 .* sin.(θe);
                      color=:black, linewidth=1.6)

    # L / R
    CairoMakie.text!(ax, -1.18, 0.0; text="L", align=(:center, :center),
                     fontsize=13, color=:black)
    CairoMakie.text!(ax,  1.18, 0.0; text="R", align=(:center, :center),
                     fontsize=13, color=:black)

    # Electrodos + etiquetas
    inside = findall(i -> xv[i]^2 + yv[i]^2 <= 1.05, eachindex(xv))
    if show_sensors && !isempty(inside)
        CairoMakie.scatter!(ax, xv[inside], yv[inside];
                            color=:black, markersize=5,
                            strokecolor=:white, strokewidth=1)
    end
    if show_channel_labels && !isempty(inside)
        for i in inside
            px, py = xv[i], yv[i]
            r = hypot(px, py)
            # Offset radial hacia fuera (o hacia +y si casi en el origen)
            if r < 1e-3
                ox, oy = 0.0, 0.05
            else
                δ = 0.055
                ox = px * (1 + δ / r) - px
                oy = py * (1 + δ / r) - py
            end
            CairoMakie.text!(ax, px + ox, py + oy;
                             text = names_v[i],
                             align = (:center, :center),
                             fontsize = Float64(label_fontsize),
                             color = (:black, 0.85))
        end
    end

    CairoMakie.Colorbar(fig[1, 2], hm; label=cblabel, width=14, labelsize=11, ticklabelsize=10)
    CairoMakie.colgap!(fig.layout, 8)
    return fig
end

"""Resolve colormap / colorbar label / colorrange from style + overrides."""
function _topomap_resolve_style(
    style::Symbol,
    vv::Vector{Float64},
    colormap,
    colorbar_label,
    clims::Union{Nothing,Tuple{Float64,Float64}},
)
    if style === :ica
        cmap = colormap === nothing ? CairoMakie.Reverse(:RdBu) : colormap
        cbl  = colorbar_label === nothing ? "ICA mixing weights (a.u.)" : String(colorbar_label)
        if clims === nothing
            M = maximum(abs, vv)
            M == 0.0 && (M = 1.0)
            cr = (-M, M)
        else
            cr = clims
        end
    elseif style === :power
        cmap = colormap === nothing ? :viridis : colormap
        cbl  = colorbar_label === nothing ? "μV²" : String(colorbar_label)
        cr   = clims === nothing ? (0.0, maximum(vv)) : clims
        cr[2] == cr[1] && (cr = (cr[1], cr[1] + 1.0))
    else
        # :auto — comportamiento legacy-friendly
        cmap = colormap === nothing ? :viridis : colormap
        cbl  = colorbar_label === nothing ? "μV²" : String(colorbar_label)
        cr   = clims === nothing ? (minimum(vv), maximum(vv)) : clims
        if cr[1] == cr[2]
            pad = max(abs(cr[1]) * 0.05, 1e-6)
            cr = (cr[1] - pad, cr[2] + pad)
        end
    end
    return cmap, cbl, cr
end
