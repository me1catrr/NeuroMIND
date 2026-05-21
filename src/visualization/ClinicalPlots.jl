# NeuroMIND/src/visualization/ClinicalPlots.jl
# Correlaciones entre métricas EEG y variables clínicas:
# EDSS, fatiga, cognición, carga lesional.

using CairoMakie

"""
    plot_clinical_correlation(eeg_vals, clin_vals, clin_name, eeg_name, band;
                              sr=nothing, save_path=nothing) -> Figure

Scatterplot de correlación entre una métrica EEG y una variable clínica.
Opcionalmente superpone la recta de regresión y muestra el estadístico de Spearman.
"""
function plot_clinical_correlation(
    eeg_vals::Vector{Float64},
    clin_vals::Vector{Float64},
    clin_name::String,
    eeg_name::String,
    band::String;
    sr::Union{Nothing,StatResult} = nothing,
    group::String = "MS",
    save_path::Union{Nothing,String} = nothing
)::CairoMakie.Figure

    valid = findall(i -> isfinite(eeg_vals[i]) && isfinite(clin_vals[i]), eachindex(eeg_vals))
    xv = clin_vals[valid]
    yv = eeg_vals[valid]

    fig = CairoMakie.Figure(size=(520, 440))
    ax  = CairoMakie.Axis(fig[1, 1];
        title  = "$(eeg_name) vs $(clin_name) [$(band)] — Grupo $(group)",
        xlabel = clin_name,
        ylabel = "$(eeg_name) (wPLI)",
        titlesize = 11,
    )

    CairoMakie.scatter!(ax, xv, yv; color=:steelblue, markersize=8,
                        strokewidth=0.8, strokecolor=:navy)

    # Recta de regresión lineal
    if length(xv) >= 3
        β1 = cov(xv, yv) / (var(xv) + 1e-12)
        β0 = mean(yv) - β1 * mean(xv)
        xl = [minimum(xv), maximum(xv)]
        CairoMakie.lines!(ax, xl, β0 .+ β1 .* xl; color=:tomato, linewidth=2)
    end

    # Texto con estadístico
    if sr !== nothing
        ann = "ρ = $(round(sr.statistic, digits=3))\np = $(round(sr.p_value, digits=4))"
        CairoMakie.text!(ax, ann; position=(minimum(xv), maximum(yv)),
                         align=(:left, :top), fontsize=10,
                         color=sr.significant ? :darkred : :gray40)
    end

    save_path !== nothing && CairoMakie.save(save_path, fig; px_per_unit=2)
    return fig
end

"""
    plot_clinical_correlation_grid(subjects, metric_fn, clin_field, cfg;
                                   band, condition, save_path) -> Figure

Panel de correlaciones para todas las variables clínicas disponibles.
"""
function plot_clinical_correlation_grid(
    subjects::Vector{Subject},
    metric_fn::Function,      # fn(subj) -> Float64
    clin_fields::Vector{String},
    cfg::PipelineConfig;
    band::String      = "ALPHA",
    condition::String = "EC",
    save_path::Union{Nothing,String} = nothing
)::CairoMakie.Figure

    ms_subjs = filter(s -> s.group == "MS", subjects)
    eeg_vals = [metric_fn(s) for s in ms_subjs]

    n_clin = length(clin_fields)
    cols   = min(3, n_clin)
    rows   = ceil(Int, n_clin / cols)

    fig = CairoMakie.Figure(size=(480 * cols, 400 * rows))
    CairoMakie.Label(fig[0, :],
        "Correlaciones clínicas — wPLI $(band) / $(condition)";
        fontsize=13, font=:bold)

    for (k, field) in enumerate(clin_fields)
        r = ((k - 1) ÷ cols) + 1
        c = ((k - 1) % cols) + 1

        clin_vals = [_get_clinical(s.clinical, field) for s in ms_subjs]
        valid     = findall(i -> isfinite(eeg_vals[i]) && isfinite(clin_vals[i]), eachindex(eeg_vals))
        length(valid) < 3 && continue

        xv = clin_vals[valid]
        yv = eeg_vals[valid]
        sr = spearman_correlation(xv, yv)

        ax = CairoMakie.Axis(fig[r, c];
            title  = field,
            xlabel = field,
            ylabel = "wPLI $(band)",
            titlesize = 10,
        )
        CairoMakie.scatter!(ax, xv, yv; color=(:steelblue, 0.7), markersize=7)
        if length(xv) >= 3
            β1 = cov(xv, yv) / (var(xv) + 1e-12)
            β0 = mean(yv) - β1 * mean(xv)
            xl = [minimum(xv), maximum(xv)]
            CairoMakie.lines!(ax, xl, β0 .+ β1 .* xl; color=:tomato, linewidth=1.8)
        end
        ann = "ρ=$(round(sr.statistic,digits=2)) p=$(round(sr.p_value,digits=3))"
        CairoMakie.text!(ax, ann; position=(minimum(xv), maximum(yv)),
                         align=(:left,:top), fontsize=8,
                         color=sr.significant ? :darkred : :gray50)
    end

    save_path !== nothing && CairoMakie.save(save_path, fig; px_per_unit=2)
    return fig
end

function _get_clinical(c::ClinicalData, field::String)::Float64
    v = getfield(c, Symbol(field))
    ismissing(v) ? NaN : Float64(v)
end
