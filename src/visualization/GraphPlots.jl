# NeuroMIND/src/visualization/GraphPlots.jl
# Guardado de figuras del pipeline por sujeto.

using CairoMakie

"""
    save_figure(fig, path; dpi=300)

Guarda una figura CairoMakie con resolución configurada.
"""
function save_figure(fig::CairoMakie.Figure, path::String; dpi::Int=300)
    px_per_unit = dpi / 96
    CairoMakie.save(path, fig; px_per_unit)
end
