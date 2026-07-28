# NeuroMIND/src/visualization/Spectra.jl
# Visualización de espectros de potencia.

using CairoMakie

"""
    plot_spectrum_grid(result::SpectralResult; xmax=50.0, cols=6) -> Figure

Grid de espectros para todos los canales (estilo BrainVision multichannel).
"""
function plot_spectrum_grid(
    result::SpectralResult;
    xmax::Real = 50.0,
    cols::Int = 6
)::CairoMakie.Figure

    n_ch  = result.meta.n_channels
    rows  = ceil(Int, n_ch / cols)
    freqs = result.freqs
    mask  = freqs .<= xmax

    fig = CairoMakie.Figure(size=(200 * cols, 160 * rows))

    for ch in 1:n_ch
        r = ((ch - 1) ÷ cols) + 1
        c = ((ch - 1) % cols) + 1
        ax = CairoMakie.Axis(fig[r, c];
            title = result.meta.channel_names[ch],
            titlesize = 9,
            xticklabelsize = 7,
            yticklabelsize = 7,
        )
        psd_ch = result.psd[ch, mask]
        CairoMakie.lines!(ax, freqs[mask], psd_ch; color=:steelblue, linewidth=1.5)
    end
    return fig
end
