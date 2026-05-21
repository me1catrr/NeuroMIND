# NeuroMIND/src/visualization/Spectra.jl
# Visualización de espectros de potencia.

using CairoMakie

const BAND_COLORS = Dict(
    "DELTA"     => :orange,
    "THETA"     => :yellow,
    "ALPHA"     => :green,
    "BETA_LOW"  => :steelblue,
    "BETA_MID"  => :royalblue,
    "BETA_HIGH" => :navy,
    "GAMMA"     => :black,
)

"""
    plot_spectrum(result::SpectralResult, channel::String;
                 xmax=50.0, highlight_band="ALPHA") -> Figure

Espectro de potencia de un canal con relleno de bandas (estilo BrainVision).
"""
function plot_spectrum(
    result::SpectralResult,
    channel::String;
    xmax::Real = 50.0,
    highlight_band::String = "ALPHA",
    bands::Dict{String,Tuple{Float64,Float64}} = Dict{String,Tuple{Float64,Float64}}()
)::CairoMakie.Figure

    ch_idx = findfirst(==(channel), result.meta.channel_names)
    ch_idx === nothing && error("Canal '$channel' no encontrado")

    freqs = result.freqs
    psd   = result.psd[ch_idx, :]
    mask  = freqs .<= xmax

    fig = CairoMakie.Figure(size=(900, 500))
    ax  = CairoMakie.Axis(fig[1, 1];
        title  = "PSD — $(channel)",
        xlabel = "Frecuencia (Hz)",
        ylabel = "Potencia (μV²)",
        limits = (0.0, xmax, 0.0, maximum(psd[mask]) * 1.15),
    )

    # Relleno por bandas
    for (bname, (f1, f2)) in bands
        bidx = findall(f -> f >= f1 && f < f2 && f <= xmax, freqs)
        isempty(bidx) && continue
        col = get(BAND_COLORS, bname, :gray)
        CairoMakie.band!(ax, freqs[bidx], zeros(length(bidx)), psd[bidx];
                          color=(col, 0.4))
    end

    # Línea del espectro
    CairoMakie.lines!(ax, freqs[mask], psd[mask]; color=:black, linewidth=2)

    # Resaltar banda seleccionada
    if haskey(bands, highlight_band)
        f1, f2 = bands[highlight_band]
        CairoMakie.vspan!(ax, f1, min(f2, xmax); color=(:mediumpurple, 0.2))
    end

    return fig
end

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
