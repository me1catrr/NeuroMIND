# NeuroMIND/src/spectral/PowerSpectrum.jl
# Análisis espectral FFT estilo BrainVision Analyzer.
# Lógica científica portada y refactorizada desde EEG_Julia/src/Spectral/FFT.jl.

"""
    compute_psd(epochs::EpochSet, cfg::PipelineConfig) -> SpectralResult

Calcula el PSD promedio (sobre epochs) usando FFT con ventana Hamming-taper
y corrección de varianza (estilo BrainVision Analyzer).
"""
function compute_psd(epochs::EpochSet, cfg::PipelineConfig)::SpectralResult
    sp   = cfg.spectral
    pct  = get(sp, "window_pct", 10.0) / 100.0
    fs   = epochs.meta.fs

    n_ch, n_samp, n_ep = size(epochs.data)

    # nfft must be >= n_samp (zero-padding allowed, truncation not)
    nfft = max(get(sp, "nfft", 1024), n_samp)

    # 1. Ventana Hamming-taper (estilo BrainVision)
    win, mw2 = _hamming_taper(n_samp, pct)

    # 2. rFFT + potencia por epoch
    n_bins = nfft ÷ 2 + 1
    P_all  = zeros(Float64, n_ch, n_bins, n_ep)

    @inbounds for ep in 1:n_ep, ch in 1:n_ch
        x       = epochs.data[ch, :, ep] .- mean(@view epochs.data[ch, :, ep])
        xw      = x .* win
        xpad    = vcat(xw, zeros(nfft - n_samp))
        X       = rfft(xpad)
        Pseg    = abs2.(X) ./ mw2
        Pseg[2:end-1] .*= 2.0       # folding: Use Full Spectrum
        Pseg  ./= nfft^2
        P_all[ch, :, ep] = Pseg
    end

    # 3. Promedio sobre epochs
    psd_mean = dropdims(mean(P_all, dims=3), dims=3)  # (ch × bins)

    # 4. Eje de frecuencias
    freqs = collect(range(0.0, fs/2; length=n_bins))

    # 5. Potencia por banda
    bp = _band_power_from_psd(psd_mean, freqs, cfg.bands)

    params = Dict{String,Any}(
        "nfft"       => nfft,
        "window"     => "hamming_taper",
        "window_pct" => pct * 100,
        "fs"         => fs
    )

    return SpectralResult(epochs.meta, psd_mean, freqs, bp, n_ep, params)
end

"""
    band_power(result::SpectralResult, band::String) -> Vector{Float64}

Devuelve la potencia media de la banda `band` para cada canal.
"""
function band_power(result::SpectralResult, band::String)::Vector{Float64}
    haskey(result.band_power, band) || error("Banda no encontrada: $band")
    return result.band_power[band]
end

# ─── Helpers privados ─────────────────────────────────────────

function _hamming_taper(n::Int, pct::Real)
    win    = ones(Float64, n)
    n_tap  = round(Int, n * pct / 2)
    α, β   = 0.54, 0.46

    if n_tap >= 2
        t_left  = 0:(n_tap - 1)
        t_right = (n - n_tap):(n - 1)
        win[1:n_tap]          .= α .- β .* cos.(2π .* t_left  ./ (n * pct))
        win[(n-n_tap+1):end]  .= α .- β .* cos.(2π .* (1.0 .- t_right ./ n) ./ pct)
    end

    mw2 = mean(win .^ 2)
    return win, mw2
end

function _band_power_from_psd(
    psd::Matrix{Float64},
    freqs::Vector{Float64},
    bands::Dict{String,Tuple{Float64,Float64}}
)::Dict{String,Vector{Float64}}

    n_ch = size(psd, 1)
    bp   = Dict{String,Vector{Float64}}()

    for (name, (flo, fhi)) in bands
        idx = findall(f -> f >= flo && f < fhi, freqs)
        if isempty(idx)
            bp[name] = fill(NaN, n_ch)
        else
            bp[name] = vec(mean(psd[:, idx], dims=2))
        end
    end
    return bp
end
