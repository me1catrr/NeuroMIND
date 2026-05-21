# NeuroMIND/src/connectivity/wPLI.jl
# Conectividad funcional mediante wPLI (Weighted Phase Lag Index).
# Portado y refactorizado desde EEG_Julia/src/Connectivity/wPLI.jl.
# Método: agregación "across segments" estilo BrainVision Analyzer.

"""
    compute_wpli(epochs::EpochSet, cfg::PipelineConfig) -> ConnectivityMatrix

Calcula matrices wPLI para cada banda de frecuencia definida en la configuración.
Usa Butterworth bandpass (orden 8) + Hilbert + agregación across segments.
"""
function compute_wpli(epochs::EpochSet, cfg::PipelineConfig)::ConnectivityMatrix
    conn    = cfg.connectivity
    order   = get(conn, "filter_order", 8)
    bands   = cfg.bands

    n_ch, n_samp, n_seg = size(epochs.data)
    fs      = epochs.meta.fs
    ch_names = epochs.meta.channel_names

    matrices = Dict{String,Matrix{Float64}}()

    for (band_name, (f1, f2)) in bands
        W = _compute_wpli_band(epochs.data, fs, f1, f2, order, n_ch, n_samp, n_seg)
        matrices[band_name] = W
    end

    params = Dict{String,Any}(
        "filter_order" => order,
        "method"       => "across_segments",
        "fs"           => fs
    )

    return ConnectivityMatrix(
        epochs.meta, "wpli", matrices, ch_names,
        get(conn, "use_csd", false) ? "CSD" : "sensor",
        n_seg,
        params
    )
end

# ─── Kernel privado ───────────────────────────────────────────

function _compute_wpli_band(
    data::Array{Float64,3},
    fs::Real,
    f1::Real,
    f2::Real,
    order::Int,
    n_ch::Int,
    n_samp::Int,
    n_seg::Int
)::Matrix{Float64}

    # Precomputar señales analíticas: z[ch] = (samples × segments)
    z = Vector{Matrix{ComplexF64}}(undef, n_ch)

    @inbounds for c in 1:n_ch
        Zc = Matrix{ComplexF64}(undef, n_samp, n_seg)
        for s in 1:n_seg
            x  = @view data[c, :, s]
            xf = _bandpass_filtfilt(x, fs, f1, f2; order)
            Zc[:, s] = _analytic_signal(xf)
        end
        z[c] = Zc
    end

    # Matriz wPLI simétrica con diagonal = 0
    W = zeros(Float64, n_ch, n_ch)
    @inbounds for i in 1:n_ch, j in (i+1):n_ch
        w = _wpli_agg(z[i], z[j])
        W[i, j] = w
        W[j, i] = w
    end
    return W
end

function _bandpass_filtfilt(x::AbstractVector, fs::Real, f1::Real, f2::Real; order::Int=8)
    nyq  = fs / 2.0
    bp   = digitalfilter(Bandpass(f1/nyq, f2/nyq), Butterworth(order))
    return filtfilt(bp, Float64.(x))
end

function _analytic_signal(x::AbstractVector{<:Real})::Vector{ComplexF64}
    N = length(x)
    X = fft(Float64.(x))
    h = zeros(Float64, N)
    if iseven(N)
        h[1] = 1.0; h[N÷2+1] = 1.0; h[2:N÷2] .= 2.0
    else
        h[1] = 1.0; h[2:(N+1)÷2] .= 2.0
    end
    return ifft(X .* h)
end

function _wpli_agg(zi::Matrix{ComplexF64}, zj::Matrix{ComplexF64})::Float64
    imv = imag.(zi .* conj.(zj))
    num = abs(sum(imv))
    den = sum(abs.(imv)) + eps()
    return num / den
end
