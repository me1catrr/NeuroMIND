# NeuroMIND/src/connectivity/wPLI.jl
# Conectividad funcional wPLI/dwPLI — arquitectura multi-método configurable.
#
# MÉTODOS (configurar via wpli_method en [connectivity] del .toml):
#
#   "hilbert"      — Filtrado Butterworth pasa-banda + señal analítica Hilbert.
#                    Método original. Rápido, válido para épocas suficientemente largas.
#                    Sensible a artefactos de borde de filtro.
#
#   "fourier_csd"  — Espectro cruzado FFT por época con ventana Hann/Hamming.
#                    Definición espectral directa. Conceptualmente más limpio.
#                    Promedia Im(Sxy) across épocas y frecuencias dentro de banda.
#
#   "multitaper"   — DPSS multitaper + espectro cruzado promediado entre tapers.
#                    Mejor compromiso varianza-sesgo. Estándar en MNE spectral_connectivity.
#                    Usa DSP.dpss (ya disponible en el proyecto).
#
# ESTIMADORES (use_dwpli en [connectivity]):
#
#   wPLI  (Vinck 2011):  |E[Im(Sxy)]| / E[|Im(Sxy)|]          rango [0, 1]
#   dwPLI (Vinck 2011):  (E[Im]²−E[Im²]) / (E[|Im|]²−E[Im²])  rango [−1, 1], no sesgado
#     → Recomendado para comparaciones grupales (n_epochs variable entre sujetos)
#
# Referencias:
#   Vinck et al. (2011) NeuroImage — wPLI y dwPLI
#   Percival & Walden (1993) — DPSS / multitaper
#   Thomson (1982) Proc. IEEE — multitaper spectral estimation

# ═══════════════════════════════════════════════════════════════
# TIPOS DE ESTIMADOR
# ═══════════════════════════════════════════════════════════════

"""Tipo abstracto base. Todos los estimadores comparten min_cycles y use_dwpli."""
abstract type AbstractWPLIEstimator end

"""
Estimador basado en transformada de Hilbert (método original NeuroMIND).

Flujo: bandpass Butterworth (order) → señal analítica → Im(z_i·conj(z_j)) → wPLI/dwPLI.
"""
struct HilbertEstimator <: AbstractWPLIEstimator
    filter_order       :: Int      # orden del filtro Butterworth (default: 8)
    use_dwpli          :: Bool
    min_cycles         :: Float64  # ciclos mínimos por época por banda
    exclude_unreliable :: Bool     # si true, saltar banda con < min_cycles ciclos
end

"""
Estimador basado en espectro cruzado FFT (Fourier CSD).

Flujo: ventana Hann/Hamming × época → FFT → Im(Sxy[f]) para f ∈ banda →
       wPLI = |Σ Im(Sxy)| / Σ|Im(Sxy)|  acumulado across (épocas × frecuencias).
"""
struct FourierCSDEstimator <: AbstractWPLIEstimator
    nfft               :: Int      # longitud FFT (0 = usar n_samp de la época)
    window             :: Symbol   # :hann | :hamming | :rect
    use_dwpli          :: Bool
    min_cycles         :: Float64
    exclude_unreliable :: Bool
end

"""
Estimador multitaper DPSS (equivalente a MNE spectral_connectivity_epochs con method="multitaper").

Flujo: K tapers DPSS × época → FFT por taper → promedio Im(Sxy) over tapers →
       wPLI = |Σ Im(Sxy_avg)| / Σ|Im(Sxy_avg)|  across (épocas × frecuencias de banda).

nw       : time-bandwidth product (típico: 2.0–5.0; default 4.0)
n_tapers : 0 = automático (K = floor(2*nw) − 1)
low_bias : descartar tapers con concentración espectral λ < 0.9
"""
struct MultitaperEstimator <: AbstractWPLIEstimator
    nw                 :: Float64  # time-bandwidth product
    n_tapers           :: Int      # 0 = automático
    use_dwpli          :: Bool
    min_cycles         :: Float64
    exclude_unreliable :: Bool
    low_bias           :: Bool
end

# ─── Constructor desde Dict de configuración ────────────────

"""
    _build_estimator(conn_cfg) -> AbstractWPLIEstimator

Construye el estimador adecuado leyendo los parámetros de [connectivity] del TOML.
"""
function _build_estimator(conn_cfg::Dict{String,Any})::AbstractWPLIEstimator
    method   = String(get(conn_cfg, "wpli_method", "hilbert"))
    use_dw   = Bool(get(conn_cfg, "use_dwpli", false))
    min_cyc  = Float64(get(conn_cfg, "min_cycles_for_wpli", 4.0))
    excl_unr = Bool(get(conn_cfg, "exclude_unreliable_bands", false))

    if method == "fourier_csd"
        csd_cfg = get(conn_cfg, "fourier_csd", Dict{String,Any}())
        win_str = String(get(csd_cfg, "window", get(conn_cfg, "window", "hann")))
        nfft    = Int(get(csd_cfg, "nfft", 0))
        return FourierCSDEstimator(nfft, Symbol(win_str), use_dw, min_cyc, excl_unr)

    elseif method == "multitaper"
        mt_cfg = get(conn_cfg, "multitaper", Dict{String,Any}())
        nw     = Float64(get(mt_cfg, "nw", 4.0))
        n_tap  = Int(get(mt_cfg, "n_tapers", 0))
        low_b  = Bool(get(mt_cfg, "low_bias", true))
        return MultitaperEstimator(nw, n_tap, use_dw, min_cyc, excl_unr, low_b)

    else  # "hilbert" (default + backward compat)
        order = Int(get(conn_cfg, "filter_order", 8))
        return HilbertEstimator(order, use_dw, min_cyc, excl_unr)
    end
end

# ═══════════════════════════════════════════════════════════════
# API PÚBLICA
# ═══════════════════════════════════════════════════════════════

"""
    compute_wpli(epochs, cfg) -> ConnectivityMatrix

Calcula matrices wPLI (o dwPLI si `use_dwpli=true`) para todas las bandas.

El método se selecciona desde `cfg.connectivity["wpli_method"]`:
  - `"hilbert"`      : Butterworth bandpass + señal analítica Hilbert (default)
  - `"fourier_csd"`  : Espectro cruzado FFT con ventana configurable
  - `"multitaper"`   : DPSS multitaper (standard MNE-compatible)

Emite `@warn` para bandas con < `min_cycles_for_wpli` ciclos/época.
Si `exclude_unreliable_bands = true`, esas bandas se omiten del resultado.
"""
function compute_wpli(epochs::EpochSet, cfg::PipelineConfig)::ConnectivityMatrix
    conn_cfg  = cfg.connectivity
    bands     = cfg.bands
    estimator = _build_estimator(conn_cfg)

    n_ch, n_samp, n_seg = size(epochs.data)
    fs       = epochs.meta.fs
    ch_names = epochs.meta.channel_names
    epoch_s  = n_samp / fs

    # ─── Validación de fiabilidad por banda ──────────────────────
    skip_bands = Set{String}()
    for (bn, (f1, _)) in bands
        n_cyc = f1 * epoch_s
        if n_cyc < estimator.min_cycles
            msg = "wPLI [$bn]: $(round(n_cyc, digits=2)) ciclos/época " *
                  "(f_low=$(f1) Hz × $(round(epoch_s, digits=2))s) < " *
                  "mínimo $(estimator.min_cycles). Estimación poco fiable."
            if estimator.exclude_unreliable
                @warn msg * " → BANDA EXCLUIDA (exclude_unreliable_bands=true)"
                push!(skip_bands, bn)
            else
                @warn msg
            end
        end
    end

    # ─── Cálculo por banda ────────────────────────────────────────
    matrices = Dict{String,Matrix{Float64}}()
    for (band_name, (f1, f2)) in bands
        band_name ∈ skip_bands && continue
        matrices[band_name] = _compute_band_matrix(
            estimator, epochs.data, fs, f1, f2, n_ch, n_samp, n_seg
        )
    end

    # ─── Metadata de la corrida ────────────────────────────────────
    method_tag     = _method_tag(estimator)
    estimator_name = estimator.use_dwpli ? "dwpli" : "wpli"

    params = Dict{String,Any}(
        "wpli_method"    => method_tag,
        "estimator"      => estimator_name,
        "use_dwpli"      => estimator.use_dwpli,
        "fs"             => fs,
        "method"         => "across_segments",
        "skipped_bands"  => collect(skip_bands),
    )
    _add_method_params!(params, estimator)

    return ConnectivityMatrix(
        epochs.meta, estimator_name, matrices, ch_names,
        Bool(get(conn_cfg, "use_csd", false)) ? "CSD" : "sensor",
        n_seg, params
    )
end

_method_tag(::HilbertEstimator)    = "hilbert"
_method_tag(::FourierCSDEstimator)  = "fourier_csd"
_method_tag(::MultitaperEstimator)  = "multitaper"

_add_method_params!(p, e::HilbertEstimator)   = (p["filter_order"] = e.filter_order)
_add_method_params!(p, e::FourierCSDEstimator) = (p["window"] = string(e.window); p["nfft"] = e.nfft)
_add_method_params!(p, e::MultitaperEstimator) = (p["nw"] = e.nw; p["n_tapers"] = e.n_tapers; p["low_bias"] = e.low_bias)

# ═══════════════════════════════════════════════════════════════
# MÉTODO 1: HILBERT
# ═══════════════════════════════════════════════════════════════

function _compute_band_matrix(
    est::HilbertEstimator,
    data::Array{Float64,3},
    fs::Real, f1::Real, f2::Real,
    n_ch::Int, n_samp::Int, n_seg::Int
)::Matrix{Float64}

    # Señal analítica por canal × época
    z = Vector{Matrix{ComplexF64}}(undef, n_ch)
    @inbounds for c in 1:n_ch
        Zc = Matrix{ComplexF64}(undef, n_samp, n_seg)
        for s in 1:n_seg
            x  = @view data[c, :, s]
            xf = _bandpass_filtfilt(x, fs, f1, f2; order=est.filter_order)
            Zc[:, s] = _analytic_signal(xf)
        end
        z[c] = Zc
    end

    # wPLI/dwPLI por par de canales
    W = zeros(Float64, n_ch, n_ch)
    @inbounds for i in 1:n_ch, j in (i+1):n_ch
        w = est.use_dwpli ? _dwpli_agg(z[i], z[j]) : _wpli_agg(z[i], z[j])
        W[i, j] = w; W[j, i] = w
    end
    return W
end

# ═══════════════════════════════════════════════════════════════
# MÉTODO 2: FOURIER CSD
# ═══════════════════════════════════════════════════════════════

"""
wPLI via espectro cruzado FFT.

Para cada época s y frecuencia f ∈ [f1, f2]:
  Im(Sxy[f]) = Im( FFT(w·xi)[f] · conj(FFT(w·xj)[f]) )

Se acumulan los 3 momentos necesarios across (épocas × frecuencias):
  Σ Im,  Σ|Im|,  ΣIm²   →   wPLI o dwPLI al final.
"""
function _compute_band_matrix(
    est::FourierCSDEstimator,
    data::Array{Float64,3},
    fs::Real, f1::Real, f2::Real,
    n_ch::Int, n_samp::Int, n_seg::Int
)::Matrix{Float64}

    nfft     = est.nfft > 0 ? est.nfft : n_samp
    win      = _make_window(est.window, nfft)
    freqs    = _fft_freqs(nfft, fs)
    band_idx = findall(f -> f1 ≤ f ≤ f2, freqs)
    isempty(band_idx) && return zeros(Float64, n_ch, n_ch)
    n_freq   = length(band_idx)

    # Pre-computar FFTs en banda: (n_ch, n_freq, n_seg)
    Xf = Array{ComplexF64,3}(undef, n_ch, n_freq, n_seg)
    @inbounds for c in 1:n_ch, s in 1:n_seg
        x = data[c, 1:min(nfft, n_samp), s]
        xpad = length(x) < nfft ? vcat(x, zeros(nfft - length(x))) : x
        Xfull = fft(win .* xpad)
        Xf[c, :, s] = Xfull[band_idx]
    end

    # Acumular momentos across (épocas × frecuencias de banda)
    W = zeros(Float64, n_ch, n_ch)
    @inbounds for i in 1:n_ch, j in (i+1):n_ch
        im_sum = 0.0; abs_sum = 0.0; sq_sum = 0.0
        for s in 1:n_seg, k in 1:n_freq
            v = imag(Xf[i,k,s] * conj(Xf[j,k,s]))
            im_sum  += v
            abs_sum += abs(v)
            sq_sum  += v * v
        end
        w = est.use_dwpli ? _dwpli_moments(im_sum, abs_sum, sq_sum) :
                            _wpli_moments(im_sum, abs_sum)
        W[i, j] = w; W[j, i] = w
    end
    return W
end

# ═══════════════════════════════════════════════════════════════
# MÉTODO 3: MULTITAPER
# ═══════════════════════════════════════════════════════════════

"""
wPLI via estimación DPSS multitaper.

Los K tapers DPSS (Discrete Prolate Spheroidal Sequences) minimizan la
fuga espectral manteniendo control óptimo del compromiso varianza-sesgo.

Para cada taper k y época s:
  Xk[f] = FFT(taper_k · x_s)

Espectro cruzado promediado sobre tapers:
  Sxy_avg[f] = (1/K) · Σ_k  Xk_i[f] · conj(Xk_j[f])

wPLI calculado desde Im(Sxy_avg) acumulado across (épocas × frecuencias de banda).

Parámetros DSP.dpss:
  N   = n_samp (longitud de época)
  nw  = time-bandwidth product (p.ej. 4.0 → K = 7 tapers óptimos)
  K   = floor(2*nw) - 1  si n_tapers = 0

Concentración espectral (λ): fracción de energía en [-W, W].
  low_bias=true: usar solo tapers con λ ≥ 0.9 (elimina últimos tapers).
"""
function _compute_band_matrix(
    est::MultitaperEstimator,
    data::Array{Float64,3},
    fs::Real, f1::Real, f2::Real,
    n_ch::Int, n_samp::Int, n_seg::Int
)::Matrix{Float64}

    tapers = _get_dpss_tapers(n_samp, est.nw, est.n_tapers; low_bias=est.low_bias)
    K      = size(tapers, 2)

    freqs    = _fft_freqs(n_samp, fs)
    band_idx = findall(f -> f1 ≤ f ≤ f2, freqs)
    isempty(band_idx) && return zeros(Float64, n_ch, n_ch)
    n_freq   = length(band_idx)

    # Pre-computar FFTs por canal, época y taper: (n_ch, n_freq, n_seg, K)
    Xft = Array{ComplexF64,4}(undef, n_ch, n_freq, n_seg, K)
    @inbounds for c in 1:n_ch, s in 1:n_seg, k in 1:K
        x = @view data[c, :, s]
        Xfull = fft(tapers[:, k] .* x)
        Xft[c, :, s, k] = Xfull[band_idx]
    end

    # Acumular momentos: promedio de Im(Sxy) sobre tapers, luego across (épocas × freq)
    W = zeros(Float64, n_ch, n_ch)
    @inbounds for i in 1:n_ch, j in (i+1):n_ch
        im_sum = 0.0; abs_sum = 0.0; sq_sum = 0.0
        for s in 1:n_seg, f_idx in 1:n_freq
            # Im del CSD promediado sobre K tapers
            v = 0.0
            for k in 1:K
                v += imag(Xft[i,f_idx,s,k] * conj(Xft[j,f_idx,s,k]))
            end
            v /= K
            im_sum  += v
            abs_sum += abs(v)
            sq_sum  += v * v
        end
        w = est.use_dwpli ? _dwpli_moments(im_sum, abs_sum, sq_sum) :
                            _wpli_moments(im_sum, abs_sum)
        W[i, j] = w; W[j, i] = w
    end
    return W
end

"""
    _get_dpss_tapers(N, nw, K; low_bias) -> Matrix{Float64}  (N × K_actual)

Wrapper sobre DSP.dpss y DSP.dpsseig.
Filtra tapers con concentración < 0.9 si low_bias = true.
"""
function _get_dpss_tapers(N::Int, nw::Float64, K::Int; low_bias::Bool=true)
    K_req  = K > 0 ? K : max(1, floor(Int, 2*nw) - 1)
    K_req  = min(K_req, N)
    tapers = DSP.dpss(N, nw, K_req)           # (N, K_req) — ya normalizados
    λ      = DSP.dpsseig(tapers, nw)           # concentraciones ∈ [0, 1]

    if low_bias
        good = findall(λ .≥ 0.9)
        isempty(good) && (good = [length(λ)])  # garantizar al menos 1 taper
        tapers = tapers[:, good]
    end
    return tapers
end

# ═══════════════════════════════════════════════════════════════
# AGREGADORES DE MOMENTOS (shared por FourierCSD y Multitaper)
# ═══════════════════════════════════════════════════════════════

@inline function _wpli_moments(im_sum::Float64, abs_sum::Float64)::Float64
    abs(im_sum) / (abs_sum + eps())
end

@inline function _dwpli_moments(im_sum::Float64, abs_sum::Float64, sq_sum::Float64)::Float64
    num = im_sum * im_sum - sq_sum
    den = abs_sum * abs_sum - sq_sum + eps()
    num / den
end

# ═══════════════════════════════════════════════════════════════
# AGREGADORES MATRICIALES (usados por Hilbert + Surrogates.jl)
# ═══════════════════════════════════════════════════════════════

"""
    _wpli_agg(zi, zj) -> Float64

wPLI clásico (Vinck et al. 2011). Rango [0, 1].
zi, zj: matrices ComplexF64 de shape (n_samp, n_seg).
"""
function _wpli_agg(zi::Matrix{ComplexF64}, zj::Matrix{ComplexF64})::Float64
    imv = imag.(zi .* conj.(zj))
    _wpli_moments(sum(imv), sum(abs.(imv)))
end

"""
    _dwpli_agg(zi, zj) -> Float64

dwPLI debiased (Vinck et al. 2011). Estimador no sesgado de (wPLI)².
Rango [-1, 1]. Valores negativos indican ausencia de acoplamiento.
"""
function _dwpli_agg(zi::Matrix{ComplexF64}, zj::Matrix{ComplexF64})::Float64
    imv = imag.(zi .* conj.(zj))
    _dwpli_moments(sum(imv), sum(abs.(imv)), sum(imv .* imv))
end

# ═══════════════════════════════════════════════════════════════
# HELPERS DE SEÑAL
# ═══════════════════════════════════════════════════════════════

"""Filtro pasa-banda Butterworth con filtfilt (fase cero)."""
function _bandpass_filtfilt(x::AbstractVector, fs::Real, f1::Real, f2::Real; order::Int=8)
    nyq = fs / 2.0
    bp  = digitalfilter(Bandpass(f1/nyq, f2/nyq), Butterworth(order))
    return filtfilt(bp, Float64.(x))
end

"""
Señal analítica vía transformada de Hilbert (FFT one-sided).
Implementación directa sin dependencias externas.
"""
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

"""Vector de frecuencias FFT: [0, fs/N, 2fs/N, ..., (N-1)fs/N]"""
function _fft_freqs(N::Int, fs::Real)::Vector{Float64}
    [k * fs / N for k in 0:(N-1)]
end

"""Ventana de análisis: :hann, :hamming, o :rect (rectangular)."""
function _make_window(win::Symbol, N::Int)::Vector{Float64}
    if win == :hann
        return [0.5 * (1.0 - cos(2π * i / (N - 1))) for i in 0:(N-1)]
    elseif win == :hamming
        return [0.54 - 0.46 * cos(2π * i / (N - 1)) for i in 0:(N-1)]
    else  # :rect
        return ones(Float64, N)
    end
end
