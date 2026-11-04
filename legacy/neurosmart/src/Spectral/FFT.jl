#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# src/Spectral/FFT.jl
#
# FAST FOURIER TRANSFORMATION (FFT)
# ============================================
# Esta rutina realiza el análisis espectral de datos EEG segmentados mediante
# Transformada Rápida de Fourier (FFT) para estimar la potencia espectral en µV².
#
# Parámetros BrainVision Analyzer (BVA):
# - Output: Power (µV²)
# - Window: Tukey (Hamming) (Window Length: 10 %)
# - Apply Variance Correction
# - Mode: Periodic
# - Frequency Resolution: 0.976 Hz (Zero-padding a 512 puntos)
# - Use BVA Full Spectrum configuration (0 Hz - Nyquist)

# Importar librerías
using Serialization
using Statistics
using Plots
using DSP
using FFTW
using CSV
using DataFrames

# El backend GR se inicializa automáticamente al crear el primer plot

# ------------------------------------------------------------------------------------
# 1. CARGA DE DATOS DESPUÉS DEL 2ND BASELINE CORRECTION
# ------------------------------------------------------------------------------------
# Se cargan los datos EEG segmentados tras la segunda corrección de baseline.
# Formato: hipermatriz (canales × muestras × segmentos)
println("=" ^ 80)
println("📊 CARGA DE DATOS DESPUÉS DEL 2ND BASELINE CORRECTION")
println("=" ^ 80)
println()

dir_baseline = stage_dir(:baseline)
path_dict_2nd_baseline = joinpath(dir_baseline, "dict_2nd_baseline_correction.bin")
dict_2nd_baseline_info = Serialization.deserialize(path_dict_2nd_baseline)

# Extraer datos e información del diccionario
eeg_2nd_baseline_corrected = dict_2nd_baseline_info["eeg_2nd_baseline_corrected"]
channels = dict_2nd_baseline_info["channels"]
fs = dict_2nd_baseline_info["fs"]
segment_length_samples = dict_2nd_baseline_info["segment_length_samples"]
n_segments = dict_2nd_baseline_info["n_segments"]
n_channels = dict_2nd_baseline_info["n_channels"]

println("✔ Datos después del 2nd baseline correction cargados desde: $(basename(path_dict_2nd_baseline))")
println("  → Dimensiones: $(size(eeg_2nd_baseline_corrected))")
println("  → Formato: (canales × muestras × segmentos)")
println("  → Número de canales: $n_channels")
println("  → Muestras por segmento: $segment_length_samples")
println("  → Número de segmentos: $n_segments")
println("  → Frecuencia de muestreo: $(fs) Hz")
println("  → Duración de cada segmento: $(round(segment_length_samples/fs, digits=3)) s")
println()

# Definir directorio de figuras y datos FFT
dir_figures_FFT = stage_dir(:FFT; kind = :figures)
dir_fft = stage_dir(:FFT)
dir_tables_FFT = stage_dir(:FFT; kind = :tables)

# Se eliminan todos los archivos previos antes de iniciar el procesamiento
# para evitar confusiones con resultados antiguos.
println("✔Eliminando resultados previos")

# Función auxiliar para limpiar archivos de una carpeta por extensión
function limpiar_carpeta(dir::String, extension::String, descripcion::String)
    if isdir(dir)
        println("→ Limpiando $descripcion en: $(basename(dir))")
        archivos = filter(f -> endswith(f, extension), readdir(dir))
        for archivo in archivos
            rm(joinpath(dir, archivo))
            println("  ✓ Eliminado: $archivo")
        end
        isempty(archivos) && println("  → No se encontraron archivos previos")
    end
    println()
end

# Limpiar carpetas
limpiar_carpeta(dir_figures_FFT, ".png", "figuras previas")
limpiar_carpeta(dir_fft, ".bin", "datos previos")

println("✔ Limpieza completada")
println()

# ------------------------------------------------------------------------------------
# DC REMOVAL (BRAINVISION)
# ------------------------------------------------------------------------------------
# Remoción de DC (corrección de continua) al estilo BrainVision antes del ventanado.
# Se calcula la media de cada segmento completo por canal y se resta antes de aplicar la ventana.

println("📘 DC Removal (BrainVision-style)")

n_channels, n_samples, n_segments = size(eeg_2nd_baseline_corrected)

# DC por canal y segmento (media del segmento completo)
dc_values = zeros(Float64, n_channels, n_segments)

# Señal sin DC (misma forma que eeg_2nd_baseline_corrected)
eeg_dc_removed = Array{Float64}(undef, n_channels, n_samples, n_segments)

@inbounds for ch in 1:n_channels
    for seg in 1:n_segments
        # media del segmento completo (global mean / DC offset)
        dc = mean(view(eeg_2nd_baseline_corrected, ch, :, seg))
        dc_values[ch, seg] = dc

        # restar DC antes de aplicar la ventana
        eeg_dc_removed[ch, :, seg] .= view(eeg_2nd_baseline_corrected, ch, :, seg) .- dc
    end
end

println("✔ DC removal (BrainVision-style) aplicado antes de la ventana")
println("  → eeg_dc_removed size = ", size(eeg_dc_removed))
println("  → dc_values size      = ", size(dc_values))
println()

# ------------------------------------------------------------------------------------
# 2. VENTANADO (HAMMING - WINDOW LENGTH 10 %)
# ------------------------------------------------------------------------------------
# Ventana tipo taper basada en Hamming según BrainVision Analyzer:
#   - Para 0 ≤ t ≤ t₁: ω(t) = α - β cos( (2π / P) * (t / SL) )
#   - Para t₁ < t < t₂: ω(t) = 1
#   - Para t₂ ≤ t ≤ SL: ω(t) = α - β cos( (2π / P) * (1 - (t / SL)) )
# donde:
#   - t₁ = SL * (P / 2), t₂ = SL * (1 - (P / 2))
#   - P = 0.10 (Window Length 10 %)
#   - SL = longitud del segmento
#   - Hamming: α = 0.54, β = 0.46
# Objetivo: reducir spectral leakage minimizando discontinuidades en los bordes.

println("📘 Ventanado (Hamming)")

Nseg = segment_length_samples
println("→ Longitud de ventana: $Nseg muestras")

# Configuración de taper según BrainVision Analyzer
# Según el manual: t₁ = SL * (P / 2), t₂ = SL * (1 - (P / 2))
# P = 0.10 (Window Length 10 %)
P = 0.10
SL = Float64(Nseg)
t1_continuous = SL * (P / 2)  # t₁ en escala continua
t2_continuous = SL * (1 - (P / 2))  # t₂ en escala continua
n_taper = round(Int, t1_continuous)  # t₁ discretizado (número de muestras en cada taper)
println("→ P (Window Length): $(P * 100)%")
println("→ t₁ = $(round(t1_continuous, digits=2)) → $n_taper muestras")
println("→ t₂ = $(round(t2_continuous, digits=2)) → $(Nseg - n_taper) muestras")
println("→ Taper izquierdo: índices 1 a $n_taper")
println("→ Zona central plana: índices $(n_taper+1) a $(Nseg-n_taper) (valor 1.0)")
println("→ Taper derecho: índices $(Nseg-n_taper+1) a $Nseg")

# Construir ventana según fórmula exacta de BrainVision Analyzer
# Fórmula Hamming: α = 0.54, β = 0.46
# Para 0 ≤ t ≤ t₁: ω(t) = α - β cos( (2π / P) * (t / SL) )
# Para t₁ < t < t₂: ω(t) = 1
# Para t₂ ≤ t ≤ SL: ω(t) = α - β cos( (2π / P) * (1 - (t / SL)) )
# Nota: En discretización, índice i corresponde a t = i-1 (índice 1 → t=0, índice Nseg → t=SL-1)
win = ones(Float64, Nseg)
if n_taper > 0 && n_taper >= 2
    alpha = 0.54
    beta = 0.46
    
    # Taper izquierdo (0 ≤ t ≤ t₁): índices 1 a n_taper
    t_left = collect(0:(n_taper-1))  # Valores de t continuos (0 a n_taper-1)
    taper_left = alpha .- beta .* cos.((2*pi / P) .* (t_left ./ SL))
    win[1:n_taper] .= taper_left
    
    # Región central (t₁ < t < t₂): ω(t) = 1
    # Ya está inicializado en ones() para índices n_taper+1 a Nseg-n_taper
    
    # Taper derecho (t₂ ≤ t ≤ SL): índices Nseg-n_taper+1 a Nseg
    taper_right_start_idx = Nseg - n_taper + 1
    t_right = collect((taper_right_start_idx-1):(Nseg-1))  # Valores de t continuos
    taper_right = alpha .- beta .* cos.((2*pi / P) .* (1 .- (t_right ./ SL)))
    win[taper_right_start_idx:Nseg] .= taper_right
    
    # Verificación según especificaciones de BrainVision
    center_start = n_taper + 1
    center_end = Nseg - n_taper
    println("→ Verificación ventana (según BrainVision):")
    println("  → win[1] = $(round(win[1], digits=4)) (debe ser ~0.08)")
    println("  → win[$n_taper] = $(round(win[n_taper], digits=4)) (debe ser ~1.0)")
    if center_start <= center_end
        center_all_one = all(win[center_start:center_end] .== 1.0)
        println("  → win[$center_start:$center_end] todos = 1.0: $center_all_one")
    end
    println("  → win[$Nseg] = $(round(win[end], digits=4)) (debe ser ~0.08)")
end

# Calcular media de win^2 para corrección de varianza
mw2 = mean(win .^ 2)
taper_frac = P / 2  # Fracción de taper (P/2 a cada lado)
println("→ Corrección de varianza:")
println("  → mean(w²) = $(round(mw2, digits=6))")
println("  → taper_frac = $taper_frac")
println()

# Visualización diagnóstica de la ventana
println("→ Generando visualización de la ventana...")
indices = 1:Nseg
p1 = plot(indices, win, 
          linewidth=2, 
          color=:blue,
          title="Ventana de Hamming (N = $Nseg)",
          xlabel="Muestra (n)",
          ylabel="Amplitud w[n]",
          legend=false,
          grid=true,
          size=(800, 400))
p2 = plot(indices, win .^ 2,
          linewidth=2,
          color=:red,
          title="Energía de la ventana: w[n]²",
          xlabel="Muestra (n)",
          ylabel="w[n]²",
          legend=false,
          grid=true,
          size=(800, 400))

# Combinar ambos plots en una figura
p_window = plot(p1, p2, layout=(2, 1), size=(800, 600))
display(p_window)

# Guardar figura en directorio de resultados
path_fig_window = joinpath(dir_figures_FFT, "hamming_window.png")
savefig(p_window, path_fig_window)
println("✔ Figura guardada en: $(basename(path_fig_window))")
println()

# Aplicar ventana a todos los segmentos (sobre eeg_dc_removed)
println("→ Aplicando ventana a todos los segmentos...")
eeg_windowed = Array{Float64}(undef, n_channels, n_samples, n_segments)

@inbounds for ch in 1:n_channels
    for seg in 1:n_segments
        eeg_windowed[ch, :, seg] .= view(eeg_dc_removed, ch, :, seg) .* win
    end
end

println("✔ Ventana aplicada a todos los segmentos")
println("  → eeg_windowed size = ", size(eeg_windowed))
println()

# Visualización del efecto de la ventana: ejemplo de un canal y segmento
ch_example = 1
seg_example = 1
channel_name = channels[ch_example]
time_segment = collect(0:(Nseg-1)) ./ fs
segment_before = eeg_dc_removed[ch_example, :, seg_example]
segment_after = segment_before .* win  # Aplicar ventana al segmento de ejemplo

p_detail = plot(layout = (2, 1), size = (1000, 600), grid = true)

plot!(p_detail[1], time_segment, segment_before,
      label = "Antes ventanado", linewidth = 2.0, color = :blue, alpha = 0.7,
      title = "Canal: $channel_name - Segmento $seg_example",
      xlabel = "Tiempo (s)", ylabel = "Amplitud (µV)", legend = :outerright)
plot!(p_detail[1], time_segment, segment_after,
      label = "Después ventanado", linewidth = 2.5, color = :red, alpha = 0.9, linestyle = :dash)

plot!(p_detail[2], time_segment, win, label = "Ventana w[n]",
      linewidth = 2.0, color = :green, title = "Ventana aplicada",
      xlabel = "Tiempo (s)", ylabel = "Amplitud", legend = :outerright)

display(p_detail)
path_fig_detail = joinpath(dir_figures_FFT, "windowing_detail_example.png")
savefig(p_detail, path_fig_detail)
println("✔ Visualización del efecto de la ventana guardada en: $(basename(path_fig_detail))")
println()

# Guarda datos ventanados y parámetros en diccionario que se actualizará durante el proceso FFT
println("💾 Guardando datos ventanados y parámetros")
println()

dict_FFT_info = Dict(
    "eeg_windowed" => eeg_windowed,
    "eeg_2nd_baseline_corrected" => eeg_2nd_baseline_corrected,
    "window" => win,
    "window_type" => "Hamming-based (Tukey with Hamming taper)",
    "taper_fraction" => taper_frac,
    "n_taper" => n_taper,
    "window_length_samples" => Nseg,
    "mean_w_squared" => mw2,
    "variance_correction_applied_after_fft" => true,  # Corrección aplicada en dominio de frecuencia (P ./= mw2)
    "dc_values" => dc_values,  # Valores DC sustraídos antes del ventanado (para posible reinserción)
    "channels" => channels,
    "fs" => fs,
    "segment_length_samples" => segment_length_samples,
    "n_segments" => n_segments,
    "n_channels" => n_channels,
    "processing_step" => "FFT - Windowing and Variance Correction",
    "window_length_percentage" => taper_frac * 100 * 2,
)

path_dict_fft = joinpath(dir_fft, "dict_FFT.bin")
Serialization.serialize(path_dict_fft, dict_FFT_info)
println("✔ Diccionario guardado (ventanado) en: $(basename(path_dict_fft))")
println("  → Se actualizará con información de zero-padding y FFT")
println("  → Ubicación: $(abspath(path_dict_fft))")
println()

# ------------------------------------------------------------------------------------
# 5. ZERO-PADDING Y RESOLUCIÓN EN FRECUENCIA
# ------------------------------------------------------------------------------------
# ZERO-PADDING:
# El zero-padding consiste en ampliar artificialmente la longitud del segmento
# añadiendo muestras con valor cero al final de la señal ventaneada antes de calcular
# la FFT. Si un segmento contiene N muestras reales (ya ventaneadas) y se define una
# longitud de FFT mayor (Nfft > N), la señal se extiende hasta Nfft rellenando las
# posiciones adicionales con ceros. La información temporal original no se modifica
# ni se altera. El zero-padding no añade información nueva ni aumenta la resolución
# espectral real; únicamente interpola el espectro y permite una representación más
# fina del eje de frecuencias.
#
# RESOLUCIÓN ESPECTRAL:
# La resolución en frecuencia de la FFT viene determinada por la frecuencia de
# muestreo y por la longitud de la transformada utilizada. La separación entre
# frecuencias adyacentes (resolución espectral) es:
#   Δf = fs / Nfft
# En este caso:
#   fs = 500 Hz
#   Nfft = 512 muestras
#   Δf = 500 / 512 = 0.9766 Hz
#
# BINS ESPECTRALES:
# Un bin espectral es un punto discreto del espectro de Fourier que representa
# la potencia de la señal en un intervalo muy estrecho de frecuencias centrado
# en una frecuencia concreta. Cada bin k corresponde a la frecuencia:
#   f[k] = k * Δf
# donde k es un índice entero. Para una FFT real (rFFT), se generan Nfft/2 + 1
# bins que cubren desde 0 Hz (DC, bin 0) hasta la frecuencia de Nyquist
# (fs/2, bin Nfft/2). Los bins intermedios están espaciados uniformemente por Δf.

println("📘 Zero-padding")

Nseg = segment_length_samples
Nfft = 512  # potencia de 2 para FFT eficiente
pad = Nfft - Nseg

println("→ Longitud original: $Nseg puntos")
println("→ Longitud FFT: $Nfft puntos (zero-padding: $pad ceros)")

# Crear array zero-padded: copiar datos ventaneados, resto en cero
eeg_padded = zeros(Float64, n_channels, Nfft, n_segments)
eeg_padded[:, 1:Nseg, :] .= eeg_windowed

println("✔ Zero-padding aplicado ($Nseg → $Nfft puntos)")
println()

# Calcular resolución en frecuencia y ejes
println("📘 Resolución en frecuencia")
Δf = fs / Nfft
f_nyq = fs / 2
println("→ Resolución Δf = $(round(Δf, digits=4)) Hz")
println("→ Frecuencia Nyquist f_Nyq = $f_nyq Hz")
println()

# Ejes de frecuencia: rFFT devuelve frecuencias 0 … f_Nyq (half-spectrum técnico)
# Tras el folding tipo BrainVision, estas frecuencias representan el espectro completo
n_bins = Nfft ÷ 2 + 1  # Número de bins espectrales (incluye DC y Nyquist)
freqs_bva = range(0, f_nyq; length=n_bins)  # 0 … f_Nyq (257 puntos), paso exacto Δf
freqs_full = collect(0:Δf:(fs - Δf))       # 0, Δf, 2Δf, ..., fs-Δf (Nfft puntos), paso exacto Δf

println("→ Bins espectrales creados:")
println("  → Número de bins (rFFT): $n_bins")
println("  → Bin 0 (DC): 0.0 Hz")
println("  → Bin $(n_bins-1) (Nyquist): $(round(f_nyq, digits=2)) Hz")
println("  → Bins intermedios: espaciados por $(round(Δf, digits=4)) Hz")
println("  → Rango completo: 0.0 - $(round(f_nyq, digits=2)) Hz")
println("→ Bins full-spectrum: $(length(freqs_full)) frecuencias")
println()

# Actualizar diccionario FFT con información de zero-padding
dict_FFT_info["freqs_bva"] = freqs_bva
dict_FFT_info["freqs_full"] = freqs_full
dict_FFT_info["Nfft"] = Nfft
dict_FFT_info["df"] = Δf
dict_FFT_info["f_nyq"] = f_nyq
dict_FFT_info["eeg_padded"] = eeg_padded
dict_FFT_info["processing_step"] = "FFT - Windowing and Zero-padding"

# Visualización del efecto del zero-padding: ejemplo de un canal y segmento
ch_example = 1
seg_example = 1
channel_name = channels[ch_example]
time_original = collect(0:(Nseg-1)) ./ fs
time_padded = collect(0:(Nfft-1)) ./ fs
segment_original = eeg_windowed[ch_example, :, seg_example]
segment_padded = eeg_padded[ch_example, :, seg_example]

p_zeropad = plot(size = (1000, 500), grid = true)

# Señal original (antes del zero-padding)
plot!(p_zeropad, time_original, segment_original,
      label = "Antes zero-padding ($Nseg puntos)", linewidth = 2.5, color = :blue, alpha = 0.8,
      title = "Canal: $channel_name - Segmento $seg_example (Efecto del zero-padding: $Nseg → $Nfft puntos)",
      xlabel = "Tiempo (s)", ylabel = "Amplitud (µV)", legend = :outerright)

# Señal con zero-padding (después)
plot!(p_zeropad, time_padded, segment_padded,
      label = "Después zero-padding ($Nfft puntos)", linewidth = 2.0, color = :red, alpha = 0.7, linestyle = :dash)

# Marcar dónde termina la señal original y empiezan los ceros
vline!(p_zeropad, [time_original[end]], linestyle = :dash, color = :gray, alpha = 0.6, linewidth = 1.5,
      label = "Fin señal original (inicio ceros)")

display(p_zeropad)
path_fig_zeropad = joinpath(dir_figures_FFT, "zero_padding_example.png")
savefig(p_zeropad, path_fig_zeropad)
println("✔ Visualización del efecto de zero-padding guardada en: $(basename(path_fig_zeropad))")
println()

# Guardar diccionario actualizado
Serialization.serialize(path_dict_fft, dict_FFT_info)
println("✔ Diccionario FFT actualizado con información de zero-padding")
println()

# ------------------------------------------------------------------------------------
# 7. CÁLCULO DE FFT Y POTENCIA ESPECTRAL
# ------------------------------------------------------------------------------------
# Pasos separados según procesamiento tipo BrainVision Analyzer:
# Paso 1: Calcular rFFT (solo transformada, sin procesamiento)
# Paso 2: Reintroducir DC en bin 0 Hz
# Paso 3: Corrección de varianza (en dominio de amplitud)
# Paso 4: Calcular potencia, folding "Use Full Spectrum" y normalización

println("=" ^ 80)
println("📘 FFT (rFFT) + Potencia por bin")
println("=" ^ 80)

# Número de bins que devuelve rFFT: Nfft/2 + 1 (incluye DC y Nyquist)
n_bins_pos = Nfft ÷ 2 + 1
println("→ Número de bins (rFFT): $n_bins_pos (0 Hz ... Nyquist)")
println()

# Paso 1: Calcular rFFT (solo transformada)
# Entrada: eeg_padded (ch × Nfft × seg) -> señal ya DC-removed + ventaneada + zero-padded
# Salida: X (ch × (Nfft÷2+1) × seg) complejo

println("📘 Paso 1: Cálculo de rFFT")
X = Array{ComplexF64}(undef, n_channels, n_bins_pos, n_segments)

@inbounds for ch in 1:n_channels
    for seg in 1:n_segments
        x = view(eeg_padded, ch, :, seg)
        X[ch, :, seg] = rfft(x)
    end
end

println("✔ rFFT calculada en todos los canales y segmentos")
println("  → size(X) = ", size(X))
println()

# Visualización del espectro FFT: ejemplo de un canal y segmento
println("📘 Visualización del espectro FFT")
ch_example = 3  # Mismo canal que el plot de potencia
seg_example = 1
channel_name = channels[ch_example]
println("  → Usando canal índice $ch_example: $channel_name")

# Obtener frecuencias y FFT para el ejemplo
freqs_example = collect(freqs_bva)  # Frecuencias en Hz (0 ... f_nyq)
X_example = X[ch_example, :, seg_example]  # FFT compleja
X_magnitude = abs.(X_example)  # Magnitud del espectro

# Crear plot del espectro FFT (magnitud)
p_fft = plot(freqs_example, X_magnitude,
             linewidth=2,
             color=:blue,
             title="Espectro FFT (Magnitud) - $channel_name (Segmento $seg_example)",
             xlabel="Frecuencia (Hz)",
             ylabel="|X(f)|",
             legend=false,
             grid=true,
             size=(1000, 500),
             xlims=(0, f_nyq))

# Agregar línea vertical en Nyquist para referencia
vline!(p_fft, [f_nyq], 
       linewidth=1, 
       linestyle=:dash, 
       color=:red, 
       alpha=0.5,
       label="Nyquist ($(round(f_nyq, digits=1)) Hz)")

display(p_fft)

# Guardar figura
path_fig_fft = joinpath(dir_figures_FFT, "fft_spectrum_example.png")
savefig(p_fft, path_fig_fft)
println("✔ Figura del espectro FFT guardada en: $(basename(path_fig_fft))")
println("  → Canal: $channel_name")
println("  → Segmento: $seg_example")
println("  → Rango de frecuencias: 0 - $(round(f_nyq, digits=2)) Hz")
println()

# Paso 2: Reintroducir DC en bin 0 Hz (X[1])
# Paso 3: Corrección de varianza (en dominio de potencia: dividir por mw2)
# Paso 4: Power = |X|^2, folding "Use Full Spectrum" y normalización
println("📘 Pasos 2-4: Post-procesamiento (DC + Power + Varianza + Folding + Normalización)")

# Array para almacenar potencia
P = Array{Float64}(undef, n_channels, n_bins_pos, n_segments)

@inbounds for ch in 1:n_channels
    for seg in 1:n_segments
        Xseg = @view X[ch, :, seg]
        
        # Paso 2: Reintroducir DC en 0 Hz (bin 1)
        # rfft usa convención FFT no normalizada -> el DC en X[1] es suma(x)
        Xseg[1] += dc_values[ch, seg] * Nfft
        
        # Paso 3: Calcular potencia
        Pseg = abs2.(Xseg)
        
        # Paso 3 (continuación): Corrección de varianza (en dominio de potencia)
        # Dividir por mw2 para corregir la energía de la ventana
        Pseg ./= mw2
        
        # Paso 4: Folding tipo "Use Full Spectrum": dobla interiores (no DC, no Nyquist)
        if n_bins_pos > 2
            Pseg[2:end-1] .*= 2
        end
        
        # Paso 4 (continuación): Normalización (dividir por Nfft^2)
        Pseg ./= (Nfft^2)
        
        P[ch, :, seg] = Pseg
    end
end

println("✔ Paso 2: DC reintroducido en bin 0 Hz")
println("✔ Paso 3: Potencia calculada y corrección de varianza aplicada (P ./= mw2)")
println("✔ Paso 4: Folding 'Use Full Spectrum' y normalización aplicados")
println("  → size(P) = ", size(P))
println("  → Bins interiores doblados (DC y Nyquist sin doblar)")
println()

# Visualización del espectro de potencia: canal 3 y segmento 1
println("📘 Visualización del espectro de potencia")
ch_power = 3
seg_power = 1
channel_name_power = channels[ch_power]
println("  → Usando canal índice $ch_power: $channel_name_power")

# Obtener frecuencias y potencia para el ejemplo
freqs_power = collect(freqs_bva)  # Frecuencias en Hz (0 ... f_nyq)
P_power = P[ch_power, :, seg_power]  # Potencia espectral

# Crear plot del espectro de potencia
p_power = plot(freqs_power, P_power,
               linewidth=2,
               color=:blue,
               title="Espectro de Potencia - $channel_name_power (Segmento $seg_power)",
               xlabel="Frecuencia (Hz)",
               ylabel="Potencia (µV²)",
               legend=false,
               grid=true,
               size=(1000, 500),
               xlims=(0, f_nyq))

# Agregar línea vertical en Nyquist para referencia
vline!(p_power, [f_nyq], 
       linewidth=1, 
       linestyle=:dash, 
       color=:red, 
       alpha=0.5,
       label="Nyquist ($(round(f_nyq, digits=1)) Hz)")

display(p_power)

# Guardar figura
path_fig_power = joinpath(dir_figures_FFT, "power_spectrum_ch3_seg1.png")
savefig(p_power, path_fig_power)
println("✔ Figura del espectro de potencia guardada en: $(basename(path_fig_power))")
println("  → Canal: $channel_name_power (canal $ch_power)")
println("  → Segmento: $seg_power")
println("  → Rango de frecuencias: 0 - $(round(f_nyq, digits=2)) Hz")
println()

# ------------------------------------------------------------------------------------
# 8. DEFINICIÓN DE BANDAS DE FRECUENCIA
# ------------------------------------------------------------------------------------
# Define bandas espectrales estándar de EEG y calcula índices correspondientes

println("📘 Definición de bandas de frecuencia")

# Bandas EEG definidas por el protocolo experimental.
# BrainVision Analyzer no impone bandas: estas se aplican seleccionando
# bins cuya frecuencia central cae dentro del intervalo especificado.
# Usamos vector de tuplas para mantener el orden (DELTA → GAMMA)

bands_hz = [
    (:DELTA,     (0.5, 4.0)),
    (:THETA,     (4.0, 8.0)),
    (:ALPHA,     (7.8, 11.7)),
    (:BETA_LOW,  (12.0, 15.0)),
    (:BETA_MID,  (15.0, 18.0)),
    (:BETA_HIGH, (18.0, 30.0)),
    (:GAMMA,     (30.0, 50.0))
]

# Función auxiliar: obtener índices de bins cuya frecuencia central cae en [fmin, fmax)
# Usamos intervalo semiabierto para evitar solapes entre bandas
function band_indices(freqs::AbstractVector{<:Real}, fmin::Real, fmax::Real)
    return findall(f -> f >= fmin && f < fmax, freqs)
end

# Calcular índices de cada banda en freqs_bva
band_idx = Dict{Symbol, Vector{Int}}()

println("Bandas de frecuencia (selección por bins):")
for (name, (fmin, fmax)) in bands_hz
    idx = band_indices(freqs_bva, fmin, fmax)
    band_idx[name] = idx
    println("→ $name : $(length(idx)) bins en [$fmin, $fmax) Hz")
end
println()

# ------------------------------------------------------------------------------------
# 9. GUARDADO DE RESULTADOS DE POTENCIA ESPECTRAL
# ------------------------------------------------------------------------------------
# Guarda la potencia espectral calculada en un diccionario separado para análisis
# posteriores (promediado, análisis por bandas, etc.).
println("=" ^ 80)
println("💾 GUARDADO DE RESULTADOS DE POTENCIA ESPECTRAL")
println("=" ^ 80)
println()

dict_FFT_power = Dict(
    "P" => P,                          # Potencia espectral (n_channels × n_bins_pos × n_segments)
    "freqs_bva" => freqs_bva,          # Frecuencias formato BrainVision (tras folding)
    "freqs_full" => freqs_full,        # Frecuencias completas (full-spectrum)
    "n_bins_pos" => n_bins_pos,        # Número de bins (Nfft/2 + 1)
    "Nfft" => Nfft,                    # Longitud de la FFT
    "df" => Δf,                        # Resolución en frecuencia
    "f_nyq" => f_nyq,                  # Frecuencia de Nyquist
    "bands_hz" => bands_hz,            # Definición de bandas de frecuencia
    "band_idx" => band_idx,            # Índices de bandas en freqs_bva
    "channels" => channels,             # Nombres de canales
    "fs" => fs,                        # Frecuencia de muestreo
    "n_segments" => n_segments,        # Número de segmentos
    "n_channels" => n_channels,        # Número de canales
    "processing_step" => "FFT - Power Spectrum",
)

path_dict_fft_power = joinpath(dir_fft, "dict_FFT_power.bin")
Serialization.serialize(path_dict_fft_power, dict_FFT_power)
println("✔ Diccionario de potencia espectral guardado en: $(basename(path_dict_fft_power))")
println("  → Ubicación: $(abspath(path_dict_fft_power))")
println()          

# ------------------------------------------------------------------------------------
# 12. PROMEDIO SOBRE SEGMENTOS + PLOTEO (ESCALA LINEAL)
# ------------------------------------------------------------------------------------
# P: (n_channels × n_bins × n_segments)
# freqs_bva: vector de frecuencias (n_bins)
# channels: vector de nombres de canal
# Promedia la potencia sobre todos los segmentos de cada canal (BrainVision: mean over segments)

println("=" ^ 80)
println("📊 PROMEDIO SOBRE SEGMENTOS + PLOTEO (ESCALA LINEAL)")
println("=" ^ 80)
println()

# 1) Promedio sobre segmentos (BrainVision-style: mean over segments)
P_mean = mean(P, dims=3)                 # (canales × bins × 1)
P_mean_2d = dropdims(P_mean, dims=3)     # (canales × bins)

println("✔ P_mean calculada: ", size(P_mean_2d), " (canales × bins)")
println()

# 2) Parámetros de visualización
xmax = 50.0                              # Gamma hasta 50 Hz
xmask = freqs_bva .<= xmax
freqs_plot = freqs_bva[xmask]

# 3) Layout tipo grid
n_rows = 6
n_cols = 6
plots_list = Vector{Any}(undef, n_rows * n_cols)

for pos in 1:(n_rows * n_cols)
    ch = pos <= size(P_mean_2d, 1) ? pos : 0

    if ch > 0
        row = ((pos - 1) ÷ n_cols) + 1
        
        # Escala vertical independiente para cada canal
        P_ch_plot = P_mean_2d[ch, xmask]
        ymin_ch = 0.0
        ymax_ch = maximum(P_ch_plot)

        p_ch = plot(
            freqs_plot, P_ch_plot,
            title = channels[ch],
            titlefontsize = 8,
            xlabel = (row == n_rows ? "Hz" : ""),
            ylabel = "",
            xlims = (0, xmax),
            ylims = (ymin_ch, ymax_ch),
            grid = false,
            legend = false,
            framestyle = :box,
            linewidth = 1.5,
            size = (150, 105)
        )

        plots_list[pos] = p_ch
    else
        plots_list[pos] = plot(
            framestyle = :none,
            grid = false,
            showaxis = false,
            size = (150, 105)
        )
    end
end

# 4) Figura final
p_grid = plot(
    plots_list...,
    layout = (n_rows, n_cols),
    size = (1200, 820),
    plot_title = "Power Spectrum by Channel (mean over segments, linear scale)"
)

display(p_grid)

# Guardar figura
path_fig = joinpath(dir_figures_FFT, "power_spectrum_grid_mean_segments_linear.png")
savefig(p_grid, path_fig)
println("✔ Figura guardada en: $(basename(path_fig))")
println()

# Calcular la desviación estándar de la potencia espectral sobre todos los segmentos de cada canal
P_std = std(P, dims=3)
SNR = P_mean ./ P_std
snr_channel = mean(SNR, dims=2)

println("→ SNR promediado por canal: $(size(snr_channel)) (canales)")
println()

# Imprimir tabla resumen
println("Tabla resumen: [Canal, AvgSegments, SNR]")
println("-"^60)
println(rpad("Canal", 15), rpad("AvgSeg", 10), "SNR")
println("-"^60)
snr_values = dropdims(snr_channel, dims=2)
for ch in 1:n_channels
    println(rpad(channels[ch], 15), rpad(string(n_segments), 10), round(snr_values[ch], digits=3))
end
println("-"^60)
println()

# Crear tabla resumen por canal (Canal, AvgSegments, SNR)
println("→ Generando tabla resumen por canal...")
df_summary = DataFrame(
    Canal = channels,
    AvgSegments = repeat([n_segments], n_channels),
    SNR = vec(dropdims(snr_channel, dims=(2, 3)))
)

# Guardar tabla resumen en el directorio de tablas
path_table_summary = joinpath(dir_tables_FFT, "channel_summary.csv")
CSV.write(path_table_summary, df_summary)
println("✔ Tabla resumen guardada en: $(basename(path_table_summary))")
println("  → Ubicación: $(abspath(path_table_summary))")
println()

# Guardar P_mean, P_std y SNR en el diccionario
dict_FFT_power["P_mean"] = P_mean
dict_FFT_power["P_std"] = P_std
dict_FFT_power["SNR"] = SNR
dict_FFT_power["snr_channel"] = snr_channel
dict_FFT_power["processing_step"] = "FFT - Power Spectrum (averaged)"

# Serializar diccionario actualizado
Serialization.serialize(path_dict_fft_power, dict_FFT_power)
println("✔ Potencia promediada guardada en diccionario")
println("  → P_mean: $(size(P_mean)) (canales × bins positivos)")
println("  → P_std: $(size(P_std))")
println("  → SNR: $(size(SNR))")
println()

# ------------------------------------------------------------------------------------
# 14. CÁLCULO DE POTENCIA POR BANDAS DE FRECUENCIA
# ------------------------------------------------------------------------------------
# Calcula la potencia promedio por banda de frecuencia para cada canal,
# usando la potencia promediada (P_mean).
# En BrainVision Analyzer, una vez definidos los bins por la FFT, todo el cálculo
# posterior trabaja sobre bins, no sobre límites continuos en Hz.
println("=" ^ 80)
println("📊 CÁLCULO DE POTENCIA POR BANDAS DE FRECUENCIA")
println("=" ^ 80)
println()

# 1) Preparación: aplanar P_mean a 2D (canales × bins)
P_mean_2d = dropdims(P_mean, dims=3)

# 2) Orden fijo de bandas (como ya está definido en bands_hz)
band_order = [:DELTA, :THETA, :ALPHA, :BETA_LOW, :BETA_MID, :BETA_HIGH, :GAMMA]

# 3) Total 0.5–50 Hz (usando bins reales de las bandas)
# En lugar de recalcular por Hz, usamos todos los bins de las bandas
# Unión de todos los bins de bandas (0.5–50 Hz real)
idx_total = sort(unique(vcat([band_idx[b] for b in band_order]...)))

# Promedio entre líneas espectrales (bins) del rango total
power_total = vec(mean(P_mean_2d[:, idx_total], dims=2))

println("→ Potencia total (average de bins reales 0.5–50 Hz) calculada")
println("  → Basado exclusivamente en bins definidos (BrainVision-style)")
println("  → Nº bins total: $(length(idx_total))")
println()

# 4) Potencia por bandas (average entre líneas espectrales)
power_by_band = Dict{Symbol, Vector{Float64}}()

for band in band_order
    idx = band_idx[band]

    if isempty(idx)
        power_by_band[band] = fill(NaN, n_channels)
        println("→ $band: sin bins (NaN)")
        continue
    end

    # Average entre líneas espectrales (bins) de la banda
    power_by_band[band] = vec(mean(P_mean_2d[:, idx], dims=2))

    println("→ $band: $(length(idx)) bins | average calculado")
end

println()

# Crear tabla de potencia por bandas
println("→ Generando tabla de potencia por bandas...")

# Crear DataFrame: filas = TOTAL + bandas, columnas = canales
band_names = [string(b) for b in band_order]
row_names = ["TOTAL"; band_names]

# Construir DataFrame con todos los Pairs
df_power_by_band = DataFrame(
    "File" => row_names,
    [channels[ch] => [power_total[ch]; [power_by_band[b][ch] for b in band_order]] 
     for ch in 1:n_channels]...
)

# Mostrar tabla por pantalla
println("Tabla de potencia por bandas:")
display(df_power_by_band)
println()

# Guardar tabla en formato CSV
path_table_power_by_band = joinpath(dir_tables_FFT, "power_by_band.csv")
CSV.write(path_table_power_by_band, df_power_by_band)
println("✔ Tabla de potencia por bandas guardada en: $(basename(path_table_power_by_band))")
println()

# Guardar potencia por bandas en el diccionario
dict_FFT_power["power_by_band"] = power_by_band
dict_FFT_power["power_total"] = power_total
dict_FFT_power["processing_step"] = "FFT - Power Spectrum (averaged) + Band Power"
Serialization.serialize(path_dict_fft_power, dict_FFT_power)
println("✔ Potencia por bandas guardada en diccionario")
println()

# ------------------------------------------------------------------------------------
# COMPARACIÓN CON RESULTADOS DE REFERENCIA (BrainVision Analyzer)
# ------------------------------------------------------------------------------------
println("📊 Comparación con resultados de Javier (BrainVision Analyzer)")
println()

# Leer archivo de Javier
path_javier = joinpath(@__DIR__, "..", "Javier_results", "M5T2cerrados.txt")
lines = filter(x -> !isempty(strip(x)), readlines(path_javier))

# Parsear encabezado y datos
header = split(strip(lines[1]), r"\s{2,}")
channels_javier = [split(ch, "-")[1] for ch in header[2:end]]  # Extraer nombre base del canal

rows = [split(strip(line), r"\s{2,}") for line in lines[2:end]]
bands_javier = [first(row) for row in rows]

# Crear diccionario similar al nuestro
power_by_band_javier = Dict{Symbol, Vector{Float64}}()
power_total_javier = Dict{String, Float64}()

# Mapeo de nombres de bandas
band_map = Dict(
    "TOTAL" => :TOTAL,
    "DELTA (0.5-4)" => :DELTA,
    "THETA (4-8)" => :THETA,
    "ALPHA (8-12)" => :ALPHA,
    "BETA LOW (12-15)" => :BETA_LOW,
    "BETA MID (15-18)" => :BETA_MID,
    "BETA HIGH (18-30)" => :BETA_HIGH,
    "GAMMA (30-50)" => :GAMMA
)

for (i, band_name) in enumerate(bands_javier)
    band_sym = get(band_map, band_name, nothing)
    if band_sym !== nothing
        values = [parse(Float64, replace(rows[i][j], ',' => '.')) for j in 2:length(header)]
        
        if band_sym == :TOTAL
            for (ch_name, val) in zip(channels_javier, values)
                power_total_javier[ch_name] = val
            end
        else
            power_by_band_javier[band_sym] = values
        end
    end
end

println("✔ Datos de Javier cargados")
println()

# Comparación: crear DataFrame
comparison_rows = []
for (ch_javier, ch_idx_javier) in zip(channels_javier, 1:length(channels_javier))
    ch_idx_ours = findfirst(==(ch_javier), channels)
    if ch_idx_ours !== nothing
        ch_name = channels[ch_idx_ours]
        
        # Comparar TOTAL
        if haskey(power_total_javier, ch_javier)
            push!(comparison_rows, (
                Canal = ch_name,
                Banda = "TOTAL",
                Javier = power_total_javier[ch_javier],
                Rafa = power_total[ch_idx_ours],
                Diferencia = power_total[ch_idx_ours] - power_total_javier[ch_javier],
                Porcentaje = ((power_total[ch_idx_ours] - power_total_javier[ch_javier]) / power_total_javier[ch_javier]) * 100
            ))
        end
        
        # Comparar bandas
        for band in band_order
            if haskey(power_by_band_javier, band) && haskey(power_by_band, band)
                val_javier = power_by_band_javier[band][ch_idx_javier]
                val_ours = power_by_band[band][ch_idx_ours]
                push!(comparison_rows, (
                    Canal = ch_name,
                    Banda = string(band),
                    Javier = val_javier,
                    Rafa = val_ours,
                    Diferencia = val_ours - val_javier,
                    Porcentaje = ((val_ours - val_javier) / val_javier) * 100
                ))
            end
        end
    end
end

df_comparison = DataFrame(comparison_rows)

# Guardar comparación
path_comparison = joinpath(dir_tables_FFT, "comparison_javier.csv")
CSV.write(path_comparison, df_comparison)
println("✔ Comparación guardada en: $(basename(path_comparison))")
println()

# Mostrar tabla de ejemplo (primer canal disponible)
ch_example = first(unique(df_comparison.Canal))
df_example = filter(row -> row.Canal == ch_example, df_comparison)
println("📊 Ejemplo: Canal $ch_example")
display(df_example)
println()

# ------------------------------------------------------------------------------------
# VISUALIZACIÓN TIPO BRAINVISION ANALYZER (BVA-LIKE)
# ------------------------------------------------------------------------------------
# Genera visualizaciones similares a BrainVision Analyzer:
# - Espectro de potencia con bandas coloreadas (estilo BVA)
# - Topomap de banda Alpha con paleta JET (estilo BVA)
println("=" ^ 80)
println("📊 VISUALIZACIÓN TIPO BRAINVISION ANALYZER (BVA-LIKE)")
println("=" ^ 80)
println()

# El backend GR ya está inicializado automáticamente

# ---------------------------
# Helpers
# ---------------------------
function load_electrode_positions(tsv_path::String)
    df = CSV.read(tsv_path, DataFrame)
    pos = Dict{String,Tuple{Float64,Float64}}()
    for r in eachrow(df)
        if hasproperty(r, :type) && r.type == "EEG"
            name = uppercase(strip(String(r.name)))
            pos[name] = (Float64(r.x), Float64(r.y))
        end
    end
    return pos
end

function channels_xy(channels::Vector{String}, posdict::Dict{String,Tuple{Float64,Float64}})
    xs = Vector{Union{Missing,Float64}}(undef, length(channels))
    ys = Vector{Union{Missing,Float64}}(undef, length(channels))
    for (i,ch) in enumerate(channels)
        key = uppercase(strip(ch))
        if haskey(posdict, key)
            xs[i], ys[i] = posdict[key]
        else
            xs[i], ys[i] = missing, missing
        end
    end
    return xs, ys
end

function idw_grid(xs, ys, vals; grid_res=200, p=2.0, epsd=1e-6)
    xg = range(-1, 1; length=grid_res)
    yg = range(-1, 1; length=grid_res)
    Z  = fill(NaN, grid_res, grid_res)

    idx = findall(i -> !(ismissing(xs[i]) || ismissing(ys[i])), eachindex(xs))
    x = Float64[xs[i] for i in idx]
    y = Float64[ys[i] for i in idx]
    v = Float64[vals[i] for i in idx]

    @inbounds for j in 1:grid_res, i in 1:grid_res
        xx = xg[i]; yy = yg[j]
        if xx^2 + yy^2 <= 1.0
            d2 = (xx .- x).^2 .+ (yy .- y).^2
            w  = 1.0 ./ ((d2 .+ epsd).^(p/2))
            Z[j,i] = sum(w .* v) / sum(w)
        end
    end
    return xg, yg, Z
end

function draw_head!(p)
    θ = range(0, 2π; length=200)
    plot!(p, sin.(θ), cos.(θ), color=:black, linewidth=2, label="")

    θn = range(-π/6, π/6; length=40)
    plot!(p, 0.12 .* sin.(θn), 1.0 .+ 0.10 .* cos.(θn), color=:black, linewidth=2, label="")

    θe = range(-π/3, π/3; length=60)
    plot!(p, -1.0 .+ 0.08 .* cos.(θe), 0.15 .* sin.(θe), color=:black, linewidth=2, label="")
    plot!(p,  1.0 .- 0.08 .* cos.(θe), 0.15 .* sin.(θe), color=:black, linewidth=2, label="")
    return p
end

# ---------------------------
# Espectro BVA-like (como captura)
# ---------------------------
function plot_bva_like_spectrum(
    freqs::AbstractVector{<:Real},
    spectrum::AbstractVector{<:Real},
    band_idx::Dict{Symbol,Vector{Int}},
    band_order::Vector{Symbol};
    channel_label::String="Oz",
    xlims=(0.0, 120.0),
    alpha_band::Symbol=:ALPHA
)
    # Colores tipo captura (legend inferior en BVA)
    band_colors = Dict(
        :DELTA     => :orange,
        :THETA     => :yellow,
        :ALPHA     => :green,
        :BETA_LOW  => :blue,
        :BETA_MID  => :blue,
        :BETA_HIGH => :blue,
        :GAMMA     => :black
    )

    p = plot(
        xlabel="Hz",
        ylabel="µV²",
        title=channel_label,
        legend=false,
        grid=false,
        xlims=xlims,
        ylims=(0.0, maximum(spectrum[freqs .<= xlims[2]]) * 1.10),
        framestyle=:box,
        size=(900, 550)
    )

    # Relleno por bandas (área bajo la curva)
    for b in band_order
        idx = get(band_idx, b, Int[])
        isempty(idx) && continue

        # recorta a xlims para que no pinte fuera
        idx2 = [i for i in idx if freqs[i] >= xlims[1] && freqs[i] <= xlims[2]]
        isempty(idx2) && continue

        plot!(p, freqs[idx2], spectrum[idx2],
              fillrange=0,
              fillalpha=0.9,
              linealpha=0.0,
              color=band_colors[b],
              label="")
    end

    # línea negra del espectro completo
    mask = (freqs .>= xlims[1]) .& (freqs .<= xlims[2])
    plot!(p, freqs[mask], spectrum[mask], color=:black, linewidth=2.0, alpha=0.95, label="")

    # banda vertical resaltada (selección alpha) — similar a la captura
    if haskey(band_idx, alpha_band) && !isempty(band_idx[alpha_band])
        idxa = [i for i in band_idx[alpha_band] if freqs[i] >= xlims[1] && freqs[i] <= xlims[2]]
        if !isempty(idxa)
            fmin = minimum(freqs[idxa])
            fmax = maximum(freqs[idxa])
            vspan!(p, (fmin, fmax), color=:mediumpurple, alpha=0.25, label="")
        end
    end

    return p
end

# ---------------------------
# Topomap BVA-like (paleta EXACTA tipo BVA = JET)
# ---------------------------
function plot_bva_like_topomap_alpha(
    channels::Vector{String},
    alpha_vals::Vector{Float64},
    tsv_path::String;
    subtitle_txt="",
    grid_res=220
)
    pos = load_electrode_positions(tsv_path)
    xs, ys = channels_xy(channels, pos)

    xg, yg, Z = idw_grid(xs, ys, alpha_vals; grid_res=grid_res, p=2.0)

    zmin = minimum(alpha_vals)
    zmax = maximum(alpha_vals)

    # Paleta tipo BrainVision (JET)
    cmap = cgrad(:jet)

    # Filtrar electrodos que están dentro del círculo de la cabeza (radio <= 1.0)
    idx_ok = findall(i -> !(ismissing(xs[i]) || ismissing(ys[i])), eachindex(xs))
    idx_inside = [i for i in idx_ok if xs[i]^2 + ys[i]^2 <= 1.0]
    
    x_inside = Float64[xs[i] for i in idx_inside]
    y_inside = Float64[ys[i] for i in idx_inside]

    # Crear el plot base con colorbar horizontal debajo
    p = contour(
        xg, yg, Z,
        fill=true,
        levels=40,
        linewidth=0.0,
        color=cmap,
        colorbar=true,
        colorbar_title="µV²",
        colorbar_orientation=:horizontal,
        colorbar_position=:bottom,
        clims=(zmin, zmax),
        aspect_ratio=1,
        xlims=(-1.15, 1.15),
        ylims=(-1.15, 1.15),
        xlabel="",
        ylabel="",
        framestyle=:box,
        grid=false,
        showaxis=true,
        size=(550, 600),
        title=isempty(subtitle_txt) ? "Mapping View" : "Mapping View\n$subtitle_txt",
    )

    draw_head!(p)

    # Dibujar solo los electrodos dentro del círculo
    scatter!(p, x_inside, y_inside,
             markersize=4, color=:black, markerstrokecolor=:white, markerstrokewidth=1,
             label="")

    return p
end

# ---------------------------
# USO: Oz + Alpha
# ---------------------------
P_mean_2d = dropdims(P_mean, dims=3)

oz_idx = findfirst(==("Oz"), channels)
oz_idx === nothing && error("Canal Oz no encontrado en channels")
oz_spectrum = P_mean_2d[oz_idx, :]

band_order = [:DELTA, :THETA, :ALPHA, :BETA_LOW, :BETA_MID, :BETA_HIGH, :GAMMA]

p_spec = plot_bva_like_spectrum(
    freqs_bva, oz_spectrum, band_idx, band_order;
    channel_label="Oz",
    xlims=(0.0, 120.0),
    alpha_band=:ALPHA
)

alpha_bins = band_idx[:ALPHA]
alpha_vals = vec(mean(P_mean_2d[:, alpha_bins], dims=2))
alpha_fmin = round(minimum(freqs_bva[alpha_bins]); digits=1)
alpha_fmax = round(maximum(freqs_bva[alpha_bins]); digits=1)

tsv_path = joinpath(electrodes_dir(), "sub-M05_ses-T2_electrodes.tsv")
p_map = plot_bva_like_topomap_alpha(
    channels, alpha_vals, tsv_path;
    subtitle_txt="$(alpha_fmin) Hz - $(alpha_fmax) Hz"
)

p_combined = plot(
    p_spec, p_map,
    layout=@layout([a{0.62w} b{0.38w}]),
    size=(1450, 600)
)

display(p_combined)

out_path = joinpath(dir_figures_FFT, "bva_like_Oz_spectrum_and_alpha_map_JET.png")
savefig(p_combined, out_path)
println("✔ Guardado: $(basename(out_path))")
println()

