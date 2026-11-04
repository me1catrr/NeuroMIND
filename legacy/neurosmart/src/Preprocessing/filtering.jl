#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# src/filtering.jl
#
# FILTRADO DE SEÑALES EEG
# =======================
# Esta rutina aplica una cascada de filtros digitales a las señales EEG cargadas
# desde el paso de IO, con el objetivo de eliminar ruido de red, artefactos de
# hardware y limitar el ancho de banda al rango de interés.
#
# PROCESO:
#   1. Carga del diccionario EEG desde data/IO (dict_EEG.bin)
#   2. Filtro Notch (50 Hz): elimina interferencia de la red eléctrica
#   3. Filtro Bandreject (100 Hz): elimina posible interferencia de hardware
#   4. Filtro Highpass (0.5 Hz): elimina deriva (drift) preservando actividad lenta
#   5. Filtro Lowpass (150 Hz): limita alias y ruido de alta frecuencia
#
# En cada paso se guarda el resultado en data/filtering/, se visualiza la respuesta
# del filtro y se compara la PSD promedio antes/después.
#
# NOTA: Los filtros pasa-altos y pasa-bajos usan filtfilt (zero phase) por defecto.
# Las frecuencias están en Hz; la frecuencia de muestreo se asume 500 Hz.

# Listado de librerías
using DSP
using Plots
# El backend GR se inicializa automáticamente al crear el primer plot
using Serialization
using InlineStrings

# Si se ejecuta este script directamente (fuera del módulo EEG_Julia),
# cargamos utilidades de rutas para disponer de `stage_dir`.
if !@isdefined(stage_dir)
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end

# ------------------------------------------------------------------------------------
# 1. CARGA DE DATOS Y CONFIGURACIÓN
# ------------------------------------------------------------------------------------
# Se carga el diccionario EEG desde el paso de IO (dict_EEG.bin) y se configuran
# las rutas de salida para los datos filtrados.

dir_io = stage_dir(:IO)
path_dict = joinpath(dir_io, "dict_EEG.bin")
if !isfile(path_dict)
    error("No se encontró $(abspath(path_dict)). Ejecuta antes src/Setup/IO.jl para generar dict_EEG.bin.")
end
dict_EEG = Serialization.deserialize(path_dict)

# Directorio base para datos filtrados
dir_filtering = stage_dir(:filtering)

# ------------------------------------------------------------------------------------
# 2. INFORMACIÓN TEMPORAL DEL REGISTRO
# ------------------------------------------------------------------------------------
fs = 500                                        # Frecuencia de muestreo (Hz)
n_muestras = length(dict_EEG[first(keys(dict_EEG))])  # Muestras por canal
tiempo_seg = collect(0:(n_muestras-1)) ./ fs    # En segundos
duracion_total = tiempo_seg[end]

println("⏱️  INFORMACIÓN TEMPORAL")
println("-" ^ 80)
println("  Frecuencia de muestreo: $(fs) Hz")
println("  Número de muestras: $(n_muestras)")
println("  Duración total: $(round(duracion_total, digits=2)) segundos ($(round(duracion_total/60, digits=2)) minutos)")
println()

# Mostramos los datos de los canales
println("📊 Datos de los canales")
println("-" ^ 80)
display(dict_EEG)
println()

# Filtrado de las señales
println("📊 Pasos a seguir en el filtrado de las señales:")
println("-" ^ 50)
println("STEP Nº1. Notch Filter")
println("STEP Nº2. Bandreject filter")
println("STEP Nº3. Highpass filter")
println("STEP Nº4. Lowpass filter")
println()

# ------------------------------------------------------------------------------------
# 3. DEFINICIÓN DE FUNCIONES DE FILTRADO Y UTILIDADES
# ------------------------------------------------------------------------------------
# Funciones para filtros Butterworth (notch, bandreject, highpass, lowpass),
# visualización de respuesta en frecuencia, PSD y guardado de señales filtradas.

"""Aplica un filtro Butterworth notch_filter a un vector x"""
function notch_filter(x::Vector{<:Real}, freq::Real, fs::Real, order::Int = 5, width::Real = 1.0)  
    # Calculamos las frecuencias inferior y superior del notch
    freq_low = freq - width/2
    freq_high = freq + width/2
    # Normalizamos las frecuencias (0–1)
    wn_low = freq_low / (fs/2)
    wn_high = freq_high / (fs/2)
    responsetype = Bandstop(wn_low, wn_high)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    y = filt(flt, x)                # señal filtrada
    return y, flt                   # devolvemos también el filtro para analizarlo
end

"""Aplica un filtro Butterworth bandreject (bandstop) a un vector x"""
function bandreject_filter(x::Vector{<:Real}, freq::Real, fs::Real, order::Int = 4, bandwidth::Real = 1.0)  
    # Calculamos las frecuencias inferior y superior del bandreject
    freq_low = freq - bandwidth/2
    freq_high = freq + bandwidth/2
    # Normalizamos las frecuencias (0–1)
    wn_low = freq_low / (fs/2)
    wn_high = freq_high / (fs/2)
    responsetype = Bandstop(wn_low, wn_high)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    y = filt(flt, x)                # señal filtrada
    return y, flt                   # devolvemos también el filtro para analizarlo
end

"""Aplica un filtro Butterworth pasa-bajos a un vector x con zero phase shift (filtfilt)
   
   Nota: filtfilt aplica el filtro dos veces (adelante y atrás), duplicando el orden efectivo.
   Para obtener orden N efectivo, se debe diseñar un filtro de orden N/2.
   
   Parámetros:
   - x: señal de entrada
   - cutoff: frecuencia de corte en Hz
   - fs: frecuencia de muestreo en Hz
   - order: orden del diseño del filtro (el orden efectivo será order × 2 con filtfilt)
   - zero_phase: si es true (por defecto), usa filtfilt para zero phase shift
"""
function lowpass_filter(x::Vector{<:Real}, cutoff::Real, fs::Real, order::Int = 5; zero_phase::Bool = true)
    wn = cutoff / (fs/2)            # Frecuencia de corte normalizada (0–1)
    responsetype = Lowpass(wn)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    if zero_phase
        y = filtfilt(flt, x)        # señal filtrada con zero phase (orden efectivo = order × 2)
    else
        y = filt(flt, x)            # señal filtrada con desplazamiento de fase (orden efectivo = order)
    end
    return y, flt                   # devolvemos también el filtro para analizarlo
end

"""Aplica un filtro Butterworth pasa-altos a un vector x con zero phase shift (filtfilt)
   
   Nota: filtfilt aplica el filtro dos veces (adelante y atrás), duplicando el orden efectivo.
   Para obtener orden N efectivo, se debe diseñar un filtro de orden N/2.
   
   Parámetros:
   - x: señal de entrada
   - cutoff: frecuencia de corte en Hz
   - fs: frecuencia de muestreo en Hz
   - order: orden del diseño del filtro (el orden efectivo será order × 2 con filtfilt)
   - zero_phase: si es true (por defecto), usa filtfilt para zero phase shift
"""
function highpass_filter(x::Vector{<:Real}, cutoff::Real, fs::Real, order::Int = 5; zero_phase::Bool = true)
    wn = cutoff / (fs/2)            # Frecuencia de corte normalizada (0–1)
    responsetype = Highpass(wn)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    if zero_phase
        y = filtfilt(flt, x)        # señal filtrada con zero phase (orden efectivo = order × 2)
    else
        y = filt(flt, x)            # señal filtrada con desplazamiento de fase (orden efectivo = order)
    end
    return y, flt                   # devolvemos también el filtro para analizarlo
end 

"""Plotea la respuesta en frecuencia de un filtro (magnitud y fase)"""
function plot_filter_response(flt, fs::Real; title::String = "Respuesta del Filtro", n_points::Int = 1024)
    # Frecuencias en rad/muestra donde evaluar (0..π)
    ω = range(0, stop = π, length = n_points)
    # Respuesta compleja del filtro en esas frecuencias
    H = freqresp(flt, ω)
    # Pasamos a Hz reales
    f = ω .* fs ./ (2π)
    # Magnitud en dB y fase (convertida a grados sexagesimales)
    mag_db   = 20 .* log10.(abs.(H))
    fase_rad = angle.(H)
    fase_deg = fase_rad .* (180 / π)
    # Aseguramos que la fase esté en el rango [-180°, 180°]
    fase_deg = mod.(fase_deg .+ 180, 360) .- 180
    # Plot de magnitud
    p_mag = plot(
        f, mag_db,
        xlabel = "Frecuencia [Hz]",
        ylabel = "Magnitud [dB]",
        title  = "$(title) - Magnitud",
        grid   = true,
        linewidth = 2,
        legend = false,
    )
    # Plot de fase
    p_phase = plot(
        f, fase_deg,
        xlabel = "Frecuencia [Hz]",
        ylabel = "Fase [degrees]",
        title  = "$(title) - Fase",
        grid   = true,
        linewidth = 2,
        ylim = (-180, 180),
        yticks = -180:90:180,
        legend = false,
    )
    # Combinamos ambos plots
    p_combined = plot(p_mag, p_phase, layout = (2, 1), size = (800, 600))
    
    return p_combined, p_mag, p_phase
end

"""Calcula la PSD promedio de un diccionario de señales EEG (sin plotear)"""
function calculate_PSD_average(dict_EEG_data::Dict, fs::Real)
    channels = collect(keys(dict_EEG_data))
    n_muestras = length(dict_EEG_data[first(channels)])
    
    # Calcular PSD para cada canal
    PSD = Dict(channel => begin
        p = welch_pgram(dict_EEG_data[channel]; fs=fs, window=hamming, nfft=n_muestras)
        (; freq = DSP.freq(p), power = DSP.power(p))
    end for channel in channels)
    
    # Calcular PSD promedio
    freqs_avg = PSD[first(channels)].freq
    avg_power = mapreduce(ch -> PSD[ch].power, +, channels) ./ length(channels)
    
    return freqs_avg, avg_power, PSD
end

"""Calcula y visualiza la PSD (Power Spectral Density) para un diccionario de señales EEG"""
function calculate_and_plot_PSD(dict_EEG_data::Dict, fs::Real; title_prefix::String = "PSD")
    println("📈 CÁLCULO DE PSD POR CANAL")
    println("-" ^ 80)
    println()
    
    # Calcular PSD promedio
    freqs_avg, avg_power, PSD = calculate_PSD_average(dict_EEG_data, fs)
    
    println("  → PSD promedio de todos los canales")
    
    power_min, power_max = extrema(avg_power)
    y_min_log = floor(log10(power_min))
    y_max_log = ceil(log10(power_max))
    ylim_log = (10.0^y_min_log, 10.0^y_max_log)
    
    p_psd_avg = plot(
        xlabel = "Frecuencia (Hz)",
        ylabel = "Potencia (µV²/Hz)",
        title = "$(title_prefix) - Promedio de todos los canales",
        legend = false,
        xlim = (0, fs / 2),
        yscale = :log10,
        ylim = ylim_log,
    )
    
    plot!(p_psd_avg, freqs_avg, avg_power; label = "Promedio", lw = 2, color = :black)
    
    display(p_psd_avg)
    println()
    
    return PSD, p_psd_avg
end

"""Compara y grafica superpuestas las PSD promedio de dos diccionarios de señales EEG"""
function compare_PSD_averages(dict_EEG_1::Dict, dict_EEG_2::Dict, fs::Real;
                              label_1::String = "Original",
                              label_2::String = "Filtrado",
                              title::String = "Comparación PSD Promedio")
    println("📈 COMPARACIÓN DE PSD PROMEDIO")
    println()
    
    # Calcular PSD promedio de ambos diccionarios
    freqs_1, avg_power_1, _ = calculate_PSD_average(dict_EEG_1, fs)
    freqs_2, avg_power_2, _ = calculate_PSD_average(dict_EEG_2, fs)
    
    # Calcular límites Y para ambos conjuntos de datos
    all_powers = vcat(avg_power_1, avg_power_2)
    power_min, power_max = extrema(all_powers)
    y_min_log = floor(log10(power_min))
    y_max_log = ceil(log10(power_max))
    ylim_log = (10.0^y_min_log, 10.0^y_max_log)
    
    # Crear plot con ambos promedios superpuestos
    p_comparison = plot(
        xlabel = "Frecuencia (Hz)",
        ylabel = "Potencia (µV²/Hz)",
        title = title,
        legend = :topright,
        xlim = (0, fs / 2),
        yscale = :log10,
        ylim = ylim_log,
    )
    
    plot!(p_comparison, freqs_1, avg_power_1; label = label_1, lw = 2, color = :blue)
    plot!(p_comparison, freqs_2, avg_power_2; label = label_2, lw = 2, color = :red)
    
    display(p_comparison)
    println()
    
    return p_comparison
end

"""Guarda un diccionario de señales filtradas en un archivo binario"""
function save_filtered_signals(dict_EEG_filtered::Dict{String, Vector{Float64}}, filename::String, base_dir::String = dir_filtering)
    # Crear el directorio si no existe
    if !isdir(base_dir)
        mkpath(base_dir)
        println("📁 Directorio creado: $(abspath(base_dir))")
    end
    # Construir el path completo
    save_path = joinpath(base_dir, filename)
    
    # Guardar el archivo
    Serialization.serialize(save_path, dict_EEG_filtered)
    println("✓ Diccionario EEG filtrado guardado en: $(abspath(save_path))")
    println()
end

# ------------------------------------------------------------------------------------
# 4. FILTRO NOTCH (50 Hz)
# ------------------------------------------------------------------------------------
# Elimina la interferencia de la red eléctrica (50 Hz en Europa).

println("📊 Filtro Notch (50 Hz)")
println("-" ^ 40)

Notch_cutoff = 50
Notch_order = 4
Notch_width = 1.0

println("Frecuencia de corte: $(Notch_cutoff) Hz")
println("Orden del filtro: $(Notch_order)")
println("Ancho de banda: $(Notch_width) Hz")
println()

# Inicializamos el diccionario para las señales filtradas
dict_EEG_Notch = Dict{String, Vector{Float64}}()

# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_notch = notch_filter(dict_EEG[first(keys(dict_EEG))], Notch_cutoff, fs, Notch_order, Notch_width)

# Aplicamos el filtro Notch a todos los canales
for (channel, signal) in dict_EEG
    y_notch, _ = notch_filter(signal, Notch_cutoff, fs, Notch_order, Notch_width)
    dict_EEG_Notch[channel] = y_notch
end

# Ploteamos la respuesta del filtro Notch
p_notch, _, _ = plot_filter_response(flt_notch, fs; title = "Filtro Notch (50 Hz)")
display(p_notch)

save_filtered_signals(dict_EEG_Notch, "dict_EEG_Notch.bin")

# Comparar PSD promedio del original vs Notch filtrado (superpuestos)
p_psd_comparison = compare_PSD_averages(
    dict_EEG, 
    dict_EEG_Notch, 
    fs;
    label_1 = "Original",
    label_2 = "Notch (50 Hz)",
    title = "Comparación PSD Promedio: Original vs Notch"
)

# ------------------------------------------------------------------------------------
# 5. FILTRO BANDREJECT (100 Hz)
# ------------------------------------------------------------------------------------
# Elimina posible interferencia específica en torno a 100 Hz (hardware).

println("📊 Filtro Bandreject (100 Hz)")
println("-" ^ 40)

Bandreject_freq = 100
Bandreject_order = 4
Bandreject_bandwidth = 1.0

println("Frecuencia central: $(Bandreject_freq) Hz")
println("Orden del filtro: $(Bandreject_order)")
println("Ancho de banda: $(Bandreject_bandwidth) Hz")
println()


println("Justificación Técnica:")
println()
println("Este filtro se incluyó de forma conservadora para comprobar si existía")
println("interferencia específica en torno a 100 Hz, posiblemente asociada al hardware.")
println("Su ancho es muy estrecho (1 Hz) para evitar una eliminación generalizada")
println("de componentes de alta frecuencia.")
println()

# Inicializamos el diccionario para las señales filtradas
dict_EEG_Bandreject = Dict{String, Vector{Float64}}()

# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_bandreject = bandreject_filter(dict_EEG_Notch[first(keys(dict_EEG_Notch))], Bandreject_freq, fs, Bandreject_order, Bandreject_bandwidth)

# Aplicamos el filtro Bandreject a todos los canales (usando dict_EEG_Notch como entrada)
println("Aplicando filtro Bandreject a todos los canales...")
for (channel, signal) in dict_EEG_Notch
    y_bandreject, _ = bandreject_filter(signal, Bandreject_freq, fs, Bandreject_order, Bandreject_bandwidth)
    dict_EEG_Bandreject[channel] = y_bandreject
    println("  ✓ Canal $(channel) filtrado")
end
println("✓ Filtro Bandreject aplicado a $(length(dict_EEG_Bandreject)) canales")
println()

# Ploteamos la respuesta del filtro Bandreject
p_bandreject, _, _ = plot_filter_response(flt_bandreject, fs; title = "Filtro Bandreject (100 Hz)")
display(p_bandreject)

save_filtered_signals(dict_EEG_Bandreject, "dict_EEG_Bandreject.bin")

# Comparar PSD promedio del Notch vs Bandreject filtrado (superpuestos)
p_psd_comparison_br = compare_PSD_averages(
    dict_EEG_Notch, 
    dict_EEG_Bandreject, 
    fs;
    label_1 = "Notch (50 Hz)",
    label_2 = "Notch + Bandreject (100 Hz)",
    title = "Comparación PSD Promedio: Notch vs Notch+Bandreject"
)

# ------------------------------------------------------------------------------------
# 6. FILTRO HIGHPASS (0.5 Hz)
# ------------------------------------------------------------------------------------
# Elimina deriva (drift) preservando actividad infra-lenta de interés clínico.

println("📊 Filtro Highpass (0.5 Hz)")
println("-" ^ 40)
Highpass_cutoff = 0.5
Highpass_order_design = 4          # Orden del diseño del filtro
Highpass_order_effective = 8       # Orden efectivo (4 × 2 con filtfilt)
println("Frecuencia de corte: $(Highpass_cutoff) Hz")
println("Orden del diseño del filtro: $(Highpass_order_design)")
println("Orden efectivo (con filtfilt): $(Highpass_order_effective)")
println()

println("Justificación Técnica:")
println()
println("Se eligió un HPF a 0.5 Hz como compromiso entre preservar señal fisiológica")
println("lenta y evitar deriva (drift) del registro. En patologías como Esclerosis")
println("Múltiple (EM), la actividad infra-lenta (<1 Hz) puede tener relevancia clínica")
println("(procesos inflamatorios / reorganización cortical). Si bien se sabe que ICA")
println("mejora con HPF ≥1 Hz, se optó por 0.5 Hz para no eliminar actividad lenta")
println("posiblemente informativa.")
println()

# Inicializamos el diccionario para las señales filtradas
dict_EEG_Highpass = Dict{String, Vector{Float64}}()

# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_highpass = highpass_filter(dict_EEG_Bandreject[first(keys(dict_EEG_Bandreject))], Highpass_cutoff, fs, Highpass_order_design)

# Aplicamos el filtro Highpass a todos los canales (usando dict_EEG_Bandreject como entrada)
println("Aplicando filtro Highpass a todos los canales...")
for (channel, signal) in dict_EEG_Bandreject
    # zero_phase = true por defecto, así que no necesitamos especificarlo explícitamente
    y_highpass, _ = highpass_filter(signal, Highpass_cutoff, fs, Highpass_order_design)
    dict_EEG_Highpass[channel] = y_highpass
    println("  ✓ Canal $(channel) filtrado")
end
println("✓ Filtro Highpass aplicado a $(length(dict_EEG_Highpass)) canales")
println()

# Ploteamos la respuesta del filtro Highpass
p_highpass, _, _ = plot_filter_response(flt_highpass, fs; title = "Filtro Highpass (0.5 Hz)")
display(p_highpass)

save_filtered_signals(dict_EEG_Highpass, "dict_EEG_Highpass.bin")

# Comparar PSD promedio del Bandreject vs Highpass filtrado (superpuestos)
p_psd_comparison_hp = compare_PSD_averages(
    dict_EEG_Bandreject, 
    dict_EEG_Highpass, 
    fs;
    label_1 = "Notch + Bandreject",
    label_2 = "Notch + Bandreject + Highpass (0.5 Hz)",
    title = "Comparación PSD Promedio: Bandreject vs Bandreject+Highpass"
)

# ------------------------------------------------------------------------------------
# 7. FILTRO LOWPASS (150 Hz)
# ------------------------------------------------------------------------------------
# Limita el ancho de banda superior; preserva gamma para análisis posteriores.

println("📊 Filtro Lowpass (150 Hz)")
println("-" ^ 40)
Lowpass_cutoff = 150
Lowpass_order_design = 4          # Orden del diseño del filtro
Lowpass_order_effective = 8       # Orden efectivo (4 × 2 con filtfilt)
println("Frecuencia de corte: $(Lowpass_cutoff) Hz")
println("Orden del diseño del filtro: $(Lowpass_order_design)")
println("Orden efectivo (con filtfilt): $(Lowpass_order_effective)")
println()

println("Justificación Técnica:")
println()
println("Se mantuvo un límite superior amplio (150 Hz) para preservar información antes")
println("de explorarla. Un corte inicial a 40–70 Hz habría impedido detectar posibles")
println("componentes de gamma alta, que en algunos estudios sobre EM se han asociado a")
println("mecanismos de plasticidad cortical. Además, un LPF demasiado cercano a la banda")
println("de interés puede crear artefactos de \"ringing\" en los bordes del filtro. La señal")
println("puede recortarse posteriormente sin problema si el análisis se enfoca en <50 Hz.")
println()

# Inicializamos el diccionario para las señales filtradas
dict_EEG_Lowpass = Dict{String, Vector{Float64}}()

# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_lowpass = lowpass_filter(dict_EEG_Highpass[first(keys(dict_EEG_Highpass))], Lowpass_cutoff, fs, Lowpass_order_design)

# Aplicamos el filtro Lowpass a todos los canales (usando dict_EEG_Highpass como entrada)
println("Aplicando filtro Lowpass a todos los canales...")
for (channel, signal) in dict_EEG_Highpass
    # zero_phase = true por defecto, así que no necesitamos especificarlo explícitamente
    y_lowpass, _ = lowpass_filter(signal, Lowpass_cutoff, fs, Lowpass_order_design)
    dict_EEG_Lowpass[channel] = y_lowpass
    println("  ✓ Canal $(channel) filtrado")
end
println("✓ Filtro Lowpass aplicado a $(length(dict_EEG_Lowpass)) canales")
println()

# Ploteamos la respuesta del filtro Lowpass
p_lowpass, _, _ = plot_filter_response(flt_lowpass, fs; title = "Filtro Lowpass (150 Hz)")
display(p_lowpass)

save_filtered_signals(dict_EEG_Lowpass, "dict_EEG_Lowpass.bin")

# Comparar PSD promedio del Highpass vs Lowpass filtrado (superpuestos)
p_psd_comparison_lp = compare_PSD_averages(
    dict_EEG_Highpass, 
    dict_EEG_Lowpass, 
    fs;
    label_1 = "Notch + Bandreject + Highpass",
    label_2 = "Notch + Bandreject + Highpass + Lowpass (150 Hz)",
    title = "Comparación PSD Promedio: Highpass vs Highpass+Lowpass"
)

