#!/usr/bin/env julia

# -*- coding: utf-8 -*-
# src/Setup/IO.jl
#
# CARGA Y PREPROCESAMIENTO INICIAL DE DATOS EEG
# ==============================================
# Esta rutina realiza la carga inicial de datos EEG desde archivos TSV y prepara
# los datos para el procesamiento posterior en el pipeline de análisis.
#
# PROCESO:
#   1. Carga datos raw desde archivo TSV (formato BrainVision/EEGLAB compatible)
#   2. Extrae información temporal (frecuencia de muestreo, duración)
#   3. Organiza datos en estructura de diccionario (canal → señal)
#   4. Guarda datos organizados en formato binario para pasos posteriores
#   5. Genera visualizaciones exploratorias en el dominio temporal:
#      - Señales temporales (canal individual y apiladas)
#   6. Análisis espectral: densidad de potencia espectral (PSD)
#      6.1. Visualización de PSD: canales superpuestos
#      6.2. Definición de bandas de frecuencia EEG
#      6.3. Visualización de PSD promedio con bandas de frecuencia
#      6.4. Cálculo de potencia integrada por banda
#   7. Comparación con resultados de referencia (BrainVision Analyzer)
#   8. Análisis estadístico de calidad de canales:
#      8.1. Visualizaciones de calidad de canales
#      8.2. Detección de canales sospechosos
#
# NOTA: Este es el primer paso del pipeline y proporciona una visión general
# de la calidad y características de los datos antes de aplicar filtros.

using CSV
using DataFrames
using Plots         # El backend GR se inicializa automáticamente al crear el primer plot
using Dates 
using DSP
using Serialization
using Statistics
using StatsBase
using StatsPlots
using FFTW

# Si este archivo se ejecuta directamente en Main, recarga utilidades de rutas
# para evitar definiciones antiguas en la sesión REPL.
if (@__MODULE__) == Main
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end

# ------------------------------------------------------------------------------------
# 1. CARGA DE DATOS RAW DESDE ARCHIVO TSV
# ------------------------------------------------------------------------------------
# Se carga el archivo TSV que contiene los datos EEG raw. El formato esperado es:
# - Primera columna: nombres de canales (Channel)
# - Columnas siguientes: muestras temporales (una columna por punto de tiempo)
# - Los datos están en unidades de microvoltios (µV)
#
# Este archivo es la entrada principal del pipeline y contiene los datos sin
# ningún procesamiento previo.

println("=" ^ 80)
println("📊 CARGA DE DATOS EEG")
println("=" ^ 80)

# Ruta al archivo de datos raw
raw_filename = "sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv"
cfg_local = @isdefined(cfg) ? cfg : nothing
dir_raw = joinpath(bids_root(cfg_local), "raw", raw_filename)

if !isfile(dir_raw)
    error("No se encontró archivo EEG raw: $(abspath(dir_raw)). Revisa cfg.data_dir y la estructura BIDS (raw_dir).")
end

# Leer el archivo TSV como DataFrame
# CSV.read carga el archivo y lo convierte en una estructura tabular
data_raw = CSV.read(dir_raw, DataFrame)   

println("✓ Archivo cargado: $(basename(dir_raw))")
println("✓ Dimensiones: $(size(data_raw, 1)) canales × $(size(data_raw, 2) - 1) muestras")
println("✓ Tabla de datos:")
display(data_raw)
println()

# ------------------------------------------------------------------------------------
# 2. INFORMACIÓN TEMPORAL DEL REGISTRO
# ------------------------------------------------------------------------------------
# Se extrae información sobre las características temporales del registro:
# - Frecuencia de muestreo (fs): número de muestras por segundo
# - Número de muestras: longitud temporal del registro
# - Duración total: tiempo total de registro en segundos y minutos
#
# Esta información es esencial para:
# - Cálculo de frecuencias en análisis espectral
# - Aplicación de filtros (requieren fs)
# - Segmentación temporal
# - Visualización de señales en el dominio temporal

println("⏱️  INFORMACIÓN TEMPORAL")
println("-" ^ 80)

# Frecuencia de muestreo (Hz)
# Define la resolución temporal: cuántas muestras se toman por segundo
fs = (cfg_local !== nothing && hasproperty(cfg_local, :filter) && hasproperty(cfg_local.filter, :fs)) ?
    cfg_local.filter.fs : 500

# Número de muestras temporales
# Se resta 1 porque la primera columna es "Channel" (nombres de canales)
n_muestras = size(data_raw, 2) - 1

# Vector de tiempo en segundos
# Se genera desde 0 hasta (n_muestras-1) y se divide por fs para obtener segundos
tiempo_seg = collect(0:(n_muestras-1)) ./ fs

# Duración total del registro
duracion_total = tiempo_seg[end]

println("  Frecuencia de muestreo: $(fs) Hz")
println("  Número de muestras: $(n_muestras)")
println("  Duración total: $(round(duracion_total, digits=2)) segundos ($(round(duracion_total/60, digits=2)) minutos)")
println()

# ------------------------------------------------------------------------------------
# 3. EXTRACCIÓN Y ORGANIZACIÓN DE CANALES
# ------------------------------------------------------------------------------------
# Se extraen los nombres de canales y se organizan los datos en estructuras
# convenientes para el procesamiento:
# - channels: vector con nombres de canales
# - EEG_matrix: matriz (canales × muestras) con todos los datos
# - dict_EEG: diccionario que mapea nombre de canal → señal (vector)
#
# El diccionario es especialmente útil porque permite acceso por nombre de canal
# y facilita el procesamiento posterior en el pipeline.

println("📡 INFORMACIÓN DE CANALES")
println("-" ^ 80)

# Nombres de canales (primera columna del DataFrame)
channels = data_raw.Channel

# Matriz de datos EEG (canales × muestras)
# Se excluye la columna "Channel" usando Not(:Channel)
# Cada fila corresponde a un canal, cada columna a una muestra temporal
EEG_matrix = Matrix(data_raw[:, Not(:Channel)])

println("  Número de canales: $(length(channels))")
println("  Canales: $(join(channels, ", "))")
println()

# Diccionario: canal → señal (vector de muestras)
# Esta estructura facilita el acceso por nombre de canal y es compatible
# con el resto del pipeline que espera datos organizados de esta manera
dict_EEG = Dict(data_raw.Channel[i] => EEG_matrix[i, :]
                for i in 1:size(data_raw, 1))

# ------------------------------------------------------------------------------------
# 4. GUARDADO DE DATOS ORGANIZADOS
# ------------------------------------------------------------------------------------
# Se guarda el diccionario de canales en formato binario nativo de Julia.
# Este formato es eficiente para matrices grandes y permite serialización
# rápida de estructuras complejas.
#
# El archivo guardado será la entrada para el siguiente paso del pipeline
# (filtrado, corrección de línea base, etc.)

# Directorio de salida para datos IO
dir_io = stage_dir(:IO; cfg = cfg_local)
path_dict = joinpath(dir_io, "dict_EEG.bin")

# Asegurar que el directorio existe
isdir(dir_io) || mkpath(dir_io)

# Serializar y guardar el diccionario
Serialization.serialize(path_dict, dict_EEG)

println("💾 Diccionario EEG raw guardado en: $(abspath(path_dict))")

# ------------------------------------------------------------------------------------
# 5. VISUALIZACIONES EXPLORATORIAS EN EL DOMINIO TEMPORAL
# ------------------------------------------------------------------------------------
# Se generan gráficos para inspección visual de los datos raw:
# - Gráfico de un canal individual (ejemplo: Cz)
# - Gráfico de todos los canales apilados (permite comparar señales)
#
# Estas visualizaciones ayudan a identificar:
# - Artefactos obvios (parpadeos, movimiento, ruido de línea)
# - Diferencias entre canales
# - Calidad general del registro

println("📊 GENERANDO VISUALIZACIONES")
println("-" ^ 80)

# Gráfico de un canal individual como ejemplo
# El canal Cz (central) suele ser representativo de la actividad EEG general
println("  → Gráfico del canal Cz (ejemplo)")
p_cz = plot(tiempo_seg, dict_EEG["Cz"], 
     xlabel = "Tiempo (s)", ylabel= "Amplitud (µV)", title="EEG - Canal Cz", legend = false)
display(p_cz)
println()

# Representación de todos los canales apilados
# Cada canal se desplaza verticalmente para facilitar la comparación visual
# Esta visualización permite detectar artefactos que afectan a múltiples canales
println("  → Gráfico de todos los canales apilados")

# Separación vertical entre canales (en µV)
# Un valor mayor aumenta la separación visual pero puede hacer que los detalles
# de amplitud sean menos visibles
sep = 80.0

# Calcular offsets (posiciones Y) para cada canal
# Los canales se apilan de arriba a abajo, con el primer canal en la parte superior
offsets = [(length(channels) - i) * sep for i in 1:length(channels)]

# Gráfico apilado
# Se transpone EEG_matrix para que cada columna sea un canal
# Se suman los offsets para desplazar verticalmente cada señal
p = plot(
    tiempo_seg,
    EEG_matrix' .+ offsets',  # transpuesta + offsets para apilar
    xlabel = "Tiempo (s)",
    ylabel = "",
    yticks = (offsets, channels),  # etiquetas Y con nombres de canales
    legend = false,
    grid = false,
    size = (1000, 600)
)

display(p)
println()

# ------------------------------------------------------------------------------------
# 6. ANÁLISIS ESPECTRAL: DENSIDAD DE POTENCIA ESPECTRAL (PSD)
# ------------------------------------------------------------------------------------
# Se calcula la PSD usando el método de Welch, que divide la señal en segmentos
# superpuestos, aplica ventanas (Hamming) y promedia los periodogramas.
#
# La PSD muestra cómo se distribuye la potencia de la señal en el dominio
# frecuencial, lo cual es esencial para:
# - Identificar bandas de frecuencia dominantes (Alpha, Beta, etc.)
# - Detectar ruido de línea (50/60 Hz)
# - Evaluar la calidad espectral de cada canal
# - Comparar características espectrales entre canales

println("📈 CÁLCULO DE PSD POR CANAL")
println("-" ^ 80)
println()

# Calcular PSD para cada canal usando método de Welch
# welch_pgram: divide la señal en segmentos, aplica ventana Hamming,
# calcula periodograma de cada segmento y promedia
# - fs: frecuencia de muestreo (necesaria para escalar frecuencias)
# - window: tipo de ventana (Hamming reduce leakage espectral)
# - nfft: tamaño de FFT (usar n_muestras para máxima resolución)
PSD = Dict(channel => begin
    p = welch_pgram(dict_EEG[channel]; fs=fs, window=hamming, nfft=n_muestras)
    # Extraer frecuencias y potencias del periodograma
    (; freq = DSP.freq(p), power = DSP.power(p))
end for channel in channels)

display(PSD)
println()

# ------------------------------------------------------------------------------------
# 6.1. VISUALIZACIÓN DE PSD: CANALES SUPERPUESTOS
# ------------------------------------------------------------------------------------
# Se superponen todas las PSD para comparar características espectrales entre canales.
# La escala logarítmica facilita visualizar rangos amplios de potencia.

println("  → PSD de todos los canales superpuestas")

# Calcular límites Y amigables para escala logarítmica
# Se buscan potencias de 10 que enmarquen todos los datos
all_powers = vcat([PSD[ch].power for ch in channels]...)
all_powers_pos = filter(>(0), all_powers)
if isempty(all_powers_pos)
    error("No hay potencias positivas para graficar en escala log10.")
end
power_min_all, power_max_all = extrema(all_powers_pos)
if power_max_all <= power_min_all
    power_max_all = power_min_all * 10
end
y_min_log_all = floor(log10(power_min_all))
y_max_log_all = ceil(log10(power_max_all))
ylim_log_all = (10.0^y_min_log_all, 10.0^y_max_log_all)

# Gráfico de PSD superpuestas
p_psd = plot(
    xlabel = "Frecuencia (Hz)",
    ylabel = "Potencia (µV²/Hz)",
    title = "PSD por canal (Welch)",
    legend = :right,
    xlim = (0, fs / 2),  # Frecuencia de Nyquist = fs/2
    yscale = :log10,     # Escala logarítmica para mejor visualización
    ylim = ylim_log_all,
)

# Añadir cada canal al gráfico
for channel in channels
    plot!(p_psd, PSD[channel].freq, PSD[channel].power; label = channel, lw = 1)
end

display(p_psd)
println()

# ------------------------------------------------------------------------------------
# 6.2. DEFINICIÓN DE BANDAS DE FRECUENCIA EEG
# ------------------------------------------------------------------------------------
# Se definen las bandas de frecuencia estándar en EEG:
# - Delta (Δ): 0.5-4 Hz - sueño profundo, patologías
# - Theta (δ): 4-8 Hz - somnolencia, meditación
# - Alpha (α): 8-12 Hz - relajación, ojos cerrados
# - Beta Low: 12-15 Hz - actividad normal despierto
# - Beta Medium: 15-18 Hz - actividad normal despierto
# - Beta High: 18-30 Hz - concentración, ansiedad
# - Gamma (γ): 30-50 Hz - procesamiento cognitivo, artefactos musculares
#
# Estas bandas se usan para cuantificar la potencia espectral en rangos
# fisiológicamente relevantes.

# Nomenclatura de bandas alineada con Pluto/BIDS/BIDS.jl
Δ = (0.5, 4.0)
δ = (4.0, 8.0)
α = (8.0, 12.0)
β_low = (12.0, 15.0)
β_medium = (15.0, 18.0)
β_high = (18.0, 30.0)
γ = (30.0, 50.0)

band_limits = Dict(:Δ => Δ, :δ => δ, :α => α, :β_low => β_low, :β_medium => β_medium, :β_high => β_high, :γ => γ)
ordered_bands = [:Δ, :δ, :α, :β_low, :β_medium, :β_high, :γ]
band_labels = Dict(
    :Δ => "DELTA (0.5-4)", :δ => "THETA (4-8)", :α => "ALPHA (8-12)",
    :β_low => "BETA LOW (12-15)", :β_medium => "BETA MID (15-18)",
    :β_high => "BETA HIGH (18-30)", :γ => "GAMMA (30-50)"
)

# Colores para visualización de bandas en gráficos
band_colors = [:lightcyan, :lavender, :lightgoldenrod, :lightgreen, :lightsalmon, :lightpink, :lightgray]

# ------------------------------------------------------------------------------------
# 6.3. VISUALIZACIÓN DE PSD PROMEDIO CON BANDAS DE FRECUENCIA
# ------------------------------------------------------------------------------------
# Se calcula la PSD promedio de todos los canales y se visualiza con las bandas
# de frecuencia marcadas. Esto proporciona una visión general del contenido
# espectral del registro y ayuda a identificar qué bandas dominan.

println("  → PSD promedio de todos los canales")

# Crear gráfico base
p_psd_avg = plot(
    xlabel = "Frecuencia (Hz)",
    ylabel = "Potencia (µV²/Hz)",
    title = "PSD promedio de todos los canales",
    legend = :outerright,
    xlim = (0, fs / 2),
    yscale = :log10,
)

# Calcular PSD promedio
# Se asume que todos los canales tienen las mismas frecuencias (mismo nfft)
freqs_avg = PSD[first(channels)].freq

# Promediar potencias de todos los canales
# mapreduce suma todas las potencias y luego divide por el número de canales
avg_power = mapreduce(ch -> PSD[ch].power, +, channels) ./ length(channels)

# Calcular límites de potencia para escala Y
avg_power_pos = filter(>(0), avg_power)
if isempty(avg_power_pos)
    error("La PSD promedio no contiene valores positivos para escala log10.")
end
power_min, power_max = extrema(avg_power_pos)
if power_max <= power_min
    power_max = power_min * 10
end

# Calcular límites Y amigables para escala logarítmica (potencias de 10)
y_min_log = floor(log10(power_min))
y_max_log = ceil(log10(power_max))
ylim_log = (10.0^y_min_log, 10.0^y_max_log)

# Actualizar el gráfico con límites Y explícitos
plot!(p_psd_avg, ylim = ylim_log)

# Añadir regiones sombreadas para cada banda de frecuencia
# fillrange crea un área sombreada entre límites positivos en cada banda
for (idx, band) in enumerate(ordered_bands)
    limits = band_limits[band]
    plot!(
        p_psd_avg,
        [limits[1], limits[2]],      # límites de frecuencia en X
        [power_max, power_max];      # altura superior del sombreado
        fillrange = power_min,       # altura inferior del sombreado
        fillcolor = band_colors[idx], # color de la banda
        fillalpha = 0.25,            # transparencia
        linecolor = :transparent,    # sin borde
        label = band_labels[band],   # etiqueta para leyenda
    )
end

# Añadir línea de PSD promedio
plot!(p_psd_avg, freqs_avg, avg_power; label = "Promedio", lw = 2, color = :black)

display(p_psd_avg)
println()

# ------------------------------------------------------------------------------------
# 6.4. CÁLCULO DE POTENCIA INTEGRADA POR BANDA
# ------------------------------------------------------------------------------------
# Se calcula la potencia total en cada banda de frecuencia para cada canal.
# La potencia se integra sumando las densidades espectrales dentro de cada banda
# y multiplicando por la resolución espectral (df).
#
# Esto permite cuantificar:
# - Qué bandas dominan en cada canal
# - Distribución relativa de potencia entre bandas
# - Comparaciones entre canales

println("  → Potencia integrada por banda y canal")

# Calcular potencia integrada por banda para cada canal
# La integración se hace sumando las densidades de potencia en cada banda
# y multiplicando por la resolución espectral (df = diferencia entre frecuencias)
band_powers = Dict(channel => begin
    freqs = PSD[channel].freq
    powers = PSD[channel].power
    
    # Resolución espectral: diferencia entre frecuencias adyacentes
    # Esto es necesario para convertir densidad (µV²/Hz) en potencia total (µV²)
    df = freqs[2] - freqs[1]
    
    # Para cada banda, sumar las potencias dentro del rango de frecuencias
    Dict(
        band => begin
            # Índices de frecuencias dentro de la banda [f_min, f_max)
            idx = (freqs .>= limits[1]) .& (freqs .< limits[2])
            # Potencia integrada = suma de densidades × resolución
            sum(powers[idx]) * df
        end
        for (band, limits) in band_limits
    )
end for channel in channels)

# Crear DataFrame con potencias por banda y canal
band_powers_df = DataFrame(Band = [band_labels[band] for band in ordered_bands])
channel_totals = Dict{String, Float64}()

# Añadir columnas para cada canal (potencia absoluta y porcentaje)
for channel in channels
    values = [band_powers[channel][band] for band in ordered_bands]
    total_power = sum(values)
    channel_totals[channel] = total_power
    
    # Potencia absoluta en µV²
    band_powers_df[!, Symbol(channel)] = values
    
    # Potencia relativa (porcentaje del total)
    # Evita división por cero si total_power es muy pequeño
    band_powers_df[!, Symbol(channel * "_pct")] = total_power ≈ 0 ? zeros(length(values)) : (values ./ total_power)
end

# Añadir fila de totales
total_row_pairs = Any[:Band => "TOTAL POWER"]
for channel in channels
    push!(total_row_pairs, Symbol(channel) => channel_totals[channel])
    push!(total_row_pairs, Symbol(channel * "_pct") => 1.0)
end
total_row = (; total_row_pairs...)

# Combinar fila de totales con el resto del DataFrame
band_powers_df = vcat(DataFrame([total_row]), band_powers_df; cols = :union)

display(band_powers_df)
println()

# ------------------------------------------------------------------------------------
# 7. COMPARACIÓN CON RESULTADOS DE REFERENCIA (BrainVision Analyzer)
# ------------------------------------------------------------------------------------
# Se cargan resultados espectrales obtenidos con BrainVision Analyzer (software
# comercial de análisis EEG) para comparación y validación.
#
# NOTA: Los valores de referencia pueden representar:
# - Mean Activity: densidad media de potencia en la banda (µV²/Hz)
# - Area: potencia integrada en la banda (µV²)
#
# La comparación permite verificar que el procesamiento en Julia produce
# resultados consistentes con herramientas estándar del campo.
# "real" (µV²/Hz) o un simple "Area" (µV²) para la potencia total.

println("📊 Lectura datos espectrales de JAVIER")
println("-" ^ 80)

# Ruta al archivo de resultados de referencia
path_candidates = [
    joinpath(project_root(), "Javier_results", "M5T2cerrados.txt"),      # layout actual (raíz del proyecto)
    joinpath(@__DIR__, "..", "Javier_results", "M5T2cerrados.txt"),      # compatibilidad legacy
]
existing_paths = filter(isfile, path_candidates)
path_javier = isempty(existing_paths) ? "" : first(existing_paths)
if isempty(path_javier)
    error("No se encontró M5T2cerrados.txt en: $(join(path_candidates, " | "))")
end

# Leer líneas del archivo y filtrar líneas vacías
lines_javier = filter(x -> !isempty(strip(x)), readlines(path_javier))

# Parsear encabezado (primera línea)
# Se separa por espacios múltiples (2 o más)
header = split(strip(lines_javier[1]), r"\s{2,}")
channels_javier = header[2:end]  # Primera columna es "Band", resto son canales

# Parsear filas de datos
rows_javier = [split(strip(line), r"\s{2,}") for line in lines_javier[2:end]]
bands_javier = first.(rows_javier)  # Primera columna de cada fila es el nombre de banda

# Parsear números (convertir comas a puntos para formato decimal)
numbers = [
    parse(Float64, replace(rows_javier[i][j], ',' => '.'))
    for i in eachindex(rows_javier), j in 2:length(header)
]

# Construir DataFrame con resultados de referencia
Javier_results = DataFrame(Band = bands_javier)
for (col_idx, channel) in enumerate(channels_javier)
    Javier_results[!, Symbol(channel)] = numbers[:, col_idx]
end

println("  → Datos espectrales de JAVIER")
display(Javier_results)
println()

# ------------------------------------------------------------------------------------
# 8. ANÁLISIS ESTADÍSTICO DE CALIDAD DE CANALES
# ------------------------------------------------------------------------------------
# Se calculan estadísticas descriptivas para cada canal que permiten evaluar
# la calidad de los datos y detectar canales problemáticos:
# - Mean: valor medio (indica offset DC)
# - Min/Max: rango de amplitudes
# - Std: desviación estándar (variabilidad, ruido)
# - RMS: raíz cuadrada de la media de cuadrados (potencia promedio)
# - Kurtosis: medida de "picos" (valores altos indican outliers/artefactos)
# - Skewness: asimetría de la distribución
#
# Canales con kurtosis muy alta (>6) o desviación estándar muy baja (<5 µV)
# pueden indicar problemas (artefactos, canales muertos, etc.)

println("📈 CÁLCULO DE ESTADÍSTICAS POR CANAL")
println("-" ^ 80)
println()

# Calcular estadísticas para cada canal
# eachrow itera sobre cada fila de la matriz (cada canal)
stats_channels = DataFrame(
    :Channel                 => channels,
    Symbol("Mean (µV)")     => mean.(eachrow(EEG_matrix)),      # Media
    Symbol("Min (µV)")      => minimum.(eachrow(EEG_matrix)),   # Mínimo
    Symbol("Max (µV)")      => maximum.(eachrow(EEG_matrix)),  # Máximo
    Symbol("Std (µV)")      => std.(eachrow(EEG_matrix)),       # Desviación estándar
    Symbol("RMS (µV)")      => rms.(eachrow(EEG_matrix)),       # RMS (root mean square)
    :Kurtosis                => kurtosis.(eachrow(EEG_matrix)), # Curtosis (medida de picos)
    :Skewness                => skewness.(eachrow(EEG_matrix)), # Asimetría
)

# ------------------------------------------------------------------------------------
# 8.1. VISUALIZACIONES DE CALIDAD DE CANALES
# ------------------------------------------------------------------------------------
# Se generan histogramas y gráficos de dispersión para identificar canales
# con características anómalas que puedan requerir atención especial.

# Histograma de kurtosis
# La kurtosis mide qué tan "picuda" es la distribución. Valores altos (>6)
# pueden indicar presencia de artefactos (picos, spikes) o canales problemáticos
println("  → Histograma de kurtosis por canal")
histogram(
    stats_channels.Kurtosis,
    bins = 20,
    xlabel = "Kurtosis",
    ylabel = "Número de canales",
    legend = false,
    title  = "Distribución de kurtosis por canal"
)

# Histograma de desviación estándar
# La desviación estándar mide la variabilidad. Valores muy bajos (<5 µV)
# pueden indicar canales muertos o con muy poca señal
println("  → Histograma de desviación estándar por canal")
histogram(
    stats_channels[!, Symbol("Std (µV)")],
    bins   = 20,
    xlabel = "Std (µV)",
    ylabel = "Número de canales",
    legend = false,
    title  = "Distribución de desviación típica por canal"
)   

# ------------------------------------------------------------------------------------
# 8.2. DETECCIÓN DE CANALES SOSPECHOSOS
# ------------------------------------------------------------------------------------
# Se usa un gráfico de dispersión (Std vs Kurtosis) para identificar canales
# que se desvían de la distribución normal, lo cual puede indicar problemas.

# Extraer variables para el scatter plot
stds   = stats_channels[!, Symbol("Std (µV)")]
kurts  = stats_channels.Kurtosis
names_ = stats_channels.Channel

# Criterio de canales sospechosos:
# - Kurtosis > 6: distribución muy picuda (posibles artefactos)
# - Std < 5 µV: variabilidad muy baja (posible canal muerto)
idx_raros = (kurts .> 6) .| (stds .< 5)

# Colores: rojo para sospechosos, azul para normales
colors = ifelse.(idx_raros, :red, :blue)

# Crear scatter plot
p_scatter = scatter(
    stds, kurts;
    color  = colors,
    xlabel = "Std (µV)",
    ylabel = "Kurtosis",
    title  = "Calidad de canales EEG (Std vs Kurtosis)",
    legend = false,
)

# Añadir etiquetas SOLO a los canales sospechosos
# Se desplazan ligeramente para mejor legibilidad
offsetx = 1.5   # desplazamiento en X
offsety = 0.2   # desplazamiento en Y

for i in eachindex(names_)
    if idx_raros[i]
        annotate!(
            p,
            stds[i] + offsetx,
            kurts[i] + offsety,
            Plots.text(string(names_[i]), 6, :red)
        )
    end
end

println("  → Scatter plot: Std vs Kurtosis (calidad de canales)")
display(p_scatter)

return path_dict