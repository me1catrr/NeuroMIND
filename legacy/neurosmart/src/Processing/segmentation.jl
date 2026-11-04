#!/usr/bin/env julia
# -*- coding: utf-8 -*- #shebang        
# src/segmentation.jl

# SEGMENTACIÓN DE SEÑALES EEG
# ============================
# Esta rutina divide la señal EEG continua en segmentos temporales de longitud fija.
# Los datos ya han sido limpiados mediante ICA y están listos para ser segmentados.
#
# PROCESO:
#   1. Carga datos continuos después de ICA cleaning
#   2. Configura parámetros de segmentación (longitud, solapamiento)
#   3. Divide la señal en segmentos consecutivos de longitud fija
#   4. Crea una hipermatriz 3D: (canales × muestras × segmentos)
#   5. Calcula estadísticas por segmento
#   6. Genera visualizaciones de ejemplo
#   7. Guarda hipermatriz segmentada y metadatos
#
# CONFIGURACIÓN:
#   - Longitud de segmento: 1.00 s (500 muestras a 500 Hz)
#   - Solapamiento: 0.00 s (segmentos consecutivos sin solapamiento)
#   - Skip bad intervals: No (no se implementa actualmente)
#
# Justificación Técnica:
# Se usaron segmentos de 1 s para favorecer la resolución espectral posterior
# (resolución ~1 Hz = fs/longitud_segmento). No se aplicó solapamiento para
# facilitar el cómputo de FFT y SNR por segmento, y para mantener segmentos
# estadísticamente independientes. En resting state con análisis espectral,
# este enfoque es válido y compatible con el método clásico de Welch si
# posteriormente se aplica promediado espectral. La hipermatriz resultante
# permite análisis por segmento (p.ej., artifact rejection) y promediado posterior.

# Listado de librerías
using Serialization
using CSV, DataFrames
using LinearAlgebra
using Statistics
using StatsBase
using Plots
# El backend GR se inicializa automáticamente al crear el primer plot

# Si se ejecuta este script directamente (fuera del módulo EEG_Julia),
# cargamos utilidades de rutas para disponer de `stage_dir`.
if !@isdefined(stage_dir)
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end

# -----------------------------------------------------------------------------
# 1. Carga de datos después de ICA cleaning
# -----------------------------------------------------------------------------
# Se cargan los datos EEG que ya han sido limpiados mediante ICA.
# Los datos están en formato de diccionario (un vector por canal) y representan
# señales continuas que serán divididas en segmentos temporales.

println("=" ^ 80)
println("📊 SEGMENTACIÓN DE SEÑALES EEG")
println("=" ^ 80)
println()

# Cargar datos desde el paso anterior (ICA cleaning)
# Estos datos ya están libres de artefactos y listos para segmentar
dir_ica_cleaning = stage_dir(:ICA)
path_dict_ica_cleaning = joinpath(dir_ica_cleaning, "dict_EEG_ICA_clean.bin")
dict_EEG_ICA_clean = Serialization.deserialize(path_dict_ica_cleaning)

println("✔ Datos de ICA cleaning cargados desde: $(basename(path_dict_ica_cleaning))")

# Obtener canales ordenados (orden alfabético para consistencia)
# Todos los canales deben tener la misma longitud (n_samples_total)
channels = sort(collect(keys(dict_EEG_ICA_clean)))
n_channels = length(channels)
n_samples_total = length(dict_EEG_ICA_clean[channels[1]])  # longitud de la señal continua

println("  → Número de canales: $n_channels")
println("  → Número total de muestras: $n_samples_total")
println()

# -----------------------------------------------------------------------------
# 2. Configuración de segmentación
# -----------------------------------------------------------------------------
# Se configuran los parámetros para dividir la señal continua en segmentos:
#   - Longitud del segmento: determina la duración de cada época
#   - Solapamiento: permite que segmentos se superpongan (0 = sin solapamiento)
#   - Paso entre segmentos: distancia entre el inicio de segmentos consecutivos
#
# CONFIGURACIÓN ACTUAL:
#   - Segmentos de 1 segundo (500 muestras a 500 Hz)
#   - Sin solapamiento (segmentos consecutivos e independientes)
#   - Resolución espectral resultante: ~1 Hz (fs / longitud_segmento)

fs = 500.0                    # Frecuencia de muestreo (Hz)
segment_length_s = 1.00       # Longitud del segmento en segundos
                               # Segmentos de 1 s son comunes para análisis espectral
segment_overlap_s = 0.00      # Solapamiento en segundos
                               # 0 = sin solapamiento (segmentos independientes)
skip_bad_intervals = false    # No saltar intervalos malos
                               # NOTA: Actualmente no se implementa esta funcionalidad

# Convertir tiempos a muestras
segment_length_samples = Int(round(segment_length_s * fs))  # 500 muestras
segment_overlap_samples = Int(round(segment_overlap_s * fs))  # 0 muestras
segment_step_samples = segment_length_samples - segment_overlap_samples  # 500 muestras
                               # Paso entre segmentos (inicio de uno al siguiente)

println("⚙️  CONFIGURACIÓN DE SEGMENTACIÓN")
println("-" ^ 80)
println("  Frecuencia de muestreo: $(fs) Hz")
println("  Longitud del segmento: $(segment_length_s) s ($(segment_length_samples) muestras)")
println("  Solapamiento: $(segment_overlap_s) s ($(segment_overlap_samples) muestras)")
println("  Paso entre segmentos: $(segment_step_samples) muestras")
println("  Saltar intervalos malos: $(skip_bad_intervals)")
println()

# Calcular número de segmentos completos que caben en la señal
# Fórmula: número de segmentos = floor((total - longitud) / paso) + 1
# Esto asegura que todos los segmentos tengan la longitud completa
n_segments = Int(floor((n_samples_total - segment_length_samples) / segment_step_samples)) + 1

println("📐 CÁLCULO DE SEGMENTOS")
println("-" ^ 80)
println("  Número de segmentos completos: $n_segments")
println("  Muestras por segmento: $segment_length_samples")
println("  Muestras totales utilizadas: $(n_segments * segment_step_samples)")
# Las muestras finales que no forman un segmento completo se descartan
if n_segments * segment_step_samples < n_samples_total
    muestras_no_usadas = n_samples_total - (n_segments * segment_step_samples)
    println("  Muestras no utilizadas (final): $muestras_no_usadas")
end
println()

# -----------------------------------------------------------------------------
# 3. Segmentación: Crear hipermatriz canales × muestras × segmentos
# -----------------------------------------------------------------------------
# Esta sección divide cada canal en segmentos consecutivos y los organiza en
# una hipermatriz 3D. Cada "rebanada" de la tercera dimensión es un segmento
# temporal completo con todos los canales.
#
# PROCESO:
#   Para cada canal:
#     1. Extrae la señal continua del diccionario
#     2. Calcula los índices de inicio y fin de cada segmento
#     3. Copia cada segmento a la hipermatriz
#     4. Los segmentos son consecutivos (sin solapamiento)
#
# RESULTADO:
#   Hipermatriz de tamaño (n_channels × segment_length_samples × n_segments)
#   - Primera dimensión: canales
#   - Segunda dimensión: muestras temporales dentro del segmento
#   - Tercera dimensión: segmentos (épocas)

println("🔄 PROCESANDO SEGMENTACIÓN")
println("-" ^ 80)

# Inicializar hipermatriz: canales × muestras × segmentos
# Esta estructura permite acceso eficiente a segmentos individuales
eeg_segmented = zeros(Float64, n_channels, segment_length_samples, n_segments)

# Segmentar cada canal independientemente
for (ch_idx, channel) in enumerate(channels)
    signal = dict_EEG_ICA_clean[channel]  # señal continua del canal
    
    # Crear segmentos para este canal
    for seg_idx in 1:n_segments
        # Calcular índices del segmento en la señal continua
        # start_idx: inicio del segmento (1-indexed)
        # end_idx: fin del segmento
        start_idx = (seg_idx - 1) * segment_step_samples + 1
        end_idx = start_idx + segment_length_samples - 1
        
        # Verificar que no excedemos el límite de la señal
        # Normalmente no debería pasar porque n_segments se calculó correctamente
        if end_idx <= length(signal)
            # Segmento completo: copiar directamente
            eeg_segmented[ch_idx, :, seg_idx] = signal[start_idx:end_idx]
        else
            # Si el último segmento excede (caso raro), truncar o rellenar
            # En este caso, como calculamos n_segments correctamente, no debería pasar
            available_samples = length(signal) - start_idx + 1
            if available_samples > 0
                # Copiar solo las muestras disponibles (truncar)
                eeg_segmented[ch_idx, 1:available_samples, seg_idx] = signal[start_idx:length(signal)]
                # El resto queda en cero (ya inicializado)
            end
        end
    end
    
    # Mostrar progreso cada 10 canales o al final
    if ch_idx % 10 == 0 || ch_idx == n_channels
        println("  ✓ Canal $ch_idx/$n_channels procesado: $channel")
    end
end

println()
println("✔ Segmentación completada")
println("  → Dimensiones de la hipermatriz: $(size(eeg_segmented))")
println("  → Formato: (canales × muestras × segmentos)")
println()

# -----------------------------------------------------------------------------
# 4. Verificación y estadísticas de segmentos
# -----------------------------------------------------------------------------
# Esta sección calcula estadísticas descriptivas para cada segmento.
# Estas estadísticas ayudan a verificar que la segmentación se realizó correctamente
# y proporcionan información sobre la variabilidad entre segmentos.
#
# ESTADÍSTICAS CALCULADAS:
#   - Media: valor promedio de amplitud en el segmento
#   - Std: desviación estándar (variabilidad)
#   - Min/Max: valores extremos de amplitud
#   - RMS: raíz cuadrada de la media de los cuadrados (medida de potencia)

println("📊 ESTADÍSTICAS DE SEGMENTACIÓN")
println("-" ^ 80)

# Calcular estadísticas por segmento
# Se calculan sobre todos los canales y muestras de cada segmento
segment_stats = DataFrame(
    Segmento = 1:n_segments,
    Media_µV = [mean(eeg_segmented[:, :, seg]) for seg in 1:n_segments],  # media global del segmento
    Std_µV = [std(vec(eeg_segmented[:, :, seg])) for seg in 1:n_segments],  # desviación estándar
    Min_µV = [minimum(eeg_segmented[:, :, seg]) for seg in 1:n_segments],    # valor mínimo
    Max_µV = [maximum(eeg_segmented[:, :, seg]) for seg in 1:n_segments],    # valor máximo
    RMS_µV = [sqrt(mean(eeg_segmented[:, :, seg].^2)) for seg in 1:n_segments]  # RMS (potencia)
)

println("  → Estadísticas por segmento (primeros 5):")
display(first(segment_stats, 5))
println()

println("  → Resumen estadístico:")
println("    Media global: $(round(mean(segment_stats.Media_µV), digits=3)) µV")
println("    Std global: $(round(mean(segment_stats.Std_µV), digits=3)) µV")
println("    RMS global: $(round(mean(segment_stats.RMS_µV), digits=3)) µV")
println()

# -----------------------------------------------------------------------------
# 5. Visualización de ejemplo: primer canal, primeros segmentos
# -----------------------------------------------------------------------------
# Esta sección genera visualizaciones para inspeccionar los segmentos creados:
#   1. Gráfico superpuesto: primeros segmentos del primer canal
#      - Permite ver la forma de los segmentos y verificar que son consecutivos
#   2. Gráfico apilado: todos los segmentos del primer canal
#      - Muestra la evolución temporal a lo largo de toda la señal
#      - Útil para detectar patrones o artefactos que se repiten

println("📈 VISUALIZACIÓN DE EJEMPLO")
println("-" ^ 80)

# Visualizar primeros 5 segmentos del primer canal como ejemplo
n_segments_plot = min(5, n_segments)  # mostrar máximo 5 segmentos
ch_example = 1                        # primer canal
channel_name = channels[ch_example]

# Crear vector de tiempo para el eje X (en segundos para un segmento)
time_segment = collect(0:(segment_length_samples-1)) ./ fs

p_segments = plot(
    xlabel = "Tiempo (s)",
    ylabel = "Amplitud (µV)",
    title = "Segmentos de ejemplo - Canal $channel_name",
    legend = :outerright,
    size = (1000, 400)
)

for seg_idx in 1:n_segments_plot
    segment_data = eeg_segmented[ch_example, :, seg_idx]
    plot!(p_segments, time_segment, segment_data, 
          label = "Segmento $seg_idx", 
          lw = 1.5, alpha = 0.7)
end

display(p_segments)
println()

# Visualización de todos los segmentos apilados (primer canal)
println("  → Visualización de todos los segmentos apilados (canal $channel_name)")
sep = 50.0  # Separación vertical entre segmentos
offsets_segments = [(n_segments - i) * sep for i in 1:n_segments]

p_stacked = plot(
    xlabel = "Tiempo (s)",
    ylabel = "",
    title = "Todos los segmentos apilados - Canal $channel_name",
    legend = false,
    grid = false,
    size = (1000, 600)
)

for seg_idx in 1:n_segments
    segment_data = eeg_segmented[ch_example, :, seg_idx]
    plot!(p_stacked, time_segment, segment_data .+ offsets_segments[seg_idx], 
          lw = 1, alpha = 0.6)
end

# Añadir etiquetas en el eje Y
yticks_positions = offsets_segments[1:max(1, div(n_segments, 10)):end]
yticks_labels = ["Seg $(i*max(1, div(n_segments, 10)))" for i in 1:length(yticks_positions)]
plot!(p_stacked, yticks = (yticks_positions, yticks_labels))

display(p_stacked)
println()

# -----------------------------------------------------------------------------
# 6. Guardar resultados
# -----------------------------------------------------------------------------
# Esta sección guarda todos los resultados del proceso de segmentación:
#   1. Hipermatriz segmentada: datos EEG organizados en segmentos
#   2. Diccionario de información: metadatos y parámetros de segmentación
#   3. Estadísticas en CSV: tabla con estadísticas por segmento

println("💾 GUARDANDO RESULTADOS")
println("-" ^ 80)

# Crear directorio de segmentación si no existe
dir_segmentation = stage_dir(:segmentation)
if !isdir(dir_segmentation)
    mkpath(dir_segmentation)
    println("  ✓ Directorio creado: $dir_segmentation")
end

# Guardar hipermatriz segmentada
# Formato: (canales × muestras × segmentos)
# Esta es la estructura principal que será usada en pasos posteriores
path_eeg_segmented = joinpath(dir_segmentation, "eeg_segmented.bin")
Serialization.serialize(path_eeg_segmented, eeg_segmented)
println("  ✓ Hipermatriz segmentada guardada en: $(basename(path_eeg_segmented))")

# Guardar información completa de segmentación
# Incluye la hipermatriz, parámetros usados, y toda la información necesaria
# para reproducir o analizar el proceso de segmentación
dict_segmentation_info = Dict(
    "eeg_segmented" => eeg_segmented,              # hipermatriz segmentada
    "channels" => channels,                         # nombres de canales
    "fs" => fs,                                     # frecuencia de muestreo
    "segment_length_s" => segment_length_s,        # longitud en segundos
    "segment_length_samples" => segment_length_samples,  # longitud en muestras
    "segment_overlap_s" => segment_overlap_s,      # solapamiento en segundos
    "segment_overlap_samples" => segment_overlap_samples,  # solapamiento en muestras
    "segment_step_samples" => segment_step_samples,  # paso entre segmentos
    "n_segments" => n_segments,                    # número de segmentos creados
    "n_channels" => n_channels,                     # número de canales
    "n_samples_total" => n_samples_total,           # muestras totales originales
    "skip_bad_intervals" => skip_bad_intervals      # parámetro (no usado actualmente)
)

path_dict_segmentation = joinpath(dir_segmentation, "dict_segmentation_info.bin")
Serialization.serialize(path_dict_segmentation, dict_segmentation_info)
println("  ✓ Información de segmentación guardada en: $(basename(path_dict_segmentation))")

# Guardar estadísticas en CSV
path_stats_csv = joinpath(dir_segmentation, "segment_statistics.csv")
CSV.write(path_stats_csv, segment_stats)
println("  ✓ Estadísticas de segmentos guardadas en: $(basename(path_stats_csv))")
println()

# -----------------------------------------------------------------------------
# 7. Resumen final
# -----------------------------------------------------------------------------
println("=" ^ 80)
println("✨ SEGMENTACIÓN COMPLETADA")
println("=" ^ 80)
println()
println("📋 RESUMEN:")
println("  • Hipermatriz creada: $(size(eeg_segmented))")
println("  • Formato: (canales × muestras × segmentos)")
println("  • Canales: $n_channels")
println("  • Muestras por segmento: $segment_length_samples ($(segment_length_s) s)")
println("  • Número de segmentos: $n_segments")
println("  • Solapamiento: $(segment_overlap_s) s")
println("  • Resolución espectral esperada: ~$(round(fs/segment_length_samples, digits=3)) Hz")
println()
println("💾 Archivos guardados en: $dir_segmentation")
println()
