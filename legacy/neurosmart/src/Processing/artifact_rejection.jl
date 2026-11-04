#!/usr/bin/env julia
# -*- coding: utf-8 -*- #shebang        
# src/artifact_rejection.jl

# RECHAZO DE ARTEFACTOS (Artifact Rejection)
# ===========================================
# Esta rutina realiza un rechazo automático de artefactos basado en umbrales de amplitud.
# Se aplica sobre datos ya segmentados y corregidos de baseline.
#
# PROCESO:
#   1. Carga datos segmentados después de baseline correction
#   2. Configura umbrales de amplitud (±70 µV por defecto)
#   3. Detecta artefactos por segmento:
#      - Verifica si algún valor en los canales usados excede los umbrales
#      - Marca segmentos completos como "buenos" o "malos"
#   4. Filtra y mantiene solo los segmentos buenos
#   5. Genera estadísticas y visualizaciones
#   6. Guarda datos filtrados y metadatos
#
# CRITERIOS DE RECHAZO:
#   - Umbral de amplitud: ±70 µV (configurable)
#   - Si cualquier valor en un segmento excede el umbral, todo el segmento se rechaza
#   - Se analizan los primeros 30 canales por defecto (configurable)
#
# NOTA: Los parámetros "before_event_ms" y "after_event_ms" están configurados
# pero no se usan en la detección actual (se guardan para referencia futura).
# La detección actual es por umbral de amplitud en todo el segmento.
#
# Justificación Técnica:
# El umbral de ±70 µV se eligió de forma conservadora para detectar artefactos
# graves (saturaciones, picos musculares extremos) sin eliminar señal neuronal válida.
# Tras la limpieza ICA, típicamente no se eliminan segmentos, lo cual indica
# que los artefactos ya fueron removidos en pasos anteriores. El criterio es
# adecuado para datos de resting state donde el movimiento voluntario es escaso.

# Listado de librerías
using Serialization
using CSV, DataFrames
using LinearAlgebra
using Statistics
using StatsBase
using Plots
# El backend GR se inicializa automáticamente al crear el primer plot

# -----------------------------------------------------------------------------
# 1. Carga de datos después de baseline correction
# -----------------------------------------------------------------------------
# Se cargan los datos EEG que ya han sido:
#   - Segmentados en épocas/segmentos de tiempo
#   - Corregidos de baseline (línea base removida)
#
# Los datos están en formato de hipermatriz 3D: (canales × muestras × segmentos)
# Cada segmento es una época temporal que será evaluada independientemente.

println("=" ^ 80)
println("📊 RECHAZO DE ARTEFACTOS (Semi-automatic Inspection)")
println("=" ^ 80)
println()

# Cargar datos desde el paso anterior (baseline correction)
dir_baseline = stage_dir(:baseline)
path_dict_baseline = joinpath(dir_baseline, "dict_1st_baseline_correction.bin")
dict_baseline_info = Serialization.deserialize(path_dict_baseline)

# Extraer información de baseline correction
# La hipermatriz contiene todos los segmentos que serán evaluados
eeg_baseline_corrected = dict_baseline_info["eeg_baseline_corrected"]  # (canales × muestras × segmentos)
channels = dict_baseline_info["channels"]                                # nombres de canales
fs = dict_baseline_info["fs"]                                             # frecuencia de muestreo
segment_length_samples = dict_baseline_info["segment_length_samples"]    # longitud de cada segmento
n_segments = dict_baseline_info["n_segments"]                             # número total de segmentos
n_channels = dict_baseline_info["n_channels"]                             # número total de canales

println("✔ Datos con baseline correction cargados desde: $(basename(path_dict_baseline))")
println("  → Dimensiones: $(size(eeg_baseline_corrected))")
println("  → Formato: (canales × muestras × segmentos)")
println("  → Número de canales: $n_channels")
println("  → Muestras por segmento: $segment_length_samples")
println("  → Número de segmentos: $n_segments")
println("  → Frecuencia de muestreo: $(fs) Hz")
println()

# -----------------------------------------------------------------------------
# 2. Configuración de artifact rejection
# -----------------------------------------------------------------------------
# Se configuran los parámetros para la detección de artefactos:
#   - Umbrales de amplitud: valores fuera de este rango se consideran artefactos
#   - Canales a usar: solo se evalúan los primeros N canales (por defecto 30)
#   - Ventanas temporales: parámetros para rechazo alrededor de eventos
#     (actualmente no se usan en la detección, solo se guardan para referencia)

n_channels_used = 30                    # Número de canales a usar para artifact rejection
                                         # Solo se evalúan estos canales (típicamente los primeros)
min_amplitude_µV = -70.0                # Amplitud mínima permitida (µV)
                                         # Valores menores se consideran artefactos (saturaciones negativas)
max_amplitude_µV = +70.0                # Amplitud máxima permitida (µV)
                                         # Valores mayores se consideran artefactos (saturaciones positivas, picos)
before_event_ms = 200                   # Marcar como malo: antes del evento (ms)
                                         # NOTA: Actualmente no se usa en la detección
after_event_ms = 300                    # Marcar como malo: después del evento (ms)
                                         # NOTA: Actualmente no se usa en la detección

# Convertir tiempos a muestras (para uso futuro)
# Estos valores se guardan pero no se usan en la detección actual
before_event_samples = Int(round(before_event_ms / 1000.0 * fs))
after_event_samples = Int(round(after_event_ms / 1000.0 * fs))

println("⚙️  CONFIGURACIÓN DE ARTIFACT REJECTION")
println("-" ^ 80)
println("  Canales usados: $n_channels_used")
println("  Amplitud mínima permitida: $(min_amplitude_µV) µV")
println("  Amplitud máxima permitida: $(max_amplitude_µV) µV")
println("  Marcar como malo antes del evento: $(before_event_ms) ms ($(before_event_samples) muestras)")
println("  Marcar como malo después del evento: $(after_event_ms) ms ($(after_event_samples) muestras)")
println()

# Verificar que tenemos al menos 30 canales
if n_channels < n_channels_used
    println("⚠️  ADVERTENCIA: Solo hay $n_channels canales disponibles, pero se requieren $n_channels_used")
    n_channels_used = n_channels
end

# Usar los primeros n_channels_used canales
channels_used = channels[1:n_channels_used]
println("  → Canales utilizados para artifact rejection: $(join(channels_used, ", "))")
println()

# -----------------------------------------------------------------------------
# 3. Detección de artefactos por segmento
# -----------------------------------------------------------------------------
# Esta sección evalúa cada segmento independientemente para detectar artefactos.
# 
# CRITERIO DE RECHAZO:
#   Un segmento se marca como "malo" si CUALQUIER valor en los canales evaluados
#   excede los umbrales de amplitud (±70 µV). Esto es un criterio conservador:
#   si hay un artefacto en cualquier canal, se rechaza todo el segmento.
#
# PROCESO:
#   1. Para cada segmento, extrae los datos de los canales usados
#   2. Calcula el mínimo y máximo de amplitud en esos canales
#   3. Compara con los umbrales configurados
#   4. Si excede, marca el segmento como "malo" e identifica qué canales violaron
#   5. Almacena estadísticas para análisis posterior

println("🔄 DETECTANDO ARTEFACTOS")
println("-" ^ 80)

# Vector para marcar segmentos como buenos (true) o malos (false)
# Inicialmente todos se asumen buenos, se marcan como malos si violan umbrales
segment_is_good = fill(true, n_segments)

# Estadísticas por segmento (para análisis y reportes)
segment_min_amplitude = zeros(n_segments)        # amplitud mínima de cada segmento
segment_max_amplitude = zeros(n_segments)        # amplitud máxima de cada segmento
segment_exceeds_threshold = fill(false, n_segments)  # si el segmento excede umbrales
segment_violation_channels = Vector{Vector{String}}(undef, n_segments)  # canales que violan

# Evaluar cada segmento
for seg_idx in 1:n_segments
    # Obtener el segmento completo (todos los canales, todas las muestras)
    segment_data = eeg_baseline_corrected[:, :, seg_idx]  # (canales × muestras)
    
    # Verificar solo los canales usados para artifact rejection
    # Esto permite enfocarse en canales específicos (p.ej., los primeros 30)
    segment_used = segment_data[1:n_channels_used, :]
    
    # Calcular min y max del segmento (en los canales usados)
    # Estos valores determinan si el segmento tiene artefactos
    seg_min = minimum(segment_used)  # valor mínimo en todo el segmento
    seg_max = maximum(segment_used)   # valor máximo en todo el segmento
    
    segment_min_amplitude[seg_idx] = seg_min
    segment_max_amplitude[seg_idx] = seg_max
    
    # Verificar si algún valor excede los umbrales
    # Si el mínimo es menor que el umbral inferior O el máximo es mayor que el umbral superior,
    # el segmento contiene artefactos
    exceeds_min = seg_min < min_amplitude_µV  # ¿hay valores por debajo del umbral?
    exceeds_max = seg_max > max_amplitude_µV  # ¿hay valores por encima del umbral?
    
    if exceeds_min || exceeds_max
        # Segmento contiene artefactos → marcar como malo
        segment_exceeds_threshold[seg_idx] = true
        segment_is_good[seg_idx] = false
        
        # Identificar qué canales específicos violan el umbral
        # Esto es útil para diagnóstico y reportes
        violating_channels = String[]
        for ch_idx in 1:n_channels_used
            ch_data = segment_used[ch_idx, :]  # datos del canal en este segmento
            # Verificar si este canal específico viola los umbrales
            if minimum(ch_data) < min_amplitude_µV || maximum(ch_data) > max_amplitude_µV
                push!(violating_channels, channels_used[ch_idx])
            end
        end
        segment_violation_channels[seg_idx] = violating_channels
    else
        # Segmento está dentro de los umbrales → es bueno
        segment_violation_channels[seg_idx] = String[]  # ningún canal viola
    end
    
    # Mostrar progreso cada 20 segmentos o al final
    if seg_idx % 20 == 0 || seg_idx == n_segments
        status = segment_is_good[seg_idx] ? "✓" : "✗"
        println("  $status Segmento $seg_idx/$n_segments: min=$(round(seg_min, digits=2)) µV, max=$(round(seg_max, digits=2)) µV")
    end
end

println()
println("✔ Detección de artefactos completada")
println()

# -----------------------------------------------------------------------------
# 4. Estadísticas de artifact rejection
# -----------------------------------------------------------------------------
# Esta sección calcula y muestra estadísticas sobre el proceso de rechazo:
#   - Número de segmentos mantenidos vs eliminados
#   - Distribución de amplitudes
#   - Detalles de segmentos eliminados
#   - Creación de DataFrame con estadísticas completas

println("📊 ESTADÍSTICAS DE ARTIFACT REJECTION")
println("-" ^ 80)

# Contar segmentos mantenidos y eliminados
n_segments_kept = sum(segment_is_good)      # segmentos que pasaron el filtro
n_segments_removed = sum(.!segment_is_good)  # segmentos rechazados

println("  → Segmentos analizados: $n_segments")
println("  → Segmentos mantenidos: $n_segments_kept")
println("  → Segmentos eliminados: $n_segments_removed")
println("  → Porcentaje mantenido: $(round(100 * n_segments_kept / n_segments, digits=2))%")
println()

if n_segments_removed > 0
    println("  → Segmentos eliminados (detalles):")
    removed_segments = findall(.!segment_is_good)
    for seg_idx in removed_segments
        println("    • Segmento $seg_idx: min=$(round(segment_min_amplitude[seg_idx], digits=2)) µV, max=$(round(segment_max_amplitude[seg_idx], digits=2)) µV")
        if !isempty(segment_violation_channels[seg_idx])
            println("      Canales con violación: $(join(segment_violation_channels[seg_idx], ", "))")
        end
    end
    println()
end

# Crear DataFrame con estadísticas
artifact_stats = DataFrame(
    Segmento = 1:n_segments,
    Min_Amplitud_µV = segment_min_amplitude,
    Max_Amplitud_µV = segment_max_amplitude,
    Excede_Umbral = segment_exceeds_threshold,
    Es_Bueno = segment_is_good,
    Canales_Violacion = [join(segment_violation_channels[i], ", ") for i in 1:n_segments]
)

println("  → Estadísticas por segmento (primeros 5):")
display(first(artifact_stats, 5))
println()

println("  → Resumen estadístico de amplitudes:")
println("    Min global: $(round(minimum(segment_min_amplitude), digits=2)) µV")
println("    Max global: $(round(maximum(segment_max_amplitude), digits=2)) µV")
println("    Media de mínimos: $(round(mean(segment_min_amplitude), digits=2)) µV")
println("    Media de máximos: $(round(mean(segment_max_amplitude), digits=2)) µV")
println()

# -----------------------------------------------------------------------------
# 5. Filtrar segmentos: mantener solo los buenos
# -----------------------------------------------------------------------------
# Esta sección crea una nueva hipermatriz que contiene únicamente los segmentos
# que pasaron el filtro de artifact rejection. Los segmentos marcados como "malos"
# se eliminan completamente de los datos.
#
# El resultado es una hipermatriz con las mismas dimensiones en canales y muestras,
# pero con menos segmentos (solo los buenos).

println("🔄 FILTRANDO SEGMENTOS")
println("-" ^ 80)

# Obtener índices de segmentos buenos
# Estos son los segmentos que pasaron el filtro y se mantendrán
good_segment_indices = findall(segment_is_good)

# Crear nueva hipermatriz solo con segmentos buenos
# Se seleccionan solo las "rebanadas" (segmentos) que son buenas
# Formato: (canales × muestras × segmentos_buenos)
eeg_artifact_rejected = eeg_baseline_corrected[:, :, good_segment_indices]
n_segments_kept_final = length(good_segment_indices)

println("  → Segmentos originales: $n_segments")
println("  → Segmentos después del filtrado: $n_segments_kept_final")
println("  → Dimensiones de la hipermatriz filtrada: $(size(eeg_artifact_rejected))")
println()

# -----------------------------------------------------------------------------
# 6. Visualización de ejemplo
# -----------------------------------------------------------------------------
# Esta sección genera visualizaciones para inspeccionar los resultados:
#   1. Gráfico temporal: muestra algunos segmentos con los umbrales marcados
#      - Segmentos buenos en azul
#      - Segmentos malos en rojo
#   2. Histogramas: distribución de amplitudes máximas y mínimas por segmento
#      - Permite ver si los umbrales son apropiados
#      - Muestra dónde están los umbrales en relación a la distribución

println("📈 VISUALIZACIÓN DE EJEMPLO")
println("-" ^ 80)

# Visualizar algunos segmentos (buenos y malos si hay)
# Se muestran los primeros 5 segmentos como ejemplo
n_segments_plot = min(5, n_segments)
ch_example = 1  # canal a visualizar (primer canal)
channel_name = channels[ch_example]

# Crear vector de tiempo para el eje X (en segundos)
time_segment = collect(0:(segment_length_samples-1)) ./ fs

p_artifact = plot(
    xlabel = "Tiempo (s)",
    ylabel = "Amplitud (µV)",
    title = "Artifact Rejection - Canal $channel_name (ejemplo)",
    legend = :outerright,
    size = (1200, 500)
)

# Añadir líneas de umbral
hline!(p_artifact, [min_amplitude_µV], 
       label = "Umbral mínimo ($(min_amplitude_µV) µV)", 
       linestyle = :dash, color = :red, alpha = 0.5)
hline!(p_artifact, [max_amplitude_µV], 
       label = "Umbral máximo ($(max_amplitude_µV) µV)", 
       linestyle = :dash, color = :red, alpha = 0.5)

# Plotear segmentos de ejemplo
for seg_idx in 1:n_segments_plot
    segment_data = eeg_baseline_corrected[ch_example, :, seg_idx]
    is_good = segment_is_good[seg_idx]
    color_seg = is_good ? :blue : :red
    alpha_seg = is_good ? 0.7 : 0.4
    label_seg = is_good ? "Seg $seg_idx (bueno)" : "Seg $seg_idx (malo)"
    
    plot!(p_artifact, time_segment, segment_data, 
          label = label_seg, 
          lw = 1.5, alpha = alpha_seg, color = color_seg)
end

display(p_artifact)
println()

# Histograma de amplitudes máximas y mínimas
p_hist = plot(
    layout = (2, 1),
    size = (800, 600)
)

histogram!(p_hist[1], segment_max_amplitude, 
          bins = 30,
          xlabel = "Amplitud máxima (µV)",
          ylabel = "Número de segmentos",
          title = "Distribución de amplitudes máximas por segmento",
          legend = false,
          color = :blue,
          alpha = 0.6)
vline!(p_hist[1], [max_amplitude_µV], 
       linestyle = :dash, color = :red, linewidth = 2,
       label = "Umbral ($(max_amplitude_µV) µV)")

histogram!(p_hist[2], segment_min_amplitude, 
          bins = 30,
          xlabel = "Amplitud mínima (µV)",
          ylabel = "Número de segmentos",
          title = "Distribución de amplitudes mínimas por segmento",
          legend = false,
          color = :blue,
          alpha = 0.6)
vline!(p_hist[2], [min_amplitude_µV], 
       linestyle = :dash, color = :red, linewidth = 2,
       label = "Umbral ($(min_amplitude_µV) µV)")

display(p_hist)
println()

# -----------------------------------------------------------------------------
# 7. Guardar resultados
# -----------------------------------------------------------------------------
# Esta sección guarda todos los resultados del proceso de artifact rejection:
#   1. Hipermatriz filtrada: datos EEG con segmentos malos eliminados
#   2. Diccionario de información: metadatos, parámetros, índices de segmentos
#   3. Estadísticas en CSV: tabla con estadísticas por segmento

println("💾 GUARDANDO RESULTADOS")
println("-" ^ 80)

# Crear directorio de artifact rejection si no existe
dir_artifact_rejection = stage_dir(:artifact_rejection)
if !isdir(dir_artifact_rejection)
    mkpath(dir_artifact_rejection)
    println("  ✓ Directorio creado: $dir_artifact_rejection")
end

# Guardar hipermatriz con artifact rejection aplicado
# Esta es la hipermatriz final que contiene solo los segmentos buenos
# Formato: (canales × muestras × segmentos_buenos)
path_eeg_artifact_rejected = joinpath(dir_artifact_rejection, "eeg_artifact_rejected.bin")
Serialization.serialize(path_eeg_artifact_rejected, eeg_artifact_rejected)
println("  ✓ Hipermatriz con artifact rejection guardada en: $(basename(path_eeg_artifact_rejected))")

# Guardar información completa de artifact rejection
# Incluye tanto los datos filtrados como los originales, parámetros usados,
# y toda la información necesaria para reproducir o analizar el proceso
dict_artifact_rejection_info = Dict(
    "eeg_artifact_rejected" => eeg_artifact_rejected,      # datos filtrados (solo segmentos buenos)
    "eeg_baseline_corrected" => eeg_baseline_corrected,    # datos originales (todos los segmentos)
    "channels" => channels,                                 # nombres de todos los canales
    "channels_used" => channels_used,                      # canales usados para artifact rejection
    "fs" => fs,                                             # frecuencia de muestreo
    "segment_length_samples" => segment_length_samples,   # longitud de cada segmento
    "n_segments_original" => n_segments,                  # número original de segmentos
    "n_segments_kept" => n_segments_kept_final,            # segmentos mantenidos
    "n_segments_removed" => n_segments_removed,            # segmentos eliminados
    "good_segment_indices" => good_segment_indices,        # índices de segmentos buenos
    "segment_is_good" => segment_is_good,                  # vector booleano (true=bueno, false=malo)
    "n_channels" => n_channels,                             # número total de canales
    "n_channels_used" => n_channels_used,                  # canales usados para evaluación
    "min_amplitude_µV" => min_amplitude_µV,               # umbral mínimo configurado
    "max_amplitude_µV" => max_amplitude_µV,                # umbral máximo configurado
    "before_event_ms" => before_event_ms,                  # parámetro (no usado actualmente)
    "after_event_ms" => after_event_ms,                    # parámetro (no usado actualmente)
    "before_event_samples" => before_event_samples,       # parámetro (no usado actualmente)
    "after_event_samples" => after_event_samples           # parámetro (no usado actualmente)
)

path_dict_artifact_rejection = joinpath(dir_artifact_rejection, "dict_artifact_rejection.bin")
Serialization.serialize(path_dict_artifact_rejection, dict_artifact_rejection_info)
println("  ✓ Información de artifact rejection guardada en: $(basename(path_dict_artifact_rejection))")

# Guardar estadísticas en CSV
path_stats_csv = joinpath(dir_artifact_rejection, "artifact_rejection_statistics.csv")
CSV.write(path_stats_csv, artifact_stats)
println("  ✓ Estadísticas de artifact rejection guardadas en: $(basename(path_stats_csv))")
println()

# -----------------------------------------------------------------------------
# 8. Resumen final
# -----------------------------------------------------------------------------
println("=" ^ 80)
println("✨ ARTIFACT REJECTION COMPLETADO")
println("=" ^ 80)
println()
println("📋 RESUMEN:")
println("  • Hipermatriz original: $(size(eeg_baseline_corrected))")
println("  • Hipermatriz filtrada: $(size(eeg_artifact_rejected))")
println("  • Formato: (canales × muestras × segmentos)")
println("  • Canales totales: $n_channels")
println("  • Canales usados para artifact rejection: $n_channels_used")
println("  • Muestras por segmento: $segment_length_samples")
println("  • Segmentos originales: $n_segments")
println("  • Segmentos mantenidos: $n_segments_kept_final")
println("  • Segmentos eliminados: $n_segments_removed")
println("  • Umbral de amplitud: $(min_amplitude_µV) µV a $(max_amplitude_µV) µV")
println("  • Marcar como malo antes del evento: $(before_event_ms) ms")
println("  • Marcar como malo después del evento: $(after_event_ms) ms")
println()
println("💾 Archivos guardados en: $dir_artifact_rejection")
println("  → eeg_artifact_rejected.bin")
println("  → dict_artifact_rejection.bin")
println("  → artifact_rejection_statistics.csv")
println()
