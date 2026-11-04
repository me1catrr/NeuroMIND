#!/usr/bin/env julia
# -*- coding: utf-8 -*- #shebang        
# src/baseline.jl

# CORRECCIÓN DE BASELINE (Baseline Correction)
# ============================================
# Esta rutina aplica corrección de baseline a datos EEG segmentados.
# Se resta la media del intervalo de baseline de cada segmento, normalizando
# la señal para que el intervalo de baseline tenga media cero.
#
# PROCESO:
#   1. Carga datos segmentados (hipermatriz: canales × muestras × segmentos)
#   2. Configura intervalo de baseline (0.00 s - 0.10 s = primeros 100 ms)
#   3. Para cada segmento y canal:
#      - Calcula la media del intervalo de baseline
#      - Resta esta media de todo el segmento
#   4. Genera estadísticas y visualizaciones de comparación
#   5. Guarda datos corregidos como "1st_baseline_correction"
#
# CONFIGURACIÓN:
#   - Tipo: Segment-based (baseline calculado por segmento)
#   - Intervalo: 0.00 s - 0.10 s (primeros 100 ms de cada segmento)
#   - Método: Resta de la media del intervalo de baseline
#
# Justificación Técnica:
# Esta corrección se mantuvo para mantener homogeneidad con proyectos previos.
# En resting state sin estímulos, no es un "baseline fisiológico" en el sentido
# tradicional, sino más bien una normalización que elimina el offset DC del
# comienzo de cada segmento. Puede entenderse como una resta de la media del
# comienzo del segmento, cumpliendo función de normalización. Alternativamente,
# podría haberse empleado demean del segmento completo, opción válida para
# estudios puramente espectrales. Esta es la primera corrección de baseline
# aplicada, por lo que se guarda como "1st_baseline_correction".

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
# 1. Carga de datos segmentados
# -----------------------------------------------------------------------------
# Se cargan los datos EEG que ya han sido segmentados en épocas temporales.
# Los datos están en formato de hipermatriz 3D: (canales × muestras × segmentos).
# Cada segmento será corregido de baseline independientemente.

println("=" ^ 80)
println("📊 CORRECCIÓN DE BASELINE")
println("=" ^ 80)
println()

# Cargar datos desde el paso anterior (segmentación)
# La hipermatriz contiene todos los segmentos que serán corregidos
dir_segmentation = stage_dir(:segmentation)
path_dict_segmentation = joinpath(dir_segmentation, "dict_segmentation_info.bin")
dict_segmentation_info = Serialization.deserialize(path_dict_segmentation)

# Extraer información de segmentación
# Todos estos parámetros son necesarios para aplicar la corrección correctamente
eeg_segmented = dict_segmentation_info["eeg_segmented"]        # hipermatriz segmentada
channels = dict_segmentation_info["channels"]                 # nombres de canales
fs = dict_segmentation_info["fs"]                             # frecuencia de muestreo
segment_length_samples = dict_segmentation_info["segment_length_samples"]  # longitud de cada segmento
n_segments = dict_segmentation_info["n_segments"]              # número de segmentos
n_channels = dict_segmentation_info["n_channels"]              # número de canales

println("✔ Datos segmentados cargados desde: $(basename(path_dict_segmentation))")
println("  → Dimensiones: $(size(eeg_segmented))")
println("  → Formato: (canales × muestras × segmentos)")
println("  → Número de canales: $n_channels")
println("  → Muestras por segmento: $segment_length_samples")
println("  → Número de segmentos: $n_segments")
println("  → Frecuencia de muestreo: $(fs) Hz")
println()

# -----------------------------------------------------------------------------
# 2. Configuración de corrección de baseline
# -----------------------------------------------------------------------------
# Se configura el intervalo de baseline que se usará para calcular la media.
# Este intervalo típicamente corresponde al inicio de cada segmento, antes de
# cualquier actividad de interés. La media de este intervalo se resta de todo
# el segmento para eliminar el offset DC.
#
# CONFIGURACIÓN ACTUAL:
#   - Intervalo: 0.00 s - 0.10 s (primeros 100 ms de cada segmento)
#   - A 500 Hz: muestras 1-50 (primeras 50 muestras)

baseline_start_s = 0.00    # Inicio del intervalo de baseline (s)
                            # Comienza al inicio del segmento (t=0)
baseline_end_s = 0.10       # Fin del intervalo de baseline (s) = 100 ms
                            # Primeros 100 ms del segmento

# Convertir tiempos a índices de muestras (1-indexed en Julia)
baseline_start_sample = Int(round(baseline_start_s * fs)) + 1  # Muestra 1 (índice base 1)
baseline_end_sample = Int(round(baseline_end_s * fs))          # Muestra 50 (a 500 Hz)

println("⚙️  CONFIGURACIÓN DE BASELINE")
println("-" ^ 80)
println("  Intervalo de baseline: $(baseline_start_s) s - $(baseline_end_s) s")
println("  Muestras de baseline: $(baseline_start_sample) - $(baseline_end_sample)")
println("  Número de muestras en baseline: $(baseline_end_sample - baseline_start_sample + 1)")
println()

# -----------------------------------------------------------------------------
# 3. Aplicar corrección de baseline
# -----------------------------------------------------------------------------
# Esta sección aplica la corrección de baseline a cada segmento de cada canal.
#
# PROCESO:
#   Para cada canal y cada segmento:
#     1. Extrae el segmento completo
#     2. Identifica el intervalo de baseline (primeros 100 ms)
#     3. Calcula la media de ese intervalo
#     4. Resta la media de todo el segmento
#
# RESULTADO:
#   Después de la corrección, el intervalo de baseline tiene media ≈ 0 en cada segmento.
#   Esto elimina el offset DC y normaliza la señal, facilitando análisis posteriores.

println("🔄 APLICANDO CORRECCIÓN DE BASELINE")
println("-" ^ 80)

# Crear copia de los datos segmentados
# Se trabaja sobre una copia para preservar los datos originales
eeg_baseline_corrected = copy(eeg_segmented)

# Aplicar corrección de baseline a cada segmento
# Se procesa cada canal y cada segmento independientemente
for ch_idx in 1:n_channels
    for seg_idx in 1:n_segments
        # Obtener el segmento completo (todas las muestras del segmento)
        segment = eeg_baseline_corrected[ch_idx, :, seg_idx]
        
        # Calcular la media del intervalo de baseline
        # Este es el valor que se restará de todo el segmento
        baseline_interval = segment[baseline_start_sample:baseline_end_sample]
        baseline_mean = mean(baseline_interval)
        
        # Restar la media del baseline de todo el segmento
        # Esto hace que el intervalo de baseline tenga media ≈ 0
        # y elimina el offset DC del segmento completo
        eeg_baseline_corrected[ch_idx, :, seg_idx] = segment .- baseline_mean
    end
    
    # Mostrar progreso cada 10 canales o al final
    if ch_idx % 10 == 0 || ch_idx == n_channels
        println("  ✓ Canal $ch_idx/$n_channels procesado: $(channels[ch_idx])")
    end
end

println()
println("✔ Corrección de baseline completada")
println("  → Dimensiones de la hipermatriz corregida: $(size(eeg_baseline_corrected))")
println()

# -----------------------------------------------------------------------------
# 4. Verificación y estadísticas
# -----------------------------------------------------------------------------
# Esta sección calcula estadísticas antes y después de la corrección para verificar
# que la corrección se aplicó correctamente. La verificación principal es que la
# media del intervalo de baseline después de la corrección debe ser ≈ 0.
#
# ESTADÍSTICAS CALCULADAS:
#   - Media global del segmento: antes y después
#   - Desviación estándar: para verificar que no cambió (solo se resta una constante)
#   - Media del intervalo de baseline: debe ser ≈ 0 después de la corrección

println("📊 ESTADÍSTICAS DE CORRECCIÓN DE BASELINE")
println("-" ^ 80)

# Calcular estadísticas antes y después de la corrección
# Se calculan sobre todos los canales de cada segmento
stats_before = DataFrame(
    Segmento = 1:n_segments,
    Media_µV = [mean(eeg_segmented[:, :, seg]) for seg in 1:n_segments],  # media global antes
    Std_µV = [std(vec(eeg_segmented[:, :, seg])) for seg in 1:n_segments],  # desviación estándar antes
    Media_Baseline_µV = [mean(eeg_segmented[:, baseline_start_sample:baseline_end_sample, seg]) for seg in 1:n_segments]  # media del baseline antes
)

stats_after = DataFrame(
    Segmento = 1:n_segments,
    Media_µV = [mean(eeg_baseline_corrected[:, :, seg]) for seg in 1:n_segments],  # media global después
    Std_µV = [std(vec(eeg_baseline_corrected[:, :, seg])) for seg in 1:n_segments],  # desviación estándar después (debe ser igual)
    Media_Baseline_µV = [mean(eeg_baseline_corrected[:, baseline_start_sample:baseline_end_sample, seg]) for seg in 1:n_segments]  # media del baseline después (debe ser ≈ 0)
)

println("  → Estadísticas ANTES de la corrección (primeros 5 segmentos):")
display(first(stats_before, 5))
println()

println("  → Estadísticas DESPUÉS de la corrección (primeros 5 segmentos):")
display(first(stats_after, 5))
println()

println("  → Verificación:")
println("    Media del baseline ANTES: $(round(mean(stats_before.Media_Baseline_µV), digits=6)) µV")
println("    Media del baseline DESPUÉS: $(round(mean(stats_after.Media_Baseline_µV), digits=6)) µV")
println("    Media global ANTES: $(round(mean(stats_before.Media_µV), digits=6)) µV")
println("    Media global DESPUÉS: $(round(mean(stats_after.Media_µV), digits=6)) µV")
println()

# -----------------------------------------------------------------------------
# 5. Visualización de ejemplo: antes vs después
# -----------------------------------------------------------------------------
# Esta sección genera visualizaciones para inspeccionar el efecto de la corrección:
#   - Gráfico comparativo: muestra segmentos antes y después superpuestos
#   - Marca el intervalo de baseline para visualizar qué se corrigió
#   - Permite verificar visualmente que la corrección funcionó correctamente

println("📈 VISUALIZACIÓN DE EJEMPLO")
println("-" ^ 80)

# Visualizar primeros 3 segmentos del primer canal como ejemplo
n_segments_plot = min(3, n_segments)  # mostrar máximo 3 segmentos
ch_example = 1                        # primer canal
channel_name = channels[ch_example]

# Crear vector de tiempo para el eje X (en segundos para un segmento)
time_segment = collect(0:(segment_length_samples-1)) ./ fs

p_comparison = plot(
    xlabel = "Tiempo (s)",
    ylabel = "Amplitud (µV)",
    title = "Corrección de Baseline - Canal $channel_name",
    legend = :outerright,
    size = (1200, 500)
)

for seg_idx in 1:n_segments_plot
    segment_before = eeg_segmented[ch_example, :, seg_idx]
    segment_after = eeg_baseline_corrected[ch_example, :, seg_idx]
    
    plot!(p_comparison, time_segment, segment_before, 
          label = "Antes - Seg $seg_idx", 
          lw = 1.5, alpha = 0.6, linestyle = :solid)
    plot!(p_comparison, time_segment, segment_after, 
          label = "Después - Seg $seg_idx", 
          lw = 1.5, alpha = 0.8, linestyle = :dash)
    
    # Marcar el intervalo de baseline
    baseline_time = collect((baseline_start_sample-1):(baseline_end_sample-1)) ./ fs
    baseline_vals = segment_before[baseline_start_sample:baseline_end_sample]
    plot!(p_comparison, baseline_time, baseline_vals, 
          label = nothing, 
          lw = 3, alpha = 0.3, color = :red, linestyle = :dot)
end

# Añadir línea vertical para marcar el fin del intervalo de baseline
vline!(p_comparison, [baseline_end_s], 
       label = "Fin baseline ($(baseline_end_s*1000) ms)", 
       linestyle = :dashdot, color = :red, alpha = 0.5)

display(p_comparison)
println()

# -----------------------------------------------------------------------------
# 6. Guardar resultados
# -----------------------------------------------------------------------------
# Esta sección guarda todos los resultados del proceso de corrección de baseline:
#   1. Hipermatriz corregida: datos EEG con baseline corregido
#   2. Diccionario de información: metadatos y parámetros de corrección
#   3. Estadísticas en CSV: tabla comparativa antes/después

println("💾 GUARDANDO RESULTADOS")
println("-" ^ 80)

# Crear directorio de baseline si no existe
dir_baseline = stage_dir(:baseline)
if !isdir(dir_baseline)
    mkpath(dir_baseline)
    println("  ✓ Directorio creado: $dir_baseline")
end

# Guardar hipermatriz con corrección de baseline
# Formato: (canales × muestras × segmentos)
# Esta es la primera corrección de baseline, por eso se guarda como "1st_baseline_correction"
path_eeg_baseline = joinpath(dir_baseline, "eeg_1st_baseline_correction.bin")
Serialization.serialize(path_eeg_baseline, eeg_baseline_corrected)
println("  ✓ Hipermatriz con corrección de baseline guardada en: $(basename(path_eeg_baseline))")

# Guardar información completa de corrección de baseline
# Incluye tanto los datos corregidos como los originales, parámetros usados,
# y toda la información necesaria para reproducir o analizar el proceso
dict_baseline_info = Dict(
    "eeg_baseline_corrected" => eeg_baseline_corrected,      # datos corregidos
    "eeg_segmented_original" => eeg_segmented,               # datos originales (para comparación)
    "channels" => channels,                                   # nombres de canales
    "fs" => fs,                                               # frecuencia de muestreo
    "segment_length_samples" => segment_length_samples,      # longitud de cada segmento
    "n_segments" => n_segments,                               # número de segmentos
    "n_channels" => n_channels,                               # número de canales
    "baseline_start_s" => baseline_start_s,                  # inicio del baseline (s)
    "baseline_end_s" => baseline_end_s,                      # fin del baseline (s)
    "baseline_start_sample" => baseline_start_sample,        # inicio del baseline (muestras)
    "baseline_end_sample" => baseline_end_sample,             # fin del baseline (muestras)
    "baseline_type" => "segment-based",                      # tipo de baseline (por segmento)
    "correction_method" => "subtract_baseline_mean"         # método usado (resta de media)
)

path_dict_baseline = joinpath(dir_baseline, "dict_1st_baseline_correction.bin")
Serialization.serialize(path_dict_baseline, dict_baseline_info)
println("  ✓ Información de corrección de baseline guardada en: $(basename(path_dict_baseline))")

# Guardar estadísticas en CSV
stats_comparison = DataFrame(
    Segmento = 1:n_segments,
    Media_Antes_µV = stats_before.Media_µV,
    Media_Despues_µV = stats_after.Media_µV,
    Std_Antes_µV = stats_before.Std_µV,
    Std_Despues_µV = stats_after.Std_µV,
    Media_Baseline_Antes_µV = stats_before.Media_Baseline_µV,
    Media_Baseline_Despues_µV = stats_after.Media_Baseline_µV,
    Diferencia_Media_µV = stats_after.Media_µV .- stats_before.Media_µV
)

path_stats_csv = joinpath(dir_baseline, "baseline_correction_statistics.csv")
CSV.write(path_stats_csv, stats_comparison)
println("  ✓ Estadísticas de corrección guardadas en: $(basename(path_stats_csv))")
println()

# -----------------------------------------------------------------------------
# 7. Resumen final
# -----------------------------------------------------------------------------
println("=" ^ 80)
println("✨ CORRECCIÓN DE BASELINE COMPLETADA")
println("=" ^ 80)
println()
println("📋 RESUMEN:")
println("  • Hipermatriz corregida: $(size(eeg_baseline_corrected))")
println("  • Formato: (canales × muestras × segmentos)")
println("  • Canales: $n_channels")
println("  • Muestras por segmento: $segment_length_samples")
println("  • Número de segmentos: $n_segments")
println("  • Intervalo de baseline: $(baseline_start_s) s - $(baseline_end_s) s ($(baseline_start_s*1000) ms - $(baseline_end_s*1000) ms)")
println("  • Muestras de baseline: $(baseline_start_sample) - $(baseline_end_sample)")
println("  • Método: Resta de la media del intervalo de baseline")
println()
println("💾 Archivos guardados en: $dir_baseline")
println("  → eeg_1st_baseline_correction.bin")
println("  → dict_1st_baseline_correction.bin")
println("  → baseline_correction_statistics.csv")
println()
