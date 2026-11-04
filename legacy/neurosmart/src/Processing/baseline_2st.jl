#!/usr/bin/env julia
# -*- coding: utf-8 -*- #shebang        
# src/baseline_2st.jl

# STEP Nº10. 2nd Baseline Correction Transformation
# --------------------------------------------------
# Baseline type: Segment-based
# Baseline interval: 0.00 s - 0.10 s (0 ms - 100 ms)
#
# Justificación Técnica:
# Se repite por coherencia metodológica, pero no implica un cambio fisiológico.
# Equivale a una segunda normalización estadística por segmento.

# Listado de librerías
using Serialization
using CSV, DataFrames
using LinearAlgebra
using Statistics
using StatsBase
using Plots
# El backend GR se inicializa automáticamente al crear el primer plot

# -----------------------------------------------------------------------------
# 1. Carga de datos después de artifact rejection
# -----------------------------------------------------------------------------
println("=" ^ 80)
println("📊 2ª CORRECCIÓN DE BASELINE")
println("=" ^ 80)
println()

dir_artifact_rejection = stage_dir(:artifact_rejection)
path_dict_artifact_rejection = joinpath(dir_artifact_rejection, "dict_artifact_rejection.bin")
dict_artifact_rejection_info = Serialization.deserialize(path_dict_artifact_rejection)

# Extraer información de artifact rejection
eeg_artifact_rejected = dict_artifact_rejection_info["eeg_artifact_rejected"]
channels = dict_artifact_rejection_info["channels"]
fs = dict_artifact_rejection_info["fs"]
segment_length_samples = dict_artifact_rejection_info["segment_length_samples"]
n_segments = dict_artifact_rejection_info["n_segments_kept"]
n_channels = dict_artifact_rejection_info["n_channels"]

println("✔ Datos con artifact rejection cargados desde: $(basename(path_dict_artifact_rejection))")
println("  → Dimensiones: $(size(eeg_artifact_rejected))")
println("  → Formato: (canales × muestras × segmentos)")
println("  → Número de canales: $n_channels")
println("  → Muestras por segmento: $segment_length_samples")
println("  → Número de segmentos: $n_segments")
println("  → Frecuencia de muestreo: $(fs) Hz")
println()

# -----------------------------------------------------------------------------
# 2. Configuración de corrección de baseline
# -----------------------------------------------------------------------------
baseline_start_s = 0.00    # Inicio del intervalo de baseline (s)
baseline_end_s = 0.10      # Fin del intervalo de baseline (s) = 100 ms

baseline_start_sample = Int(round(baseline_start_s * fs)) + 1  # Muestra 1 (índice base 1)
baseline_end_sample = Int(round(baseline_end_s * fs))          # Muestra 50

println("⚙️  CONFIGURACIÓN DE BASELINE")
println("-" ^ 80)
println("  Intervalo de baseline: $(baseline_start_s) s - $(baseline_end_s) s")
println("  Muestras de baseline: $(baseline_start_sample) - $(baseline_end_sample)")
println("  Número de muestras en baseline: $(baseline_end_sample - baseline_start_sample + 1)")
println()

# -----------------------------------------------------------------------------
# 3. Aplicar corrección de baseline
# -----------------------------------------------------------------------------
println("🔄 APLICANDO CORRECCIÓN DE BASELINE")
println("-" ^ 80)

# Crear copia de los datos con artifact rejection
eeg_2nd_baseline_corrected = copy(eeg_artifact_rejected)

# Aplicar segunda corrección de baseline a cada segmento
for ch_idx in 1:n_channels
    for seg_idx in 1:n_segments
        # Obtener el segmento
        segment = eeg_2nd_baseline_corrected[ch_idx, :, seg_idx]
        
        # Calcular la media del intervalo de baseline
        baseline_interval = segment[baseline_start_sample:baseline_end_sample]
        baseline_mean = mean(baseline_interval)
        
        # Restar la media del baseline de todo el segmento
        eeg_2nd_baseline_corrected[ch_idx, :, seg_idx] = segment .- baseline_mean
    end
    
    if ch_idx % 10 == 0 || ch_idx == n_channels
        println("  ✓ Canal $ch_idx/$n_channels procesado: $(channels[ch_idx])")
    end
end

println()
println("✔ Segunda corrección de baseline completada")
println("  → Dimensiones de la hipermatriz corregida: $(size(eeg_2nd_baseline_corrected))")
println()

# -----------------------------------------------------------------------------
# 4. Verificación y estadísticas
# -----------------------------------------------------------------------------
println("📊 ESTADÍSTICAS DE CORRECCIÓN DE BASELINE")
println("-" ^ 80)

# Calcular estadísticas antes y después de la corrección
stats_before = DataFrame(
    Segmento = 1:n_segments,
    Media_µV = [mean(eeg_artifact_rejected[:, :, seg]) for seg in 1:n_segments],
    Std_µV = [std(vec(eeg_artifact_rejected[:, :, seg])) for seg in 1:n_segments],
    Media_Baseline_µV = [mean(eeg_artifact_rejected[:, baseline_start_sample:baseline_end_sample, seg]) for seg in 1:n_segments]
)

stats_after = DataFrame(
    Segmento = 1:n_segments,
    Media_µV = [mean(eeg_2nd_baseline_corrected[:, :, seg]) for seg in 1:n_segments],
    Std_µV = [std(vec(eeg_2nd_baseline_corrected[:, :, seg])) for seg in 1:n_segments],
    Media_Baseline_µV = [mean(eeg_2nd_baseline_corrected[:, baseline_start_sample:baseline_end_sample, seg]) for seg in 1:n_segments]
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
println("📈 VISUALIZACIÓN DE EJEMPLO")
println("-" ^ 80)

# Visualizar primeros 3 segmentos del primer canal
n_segments_plot = min(3, n_segments)
ch_example = 1
channel_name = channels[ch_example]

time_segment = collect(0:(segment_length_samples-1)) ./ fs  # Tiempo en segundos para un segmento

p_comparison = plot(
    xlabel = "Tiempo (s)",
    ylabel = "Amplitud (µV)",
    title = "2ª Corrección de Baseline - Canal $channel_name",
    legend = :outerright,
    size = (1200, 500)
)

for seg_idx in 1:n_segments_plot
    segment_before = eeg_artifact_rejected[ch_example, :, seg_idx]
    segment_after = eeg_2nd_baseline_corrected[ch_example, :, seg_idx]
    
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
println("💾 GUARDANDO RESULTADOS")
println("-" ^ 80)

# Crear directorio de baseline si no existe
dir_baseline = stage_dir(:baseline)
if !isdir(dir_baseline)
    mkpath(dir_baseline)
    println("  ✓ Directorio creado: $dir_baseline")
end

# Guardar hipermatriz con segunda corrección de baseline
path_eeg_2nd_baseline = joinpath(dir_baseline, "eeg_2nd_baseline_correction.bin")
Serialization.serialize(path_eeg_2nd_baseline, eeg_2nd_baseline_corrected)
println("  ✓ Hipermatriz con 2ª corrección de baseline guardada en: $(basename(path_eeg_2nd_baseline))")

# Guardar información de segunda corrección de baseline
dict_2nd_baseline_info = Dict(
    "eeg_2nd_baseline_corrected" => eeg_2nd_baseline_corrected,
    "eeg_artifact_rejected_original" => eeg_artifact_rejected,  # Guardar también el original para comparación
    "channels" => channels,
    "fs" => fs,
    "segment_length_samples" => segment_length_samples,
    "n_segments" => n_segments,
    "n_channels" => n_channels,
    "baseline_start_s" => baseline_start_s,
    "baseline_end_s" => baseline_end_s,
    "baseline_start_sample" => baseline_start_sample,
    "baseline_end_sample" => baseline_end_sample,
    "baseline_type" => "segment-based",
    "correction_method" => "subtract_baseline_mean",
    "correction_number" => 2
)

path_dict_2nd_baseline = joinpath(dir_baseline, "dict_2nd_baseline_correction.bin")
Serialization.serialize(path_dict_2nd_baseline, dict_2nd_baseline_info)
println("  ✓ Información de 2ª corrección de baseline guardada en: $(basename(path_dict_2nd_baseline))")

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

path_stats_csv = joinpath(dir_baseline, "2nd_baseline_correction_statistics.csv")
CSV.write(path_stats_csv, stats_comparison)
println("  ✓ Estadísticas de 2ª corrección guardadas en: $(basename(path_stats_csv))")
println()

# -----------------------------------------------------------------------------
# 7. Resumen final
# -----------------------------------------------------------------------------
println("=" ^ 80)
println("✨ 2ª CORRECCIÓN DE BASELINE COMPLETADA")
println("=" ^ 80)
println()
println("📋 RESUMEN:")
println("  • Hipermatriz corregida: $(size(eeg_2nd_baseline_corrected))")
println("  • Formato: (canales × muestras × segmentos)")
println("  • Canales: $n_channels")
println("  • Muestras por segmento: $segment_length_samples")
println("  • Número de segmentos: $n_segments")
println("  • Intervalo de baseline: $(baseline_start_s) s - $(baseline_end_s) s ($(baseline_start_s*1000) ms - $(baseline_end_s*1000) ms)")
println("  • Muestras de baseline: $(baseline_start_sample) - $(baseline_end_sample)")
println("  • Método: Resta de la media del intervalo de baseline")
println()
println("💾 Archivos guardados en: $dir_baseline")
println("  → eeg_2nd_baseline_correction.bin")
println("  → dict_2nd_baseline_correction.bin")
println("  → 2nd_baseline_correction_statistics.csv")
println()
