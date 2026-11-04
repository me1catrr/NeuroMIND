### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ a2fca222-2c13-11f1-bfbb-ffe9fd9d5381
begin
using PlutoUI
PlutoUI.TableOfContents(title = "Contenido")
end

# ╔═╡ 4a417ce1-7265-4a51-954f-1dacb925486d
begin
include(joinpath(@__DIR__, "..", "_template_base.jl"))
using .PlutoTemplateBase
using CSV, DataFrames, Serialization, Statistics, StatsBase, DSP, Plots, Dates, InlineStrings, LinearAlgebra, Random, FFTW
end

# ╔═╡ b17f8e3c-8437-4658-baee-9fccb55991ea
begin
# Si se ejecuta este script directamente (fuera del módulo EEG_Julia),
# cargamos utilidades de rutas para disponer de `stage_dir`.
if !@isdefined(stage_dir)
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end
end

# ╔═╡ 9e06da1a-8985-43ed-a13e-4b0b4e34bbdd
md"""
**PAQUETES CARGADOS**
"""

# ╔═╡ 8c8b3d9d-a425-4c48-8fdc-ac987a3cdb39
notebook_intro("PROCESSING")

# ╔═╡ 9139e5ac-7c0a-43cd-aed4-44eb4eb3b115
md"""
## Segmentación

After ICA cleaning (Phase 2b), the continuous EEG is free of the main artefactual sources and is ready to be divided into fixed-length epochs for spectral and connectivity analysis.

Segmentation produces a 3D hypermatrix

```math
\mathbf{E} \in \mathbb{R}^{C \times L \times K}
```

channels × samples per segment × segments.

Each segment is a contiguous block of

```math
L = F_s \cdot T_{\mathrm{seg}}
```

samples (e.g. $T_{\mathrm{seg}} = 1\,\text{s}$ at $F_s = 500\,\text{Hz}$ gives $L = 500$).

The step between consecutive segment onsets is

```math
\text{step} = L - \text{overlap}
```

With zero overlap, $\text{step} = L$ and segments are statistically independent, which simplifies FFT-based spectral estimation and later averaging.

The number of full segments is

```math
K = \left\lfloor \frac{N_{\mathrm{total}} - L}{\text{step}} \right\rfloor + 1
```

where $N_{\mathrm{total}}$ is the length of the continuous signal; any trailing samples that do not form a complete segment are discarded.

The 1 s segment length balances spectral resolution ($\approx F_s/L \approx 1\,\text{Hz}$, adequate for band-based connectivity) with the assumption of weak stationarity within each segment, which is required for stable cross-spectral and wPLI estimation.

The implementation is in `src/segmentation.jl`, which loads `data/ICA/dict_EEG_ICA_clean.bin`, applies the splitting, computes per-segment statistics (mean, std, min, max, RMS), generates example plots (overlaid and stacked segments), and saves the hypermatrix and metadata to `data/segmentation/`.

Algorithm 1 and Table 1 summarise the step; the box below formalises Pipeline Phase 3.

!!! info "Algoritmo Segmentación en épocas de longitud fija (`segmentation.jl`)"

    **Entrada:**  
    `dict_EEG` (diccionario de canales con EEG limpio tras ICA),  
    `F_s`, `T_seg`, `overlap`

    1. Cálculo de parámetros de segmentación

       ```text
       L      ← floor(F_s · T_seg)
       step   ← L − overlap
       N_total ← length(dict_EEG[channels[1]])
       ```

    2. Número de segmentos

       ```text
       K ← floor((N_total − L) / step) + 1
       ```

    3. Inicialización de la hipermatriz

       ```text
       Allocate E ∈ ℝ^{C × L × K}
       ```

    4. Segmentación por canal

       ```text
       for c = 1 to C

           signal ← dict_EEG[channels[c]]

           for k = 1 to K

               start ← (k − 1) · step + 1
               end   ← start + L − 1

               if end ≤ N_total
                   E[c,:,k] ← signal[start : end]
               else
                   copy available samples
                   pad or truncate if needed
               end

           end

       end
       ```

    **Salida:**  
    `E` (hipermatriz de segmentos),  
    `K` (número de segmentos),  
    `L` (muestras por segmento),  
    `step` (desplazamiento entre segmentos)

    Opcionalmente: estadísticas por segmento  
    `(mean, std, min, max, RMS)`

!!! info "Pipeline Fase 3: Segmentación"

    **Objetivo.**
    - Dividir el EEG continuo limpio tras ICA en segmentos de longitud fija sin solapamiento.
    - Construir una hipermatriz tridimensional (canales × muestras × segmentos).
    - Generar metadatos y estadísticas necesarias para las fases posteriores de corrección de baseline y estimación de conectividad.

    **Entrada.**
    - EEG limpio tras ICA: `data/ICA/dict_EEG_ICA_clean.bin`  
      (señales continuas indexadas por canal).

    - Parámetros de segmentación:

      ```math
      F_s = 500\,\text{Hz}
      ```

      ```math
      T_{\mathrm{seg}} = 1\,\text{s}
      ```

      ```math
      \text{overlap} = 0\,\text{s}
      ```

    **Salida.**
    - `data/segmentation/eeg_segmented.bin`: hipermatriz

      ```math
      \mathbf{E} \in \mathbb{R}^{C \times L \times K}
      ```

      donde

      ```math
      L = 500
      ```

      y

      ```math
      K = \left\lfloor \frac{N_{\mathrm{total}} - L}{\text{step}} \right\rfloor + 1
      ```

    - `data/segmentation/dict_segmentation_info.bin`: diccionario con:
        - `eeg_segmented`
        - `channels`
        - `fs`
        - `segment_length_s`
        - `segment_length_samples`
        - `segment_overlap_s`
        - `segment_step_samples`
        - `n_segments`
        - `n_channels`
        - `n_samples_total`

    - `data/segmentation/segment_statistics.csv`: estadísticas por segmento  
      (media, desviación estándar, mínimo, máximo y RMS) para verificación.

    - Opcional: gráficos de control de calidad  
      (primer canal con segmentos superpuestos o segmentos apilados).

    **Pasos de procesamiento (alineados con `src/segmentation.jl`).**

    1. **Carga de datos.**  
       Cargar `dict_EEG_ICA_clean.bin` y obtener:

       - `channels`
       - número total de muestras

       ```math
       N_{\mathrm{total}}
       ```

    2. **Definición de parámetros de segmentación.**  
       Calcular la longitud del segmento:

       ```math
       L = F_s \cdot T_{\mathrm{seg}}
       ```

       Definir el paso entre segmentos (sin solapamiento):

       ```math
       \text{step} = L
       ```

       Calcular el número total de segmentos:

       ```math
       K = \left\lfloor \frac{N_{\mathrm{total}} - L}{\text{step}} \right\rfloor + 1
       ```

    3. **Construcción de la hipermatriz.**  
       Reservar memoria para

       ```math
       \mathbf{E} \in \mathbb{R}^{C \times L \times K}
       ```

       Para cada canal \(c\) y segmento \(k\), copiar el bloque correspondiente de señal:

       ```math
       \mathbf{E}[c,:,k] =
       \text{signal}[(k-1)\text{step}+1 : (k-1)\text{step}+L]
       ```

    4. **Cálculo de estadísticas.**  
       Para cada segmento calcular:

       - media  
       - desviación estándar  
       - mínimo  
       - máximo  
       - RMS  

       Guardar los resultados en `segment_statistics.csv`.

    5. **Visualización (opcional).**  
       Generar gráficos de control de calidad, por ejemplo:

       - primeros segmentos superpuestos del primer canal  
       - todos los segmentos apilados  

    6. **Guardado de resultados.**  
       - Serializar la hipermatriz en `eeg_segmented.bin`.
       - Construir el diccionario `dict_segmentation_info` con los metadatos y guardarlo en `dict_segmentation_info.bin`.

       Los datos resultantes se utilizan como entrada para la **Fase 4: corrección de baseline**.
"""

# ╔═╡ 91f7264e-9f4e-4ea1-95c0-77d9de3e6bbd
md"
## Carga de datos (ICA cleaning)

Se cargan los datos **EEG** que ya han sido limpiados mediante **ICA**. Los datos están en formato de diccionario (un vector por canal) y representan señales continuas que serán divididas en segmentos temporales.
"

# ╔═╡ 30ff2319-076f-43c2-a134-9016be56e46c
begin
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
end

# ╔═╡ c220c177-7b8e-40cd-8d99-4aae100de7fc
begin
# Obtener canales ordenados (orden alfabético para consistencia)
# Todos los canales deben tener la misma longitud (n_samples_total)
channels = sort(collect(keys(dict_EEG_ICA_clean)))
n_channels = length(channels)
n_samples_total = length(dict_EEG_ICA_clean[channels[1]])  # longitud de la señal continua

println("  → Número de canales: $n_channels")
println("  → Número total de muestras: $n_samples_total")
println()
end

# ╔═╡ 37f24271-0da7-4236-9d84-6210a1e5ac0f
md"
## Configuración de segmentación

Se configuran los parámetros para dividir la señal continua en segmentos:
   - Longitud del segmento: determina la duración de cada época
   - Solapamiento: permite que segmentos se superpongan (0 = sin solapamiento)
   - Paso entre segmentos: distancia entre el inicio de segmentos consecutivos

**Configuración actual:**

   - Segmentos de 1 segundo (500 muestras a 500 Hz)
   - Sin solapamiento (segmentos consecutivos e independientes)
   - Resolución espectral resultante: ~1 Hz (fs / longitud_segmento)
"

# ╔═╡ a83800b1-9df8-4536-9812-525bcad1902a
begin
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
println("-" ^ 50)
println("Frecuencia de muestreo: $(fs) Hz")
println("Longitud del segmento: $(segment_length_s) s ($(segment_length_samples) muestras)")
println("Solapamiento: $(segment_overlap_s) s ($(segment_overlap_samples) muestras)")
println("Paso entre segmentos: $(segment_step_samples) muestras")
println("Saltar intervalos malos: $(skip_bad_intervals)")
println()
end

# ╔═╡ 08cdbdf3-e827-4aea-a40b-6f157fe0ba5b
begin
# Calcular número de segmentos completos que caben en la señal
# Fórmula: número de segmentos = floor((total - longitud) / paso) + 1
# Esto asegura que todos los segmentos tengan la longitud completa
n_segments = Int(floor((n_samples_total - segment_length_samples) / segment_step_samples)) + 1

println("📐 CÁLCULO DE SEGMENTOS")
println("-" ^ 40)
println("Número de segmentos completos: $n_segments")
println("Muestras por segmento: $segment_length_samples")
println("Muestras totales utilizadas: $(n_segments * segment_step_samples)")
# Las muestras finales que no forman un segmento completo se descartan
if n_segments * segment_step_samples < n_samples_total
    muestras_no_usadas = n_samples_total - (n_segments * segment_step_samples)
    println("Muestras no utilizadas (final): $muestras_no_usadas")
end
println()
end

# ╔═╡ 00b7e07b-b25b-4a72-9a96-59e7f18430bb
md"
## Segmentación 

Crear hipermatriz canales × muestras × segmentos

Esta sección divide cada canal en segmentos consecutivos y los organiza en una hipermatriz 3D. Cada **rebanada** de la tercera dimensión es un segmento temporal completo con todos los canales.

PROCESO:

Para cada canal:

- Extrae la señal continua del diccionario
- Calcula los índices de inicio y fin de cada segmento
- Copia cada segmento a la hipermatriz
- Los segmentos son consecutivos (sin solapamiento)

RESULTADO:

Hipermatriz de tamaño (n_channels × segment_length_samples × n_segments)

   - Primera dimensión: canales
   - Segunda dimensión: muestras temporales dentro del segmento
   - Tercera dimensión: segmentos (épocas)
"

# ╔═╡ f6f6044b-2d6c-4116-b4c9-529cdc7c7d1f
begin
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
end

# ╔═╡ 98d2ba11-2456-4a9b-9e0a-ab20b5bd276d
md"
## Verificación de segmentos

Esta sección calcula estadísticas descriptivas para cada segmento. Estas estadísticas ayudan a verificar que la segmentación se realizó correctamente y proporcionan información sobre la variabilidad entre segmentos.

**Estadísticas calculadas**:

   - **Media**: valor promedio de amplitud en el segmento
   - **Std**: desviación estándar (variabilidad)
   - **Min/Max**: valores extremos de amplitud
   - **RMS**: raíz cuadrada de la media de los cuadrados (medida de potencia)
"

# ╔═╡ 7c15972a-dde7-45f6-8d18-b62da3587cc0
begin
println("📊 ESTADÍSTICAS DE SEGMENTACIÓN")
println("-" ^ 70)

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
end

# ╔═╡ 212b3eb6-3a86-4e41-a230-6c8b92147573
md"
## Visualización de ejemplo

Esta sección genera visualizaciones para inspeccionar los segmentos creados:

1. **Gráfico superpuesto**: primeros segmentos del primer canal
      - Permite ver la forma de los segmentos y verificar que son consecutivos
2. **Gráfico apilado**: todos los segmentos del primer canal
      - Muestra la evolución temporal a lo largo de toda la señal
      - Útil para detectar patrones o artefactos que se repiten
"

# ╔═╡ cec9cd67-27db-4736-85c5-2af2a0093364
begin
println("📈 VISUALIZACIÓN DE EJEMPLO")
println("-" ^ 60)

# Visualizar primeros 5 segmentos del primer canal como ejemplo
n_segments_plot = min(5, n_segments)  # mostrar máximo 5 segmentos
ch_example = 1                        # primer canal
channel_name = channels[ch_example]

# Crear vector de tiempo para el eje X (en segundos para un segmento)
time_segment = collect(0:(segment_length_samples-1)) ./ fs
end

# ╔═╡ 7c3becba-7fc8-4eaf-86c6-2d92e08b04fa
begin
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

p_segments
end

# ╔═╡ 9567bffe-65b2-4073-ae85-ee339577f9be
begin
# Visualización de todos los segmentos apilados (primer canal)
println("  → Visualización de todos los segmentos apilados (canal $channel_name)")
sep = 50.0  # Separación vertical entre segmentos
offsets_segments = [(n_segments - i) * sep for i in 1:n_segments]
end

# ╔═╡ a892d381-a5a7-4659-bb5b-19e65183522f
begin
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

p_stacked

end

# ╔═╡ 7f0a28b2-2cd9-4910-8312-dad488e40399
md"
## Guardar resultados

Esta sección guarda todos los resultados del proceso de segmentación:
   1. **Hipermatriz segmentada:** datos EEG organizados en segmentos
   2. **Diccionario de información:** metadatos y parámetros de segmentación
   3. **Estadísticas en CSV:** tabla con estadísticas por segmento
"

# ╔═╡ bd6dd4b5-ed90-4b30-9c0f-a04f7a8adcf7
begin
println("💾 GUARDANDO RESULTADOS")
println("-" ^ 60)

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
end

# ╔═╡ a877c6a1-25b2-43b3-9191-6da2522b65e9
md"
## Resumen segmentación
"

# ╔═╡ 56795636-501e-4928-a722-d515eedb4dd3
begin
println("=" ^ 60)
println("✨ SEGMENTACIÓN COMPLETADA")
println("=" ^ 60)
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
end

# ╔═╡ 720b30e4-bc60-48f5-99bb-d1c7fee7ae57
md"""
!!! info "Algoritmo Corrección de baseline por segmento (`baseline.jl`)"

    **Entrada:**  
    ``E ∈ ℝ^{C × L × K}`` (EEG segmentado),  
    `F_s`, `t_start`, `t_end`

    1. Cálculo de índices del intervalo de baseline

       ```text
       n_start ← floor(t_start · F_s) + 1
       n_end   ← floor(t_end   · F_s)
       ```

    2. Inicialización de la señal corregida

       ```text
       E_corr ← copy(E)
       ```

    3. Corrección por canal y segmento

       ```text
       for c = 1 to C
           for k = 1 to K

               segment ← E_corr[c,:,k]

               baseline_interval ← segment[n_start : n_end]

               μ ← mean(baseline_interval)

               E_corr[c,:,k] ← segment − μ

           end
       end
       ```

    **Salida:**  
    `E_corr` (EEG segmentado con corrección de baseline aplicada)

!!! info "Pipeline Fase 4: Corrección de baseline"

    **Objetivo.**
    - Eliminar el desplazamiento DC de cada segmento restando la media de un intervalo de baseline.
    - Utilizar como baseline los primeros 100 ms de cada segmento.
    - Generar EEG segmentado corregido y estadísticas de verificación para los análisis posteriores de conectividad y espectro.

    **Entrada.**
    - EEG segmentado: `data/segmentation/dict_segmentation_info.bin`, que contiene:
        - `eeg_segmented` (hipermatriz \(C \times L \times K\))
        - `channels`
        - `fs`
        - `segment_length_samples`
        - `n_segments`
        - `n_channels`

    - Ventana de baseline:

      ```math
      t_{\mathrm{start}} = 0\,\text{s}, \quad t_{\mathrm{end}} = 0.1\,\text{s}
      ```

      (primeros **100 ms** de cada segmento).

    **Salida.**
    - `data/baseline/eeg_1st_baseline_correction.bin`: hipermatriz corregida

      ```math
      \mathbf{E}^{\mathrm{corr}}
      ```

      (canales × muestras × segmentos).

    - `data/baseline/dict_1st_baseline_correction.bin`: diccionario con:
        - `eeg_baseline_corrected`
        - `eeg_segmented_original`
        - `channels`
        - `fs`
        - parámetros de segmentación y baseline
        - `baseline_type = "segment-based"`
        - `correction_method = "subtract_baseline_mean"`

    - `data/baseline/baseline_correction_statistics.csv`: estadísticas por segmento  
      (media, desviación estándar, media del baseline antes y después de la corrección y diferencia).

    **Pasos de procesamiento (alineados con `src/baseline.jl`).**

    1. **Carga de datos.**  
       Cargar `dict_segmentation_info.bin` y extraer:

       - `eeg_segmented`
       - `channels`
       - `fs`
       - `segment_length_samples`
       - `n_segments`
       - `n_channels`

    2. **Cálculo de índices de baseline.**  
       Determinar los índices correspondientes al intervalo de baseline:

       ```math
       n_{\mathrm{start}} = \lfloor t_{\mathrm{start}} F_s \rfloor + 1
       ```

       ```math
       n_{\mathrm{end}} = \lfloor t_{\mathrm{end}} F_s \rfloor
       ```

       Por ejemplo, con

       ```math
       F_s = 500\,\text{Hz}
       ```

       el intervalo \([0, 0.1]\) s corresponde a las muestras **1–50**.

    3. **Corrección de baseline.**  
       Para cada canal \(c\) y segmento \(k\), calcular la media del baseline:

       ```math
       \mu_{c,k} = \text{mean}(\mathbf{E}[c, n_{\mathrm{start}}:n_{\mathrm{end}}, k])
       ```

       y restarla a todo el segmento:

       ```math
       \mathbf{E}^{\mathrm{corr}}[c,:,k] = \mathbf{E}[c,:,k] - \mu_{c,k}
       ```

    4. **Verificación.**  
       Calcular por segmento:

       - media  
       - desviación estándar  
       - media del baseline antes de la corrección  
       - media del baseline después de la corrección  

       Confirmar que la media del baseline posterior sea aproximadamente

       ```math
       0
       ```

       Opcionalmente generar un gráfico comparativo (por ejemplo, primer canal y primeros segmentos) con el intervalo de baseline marcado.

    5. **Guardado de resultados.**  
       - Serializar `eeg_baseline_corrected` en `eeg_1st_baseline_correction.bin`.
       - Construir el diccionario `dict_1st_baseline_correction` y guardarlo en `dict_1st_baseline_correction.bin`.
       - Exportar las estadísticas a `baseline_correction_statistics.csv`.

       Los datos resultantes se utilizan como entrada para la **Fase 5: rechazo de artefactos**.
"""

# ╔═╡ 381398ab-95f6-4fe2-81c6-a0478915c7c0
md"
## Configuración baseline

Se configura el **intervalo de baseline** que se usará para calcular la media. Este intervalo típicamente corresponde al inicio de cada segmento, antes de cualquier actividad de interés. La media de este intervalo se resta de todo el segmento para eliminar el **offset DC**.

**Configuración**:
   - **Intervalo**: 0.00 s - 0.10 s (primeros 100 ms de cada segmento)
   - A 500 Hz: muestras 1-50 (primeras 50 muestras)
"

# ╔═╡ 7a7fd84b-36f0-4a7f-aacf-338a6a523311
begin
baseline_start_s = 0.00    # Inicio del intervalo de baseline (s)
                            # Comienza al inicio del segmento (t=0)
baseline_end_s = 0.10       # Fin del intervalo de baseline (s) = 100 ms
                            # Primeros 100 ms del segmento

# Convertir tiempos a índices de muestras (1-indexed en Julia)
baseline_start_sample = Int(round(baseline_start_s * fs)) + 1  # Muestra 1 (índice base 1)
baseline_end_sample = Int(round(baseline_end_s * fs))          # Muestra 50 (a 500 Hz)

println("⚙️  CONFIGURACIÓN DE BASELINE")
println("-" ^ 60)
println("  Intervalo de baseline: $(baseline_start_s) s - $(baseline_end_s) s")
println("  Muestras de baseline: $(baseline_start_sample) - $(baseline_end_sample)")
println("  Número de muestras en baseline: $(baseline_end_sample - baseline_start_sample + 1)")
println()
end

# ╔═╡ 584aed06-bed9-4158-8b2d-11f8fe72b693
md"
## Aplicar corrección baseline

Esta sección aplica la **corrección de baseline** a cada segmento de cada canal.

**Proceso:**

Para cada canal y cada segmento:
     1. Extrae el segmento completo
     2. Identifica el intervalo de baseline (primeros 100 ms)
     3. Calcula la media de ese intervalo
     4. Resta la media de todo el segmento

**Resultado:**

Después de la corrección, el intervalo de baseline tiene media ≈ 0 en cada segmento. Esto elimina el **offset DC** y normaliza la señal, facilitando análisis posteriores.
"

# ╔═╡ aa01d3bf-16c6-4436-9eb6-55f737e4fdba
begin
println("🔄 APLICANDO CORRECCIÓN DE BASELINE")
println("-" ^ 60)
	
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
println("  → Dimensiones de la hipermatriz corregida $(size(eeg_baseline_corrected))")
println()
end

# ╔═╡ 0be498cd-14f4-45d1-8028-40a793716b0f
md"
## Verificación y estadísticas

Esta sección calcula estadísticas antes y después de la corrección para verificar que la corrección se aplicó correctamente. La verificación principal es que la media del intervalo de baseline después de la corrección debe ser ≈ 0.

**Estadísticas calculadas:**
   - **Media global del segmento**: antes y después
   - **Desviación estándar**: para verificar que no cambió (solo se resta una constante)
   - **Media del intervalo de baseline**: debe ser ≈ 0 después de la corrección
"

# ╔═╡ d9e86f10-3244-476f-aca1-2c26fc5b592a
begin
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
end

# ╔═╡ bc39f4d8-0192-41c4-885a-d09f65de298a
md"
## Visualización antes vs después

Esta sección genera visualizaciones para inspeccionar el efecto de la corrección:
   - Gráfico comparativo: muestra segmentos antes y después superpuestos
   - Marca el intervalo de baseline para visualizar qué se corrigió
   - Permite verificar visualmente que la corrección funcionó correctamente
"

# ╔═╡ bc400aa6-b180-4b02-8abc-60da09acb937
begin
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

p_comparison
end

# ╔═╡ 76b687c1-cfc9-4206-a9bd-3859932e9c09
md"
## Guardar resultados

Esta sección guarda todos los resultados del proceso de corrección de baseline:

   1. **Hipermatriz corregida:** datos **EEG** con baseline corregido
   2. **Diccionario de información:** metadatos y parámetros de corrección
   3. **Estadísticas en CSV:** tabla comparativa antes/después
"

# ╔═╡ 47ce252b-c3db-4f75-85c1-0fda4b34da80
begin
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

CSV.write(path_stats_csv, stats_comparison)
println("  ✓ Estadísticas de corrección guardadas en: $(basename(path_stats_csv))")
println()
end

# ╔═╡ 812069cd-02cf-4e48-8aa7-94753f6735a4
md"
## Resumen correción baseline (1ª)
"

# ╔═╡ 508362b0-0a1a-4eec-a26f-1c45db18563c
begin
println("=" ^ 60)
println("✨ CORRECCIÓN DE BASELINE COMPLETADA")
println("=" ^ 60)
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
end

# ╔═╡ c07064a7-0303-4299-a080-c45f0efd5c4f
md"""
!!! info "Algoritmo Rechazo de artefactos basado en amplitud (`artifact_rejection.jl`)"

    **Entrada:**  
    ``E_bl ∈ ℝ^{C × L × K}`` (segmentos corregidos por baseline),  
    lista de canales, `N_used`, `A_min`, `A_max`

    1. Selección de canales a inspeccionar

       ```text
       C_used ← {1, … , min(N_used, C)}
       ```

    2. Inicialización

       ```text
       segment_is_good[1:K] ← true
       initialise arrays for:
           per-segment min amplitude
           per-segment max amplitude
           violating channel lists
       ```

    3. Evaluación de cada segmento

       ```text
       for k = 1 to K

           segment_used ← E_bl[C_used,:,k]

           m_k ← min(segment_used)
           M_k ← max(segment_used)

           if (m_k < A_min) or (M_k > A_max)

               segment_is_good[k] ← false

               V_k ← { c ∈ C_used :
                       min_n(E_bl[c,n,k]) < A_min
                       or
                       max_n(E_bl[c,n,k]) > A_max }

           else

               V_k ← ∅

           end

       end
       ```

    4. Identificación de segmentos válidos

       ```text
       good_indices ← { k : segment_is_good[k] = true }
       ```

    5. Construcción de la hipermatriz filtrada

       ```text
       E_AR ← E_bl[:,:,good_indices]
       ```

    **Salida:**  
    `E_AR` (EEG con segmentos artefactuales eliminados)  
    `segment_is_good` (vector booleano por segmento)  
    estadísticas por segmento y listas de canales que violan los umbrales

!!! info "Pipeline Fase 5: Rechazo de artefactos"

    **Objetivo.**
    - Aplicar un criterio conservador basado en amplitud para descartar segmentos con artefactos evidentes.
    - Detectar saturaciones o ráfagas EMG extremas que puedan contaminar el análisis.
    - Conservar únicamente los segmentos limpios de EEG en reposo para los análisis posteriores de conectividad y espectro.

    **Entrada.**
    - EEG segmentado y corregido por baseline (Fase 4):  
      `data/baseline/dict_1st_baseline_correction.bin`, que contiene:
        - `eeg_baseline_corrected` (hipermatriz \(C \times L \times K\))
        - `channels`
        - `fs`
        - `segment_length_samples`
        - `n_segments`
        - `n_channels`

    - Configuración de rechazo de artefactos:

      ```math
      N_{\mathrm{used}} = 30
      ```

      canales inspeccionados.

      ```math
      A_{\min} = -70\,\mu\text{V}, \quad A_{\max} = +70\,\mu\text{V}
      ```

      Umbrales de amplitud.

      Los parámetros `before_event_ms` y `after_event_ms` se almacenan para posible uso futuro, pero no se aplican en esta implementación.

    **Salida.**
    - `data/artifact_rejection/eeg_artifact_rejected.bin`: hipermatriz filtrada con solo segmentos válidos

      ```math
      C \times L \times K_{\mathrm{good}}
      ```

    - `data/artifact_rejection/dict_artifact_rejection.bin`: diccionario con:
        - hipermatriz original y filtrada
        - lista de canales
        - umbrales de amplitud
        - índices de segmentos válidos
        - máscara booleana de segmentos buenos
        - conteos por segmento

    - `data/artifact_rejection/artifact_rejection_statistics.csv`: estadísticas por segmento (mínimos, máximos, flags y canales que violan los umbrales).

    **Pasos de procesamiento (alineados con `src/artifact_rejection.jl`).**

    1. **Carga de datos.**  
       Cargar `dict_1st_baseline_correction.bin` y extraer:

       - `eeg_baseline_corrected`
       - `channels`
       - `fs`
       - `segment_length_samples`
       - `n_segments`
       - `n_channels`

    2. **Configuración de parámetros.**  
       Definir:

       ```math
       N_{\mathrm{used}},\; A_{\min},\; A_{\max}
       ```

       Calcular equivalentes en muestras para `before_event_ms` y `after_event_ms` (para uso futuro).  
       Seleccionar `channels_used` como los primeros

       ```math
       N_{\mathrm{used}}
       ```

       canales.

    3. **Detección de artefactos.**  
       Para cada segmento \(k\):

       - calcular

       ```math
       m_k = \min(E_{c,n,k})
       ```

       ```math
       M_k = \max(E_{c,n,k})
       ```

       considerando únicamente `channels_used`.

       El segmento se marca como inválido si

       ```math
       m_k < A_{\min} \quad \text{o} \quad M_k > A_{\max}
       ```

       Registrar los canales que violan los umbrales.

    4. **Cálculo de estadísticas.**  
       Construir un `DataFrame` con:

       - amplitud mínima por segmento  
       - amplitud máxima por segmento  
       - flag de segmento válido/inválido  
       - canales que superan los umbrales  

       Calcular también:

       - número de segmentos conservados  
       - número de segmentos rechazados  
       - porcentaje de rechazo  
       - mínimos y máximos globales.

    5. **Filtrado de segmentos.**  
       Construir `good_segment_indices` a partir de `segment_is_good`.

       Crear la hipermatriz filtrada:

       ```math
       \mathbf{E}^{\mathrm{AR}}
       ```

       seleccionando únicamente los segmentos válidos a lo largo de la tercera dimensión.

    6. **Visualización (opcional).**  
       Generar:

       - gráficos de segmentos de ejemplo con umbrales superpuestos  
       - histogramas de amplitudes máximas y mínimas por segmento con los umbrales marcados.

    7. **Guardado de resultados.**  
       - Serializar `eeg_artifact_rejected` en `eeg_artifact_rejected.bin`.
       - Construir el diccionario `dict_artifact_rejection` con todos los metadatos y guardarlo en `dict_artifact_rejection.bin`.
       - Exportar las estadísticas a `artifact_rejection_statistics.csv`.

       Los datos resultantes se utilizan como entrada para la **Fase 6: segunda corrección opcional de baseline**.
"""

# ╔═╡ 13b382df-1a09-4478-81a3-7aa9d4ab62f7
md"
## Carga de datos (Baseline 1ª)
"

# ╔═╡ 16480b91-338f-4f7f-b4bd-89f6e73bde55
begin
println("✔ Datos con baseline correction cargados desde: $(basename(path_dict_baseline))")
println("  → Dimensiones: $(size(eeg_baseline_corrected))")
println("  → Formato: (canales × muestras × segmentos)")
println("  → Número de canales: $n_channels")
println("  → Muestras por segmento: $segment_length_samples")
println("  → Número de segmentos: $n_segments")
println("  → Frecuencia de muestreo: $(fs) Hz")
println()
end

# ╔═╡ 021b17fa-84b6-4839-96e6-a5710fae3533
md"
## Configuración artifact rejection

Se configuran los parámetros para la **detección de artefactos**:
   - **Umbrales de amplitud:** valores fuera de este rango se consideran artefactos
   - **Canales a usar:** solo se evalúan los primeros N canales (por defecto 30)
   - **Ventanas temporales:** parámetros para rechazo alrededor de eventos
     (actualmente no se usan en la detección, solo se guardan para referencia)
"

# ╔═╡ e790286e-3a1d-40a0-b323-cd3d9a81ed20
begin
n_channels_used = 30    # Número de canales a usar para artifact rejection
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
end

# ╔═╡ 35f75188-5502-443e-b78a-c0dd7c91f00f
md"
## Detección artefactos por segmento

Esta sección evalúa cada segmento independientemente para **detectar artefactos**.

**Criterio de rechazo:**

Un segmento se marca como **malo** si CUALQUIER valor en los canales evaluados excede los umbrales de amplitud (±70 µV). 

Esto es un criterio conservador: si hay un artefacto en cualquier canal, se rechaza todo el segmento.

**Proceso:**
   1. Para cada segmento, extrae los datos de los canales usados
   2. Calcula el mínimo y máximo de amplitud en esos canales
   3. Compara con los umbrales configurados
   4. Si excede, marca el segmento como **malo** e identifica qué canales violaron
   5. Almacena estadísticas para análisis posterior
"

# ╔═╡ 0551b30d-150d-4050-9981-c52e6dba2630
begin
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
end

# ╔═╡ 7cd221d4-a26c-4be7-a29b-fcfbff3b8938
md"
## Estadísticas artifact rejection

Esta sección calcula y muestra estadísticas sobre el proceso de rechazo:
   - Número de segmentos mantenidos vs eliminados
   - Distribución de amplitudes
   - Detalles de segmentos eliminados
   - Creación de DataFrame con estadísticas completas
"

# ╔═╡ 2fc078ca-f662-4d67-b564-7c7ed703b629
begin
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
end

# ╔═╡ f56d927a-74d2-48b7-b7f3-efc9ed18e006
md"
## Filtrar segmentos

Esta sección crea una nueva hipermatriz que contiene únicamente los segmentos que pasaron el filtro de artifact rejection. Los segmentos marcados como **malos** se eliminan completamente de los datos.

El resultado es una hipermatriz con las mismas dimensiones en canales y muestras, pero con menos segmentos (solo los buenos).
"

# ╔═╡ 48f841e8-b704-414a-bbb1-d002a6e79755
begin
println("🔄 FILTRANDO SEGMENTOS")
println("-" ^ 60)

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
end

# ╔═╡ 487d2c06-70ff-4537-83f6-7bd07cfa274b
md"
## Visualización de ejemplo

Esta sección genera visualizaciones para inspeccionar los resultados:

   1. **Gráfico temporal:** muestra algunos segmentos con los umbrales marcados 
      - Segmentos buenos en azul
      - Segmentos malos en rojo

   2. **Histogramas:** distribución de amplitudes máximas y mínimas por segmento
      - Permite ver si los umbrales son apropiados
      - Muestra dónde están los umbrales en relación a la distribución
"

# ╔═╡ 77688f50-5b5c-4a41-89d5-1c19c0f2590e
begin
# Visualizar algunos segmentos (buenos y malos si hay)
# Se muestran los primeros 5 segmentos como ejemplo

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

p_artifact
end

# ╔═╡ 9715fadc-9d3b-45d4-a2f7-67e0a852be80
begin
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

p_hist
end

# ╔═╡ 37c074eb-6141-42b9-9a6b-dab79d781246
md"
# Corrección de baseline (2ª)

Una segunda **corrección de baseline** (opcional) puede aplicarse a los segmentos que han pasado el rechazo de artefactos (salida de la Fase 5) si el análisis lo requiere; este paso se formaliza como la Fase 6 del pipeline.

Si se utiliza, se aplica el mismo procedimiento basado en segmentos que en la Fase 4 (restar, para cada segmento, la media de una ventana de baseline) al conjunto reducido de segmentos.

El resultado se utiliza posteriormente como entrada para la estimación de conectividad (Fase 7).
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
InlineStrings = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CSV = "~0.10.16"
DSP = "~0.8.4"
DataFrames = "~1.8.1"
FFTW = "~1.10.0"
InlineStrings = "~1.4.5"
Plots = "~1.41.6"
PlutoUI = "~0.7.79"
StatsBase = "~0.34.10"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.3"
manifest_format = "2.0"
project_hash = "58c04f85bb9bb7c8ac4651905c2a38bdaaf6dbbb"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "8d8e0b0f350b8e1c91420b5e64e5de774c2f0f4d"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.16"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a21c5464519504e41e0cbc91f0188e8ca23d7440"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.3.0+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "5989debfc3b38f736e69724818210c67ffee4352"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.4"

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

    [deps.DSP.weakdeps]
    OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d8928e9169ff76c6281f39a659f9bca3a573f24c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.1"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "66381d7059b5f3f6162f28831854008040a4e905"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InlineStrings]]
git-tree-sha1 = "8f3d257792a522b4601c24a577954b0a8cd7334d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.5"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.15.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "OpenSSL_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.0+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "OpenSSL_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.3+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ff69a2b1330bcb730b9ac1ab7dd680176f5896b8"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.1010+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2025.5.20"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.29+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.4+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.44.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.12.1"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "cb20a4eacda080e517e4deb9cfb6c7c518131265"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.6"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "Setfield", "SparseArrays"]
git-tree-sha1 = "2d99b4c8a7845ab1342921733fa29366dae28b24"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.1"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"
    PolynomialsRecipesBaseExt = "RecipesBase"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "07a921781cab75691315adc645096ed5e370cb77"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "211530a7dc76ab59087f4d4d1fc3f086fbe87594"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.2.3"

    [deps.PrettyTables.extensions]
    PrettyTablesTypstryExt = "Typstry"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.12.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5acc6a41b3082920f79ca3c759acbcecf18a8d78"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "d05693d339e37d6ab134c5ab53c29fce5ee5d7d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.4"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "28145feabf717c5d65c1d5e09747ee7b1ff3ed13"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.3"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.1+2"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e015f211ebb898c8180887012b938f3851e719ac"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.55+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.64.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.7.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╠═a2fca222-2c13-11f1-bfbb-ffe9fd9d5381
# ╠═9e06da1a-8985-43ed-a13e-4b0b4e34bbdd
# ╠═4a417ce1-7265-4a51-954f-1dacb925486d
# ╠═8c8b3d9d-a425-4c48-8fdc-ac987a3cdb39
# ╠═b17f8e3c-8437-4658-baee-9fccb55991ea
# ╠═9139e5ac-7c0a-43cd-aed4-44eb4eb3b115
# ╠═91f7264e-9f4e-4ea1-95c0-77d9de3e6bbd
# ╠═30ff2319-076f-43c2-a134-9016be56e46c
# ╠═c220c177-7b8e-40cd-8d99-4aae100de7fc
# ╠═37f24271-0da7-4236-9d84-6210a1e5ac0f
# ╠═a83800b1-9df8-4536-9812-525bcad1902a
# ╠═08cdbdf3-e827-4aea-a40b-6f157fe0ba5b
# ╠═00b7e07b-b25b-4a72-9a96-59e7f18430bb
# ╠═f6f6044b-2d6c-4116-b4c9-529cdc7c7d1f
# ╠═98d2ba11-2456-4a9b-9e0a-ab20b5bd276d
# ╠═7c15972a-dde7-45f6-8d18-b62da3587cc0
# ╠═212b3eb6-3a86-4e41-a230-6c8b92147573
# ╠═cec9cd67-27db-4736-85c5-2af2a0093364
# ╠═7c3becba-7fc8-4eaf-86c6-2d92e08b04fa
# ╠═9567bffe-65b2-4073-ae85-ee339577f9be
# ╠═a892d381-a5a7-4659-bb5b-19e65183522f
# ╠═7f0a28b2-2cd9-4910-8312-dad488e40399
# ╠═bd6dd4b5-ed90-4b30-9c0f-a04f7a8adcf7
# ╠═a877c6a1-25b2-43b3-9191-6da2522b65e9
# ╠═56795636-501e-4928-a722-d515eedb4dd3
# ╠═720b30e4-bc60-48f5-99bb-d1c7fee7ae57
# ╠═381398ab-95f6-4fe2-81c6-a0478915c7c0
# ╠═7a7fd84b-36f0-4a7f-aacf-338a6a523311
# ╠═584aed06-bed9-4158-8b2d-11f8fe72b693
# ╠═aa01d3bf-16c6-4436-9eb6-55f737e4fdba
# ╠═0be498cd-14f4-45d1-8028-40a793716b0f
# ╠═d9e86f10-3244-476f-aca1-2c26fc5b592a
# ╠═bc39f4d8-0192-41c4-885a-d09f65de298a
# ╠═bc400aa6-b180-4b02-8abc-60da09acb937
# ╠═76b687c1-cfc9-4206-a9bd-3859932e9c09
# ╠═47ce252b-c3db-4f75-85c1-0fda4b34da80
# ╠═812069cd-02cf-4e48-8aa7-94753f6735a4
# ╠═508362b0-0a1a-4eec-a26f-1c45db18563c
# ╠═c07064a7-0303-4299-a080-c45f0efd5c4f
# ╠═13b382df-1a09-4478-81a3-7aa9d4ab62f7
# ╠═16480b91-338f-4f7f-b4bd-89f6e73bde55
# ╠═021b17fa-84b6-4839-96e6-a5710fae3533
# ╠═e790286e-3a1d-40a0-b323-cd3d9a81ed20
# ╠═35f75188-5502-443e-b78a-c0dd7c91f00f
# ╠═0551b30d-150d-4050-9981-c52e6dba2630
# ╠═7cd221d4-a26c-4be7-a29b-fcfbff3b8938
# ╠═2fc078ca-f662-4d67-b564-7c7ed703b629
# ╠═f56d927a-74d2-48b7-b7f3-efc9ed18e006
# ╠═48f841e8-b704-414a-bbb1-d002a6e79755
# ╠═487d2c06-70ff-4537-83f6-7bd07cfa274b
# ╠═77688f50-5b5c-4a41-89d5-1c19c0f2590e
# ╠═9715fadc-9d3b-45d4-a2f7-67e0a852be80
# ╠═37c074eb-6141-42b9-9a6b-dab79d781246
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
