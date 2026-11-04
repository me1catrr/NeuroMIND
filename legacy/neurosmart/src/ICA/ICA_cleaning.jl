#!/usr/bin/env julia
# -*- coding: utf-8 -*- #shebang        
# src/ICA/ICA_cleaning.jl

# LIMPIEZA DE DATOS MEDIANTE ICA (ICA Cleaning)
# =============================================
# Esta rutina evalúa, identifica y elimina componentes artefactuales de los datos
# descompuestos mediante ICA. Utiliza un sistema automático de evaluación basado
# en features espaciales, espectrales y estadísticas para clasificar componentes
# como artefactos o señales neuronales.
#
# PROCESO COMPLETO:
#   1. Preliminares: cálculo de métricas de referencia (datos pre-ICA)
#   2. Carga de resultados ICA (matrices S, A, W_total)
#   3. Carga de posiciones de electrodos (para mapas topográficos)
#   4. Generación automática de visualizaciones para todos los ICs
#   5. Evaluación automática mediante sistema de features y scores
#   6. Selección y eliminación de componentes artefactuales
#   7. Reconstrucción de datos limpios
#   8. Visualización comparativa antes/después
#
# EVALUACIÓN Y SELECCIÓN DE COMPONENTES ICA
# ------------------------------------------
# Preliminares:
#   - Cálculo de kurtosis y energía para datos pre-ICA (Lowpass)
#   - Métricas de referencia para comparación posterior
#
# Métricas cuantitativas por componente ICA:
#   - Kurtosis y energía de cada componente (filas de S)
#   - Identificación de componentes con características extremas
#
# Generación automática de visualizaciones para todos los ICs:
#   - Mapas topográficos (columnas de A) con proyección 2D de electrodos
#   - Time courses (filas de S) - actividad temporal de cada componente
#   - Espectro de potencia (PSD usando Welch) - características espectrales
#   - Guardado de resúmenes visuales en results/figures/ICA_cleaning/
#
# Evaluación automática de componentes mediante features y scores:
#   Features calculadas:
#     - frontal_ratio: ratio de potencia en canales frontales
#     - temporal_ratio: ratio de potencia en canales temporales
#     - blink_ratio: ratio espectral de parpadeo (0.5-4 Hz / 4-40 Hz)
#     - emg_ratio: ratio espectral de EMG (30-80 Hz / 1-30 Hz)
#     - line_ratio: ratio de ruido de línea eléctrica (48-52 Hz / total)
#     - kurtosis: medida de no-gaussianidad
#     - extreme_frac: fracción de muestras extremas (±5σ)
#     - corrEOG: correlación con canales EOG (si disponibles)
#
#   Scores calculados (combinaciones ponderadas de features):
#     - ocular_score: identifica componentes oculares (parpadeos)
#     - muscle_score: identifica componentes de EMG
#     - line_score: identifica ruido de línea eléctrica
#     - jump_score: identifica saltos y ruido bruto
#     - artifact_score: máximo de los scores anteriores (score global)
#
#   Etiquetado automático:
#     - "parpadeo": componente ocular
#     - "músculo": componente de EMG
#     - "línea": ruido de línea eléctrica
#     - "salto": saltos/ruido bruto
#     - "neuronal?": probable señal neuronal (artifact_score ≤ 1.5)
#
#   Guardado de evaluación en results/tables/ICA_components_evaluation.csv
#
# RECONSTRUCCIÓN DE DATOS LIMPIOS
# --------------------------------
# Selección automática:
#   - Componentes con artifact_score > 1.5 se marcan como artefactuales
#   - Umbral configurable (por defecto: 1.5)
#
# Eliminación de artefactos:
#   - Anulación de componentes artefactuales: S_clean[IC_art, :] = 0
#   - Solo se eliminan las componentes marcadas como artefactos
#
# Reconstrucción:
#   - X_clean = A · S_clean (sin término de media)
#   - Los datos se reconstruyen usando solo componentes "buenas"
#
# Guardado:
#   - Datos limpios en dict_EEG_ICA_clean.bin
#   - Información completa en dict_EEG_ICA_full.bin (incluye metadatos)
#
# VISUALIZACIÓN COMPARATIVA
# --------------------------
# Comparación antes vs después del ICA:
#   - Segmentos temporales de EEG por canal (superpuestos)
#   - Espectros de potencia (PSD) por canal (comparación espectral)
#   - Kurtosis y energía por canal (gráficos de barras)
#   - Permite evaluar visualmente la efectividad de la limpieza
#

# Listado de librerías
using Serialization
using CSV, DataFrames
using LinearAlgebra
using Statistics
using StatsBase  # para kurtosis
using Plots   # para las gráficas
using DSP   # asegúrate de tenerlo ya cargado
# El backend GR se inicializa automáticamente al crear el primer plot
using FFTW

# Si se ejecuta este script directamente (fuera del módulo EEG_Julia),
# cargamos utilidades de rutas para disponer de `stage_dir`.
if !@isdefined(stage_dir)
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end

# -----------------------------------------------------------------------------
# 0. Métricas cuantitativas usando datos de dict_EEG_Lowpass.bin
# -----------------------------------------------------------------------------
# Esta sección calcula métricas de referencia sobre los datos filtrados
# (antes del ICA) para poder comparar después con los datos limpios.
# Las métricas calculadas son:
#   - Kurtosis: medida de no-gaussianidad (valores altos indican picos/outliers)
#   - Energía: varianza de la señal (medida de potencia)
# Estas métricas ayudan a evaluar la efectividad de la limpieza ICA.

println("=========================================================================")
println("  Preliminares: Métricas cuantitativas (datos de dict_EEG_Lowpass.bin)")
println("=========================================================================\n")

# Cargar datos filtrados (antes del ICA) con variables locales.
# Esto evita conflictos de redefinición al copiar este bloque en Pluto.
dict_EEG_Lowpass = let
    dir_filtering_local = stage_dir(:filtering)
    path_dict_lowpass_local = joinpath(dir_filtering_local, "dict_EEG_Lowpass.bin")

    if !isfile(path_dict_lowpass_local)
        error("No se encontró $(abspath(path_dict_lowpass_local)). Ejecuta antes src/Preprocessing/filtering.jl para generar dict_EEG_Lowpass.bin (y src/Setup/IO.jl si parte desde cero).")
    end

    println("✓ Archivo: $(basename(path_dict_lowpass_local))\n")
    Serialization.deserialize(path_dict_lowpass_local)
end

# Obtener canales ordenados (orden alfabético para consistencia)
channels_lowpass = sort(collect(keys(dict_EEG_Lowpass)))

# Calcular métricas para cada canal
# Kurtosis: mide la "cola" de la distribución (valores extremos)
# Energía: varianza de la señal (potencia)
kurtosis_vals = Float64[]
energy_vals = Float64[]

for ch in channels_lowpass
    signal = dict_EEG_Lowpass[ch]
    push!(kurtosis_vals, kurtosis(signal))  # kurtosis de la señal
    push!(energy_vals, var(signal))          # varianza = energía
end

# Crear DataFrame con las métricas para visualización
df_metrics = DataFrame(
    Canal = channels_lowpass,
    Kurtosis = kurtosis_vals,
    Energía = energy_vals
)

println("Kurtosis y Energía por canal (post filtrado):")
display(df_metrics)
println()

# -----------------------------------------------------------------------------
# 1. Cargar resultados de ICA
# -----------------------------------------------------------------------------
# Se cargan los resultados del análisis ICA realizado previamente:
#   - S: componentes independientes (señales fuente)
#   - W_total: matriz de desmezcla total
#   - A: matriz de mezcla (columnas = mapas topográficos de cada IC)
#   - channels: orden de canales usado en el ICA
#
# También se calculan métricas básicas (kurtosis y energía) para cada componente,
# que servirán como referencia para la evaluación posterior.

println("====================================")
println("  ICA Visualización y Limpieza")
println("====================================\n")

# Cargar resultados del análisis ICA
dir_ica       = stage_dir(:ICA)
path_dict_ica = joinpath(dir_ica, "dict_EEG_ICA.bin")

dict_ICA = Serialization.deserialize(path_dict_ica)

# Extraer matrices y parámetros del diccionario
S        = dict_ICA["S"]          # (n_comp × n_samples) - componentes independientes
W_total  = dict_ICA["W_total"]    # (n_comp × n_channels) - matriz de desmezcla
A        = dict_ICA["A"]          # (n_channels × n_comp) - matriz de mezcla (mapas topográficos)
channels = dict_ICA["channels"]   # Vector{String} con el orden de canales
max_iter = dict_ICA["max_iter"]   # parámetros del ICA (para referencia)
tol      = dict_ICA["tol"]         # parámetros del ICA (para referencia)

println("✔ Carga de datos de ICA:")
println("  S (n_comp × n_samples):        ", size(S))
println("  W_total (n_comp × n_channels):  ", size(W_total))
println("  A (n_channels × n_comp):        ", size(A))
println("  channels: ", channels, "\n")

# Calcular límite de color para mapas topográficos
# Se usa el máximo absoluto de todos los pesos para normalizar la escala de colores
global topo_clim = maximum(abs, vec(A))

# Métricas cuantitativas para los datos de ICA
# Se calculan kurtosis y energía para cada componente (cada fila de S)
# Estas métricas ayudan a identificar componentes con características extremas
kurtosis_vals = kurtosis.(eachrow(S))  # kurtosis por componente
energy_vals   = var.(eachrow(S))        # energía (varianza) por componente

println("Kurtosis y Energía por componente ICA:")
display(DataFrame(Componente = 1:size(S, 1), Kurtosis = kurtosis_vals, Energía = energy_vals))
println()

# -----------------------------------------------------------------------------
# 2. Cargar archivo TSV con posiciones de electrodos
# -----------------------------------------------------------------------------
# Se cargan las coordenadas 3D (x, y, z) de los electrodos desde un archivo TSV.
# Estas coordenadas son necesarias para generar los mapas topográficos de las
# componentes ICA. Los mapas topográficos muestran cómo se distribuyen espacialmente
# los pesos de cada componente sobre el cuero cabelludo.

dir_electrodes = electrodes_dir()
path_elec      = joinpath(dir_electrodes, "sub-M05_ses-T2_electrodes.tsv")

# Leer archivo TSV con información de electrodos
df_elec = CSV.read(path_elec, DataFrame; delim = '\t')

println("📍 Electrodos cargados:")
display(first(df_elec, 5))
println()

# Filtramos solo electrodos de tipo EEG
# El archivo puede contener otros tipos de electrodos (EOG, EMG, etc.)
df_eeg = filter(row -> row.type == "EEG", df_elec)

println("Electrodos EEG detectados:")
display(df_eeg[:, [:name, :x, :y, :z]])

# Crear diccionario nombre → coordenadas (x, y, z)
# Facilita el acceso rápido a las coordenadas por nombre de canal
elec_pos = Dict{String, Tuple{Float64,Float64,Float64}}()

for row in eachrow(df_eeg)
    elec_pos[String(row.name)] = (Float64(row.x), Float64(row.y), Float64(row.z))
end

# Verificar que todos los canales del ICA tengan coordenadas
# Si falta alguna posición, no se podrá dibujar el mapa topográfico completo
missing_ch = String[]
for ch in channels
    if !haskey(elec_pos, ch)
        push!(missing_ch, ch)
    end
end

if !isempty(missing_ch)
    println("⚠ Ojo: faltan posiciones para estos canales: ", missing_ch)
else
    println("✅ Todos los canales tienen coordenadas.\n")
end

"""
    project_elec_xy(elec_pos::Dict{String,Tuple{Float64,Float64,Float64}}, channels)

Proyecta coordenadas 3D (x,y,z) sobre un disco 2D (top-down)
usando una proyección azimutal tipo EEGLAB.

Esta función convierte las coordenadas 3D de los electrodos en coordenadas 2D
para poder dibujar mapas topográficos en un plano. La proyección es similar
a la usada en EEGLAB: los electrodos en la parte superior de la cabeza (z alto)
se proyectan cerca del centro, y los del ecuador (z=0) se proyectan cerca del borde.

PROCESO:
  1. Normaliza las coordenadas 3D a una esfera unitaria
  2. Calcula el ángulo azimutal (θ) en el plano horizontal
  3. Convierte la altura (z) en un radio (r) en el plano 2D
  4. Proyecta usando coordenadas polares: (x, y) = (r*cos(θ), r*sin(θ))

Devuelve:
    xs2d::Vector{Float64}, ys2d::Vector{Float64}, mask_haspos::BitVector
"""
function project_elec_xy(elec_pos::Dict{String,Tuple{Float64,Float64,Float64}}, channels)
    xs2d = Float64[]
    ys2d = Float64[]
    haspos = BitVector(undef, length(channels))

    # Normalizamos a radio 1 y proyectamos
    # Cada electrodo se proyecta desde su posición 3D a coordenadas 2D
    for (i, ch) in enumerate(channels)
        if haskey(elec_pos, ch)
            (x, y, z) = elec_pos[ch]
            # Normalizar a esfera unitaria
            r3 = sqrt(x^2 + y^2 + z^2) + eps()
            x /= r3; y /= r3; z /= r3

            # Ángulo azimutal en el plano horizontal (0 a 2π)
            θ = atan(y, x)
            
            # Convertir altura (z) a radio en el plano 2D
            # z=1 (vértice) -> r=0, z=0 (ecuador) -> r≈1
            r = (π/2 - asin(z)) / (π/2)         # en [0,1]

            # Proyección a coordenadas cartesianas 2D
            push!(xs2d, r * cos(θ))
            push!(ys2d, r * sin(θ))
            haspos[i] = true
        else
            # Si no hay posición, usar NaN (no se dibujará)
            push!(xs2d, NaN)
            push!(ys2d, NaN)
            haspos[i] = false
        end
    end

    return xs2d, ys2d, haspos
end

xs2d_all, ys2d_all, haspos = project_elec_xy(elec_pos, channels)

# -----------------------------------------------------------------------------
# 3. Función para dibujar topografía de un componente ICA (A)   
# -----------------------------------------------------------------------------
"""
    make_topomap_grid(xs, ys, vals; radius=1.0, n_grid=100, sigma=0.12)

Genera una rejilla (X, Y, Z) para topomap interpolando `vals` en (xs,ys)
mediante un kernel gaussiano con ancho `sigma`.

Esta función crea una rejilla interpolada para visualizar mapas topográficos.
Interpola los valores de los electrodos (vals) en una rejilla regular usando
un kernel gaussiano. Esto permite crear mapas de contorno suaves que muestran
la distribución espacial de los pesos de cada componente ICA.

PROCESO:
  1. Crea una rejilla regular de n_grid × n_grid puntos
  2. Para cada punto de la rejilla dentro del círculo:
     - Calcula distancias a todos los electrodos
     - Aplica pesos gaussianos según la distancia
     - Interpola el valor como promedio ponderado
  3. Puntos fuera del círculo quedan como NaN (no se dibujan)

Parámetros:
  - radius: radio del círculo de la cabeza (típicamente 1.0)
  - n_grid: número de puntos en cada dimensión (mayor = más suave pero más lento)
  - sigma: ancho del kernel gaussiano (mayor = interpolación más suave)

Sólo rellena puntos dentro del círculo de radio `radius`.
"""
function make_topomap_grid(xs, ys, vals; radius=1.0, n_grid=100, sigma=0.12)
    # Crear rejilla regular de coordenadas
    xi = range(-radius, radius; length=n_grid)
    yi = range(-radius, radius; length=n_grid)

    # Matriz de valores interpolados (inicializada con NaN)
    Z = fill(NaN, n_grid, n_grid)

    # Para cada punto de la rejilla
    for (ix, xg) in enumerate(xi), (iy, yg) in enumerate(yi)
        # Verificar que el punto esté dentro del círculo
        r = hypot(xg, yg)
        r > radius && continue

        # Calcular pesos gaussianos basados en distancia a electrodos
        # Los electrodos más cercanos tienen mayor peso
        d2 = (xs .- xg).^2 .+ (ys .- yg).^2  # distancias al cuadrado
        w  = @. exp(-d2 / (2*sigma^2))        # pesos gaussianos

        sw = sum(w)
        sw == 0 && continue  # Si no hay electrodos cercanos, saltar

        # Interpolación: promedio ponderado de valores de electrodos
        Z[iy, ix] = sum(w .* vals) / sw
    end

    return xi, yi, Z
end

# 4.1. Función para dibujar topografía de un componente ICA
"""
    plot_ic_topography(ic; clim=nothing, n_grid=100)

Dibuja el mapa topográfico de un componente ICA (tipo EEGLAB).

El mapa topográfico muestra cómo se distribuyen espacialmente los pesos de una
componente ICA sobre el cuero cabelludo. Esto es crucial para identificar el
tipo de artefacto: componentes oculares tienen alta potencia frontal, componentes
de EMG tienen alta potencia temporal, etc.

ELEMENTOS DEL MAPA:
  - Mapa interpolado con isolíneas: distribución suave de los pesos espaciales
  - Círculo de cabeza: contorno del cuero cabelludo
  - Nariz y orejas: referencias anatómicas
  - Electrodos: puntos negros con etiquetas de nombres de canales
  - Escala de colores: rojo = valores positivos, azul = valores negativos

Parámetros:
  - ic: índice del componente ICA a visualizar (1 a n_comp)
  - clim: límites de color (opcional). Si es nothing, se usa auto-escala
  - n_grid: resolución de la rejilla de interpolación (mayor = más suave)

Devuelve un objeto Plot que puede ser mostrado o guardado.
"""
function plot_ic_topography(ic::Int; clim = nothing, n_grid::Int = 100)
    n_channels, n_comp = size(A)
    @assert 1 ≤ ic ≤ n_comp "ic debe estar entre 1 y $n_comp"

    topo = A[:, ic]              # pesos espaciales del IC

    # Usamos únicamente electrodos con posición conocida
    xs = xs2d_all[haspos]
    ys = ys2d_all[haspos]
    vals = topo[haspos]

    # Generar rejilla interpolada
    xi, yi, Z = make_topomap_grid(xs, ys, vals; n_grid = n_grid)

    # --- Mapa interpolado con isolíneas ---
    plt = contourf(
        xi, yi, Z;
        levels = 15,
        colorbar = true,
        aspect_ratio = 1,
        xlims = (-1, 1),
        ylims = (-1, 1),
        title = "Mapa componente $ic",
        legend = false,
        color = :jet,
        clim = clim !== nothing ? clim : :auto,
    )

    # --- Cabeza: círculo, nariz y orejas ---
    θ = range(0, 2π; length=200)
    plot!(plt, cos.(θ), sin.(θ); lw=2, color=:black)   # contorno cabeza

    # nariz
    nose_x = [0.0, 0.08, 0.0, -0.08]
    nose_y = [1.0, 1.08, 1.16, 1.08]
    plot!(plt, nose_x, nose_y; lw=2, color=:black)

    # orejas
    ear_xL = [-1.02, -1.1, -1.12, -1.1, -1.02]
    ear_yL = [ 0.2,  0.18,  0.0, -0.18, -0.2]
    ear_xR =  .*(-1, ear_xL)
    ear_yR = ear_yL
    plot!(plt, ear_xL, ear_yL; lw=2, color=:black)
    plot!(plt, ear_xR, ear_yR; lw=2, color=:black)

    # --- Electrodos: puntos + etiquetas ---
    scatter!(plt, xs, ys; ms=4, mc=:black, msw=0.5, ma=0.8)

    for (i, ch) in enumerate(channels)
        haspos[i] || continue
        x = xs2d_all[i]
        y = ys2d_all[i]
        annotate!(plt, x, y, text(ch, 7, :black))
    end

    # quitamos ejes numéricos
    plot!(plt, xticks = nothing, yticks = nothing, framestyle = :box)

    return plt
end

# -----------------------------------------------------------------------------
# 4. Función para mostrar la serie temporal de un IC
# -----------------------------------------------------------------------------
"""
    plot_ic_timecourse(ic; fs=nothing)

Visualiza la serie temporal (time course) de un componente ICA.

El time course muestra la actividad temporal de la componente a lo largo de
todo el registro. Es útil para identificar:
  - Parpadeos: picos grandes y abruptos
  - EMG: actividad de alta frecuencia
  - Saltos: valores extremos aislados
  - Actividad neuronal: oscilaciones más regulares

Parámetros:
  - ic: índice del componente ICA a visualizar
  - fs: frecuencia de muestreo (Hz). Si es nothing, se usa escala de muestras

Devuelve un objeto Plot con la serie temporal.
"""
function plot_ic_timecourse(ic::Int; fs::Union{Nothing,Float64} = nothing)
    n_comp, n_samples = size(S)
    @assert 1 ≤ ic ≤ n_comp "ic debe estar entre 1 y $n_comp"

    sig = S[ic, :]

    if fs === nothing
        t = 1:n_samples
        xl = "Muestra"
    else
        t = (0:(n_samples-1)) ./ fs
        xl = "Tiempo (s)"
    end

    plt = plot(
        t, sig,
        xlabel = xl,
        ylabel = "Amplitud",
        title = "Actividad componente $ic",
        legend = false,
    )

    return plt
end

# 5. PSD de un componente ICA
"""
    calculate_PSD_signal(signal::AbstractVector, fs::Real)

Calcula la PSD de una señal (p.ej. un componente ICA) usando Welch.
Devuelve (freqs, power).
"""
function calculate_PSD_signal(signal::AbstractVector, fs::Real)
    n_muestras = length(signal)

    # quitar DC
    sig = signal .- mean(signal)

    p = welch_pgram(sig; fs = fs, window = hamming, nfft = n_muestras)
    freqs = DSP.freq(p)
    power = DSP.power(p)   # unidades ~ µV²/Hz si la señal está en µV

    return freqs, power
end

"""
    plot_ic_spectrum(ic; fs, fmax=50.0)

Calcula y grafica el espectro de potencia (PSD) de un componente ICA.

El espectro de potencia muestra la distribución de energía en frecuencias.
Es crucial para identificar tipos de artefactos:
  - Parpadeos: alta potencia en 0.5-4 Hz (banda delta/baja)
  - EMG: alta potencia en 30-80 Hz (banda alta)
  - Línea eléctrica: pico en 50 Hz (o 60 Hz según región)
  - Actividad neuronal: potencia distribuida en bandas típicas (alpha, beta, etc.)

El método usado es Welch (promedio de periodogramas con ventanas de Hamming),
que proporciona estimaciones más suaves y estables que el periodograma simple.

Parámetros:
  - ic: índice del componente ICA a visualizar
  - fs: frecuencia de muestreo (Hz) - requerido
  - fmax: frecuencia máxima a mostrar (Hz). Por defecto 50.0

Devuelve un objeto Plot con el espectro en escala logarítmica.
"""
function plot_ic_spectrum(ic::Int; fs::Real, fmax::Real = 50.0)
    n_comp, n_samples = size(S)
    @assert 1 ≤ ic ≤ n_comp "ic debe estar entre 1 y $n_comp"

    signal = S[ic, :]

    freqs, power = calculate_PSD_signal(signal, fs)

    # Recortar a fmax
    idx = freqs .<= fmax
    freqs = freqs[idx]
    power = power[idx]

    # Escala logarítmica en potencia, igual que en tu función de canales
    power_min, power_max = extrema(power)
    y_min_log = floor(log10(power_min))
    y_max_log = ceil(log10(power_max))
    ylim_log = (10.0^y_min_log, 10.0^y_max_log)

    p_spec = plot(
        xlabel = "Frecuencia (Hz)",
        ylabel = "Potencia (µV²/Hz)",
        title  = "PSD componente $ic",
        legend = false,
        xlim   = (0, fmax),
        yscale = :log10,
        ylim   = ylim_log,
    )

    plot!(p_spec, freqs, power; lw = 2)  # opcional: label = "IC $ic"

    return p_spec
end


# Figura tipo EEGLAB
"""
    plot_ic_summary(ic; fs)

Crea una figura resumen tipo EEGLAB con tres paneles para un componente ICA.

Esta función genera una visualización completa que combina las tres vistas
principales de una componente ICA, similar a la interfaz de EEGLAB. Es la
visualización estándar para inspección manual de componentes.

PANELES:
  1. Mapa topográfico (arriba izquierda):
     - Muestra la distribución espacial de los pesos
     - Útil para identificar localización del artefacto
  2. Serie temporal (arriba derecha):
     - Muestra la actividad temporal de la componente
     - Útil para identificar patrones temporales
  3. Espectro de potencia (abajo, ancho completo):
     - Muestra la distribución espectral
     - Útil para identificar características de frecuencia

Esta visualización combinada permite una evaluación completa de cada componente
para decidir si es artefacto o señal neuronal.

Parámetros:
  - ic: índice del componente ICA a visualizar
  - fs: frecuencia de muestreo (Hz) - requerido

Devuelve un objeto Plot con el layout combinado (2 filas, 2 columnas).
"""
function plot_ic_summary(ic::Int; fs::Float64)
    p_topo     = plot_ic_topography(ic; clim = (-topo_clim, topo_clim))
    p_activity = plot_ic_timecourse(ic; fs = fs)
    p_spec     = plot_ic_spectrum(ic; fs = fs)

    l = @layout [a b; c c]

    fig = plot(p_topo, p_activity, p_spec,
               layout = l, size = (900, 700))

    display(fig)
    return fig
end

# Generar y guardar plot_ic_summary para todas las componentes ICA
# -----------------------------------------------------------------------------
# Esta sección genera automáticamente visualizaciones completas (tipo EEGLAB)
# para todas las componentes ICA. Cada visualización incluye:
#   - Mapa topográfico (distribución espacial)
#   - Serie temporal (actividad temporal)
#   - Espectro de potencia (distribución espectral)
#
# Estas visualizaciones se guardan como archivos PNG y permiten inspección
# manual de todas las componentes para verificar la evaluación automática.
# -----------------------------------------------------------------------------
println("=" ^ 60)
println("  Generando resúmenes visuales de componentes ICA")
println("=" ^ 60 * "\n")

# Crear directorio de salida si no existe
# Las figuras se guardan en results/figures/ICA_cleaning/
dir_output = stage_dir(:ICA_cleaning; kind = :figures)
isdir(dir_output) || mkpath(dir_output)

# Obtener número de componentes ICA
n_comp = size(S, 1)
fs = 500.0  # frecuencia de muestreo (Hz)

println("Generando resúmenes para $n_comp componentes ICA...")
println("  → Directorio de salida: $dir_output")
println()

# Generar visualización para cada componente
for ic in 1:n_comp
    println("  Generando resumen para componente $ic/$n_comp...")
    
    # Generar el resumen completo (topografía + time course + espectro)
    fig = plot_ic_summary(ic; fs = fs)
    
    # Guardar la figura como PNG
    # Formato de nombre: IC_001.png, IC_002.png, etc. (con padding de 3 dígitos)
    filename = joinpath(dir_output, "IC_$(lpad(ic, 3, '0')).png")
    savefig(fig, filename)
    
    println("    ✓ Guardado en: $filename")
end

println("\n✅ Todos los resúmenes de componentes ICA han sido generados y guardados.")
println("   Directorio: $dir_output\n")

# Scores para las componentes 
# ----------------------------
# Score_1. Ratio espacial-temporal
# Score_2. Potencias de banda PSD
# Score_3. Blink Ratio
# Score_4. EMG Ratio
# Score_5. Line Ratio
# Score_6. Kurtosis
# Score_7. Fracción por muestras
# Score_8. Correlación con EOG

#############################
# Funciones auxiliares básicas
#############################

"""
    zscore(v::AbstractVector)

Normalización z-score simple de un vector.
Convierte los valores a desviaciones estándar desde la media.
Si la desviación estándar es 0, devuelve un vector de ceros.

Esta función se usa para normalizar features antes de calcular scores,
permitiendo comparar features con diferentes escalas.
"""
function zscore(v::AbstractVector)
    m = mean(v)
    s = std(v)
    if s == 0
        return zeros(length(v))
    else
        return (v .- m) ./ s
    end
end

"""
    kurtosis_simple(x::AbstractVector)

Calcula kurtosis poblacional sencilla (sin restar 3).
La kurtosis mide la "cola" de la distribución: valores altos indican
distribuciones con picos pronunciados o valores extremos.

Se usa para identificar componentes con características extremas
(como artefactos de parpadeo o saltos).
"""
function kurtosis_simple(x::AbstractVector)
    m = mean(x)
    xc = x .- m  # centrar
    s2 = mean(xc.^2)  # varianza
    if s2 == 0
        return 0.0
    end
    m4 = mean(xc.^4)  # momento de cuarto orden
    return m4 / (s2^2)  # kurtosis (sin restar 3, no necesario para comparaciones)
end

"""
    bandpower(x::AbstractVector, fs::Real, f_low::Real, f_high::Real)

Calcula la potencia en una banda de frecuencias usando FFT (periodograma sencillo).

Esta función calcula cuánta energía tiene una señal en un rango de frecuencias
específico. Se usa para identificar componentes con características espectrales
típicas de artefactos:
  - Parpadeos: alta potencia en 0.5-4 Hz
  - EMG: alta potencia en 30-80 Hz
  - Línea eléctrica: alta potencia en 48-52 Hz

Parámetros:
  - x: señal temporal
  - fs: frecuencia de muestreo
  - f_low, f_high: límites de la banda de frecuencias (Hz)

Devuelve la potencia total en la banda especificada.
"""
function bandpower(x::AbstractVector, fs::Real, f_low::Real, f_high::Real)
    x = x .- mean(x)  # eliminar DC
    N = length(x)

    # FFT real -> solo frecuencias positivas (más eficiente)
    X = rfft(x)
    freqs = (0:length(X)-1) .* (fs / N)  # frecuencias correspondientes

    # PSD simple (periodograma)
    # No nos preocupa la constante exacta, solo comparaciones relativas
    psd = abs.(X).^2 ./ (fs * N)

    # Seleccionar frecuencias en la banda deseada
    idx = (freqs .>= f_low) .& (freqs .<= f_high)
    return sum(psd[idx])  # potencia total en la banda
end

#############################
# Cálculo de features por IC
#############################

"""
    compute_features(A, S, fs; frontal_idxs, temporal_idxs, eog = nothing)

Calcula ratios espaciales y espectrales + estadísticos para todos los ICs.

Esta función extrae características (features) de cada componente ICA que ayudan
a identificar si es un artefacto o señal neuronal. Las features calculadas son:

FEATURES ESPACIALES (basadas en la matriz A):
  - frontal_ratio: ratio de potencia en canales frontales vs no frontales
    (alto = componente frontal, típico de parpadeos)
  - temporal_ratio: ratio de potencia en canales temporales vs no temporales
    (alto = componente temporal, típico de EMG)

FEATURES ESPECTRALES (basadas en la señal temporal S):
  - blink_ratio: ratio potencia 0.5-4 Hz / 4-40 Hz (alto = parpadeo)
  - emg_ratio: ratio potencia 30-80 Hz / 1-30 Hz (alto = actividad muscular)
  - line_ratio: ratio potencia 48-52 Hz / total (alto = ruido de línea eléctrica)

FEATURES ESTADÍSTICAS:
  - kurtosis: medida de no-gaussianidad (alto = picos/valores extremos)
  - extreme_frac: fracción de muestras fuera de ±5 desviaciones estándar
    (alto = saltos/ruido bruto)

FEATURES DE CORRELACIÓN:
  - corrEOG: máxima correlación con canales EOG (si disponibles)
    (alto = componente relacionado con movimiento ocular)

Parámetros:
  - A: matriz canales × ICs (pesos espaciales / mapas topográficos)
  - S: matriz ICs × muestras (activaciones temporales)
  - fs: frecuencia de muestreo
  - frontal_idxs: índices de canales frontales
  - temporal_idxs: índices de canales temporales
  - eog: opcional, matriz canales_eog × muestras

Devuelve un NamedTuple con todas las features calculadas.
"""
function compute_features(
    A::AbstractMatrix,    # canales x ICs
    S::AbstractMatrix,    # ICs x muestras
    fs::Real;
    frontal_idxs::Vector{Int},
    temporal_idxs::Vector{Int},
    eog::Union{Nothing, AbstractMatrix} = nothing  # eog: canales_eog x muestras
)

    n_chan, n_ic = size(A)
    _, n_samp    = size(S)

    # Inicializar vectores para almacenar features de cada IC
    frontal_ratio  = zeros(n_ic)
    temporal_ratio = zeros(n_ic)
    blink_ratio    = zeros(n_ic)
    emg_ratio      = zeros(n_ic)
    line_ratio     = zeros(n_ic)
    kurtosis_vec   = zeros(n_ic)
    extreme_frac   = zeros(n_ic)
    corrEOG        = zeros(n_ic)

    # Calcular índices complementarios (canales no frontales, no temporales)
    # Estos se usan como denominadores en los ratios espaciales
    all_idxs = collect(1:n_chan)
    nonfrontal_idxs  = setdiff(all_idxs, frontal_idxs)
    nontemporal_idxs = setdiff(all_idxs, temporal_idxs)

    # Calcular features para cada componente ICA
    for k in 1:n_ic
        map_k = A[:, k]  # mapa topográfico del IC k (columna k de A)
        s_k   = S[k, :]  # señal temporal del IC k (fila k de S)

        # --- RATIOS ESPACIALES ---
        # Comparan la potencia del componente en regiones específicas vs el resto
        # Componentes de parpadeo suelen tener alta potencia frontal
        # Componentes de EMG suelen tener alta potencia temporal
        num_f  = mean(abs.(map_k[frontal_idxs]))      # potencia promedio en frontales
        den_f  = mean(abs.(map_k[nonfrontal_idxs])) + 1e-12  # potencia promedio en no frontales
        frontal_ratio[k] = num_f / den_f

        num_t  = mean(abs.(map_k[temporal_idxs]))     # potencia promedio en temporales
        den_t  = mean(abs.(map_k[nontemporal_idxs])) + 1e-12  # potencia promedio en no temporales
        temporal_ratio[k] = num_t / den_t

        # --- POTENCIAS DE BANDA ESPECTRALES ---
        # Calculan potencia en diferentes bandas de frecuencia
        # Los artefactos tienen características espectrales distintivas
        P_total = bandpower(s_k, fs, 0.5, 100)  # potencia total (0.5-100 Hz)
        P_0_4   = bandpower(s_k, fs, 0.5, 4)    # banda delta/baja (parpadeos)
        P_4_40  = bandpower(s_k, fs, 4, 40)      # banda típica neuronal
        P_1_30  = bandpower(s_k, fs, 1, 30)      # banda amplia (referencia)
        P_30_80 = bandpower(s_k, fs, 30, 80)    # banda alta (EMG)
        P_48_52 = bandpower(s_k, fs, 48, 52)     # banda de línea eléctrica (50 Hz)

        # Ratios espectrales: comparan potencia en bandas características
        blink_ratio[k] = P_0_4   / (P_4_40  + 1e-12)  # alto = parpadeo
        emg_ratio[k]   = P_30_80 / (P_1_30  + 1e-12)  # alto = EMG
        line_ratio[k]  = P_48_52 / (P_total + 1e-12)  # alto = ruido de línea

        # --- KURTOSIS Y MUESTRAS EXTREMAS ---
        # Identifican componentes con distribuciones no gaussianas o valores extremos
        # Típico de artefactos como saltos, parpadeos grandes, etc.
        kurtosis_vec[k] = kurtosis_simple(s_k)  # kurtosis de la señal

        # Fracción de muestras que están fuera de ±5 desviaciones estándar
        # Valores altos indican saltos o ruido bruto
        μ  = mean(s_k)
        σ  = std(s_k)
        thr = μ + 5σ   # umbral superior
        thrn = μ - 5σ  # umbral inferior
        if σ == 0
            extreme_frac[k] = 0.0
        else
            n_ext = count(t -> (t > thr) || (t < thrn), s_k)  # contar extremos
            extreme_frac[k] = n_ext / n_samp  # fracción de muestras extremas
        end

        # --- CORRELACIÓN CON EOG (SI HAY) ---
        # Si hay canales EOG disponibles, calcula la correlación máxima
        # Componentes de parpadeo/movimiento ocular tienen alta correlación con EOG
        if eog !== nothing
            # eog: canales_eog × muestras
            # calcula máxima correlación absoluta con cualquier canal EOG
            maxcorr = 0.0
            for e in axes(eog, 1)
                e_sig = eog[e, :]
                if std(e_sig) > 0 && std(s_k) > 0
                    c = cor(s_k, e_sig)  # correlación de Pearson
                    maxcorr = max(maxcorr, abs(c))  # mantener la máxima (absoluta)
                end
            end
            corrEOG[k] = maxcorr
        else
            corrEOG[k] = 0.0  # sin EOG, no hay correlación
        end
    end

    return (
        frontal_ratio  = frontal_ratio,
        temporal_ratio = temporal_ratio,
        blink_ratio    = blink_ratio,
        emg_ratio      = emg_ratio,
        line_ratio     = line_ratio,
        kurtosis       = kurtosis_vec,
        extreme_frac   = extreme_frac,
        corrEOG        = corrEOG,
    )
end

#############################
# Funciones de score
#############################
# Estas funciones combinan features normalizadas (z-scores) para crear scores
# que indican la probabilidad de que un componente sea un tipo específico
# de artefacto. Los scores se calculan como combinaciones lineales ponderadas
# de las features relevantes.

"""
    ocular_score(frontal_z, blink_z, corrEOG_z; w1=0.4, w2=0.4, w3=0.2)

Score para identificar componentes oculares (parpadeos, movimientos oculares).

Combina tres features normalizadas:
  - frontal_z: ratio espacial frontal (peso 0.4)
  - blink_z: ratio espectral de parpadeo (peso 0.4)
  - corrEOG_z: correlación con EOG (peso 0.2)

Un score alto indica que el componente probablemente es un artefacto ocular.
"""
function ocular_score(frontal_z, blink_z, corrEOG_z;
                      w1 = 0.4, w2 = 0.4, w3 = 0.2)
    return w1 .* frontal_z .+ w2 .* blink_z .+ w3 .* corrEOG_z
end

"""
    muscle_score(emg_z, temporal_z, kurtosis_z; v1=0.5, v2=0.3, v3=0.2)

Score para identificar componentes de actividad muscular (EMG).

Combina tres features normalizadas:
  - emg_z: ratio espectral EMG (peso 0.5) - más importante
  - temporal_z: ratio espacial temporal (peso 0.3)
  - kurtosis_z: kurtosis (peso 0.2)

Un score alto indica que el componente probablemente es actividad muscular.
"""
function muscle_score(emg_z, temporal_z, kurtosis_z;
                      v1 = 0.5, v2 = 0.3, v3 = 0.2)
    return v1 .* emg_z .+ v2 .* temporal_z .+ v3 .* kurtosis_z
end

"""
    line_score(line_z)

Score para identificar componentes de ruido de línea eléctrica (50/60 Hz).

Usa directamente el z-score del line_ratio (peso 1.0).
Un score alto indica que el componente contiene principalmente ruido de línea.
"""
function line_score(line_z)
    return line_z      # aquí el peso es 1 directamente
end

"""
    jump_score(kurtosis_z, extreme_z; u1=0.5, u2=0.5)

Score para identificar componentes con saltos o ruido bruto.

Combina dos features normalizadas con pesos iguales:
  - kurtosis_z: kurtosis (peso 0.5)
  - extreme_z: fracción de muestras extremas (peso 0.5)

Un score alto indica que el componente tiene valores extremos o saltos.
"""
function jump_score(kurtosis_z, extreme_z;
                    u1 = 0.5, u2 = 0.5)
    return u1 .* kurtosis_z .+ u2 .* extreme_z
end

#############################
# Función global: evalúa todos los ICs
#############################

"""
    evaluate_ics(A, S, fs; frontal_idxs, temporal_idxs, eog = nothing, artifact_thresh = 1.5)

Evalúa todos los componentes ICA y los etiqueta como artefactos o señales neuronales.

Esta función es el núcleo del sistema de evaluación automática. Realiza los siguientes pasos:

1. Calcula features crudos para cada componente (ratios espaciales, espectrales, estadísticos)
2. Normaliza las features usando z-score (para comparar en la misma escala)
3. Calcula scores por tipo de artefacto (ocular, muscular, línea, saltos)
4. Determina el score global de artefacto (máximo de los scores individuales)
5. Etiqueta cada componente según el tipo de artefacto más probable
6. Calcula un neural_score (negativo del artifact_score) para identificar componentes neuronales

PROCESO DE ETIQUETADO:
  - Si artifact_score > artifact_thresh: etiqueta como tipo de artefacto correspondiente
  - Si artifact_score ≤ artifact_thresh: etiqueta como "neuronal?" (probable señal neuronal)

Parámetros:
  - A: matriz canales × ICs (mezcla / pesos espaciales)
  - S: matriz ICs × muestras (activaciones temporales)
  - fs: frecuencia de muestreo
  - frontal_idxs, temporal_idxs: índices de canales frontales / temporales
  - eog: opcional, matriz canales_eog × muestras
  - artifact_thresh: umbral de score global para etiquetar como artefacto (default: 1.5)

Devuelve un DataFrame con features, scores y etiqueta final por IC.
"""
function evaluate_ics(
    A::AbstractMatrix,
    S::AbstractMatrix,
    fs::Real;
    frontal_idxs::Vector{Int},
    temporal_idxs::Vector{Int},
    eog::Union{Nothing, AbstractMatrix} = nothing,
    artifact_thresh::Real = 1.5
)

    n_chan, n_ic = size(A)

    # 1) FEATURES crudos
    # Calcula todas las características (ratios, kurtosis, etc.) para cada IC
    feats = compute_features(A, S, fs;
                             frontal_idxs = frontal_idxs,
                             temporal_idxs = temporal_idxs,
                             eog = eog)

    # 2) Normalización (z-score) de cada feature
    # Esto permite comparar features con diferentes escalas y unidades
    # Un z-score de 1.5 significa que el valor está 1.5 desviaciones estándar por encima de la media
    frontal_z   = zscore(feats.frontal_ratio)
    temporal_z  = zscore(feats.temporal_ratio)
    blink_z     = zscore(feats.blink_ratio)
    emg_z       = zscore(feats.emg_ratio)
    line_z      = zscore(feats.line_ratio)
    kurtosis_z  = zscore(feats.kurtosis)
    extreme_z   = zscore(feats.extreme_frac)
    corrEOG_z   = zscore(feats.corrEOG)

    # 3) Scores por tipo de artefacto
    # Cada score combina features relevantes para identificar un tipo específico
    ocular   = ocular_score(frontal_z, blink_z, corrEOG_z)  # score de parpadeo
    muscle   = muscle_score(emg_z, temporal_z, kurtosis_z)  # score de EMG
    line     = line_score(line_z)                            # score de línea eléctrica
    jump     = jump_score(kurtosis_z, extreme_z)             # score de saltos

    # 4) Score global de artefacto y etiqueta
    # El score global es el máximo de los scores individuales
    # Esto identifica el tipo de artefacto más probable
    artifact_global = similar(ocular)
    label           = Vector{String}(undef, n_ic)

    for k in 1:n_ic
        scores_k = (ocular[k], muscle[k], line[k], jump[k])
        maxscore = maximum(scores_k)  # máximo de los 4 scores
        artifact_global[k] = maxscore

        # Etiquetar según qué score fue el máximo
        if maxscore > artifact_thresh
            # Decide etiqueta según qué score fue mayor
            if maxscore == ocular[k]
                label[k] = "parpadeo"
            elseif maxscore == muscle[k]
                label[k] = "músculo"
            elseif maxscore == line[k]
                label[k] = "línea"
            else
                label[k] = "salto"
            end
        else
            # Score bajo = probablemente señal neuronal
            label[k] = "neuronal?"
        end
    end

    # 5) neural_score (simple: negativo del global)
    # Un score alto de artefacto → score bajo (negativo) de neural
    # Un score bajo de artefacto → score alto (menos negativo) de neural
    neural_score = -artifact_global

    # 6) Construimos DataFrame resumen
    # Incluye todas las features, scores y la etiqueta final
    df = DataFrame(
        IC             = collect(1:n_ic),
        frontal_ratio  = feats.frontal_ratio,
        temporal_ratio = feats.temporal_ratio,
        blink_ratio    = feats.blink_ratio,
        emg_ratio      = feats.emg_ratio,
        line_ratio     = feats.line_ratio,
        kurtosis       = feats.kurtosis,
        extreme_frac   = feats.extreme_frac,
        corrEOG        = feats.corrEOG,
        ocular_score   = ocular,
        muscle_score   = muscle,
        line_score     = line,
        jump_score     = jump,
        artifact_score = artifact_global,
        neural_score   = neural_score,
        label          = label
    )

    return df
end

# -----------------------------------------------------------------------------
# Definir índices de canales frontales y temporales
# -----------------------------------------------------------------------------
# Esta sección identifica automáticamente los canales frontales y temporales
# basándose en la nomenclatura estándar 10-20. Estos canales se usan para
# calcular los ratios espaciales que ayudan a identificar tipos de artefactos:
#   - Canales frontales: típicamente afectados por parpadeos y movimientos oculares
#   - Canales temporales: típicamente afectados por actividad muscular (EMG)
#
# La identificación se hace por nombre de canal usando prefijos estándar.

# Identificar canales frontales (F, FC, FT, FP)
# Estos canales están en la parte frontal del cuero cabelludo
frontal_idxs = Int[]
for (i, ch) in enumerate(channels)
    ch_upper = uppercase(ch)  # convertir a mayúsculas para comparación
    # Buscar canales que empiecen con F, FC, FT o FP
    if startswith(ch_upper, "F") || startswith(ch_upper, "FC") || 
       startswith(ch_upper, "FT") || startswith(ch_upper, "FP")
        push!(frontal_idxs, i)
    end
end

# Identificar canales temporales (T, TP)
# Estos canales están en las regiones temporales (lados de la cabeza)
temporal_idxs = Int[]
for (i, ch) in enumerate(channels)
    ch_upper = uppercase(ch)
    # Buscar canales que empiecen con T o TP
    if startswith(ch_upper, "T") || startswith(ch_upper, "TP")
        push!(temporal_idxs, i)
    end
end

println("📍 Canales frontales detectados: ", [channels[i] for i in frontal_idxs])
println("📍 Canales temporales detectados: ", [channels[i] for i in temporal_idxs])
println()

# Frecuencia de muestreo (necesaria para cálculos espectrales)
fs = 500.0

# Canales EOG (no disponibles en estos datos)
# Si hubiera canales EOG, se usarían para calcular correlaciones
# que ayudarían a identificar componentes oculares
eog = nothing

# Evaluar todos los ICs usando el sistema automático
# Esta función calcula features, scores y etiquetas para todas las componentes
df_all = evaluate_ics(A, S, fs;
                      frontal_idxs = frontal_idxs,
                      temporal_idxs = temporal_idxs,
                      eog = eog)
display(df_all)

# Guardar resultados de evaluación en tabla
dir_tables = stage_dir(:ICA_cleaning; kind = :tables)
isdir(dir_tables) || mkpath(dir_tables)

path_table = joinpath(dir_tables, "ICA_components_evaluation.csv")
CSV.write(path_table, df_all)
println("\n💾 Tabla de evaluación de componentes ICA guardada en: $path_table\n")

# -----------------------------------------------------------------------------
# 5. Seleccionar componentes a eliminar y reconstruir EEG limpio
# -----------------------------------------------------------------------------
# Esta sección identifica automáticamente los componentes artefactuales basándose
# en los scores calculados y luego reconstruye los datos EEG eliminando esos componentes.
#
# PROCESO:
#   1. Identificar componentes con artifact_score > umbral
#   2. Agrupar por tipo de artefacto (parpadeo, músculo, línea, salto)
#   3. Anular las componentes artefactuales en S_clean (poner a cero)
#   4. Reconstruir datos limpios: X_clean = A * S_clean
#   5. Guardar datos limpios y metadatos

println("====================================")
println("  Selección de Componentes Artefactuales")
println("====================================\n")

# Identificar componentes artefactuales basándose en los scores
# Por defecto, marcamos como artefactos aquellos con artifact_score > 1.5
# Este umbral puede ajustarse según la calidad de los datos y la sensibilidad deseada
artifact_thresh = 1.5
bad_ics = Int[]  # Lista de componentes a eliminar
artifacts_by_type = Dict{String, Vector{Int}}(
    "parpadeo" => Int[],
    "músculo" => Int[],
    "línea" => Int[],
    "salto" => Int[],
    "otros" => Int[]
)

# Recorrer el DataFrame de evaluación y seleccionar componentes artefactuales
if hasproperty(df_all, :artifact_score) && hasproperty(df_all, :label)
    for (idx, row) in enumerate(eachrow(df_all))
        # Si el score de artefacto supera el umbral y no está etiquetado como neuronal
        if row.artifact_score > artifact_thresh && row.label != "neuronal?"
            push!(bad_ics, row.IC)  # Agregar a la lista de componentes a eliminar
            label = row.label
            
            # Agrupar por tipo de artefacto para estadísticas
            if label == "parpadeo"
                push!(artifacts_by_type["parpadeo"], row.IC)
            elseif label == "músculo"
                push!(artifacts_by_type["músculo"], row.IC)
            elseif label == "línea"
                push!(artifacts_by_type["línea"], row.IC)
            elseif label == "salto"
                push!(artifacts_by_type["salto"], row.IC)
            else
                push!(artifacts_by_type["otros"], row.IC)
            end
        end
    end
end

# Mostrar resumen de componentes artefactuales
println("📊 Resumen de componentes artefactuales (umbral: $artifact_thresh):")
println("-" ^ 60)
total_artifacts = length(bad_ics)
println("  Total de componentes artefactuales: $total_artifacts")

if total_artifacts > 0
    for (tipo, ics) in artifacts_by_type
        if !isempty(ics)
            println("  • $tipo: IC$(join(ics, ", IC")) ($(length(ics)) componente$(length(ics) > 1 ? "s" : ""))")
        end
    end
    println("\n  ICs a eliminar: $(join(bad_ics, ", "))")
else
    println("  ✅ No se detectaron componentes artefactuales con el umbral actual.")
end

# Alternativamente, puedes descomentar y modificar manualmente:
#bad_ics = [1, 5, 7]   # <── MODIFÍCAME según tus inspecciones visuales

println("\n" * "=" ^ 60)
println("  Reconstrucción de EEG Limpio")
println("=" ^ 60 * "\n")

# Crear copia de S y anular las componentes artefactuales
# Al poner S_clean[ic, :] = 0, eliminamos la contribución de ese componente
S_clean = copy(S)
for ic in bad_ics
    @assert 1 ≤ ic ≤ size(S, 1) "IC $ic está fuera del rango [1, $(size(S, 1))]"
    S_clean[ic, :] .= 0.0  # Anular completamente la componente artefactual
end

# Reconstrucción: X_clean = A * S_clean
# La matriz A contiene los mapas topográficos (cómo se mezclan las componentes)
# Al multiplicar A por S_clean, reconstruimos los datos EEG sin los artefactos
# Solo las componentes "buenas" (no anuladas) contribuyen a la reconstrucción
X_clean = A * S_clean

println("✓ Componentes originales: $(size(S, 1))")
println("✓ Componentes eliminados: $(length(bad_ics))")
println("✓ Componentes restantes: $(size(S, 1) - length(bad_ics))")
println("✓ Señal reconstruida (limpia): $(size(X_clean, 1)) canales × $(size(X_clean, 2)) muestras\n")

# Guardar en formato de diccionario por canal
# Convierte la matriz X_clean de vuelta al formato original (diccionario por canal)
# Esto facilita el uso posterior de los datos limpios
dict_EEG_clean = Dict{String, Vector{Float64}}()
for (i, ch) in enumerate(channels)
    dict_EEG_clean[ch] = vec(X_clean[i, :])  # Cada fila de X_clean es un canal
end

path_dict_clean = joinpath(dir_ica, "dict_EEG_ICA_clean.bin")
Serialization.serialize(path_dict_clean, dict_EEG_clean)
println("💾 EEG limpio guardado en: $path_dict_clean")

# -----------------------------------------------------------------------------
# 6. Guardar versión completa con info ICA + cleaning
# -----------------------------------------------------------------------------
dict_ICA_full = Dict(
    "S"        => S,
    "S_clean"  => S_clean,
    "W_total"  => W_total,
    "A"        => A,
    "channels" => channels,
    "bad_ics"  => bad_ics,
    "artifacts_by_type" => artifacts_by_type,
    "artifact_thresh" => artifact_thresh,
    "df_evaluation" => df_all,  # Guardar también el DataFrame de evaluación
    "max_iter" => max_iter,
    "tol"      => tol,
)

path_dict_ica_full = joinpath(dir_ica, "dict_EEG_ICA_full.bin")
Serialization.serialize(path_dict_ica_full, dict_ICA_full)

println("💾 Info completa con bad_ics guardada en: $path_dict_ica_full\n")
println("✨ ICA cleaning finalizado correctamente.")

# -----------------------------------------------------------------------------
# 7. Preparar datos para visualización comparativa
# -----------------------------------------------------------------------------
# Esta sección prepara los datos para comparar visualmente el EEG antes y después
# de la limpieza ICA. Se cargan los datos originales (filtrados) y los datos limpios,
# y se convierten a formato de matriz para facilitar la visualización.

println("=" ^ 60)
println("  Preparación de datos para visualización")
println("=" ^ 60 * "\n")

# Cargar datos antes del ICA (Lowpass)
# Estos son los datos filtrados pero aún con artefactos.
# Recalculamos la ruta localmente para no depender de variables globales previas.
dict_EEG_Lowpass = let
    dir_filtering_local = stage_dir(:filtering)
    path_dict_lowpass_local = joinpath(dir_filtering_local, "dict_EEG_Lowpass.bin")
    Serialization.deserialize(path_dict_lowpass_local)
end

# Cargar datos después del ICA (limpio)
# Estos son los datos después de eliminar componentes artefactuales
dict_EEG_ICA_clean = Serialization.deserialize(path_dict_clean)

# Convertir diccionarios a matrices en el orden correcto de canales
# Usamos el orden de 'channels' que viene del ICA para mantener consistencia
n_channels = length(channels)
n_samples = length(dict_EEG_Lowpass[channels[1]])

# Crear matrices para comparación
# Cada fila es un canal, cada columna es una muestra temporal
eeg_orig = zeros(Float64, n_channels, n_samples)   # Antes del ICA (Lowpass)
eeg_clean = zeros(Float64, n_channels, n_samples)  # Después del ICA (limpio)

# Rellenar matrices con datos de los diccionarios
# Se mantiene el mismo orden de canales para comparación directa
for (i, ch) in enumerate(channels)
    if haskey(dict_EEG_Lowpass, ch)
        eeg_orig[i, :] = dict_EEG_Lowpass[ch]
    else
        println("⚠ Canal $ch no encontrado en dict_EEG_Lowpass")
    end
    
    if haskey(dict_EEG_ICA_clean, ch)
        eeg_clean[i, :] = dict_EEG_ICA_clean[ch]
    else
        println("⚠ Canal $ch no encontrado en dict_EEG_ICA_clean")
    end
end

println("✓ Datos antes del ICA (Lowpass): $(size(eeg_orig, 1)) canales × $(size(eeg_orig, 2)) muestras")
println("✓ Datos después del ICA (limpio): $(size(eeg_clean, 1)) canales × $(size(eeg_clean, 2)) muestras\n")

# -----------------------------------------------------------------------------
# 8. Funciones auxiliares para visualización
# -----------------------------------------------------------------------------

# Alias para compute_psd usando la función existente
compute_psd(signal::AbstractVector, fs::Real) = calculate_PSD_signal(signal, fs)

# Función para calcular kurtosis de un canal
kurtosis_channel(signal::AbstractVector) = kurtosis(signal)

"""
    plot_eeg_segment(eeg_orig, eeg_clean, fs; channels=[1], t_start=0.0, t_len=5.0, chan_labels=nothing)

Visualiza segmentos temporales de EEG comparando datos originales vs limpios.

Esta función permite inspeccionar visualmente la efectividad de la limpieza ICA
mostrando superpuestas las señales antes y después del procesamiento. Es útil para
verificar que los artefactos se han eliminado correctamente sin afectar la señal neuronal.

Parámetros:
  - eeg_orig, eeg_clean: matrices (canales × muestras) - datos antes y después
  - fs: frecuencia de muestreo
  - channels: índices de canales a mostrar (puede ser un vector)
  - t_start: inicio del segmento en segundos
  - t_len: duración de la ventana en segundos
  - chan_labels: opcional, nombres de canales para etiquetar los gráficos

Cada canal se muestra en un subplot separado con ambas señales superpuestas.
"""
function plot_eeg_segment(eeg_orig::AbstractMatrix,
                          eeg_clean::AbstractMatrix,
                          fs::Real;
                          channels = [1],
                          t_start = 0.0,
                          t_len = 5.0,
                          chan_labels = nothing)

    n_chan, n_samp = size(eeg_orig)
    @assert size(eeg_clean) == size(eeg_orig) "Dimensiones antes/después deben coincidir"

    # índices de muestra
    i_start = max(1, Int(floor(t_start * fs)) + 1)
    i_end   = min(n_samp, i_start + Int(floor(t_len * fs)) - 1)
    t = (i_start:i_end) ./ fs

    n_plots = length(channels)
    plt = plot(layout = (n_plots, 1), size=(900, 250*n_plots))

    for (i, ch) in enumerate(channels)
        y_orig  = eeg_orig[ch, i_start:i_end]
        y_clean = eeg_clean[ch, i_start:i_end]

        ch_name = isnothing(chan_labels) ? "Ch $ch" : string(chan_labels[ch])

        plot!(plt[i], t, y_orig,
              label = "Original",
              title = "Canal $ch_name (t = $(t_start)s–$(t_start+t_len)s)",
              xlabel = "Tiempo (s)",
              ylabel = "µV")

        plot!(plt[i], t, y_clean,
              label = "Limpio",
              lw = 1.5)
    end

    display(plt)
    return plt
end


# -----------------------------------------------------------------------------
# 9. Visualización comparativa: antes vs después del ICA
# -----------------------------------------------------------------------------
println("=" ^ 60)
println("  Visualización Comparativa")
println("=" ^ 60 * "\n")

# Ejemplo: visualizar segmentos de EEG de los canales 1, 10 y 20
# desde el segundo 10 hasta el segundo 14
println("📊 Visualizando segmentos de EEG (canales 1, 10, 20)...")
plot_eeg_segment(eeg_orig, eeg_clean, fs;
                 channels = [1, 10, 20],
                 t_start = 10.0,
                 t_len = 4.0,
                 chan_labels = channels)   

                 """
plot_channel_psd(eeg_orig, eeg_clean, fs;
                 chan = 1,
                 fmax = 80.0,
                 chan_labels = nothing)
"""
function plot_channel_psd(eeg_orig::AbstractMatrix,
                          eeg_clean::AbstractMatrix,
                          fs::Real;
                          chan::Int = 1,
                          fmax::Real = 80.0,
                          chan_labels = nothing)

    x_orig  = eeg_orig[chan, :]
    x_clean = eeg_clean[chan, :]

    f_orig, psd_orig  = compute_psd(x_orig, fs)
    f_clean, psd_clean = compute_psd(x_clean, fs)

    # limitar a fmax
    idx_orig  = f_orig .<= fmax
    idx_clean = f_clean .<= fmax

    ch_name = isnothing(chan_labels) ? "Ch $chan" : string(chan_labels[chan])

    plt = plot(f_orig[idx_orig], psd_orig[idx_orig],
               label = "Original",
               xlabel = "Frecuencia (Hz)",
               ylabel = "PSD (u.a.)",
               title = "PSD antes vs después – $ch_name",
               yscale = :log10)

    plot!(plt, f_clean[idx_clean], psd_clean[idx_clean],
          label = "Limpio")

    display(plt)
    return plt
end

println("\n📊 Visualizando PSD del canal 1...")
plot_channel_psd(eeg_orig, eeg_clean, fs;
                 chan = 1,
                 fmax = 80.0,
                 chan_labels = channels)

"""
    compute_kurtosis_per_channel(eeg::AbstractMatrix)

Calcula kurtosis por canal de una matriz canales × muestras.

La kurtosis mide la "cola" de la distribución. Valores altos indican señales
con picos pronunciados o valores extremos (típico de artefactos).
Después de la limpieza ICA, la kurtosis debería disminuir en canales afectados por artefactos.

Devuelve un vector con la kurtosis de cada canal.
"""
function compute_kurtosis_per_channel(eeg::AbstractMatrix)
    n_chan, _ = size(eeg)
    k = zeros(n_chan)
    for ch in 1:n_chan
        k[ch] = kurtosis_channel(eeg[ch, :])
    end
    return k
end

"""
    compute_energy_per_channel(eeg::AbstractMatrix)

Calcula energía media por canal (mean(x.^2)).

La energía es una medida de la potencia de la señal. Después de la limpieza ICA,
la energía debería disminuir ligeramente (se eliminan artefactos) pero no demasiado
(si disminuye mucho, puede indicar que se eliminaron componentes neuronales importantes).

Devuelve un vector con la energía media de cada canal.
"""
function compute_energy_per_channel(eeg::AbstractMatrix)
    n_chan, _ = size(eeg)
    e = zeros(n_chan)
    for ch in 1:n_chan
        x = eeg[ch, :]
        e[ch] = mean(x.^2)  # energía = promedio del cuadrado
    end
    return e
end

"""
plot_kurtosis_energy(eeg_orig, eeg_clean; chan_labels = nothing)

Hace dos figuras con comparación superpuesta:
1) Kurtosis antes vs después por canal (en la misma gráfica)
2) Energía antes vs después por canal (en la misma gráfica)
"""
function plot_kurtosis_energy(eeg_orig::AbstractMatrix,
                              eeg_clean::AbstractMatrix;
                              chan_labels = nothing)

    @assert size(eeg_orig) == size(eeg_clean)

    n_chan, _ = size(eeg_orig)
    chans = 1:n_chan

    labels_x = isnothing(chan_labels) ? string.(chans) : string.(chan_labels)

    # --- KURTOSIS ---
    k_orig  = compute_kurtosis_per_channel(eeg_orig)
    k_clean = compute_kurtosis_per_channel(eeg_clean)

    # Gráfico de barras superpuestas para Kurtosis
    plt_k = bar(labels_x, [k_orig k_clean],
                label = ["Original" "Limpio"],
                xlabel = "Canal",
                ylabel = "Kurtosis",
                title = "Kurtosis: Antes vs Después de ICA",
                color = [:blue :green],
                alpha = 0.7,
                xticks = (1:length(labels_x), labels_x),
                xrotation = 45,
                size = (1200, 400),
                legend = :topright)

    display(plt_k)

    # --- ENERGÍA ---
    e_orig  = compute_energy_per_channel(eeg_orig)
    e_clean = compute_energy_per_channel(eeg_clean)

    # Gráfico de barras superpuestas para Energía
    plt_e = bar(labels_x, [e_orig e_clean],
                label = ["Original" "Limpio"],
                xlabel = "Canal",
                ylabel = "Energía media (µV²)",
                title = "Energía: Antes vs Después de ICA",
                color = [:blue :green],
                alpha = 0.7,
                xticks = (1:length(labels_x), labels_x),
                xrotation = 45,
                size = (1200, 400),
                legend = :topright)

    display(plt_e)

    return plt_k, plt_e
end

println("\n📊 Visualizando Kurtosis y Energía por canal...")
plot_kurtosis_energy(eeg_orig, eeg_clean; chan_labels = channels)

println("\n✨ Visualización comparativa completada.")
