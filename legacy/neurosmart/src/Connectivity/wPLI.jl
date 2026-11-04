#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# src/Connectivity/wPLI.jl
#
# CONECTIVIDAD / wPLI (Weighted Phase Lag Index)
# ===============================================
# Esta rutina calcula el índice de desfase de fase ponderado (wPLI) por bandas de frecuencia
# a partir de datos EEG en espacio CSD (output de CSD.jl).
#
# DESCRIPCIÓN:
# El wPLI es una medida de conectividad funcional que cuantifica la sincronización de fase
# entre pares de canales EEG. A diferencia del PLI estándar, el wPLI pondera las contribuciones
# por la magnitud de la parte imaginaria del producto cruzado espectral, reduciendo la sensibilidad
# al ruido y mejorando la robustez estadística.
#
# MÉTODO (ESTILO BRAINVISION ANALYZER):
# 1. Filtrado bandpass por banda de frecuencia (Butterworth orden 8 + filtfilt para fase estable)
# 2. Transformada de Hilbert para obtener señal analítica (fase instantánea)
# 3. Cálculo de wPLI agregando sobre TODAS las muestras de TODOS los segmentos simultáneamente
#    (enfoque "across segments" de Analyzer, no promediar métricas por segmento)
# 4. Matriz de conectividad simétrica con diagonales en 0.0
#
# PROCESO:
#   1. Carga datos EEG en espacio CSD desde dict_csd.bin
#   2. Carga definición de bandas de frecuencia desde dict_FFT_power.bin (o dict_FFT.bin)
#   3. Para cada banda de frecuencia:
#      a. Filtrado bandpass (orden 8) de cada canal y segmento
#      b. Transformada de Hilbert para obtener señal analítica
#      c. Agregación de señales analíticas en matrices (samples × segments) por canal
#      d. Cálculo de wPLI agregando sobre todas las muestras y segmentos simultáneamente
#         (enfoque "across segments" estilo Analyzer)
#   4. Guarda matrices de conectividad, tablas y figuras
#
# ENTRADA ESPERADA:
#   - data/CSD/dict_csd.bin: datos EEG transformados con CSD
#   - data/FFT/dict_FFT_power.bin (o dict_FFT.bin): definición de bandas de frecuencia
#
# SALIDAS:
#   - data/wPLI/dict_wpli.bin: diccionario completo con resultados
#   - results/tables/wPLI/*.csv: matrices y listas de conexiones (edges)
#   - results/figures/wPLI/*.png: mapas de calor (heatmaps) de conectividad

# ------------------------------------------------------------------------------------
# IMPORTACIÓN DE LIBRERÍAS
# ------------------------------------------------------------------------------------
using CSV              # Lectura y escritura de archivos CSV
using DataFrames        # Manipulación de datos tabulares
using Dates             # Formateo de fechas para logs
using Serialization     # Serialización de datos binarios
using Statistics        # Funciones estadísticas básicas
using LinearAlgebra     # Operaciones de álgebra lineal
using FFTW              # Transformada rápida de Fourier (para Hilbert)
using DSP               # Procesamiento de señales digitales (filtros)
# Importación selectiva de CairoMakie para evitar ambigüedades con otros módulos
using CairoMakie: Figure, Axis, heatmap!, Colorbar, save

# Si se ejecuta este script directamente (fuera del módulo EEG_Julia),
# cargamos utilidades de rutas para disponer de `stage_dir`.
if !@isdefined(stage_dir)
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end

# ------------------------------------------------------------------------------------
# 1. CONFIGURACIÓN DE RUTAS
# ------------------------------------------------------------------------------------
# Definición de rutas relativas desde el directorio del script
# Nota: @__DIR__ apunta a src/Connectivity/, por lo que subimos dos niveles (.., ..)
#       para llegar a la raíz del proyecto

# Rutas de entrada: datos CSD
dir_csd = stage_dir(:CSD)
path_dict_csd = joinpath(dir_csd, "dict_csd.bin")

# Rutas de entrada: definición de bandas de frecuencia desde FFT
dir_fft = stage_dir(:FFT)
path_dict_fft_power = joinpath(dir_fft, "dict_FFT_power.bin")  # Intento primero con dict_FFT_power
path_dict_fft = joinpath(dir_fft, "dict_FFT.bin")              # Fallback a dict_FFT

# Rutas de salida: datos procesados con wPLI
dir_wpli_data = stage_dir(:wPLI)
out_dict_wpli = joinpath(dir_wpli_data, "dict_wpli.bin")

# Rutas de salida: logs y resultados
dir_logs    = stage_dir(:wPLI; kind = :logs)      # Archivos de log
dir_tables  = stage_dir(:wPLI; kind = :tables)    # Tablas de resultados
dir_figures = stage_dir(:wPLI; kind = :figures)   # Figuras y visualizaciones

# ------------------------------------------------------------------------------------
# 2. LIMPIEZA DE RESULTADOS PREVIOS
# ------------------------------------------------------------------------------------
# Se eliminan todos los archivos previos antes de iniciar el procesamiento
# para evitar confusiones con resultados antiguos y asegurar reproducibilidad.
println("✔ Eliminando resultados previos")

# Función auxiliar para limpiar archivos de una carpeta por extensión
# Parámetros:
#   - dir: directorio a limpiar
#   - extension: extensión de archivos a eliminar (ej: ".png", ".bin", ".csv")
#   - descripcion: descripción para mensajes informativos
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

# Limpiar carpetas de resultados previos
limpiar_carpeta(dir_figures, ".png", "figuras previas")
limpiar_carpeta(dir_wpli_data, ".bin", "datos previos")
limpiar_carpeta(dir_tables, ".csv", "tablas previas")
limpiar_carpeta(dir_logs, ".log", "logs previos")

println("✔ Limpieza completada")
println()

# Crear carpetas si no existen (necesario antes de escribir logs y resultados)
mkpath(dir_wpli_data)
mkpath(dir_logs)
mkpath(dir_tables)
mkpath(dir_figures)

# Generar identificador único para esta ejecución (timestamp)
run_id   = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
log_path = joinpath(dir_logs, "wPLI_" * run_id * ".log")

# ------------------------------------------------------------------------------------
# 3. SISTEMA DE LOGGING
# ------------------------------------------------------------------------------------
# Función auxiliar para escribir mensajes con timestamp en el archivo de log
function logmsg(io, s)
    println(io, "[$(Dates.format(now(), "HH:MM:SS"))] $s")
end

# Inicializar archivo de log con información de inicio
open(log_path, "w") do io
    logmsg(io, "Inicio wPLI.jl")
    logmsg(io, "path_dict_csd = $path_dict_csd")
    logmsg(io, "out_dict_wpli = $out_dict_wpli")
end

# ------------------------------------------------------------------------------------
# 4. CARGA DE DATOS CSD
# ------------------------------------------------------------------------------------
# Se cargan los datos EEG transformados con CSD.
# Formato esperado: hipermatriz 3D (canales × muestras × segmentos)
println("=" ^ 80)
println("📥 CARGA DE DATOS CSD")
println("=" ^ 80)
println()

# Deserializar diccionario con datos y metadatos CSD
dict_csd = Serialization.deserialize(path_dict_csd)

# Extraer datos e información del diccionario
eeg_csd   = dict_csd["eeg_csd"]          # Array 3D: (canales, muestras, segmentos)
ch_names  = dict_csd["channels"]         # Vector{String}: nombres de canales
fs        = Float64(dict_csd["fs"])      # Float64: frecuencia de muestreo (Hz)

n_ch, n_samp, n_seg = size(eeg_csd)

# Mostrar información de los datos cargados
println("✔ CSD cargado: (ch, samples, seg) = $(size(eeg_csd))")
println("  → Frecuencia de muestreo: $fs Hz")
println("  → Número de canales: $n_ch")
println("  → Muestras por segmento: $n_samp")
println("  → Número de segmentos: $n_seg")
println("  → Duración de cada segmento: $(round(n_samp/fs, digits=3)) s")
println()

# Registrar información en el log
open(log_path, "a") do io
    logmsg(io, "Datos CSD: n_ch=$n_ch, n_samp=$n_samp, n_seg=$n_seg, fs=$fs")
end

# ------------------------------------------------------------------------------------
# 5. CARGA DE DEFINICIÓN DE BANDAS DE FRECUENCIA
# ------------------------------------------------------------------------------------
# Se cargan las bandas de frecuencia definidas previamente en el análisis FFT.
# Esto asegura consistencia entre diferentes análisis del mismo dataset.
println("=" ^ 80)
println("📊 CARGA DE BANDAS DE FRECUENCIA DESDE FFT")
println("=" ^ 80)
println()

# Intentar cargar dict_FFT_power.bin primero, luego dict_FFT.bin como fallback
# Si no se encuentran bandas, usar las definiciones estándar de FFT.jl
bands_hz = nothing
source_file = ""

if isfile(path_dict_fft_power)
    dict_fft_power = Serialization.deserialize(path_dict_fft_power)
    bands_hz = get(dict_fft_power, "bands_hz", nothing)
    if bands_hz !== nothing
        source_file = basename(path_dict_fft_power)
        println("✔ Bandas cargadas desde: $source_file")
    end
end

# Si no se encontraron en dict_FFT_power, intentar dict_FFT.bin
if bands_hz === nothing && isfile(path_dict_fft)
    dict_fft = Serialization.deserialize(path_dict_fft)
    bands_hz = get(dict_fft, "bands_hz", nothing)
    if bands_hz !== nothing
        source_file = basename(path_dict_fft)
        println("✔ Bandas cargadas desde: $source_file")
    end
end

# Si aún no se encontraron, usar definiciones estándar de FFT.jl
# (mismas bandas que se definen en FFT.jl líneas 598-606)
if bands_hz === nothing
    println("⚠ No se encontraron bandas en archivos FFT. Usando definiciones estándar.")
    bands_hz = [
        (:DELTA,     (0.5, 4.0)),
        (:THETA,     (4.0, 8.0)),
        (:ALPHA,     (7.8, 11.7)),
        (:BETA_LOW,  (12.0, 15.0)),
        (:BETA_MID,  (15.0, 18.0)),
        (:BETA_HIGH, (18.0, 30.0)),
        (:GAMMA,     (30.0, 50.0))
    ]
    source_file = "definiciones estándar"
    println("✔ Usando bandas estándar (mismas que FFT.jl)")
end

# Convertir bandas a formato Dict{String, Tuple{Float64, Float64}} para compatibilidad
# bands_hz viene como Vector{Tuple{Symbol, Tuple{Float64, Float64}}}
bands = Dict{String, Tuple{Float64, Float64}}()
for (sym, (f1, f2)) in bands_hz
    bands[string(sym)] = (Float64(f1), Float64(f2))
end

# Mostrar bandas cargadas
println("  → Bandas de frecuencia cargadas:")
for (bname, (f1, f2)) in bands
    println("    • $(bname): $(f1)-$(f2) Hz")
end
println()

# Registrar información en el log
open(log_path, "a") do io
    logmsg(io, "Bandas cargadas desde: $source_file")
    logmsg(io, "Bandas cargadas: $(length(bands)) bandas")
    for (bname, (f1, f2)) in bands
        logmsg(io, "  $(bname): $(f1)-$(f2) Hz")
    end
end

# ------------------------------------------------------------------------------------
# 6. FUNCIONES AUXILIARES: FILTRADO Y TRANSFORMADA DE HILBERT
# ------------------------------------------------------------------------------------
# Implementación de filtrado bandpass y transformada de Hilbert para obtener
# la señal analítica (fase instantánea) necesaria para el cálculo de wPLI.

# ------------------------------------------------------------------------------------
# 6.1. Filtro Bandpass con filtfilt
# ------------------------------------------------------------------------------------
# Filtro Butterworth bandpass con filtrado bidireccional (filtfilt) para fase estable (0-lag).
# Esto es crítico para el cálculo de wPLI, ya que preserva la fase original de la señal.
# En BrainVision Analyzer, los filtros bandpass para transformaciones tipo ERS/ERD y
# demodulación se especifican típicamente como Butterworth zero-phase de orden 8.
# Parámetros:
#   - x: señal temporal a filtrar
#   - fs: frecuencia de muestreo (Hz)
#   - f1: frecuencia inferior de la banda (Hz)
#   - f2: frecuencia superior de la banda (Hz)
#   - order: orden del filtro Butterworth (default: 8, estilo Analyzer)
# Retorna:
#   - Vector{Float64}: señal filtrada con fase preservada
function bandpass_filtfilt(x::AbstractVector{<:Real}, fs::Real, f1::Real, f2::Real;
                           order::Int=8)
    # Normalizar frecuencias (0-1) dividiendo por la frecuencia de Nyquist (fs/2)
    wn_low = f1 / (fs/2)
    wn_high = f2 / (fs/2)
    
    # Crear filtro bandpass (sin argumento fs, ya que las frecuencias están normalizadas)
    responsetype = Bandpass(wn_low, wn_high)
    designmethod = Butterworth(order)
    bp = digitalfilter(responsetype, designmethod)
    
    # Aplicar filtrado bidireccional (filtfilt) para fase cero
    return filtfilt(bp, Float64.(x))
end

# ------------------------------------------------------------------------------------
# 6.2. Transformada de Hilbert (Señal Analítica)
# ------------------------------------------------------------------------------------
# Calcula la señal analítica mediante transformada de Hilbert usando FFT.
# La señal analítica z(t) = x(t) + i·H[x(t)] contiene información de amplitud y fase.
# Parámetros:
#   - x: señal temporal real
# Retorna:
#   - Vector{ComplexF64}: señal analítica del mismo tamaño que x
function analytic_signal(x::AbstractVector{<:Real})
    N = length(x)
    X = fft(Float64.(x))
    h = zeros(Float64, N)

    # Construir filtro de Hilbert en el dominio de frecuencia
    # Para N par: h[1]=1, h[N/2+1]=1, h[2:N/2]=2
    # Para N impar: h[1]=1, h[2:(N+1)/2]=2
    if iseven(N)
        h[1] = 1.0
        h[N÷2 + 1] = 1.0
        h[2:N÷2] .= 2.0
    else
        h[1] = 1.0
        h[2:(N+1)÷2] .= 2.0
    end

    return ifft(X .* h)
end

# ------------------------------------------------------------------------------------
# 7. FUNCIÓN CORE: CÁLCULO DE wPLI (ESTILO ANALYZER)
# ------------------------------------------------------------------------------------
# Calcula el Weighted Phase Lag Index (wPLI) agregando sobre TODAS las muestras
# de TODOS los segmentos simultáneamente (enfoque "across segments" de BrainVision Analyzer).
# 
# En Analyzer, la filosofía de conectividad en frecuencia se apoya en cantidades tipo
# cross-spectrum agregadas "across segments" (promediadas/agrupadas a través de segmentos),
# no en calcular una métrica por segmento y luego promediarla.
#
# wPLI(i,j) = |sum(Im)| / sum(|Im|)
# donde Im = imag(z_i * conj(z_j)), agregando sobre todas las muestras y segmentos.
# 
# Este enfoque es más estable con ventanas de 1s y más cercano al método de Analyzer
# de promediar cantidades espectrales a través de segmentos.
# Parámetros:
#   - zi: señal analítica del canal i (Matrix{ComplexF64} de tamaño samples × segments)
#   - zj: señal analítica del canal j (Matrix{ComplexF64} de tamaño samples × segments)
# Retorna:
#   - Float64: valor de wPLI entre 0 y 1
function wpli_from_analytic_agg(zi::AbstractMatrix{ComplexF64}, zj::AbstractMatrix{ComplexF64})
    # Parte imaginaria del producto cruzado espectral (agregado sobre todas las muestras y segmentos)
    imv = imag.(zi .* conj.(zj))
    
    # Numerador: valor absoluto de la suma (agregación "across segments")
    num = abs(sum(imv))
    
    # Denominador: suma del valor absoluto (con epsilon para evitar división por cero)
    den = sum(abs.(imv)) + eps()
    
    return num / den
end

# ------------------------------------------------------------------------------------
# 8. CÁLCULO DE wPLI POR BANDAS DE FRECUENCIA (ESTILO ANALYZER)
# ------------------------------------------------------------------------------------
# Para cada banda de frecuencia definida, se calcula la matriz de conectividad wPLI
# entre todos los pares de canales, agregando sobre TODAS las muestras de TODOS los
# segmentos simultáneamente (enfoque "across segments" de BrainVision Analyzer).
println("=" ^ 80)
println("🧠 CÁLCULO wPLI POR BANDAS DE FRECUENCIA (ESTILO ANALYZER)")
println("=" ^ 80)
println()

# Diccionario para almacenar resultados por banda
results = Dict{String, Any}()

# Iterar sobre bandas en orden establecido (reproducibilidad)
# Ordenar las bandas para asegurar orden consistente entre ejecuciones
band_names_ordered = sort(collect(keys(bands)))

for bname in band_names_ordered
    f1, f2 = bands[bname]
    println("→ Banda $(bname) ($(f1)-$(f2) Hz)")

    # Precomputar señales analíticas por canal
    # z[c] será una matriz (samples × segments) para facilitar agregación "across segments"
    z = Vector{Matrix{ComplexF64}}(undef, n_ch)

    @inbounds for c in 1:n_ch
        # Matriz para almacenar señal analítica de este canal (samples × segments)
        Zc = Matrix{ComplexF64}(undef, n_samp, n_seg)
        
        for s in 1:n_seg
            # Extraer señal del canal c en el segmento s
            x = @view eeg_csd[c, :, s]
            
            # Filtrar bandpass para la banda actual (orden 8, estilo Analyzer)
            xf = bandpass_filtfilt(x, fs, f1, f2; order=8)
            
            # Calcular señal analítica (Hilbert)
            Zc[:, s] = analytic_signal(xf)
        end
        
        z[c] = Zc
    end

    # Matriz wPLI (simétrica, diagonal 0.0 como en Analyzer)
    # W[i,j] = wPLI entre canal i y canal j
    W = zeros(Float64, n_ch, n_ch)

    @inbounds for i in 1:n_ch
        for j in (i+1):n_ch
            # Calcular wPLI agregando sobre todas las muestras y segmentos
            # (enfoque "across segments" estilo Analyzer)
            w = wpli_from_analytic_agg(z[i], z[j])
            
            # Matriz simétrica
            W[i,j] = w
            W[j,i] = w
        end
    end
    
    # Diagonales en 0.0 (Analyzer trabaja con pares, no auto-conectividad)
    # Ya están en 0.0 por zeros(), pero lo explicitamos para claridad

    # Guardar matriz wPLI como CSV (formato tabla con nombres de canales)
    mat_csv = joinpath(dir_tables, "wPLI_$(bname)_matrix.csv")
    dfW = DataFrame(W, Symbol.(ch_names))
    dfW = hcat(DataFrame(channel = ch_names), dfW)
    CSV.write(mat_csv, dfW)
    println("  ✔ Matriz guardada: $(basename(mat_csv))")

    # Guardar lista de conexiones (edges) como CSV
    # Formato: from, to, wpli (solo pares únicos, i < j)
    edges_csv = joinpath(dir_tables, "wPLI_$(bname)_edges.csv")
    rows = DataFrame(from=String[], to=String[], wpli=Float64[])
    for i in 1:n_ch, j in (i+1):n_ch
        push!(rows, (ch_names[i], ch_names[j], W[i,j]))
    end
    CSV.write(edges_csv, rows)
    println("  ✔ Edges guardados: $(basename(edges_csv))")

    # Generar figura heatmap de conectividad
    fig_png = joinpath(dir_figures, "wPLI_$(bname)_heatmap.png")
    let
        f = Figure(size=(900, 750))
        ax = Axis(f[1,1], 
                  title="wPLI - $(bname) ($(f1)-$(f2) Hz)", 
                  xlabel="Canal", 
                  ylabel="Canal")
        hm = heatmap!(ax, W)
        Colorbar(f[1,2], hm)
        save(fig_png, f)
    end
    println("  ✔ Figura guardada: $(basename(fig_png))")
    println()

    # Registrar en el log
    open(log_path, "a") do io
        logmsg(io, "Banda $bname: matrix=$mat_csv, edges=$edges_csv, fig=$fig_png")
    end

    # Almacenar resultados en diccionario
    results[bname] = Dict(
        "band" => (f1, f2),
        "W" => W,
        "matrix_csv" => mat_csv,
        "edges_csv" => edges_csv,
        "figure_png" => fig_png
    )
end

# ------------------------------------------------------------------------------------
# 9. GUARDADO DE RESULTADOS Y METADATOS
# ------------------------------------------------------------------------------------
# Se guarda un diccionario completo con todos los resultados, metadatos y referencias
# a los archivos generados, siguiendo el mismo patrón que otros scripts del pipeline.
println("=" ^ 80)
println("💾 GUARDADO DE RESULTADOS Y METADATOS")
println("=" ^ 80)
println()

# Construir diccionario completo con metadatos
dict_wpli = Dict{String, Any}()
dict_wpli["space"] = get(dict_csd, "space", "CSD")      # Espacio de los datos (CSD)
dict_wpli["input"] = "dict_csd.bin"                     # Archivo de entrada
dict_wpli["fs"] = fs                                    # Frecuencia de muestreo
dict_wpli["channels"] = ch_names                        # Nombres de canales
dict_wpli["bands"] = bands                              # Definición de bandas usadas
dict_wpli["wpli"] = results                             # Resultados por banda

# Asegurar que el directorio existe
mkpath(dir_wpli_data)

# Serializar y guardar resultados
Serialization.serialize(out_dict_wpli, dict_wpli)

println("✔ Guardado final: $(basename(out_dict_wpli))")
println("  → Ubicación: $(abspath(out_dict_wpli))")
println()

# Registrar guardado en el log
open(log_path, "a") do io
    logmsg(io, "Guardado dict_wpli: $out_dict_wpli")
    logmsg(io, "Fin OK")
end

# ------------------------------------------------------------------------------------
# 10. VERIFICACIÓN AUTOMÁTICA DE RESULTADOS
# ------------------------------------------------------------------------------------
# Chequeos rápidos de sanidad para detectar errores en la implementación
function verificar_resultados_wpli()
    all_ok = true
    
    for bname in band_names_ordered
        W = results[bname]["W"]
        
        # Chequeos: rango [0,1], diagonal=0, simetría, NaN/Inf
        minW, maxW = minimum(W), maximum(W)
        diag_max = maximum(abs.(diag(W)))
        sym_err = maximum(abs.(W .- W'))
        has_invalid = any(isnan.(W)) || any(isinf.(W))
        
        band_ok = (minW ≥ -1e-8) && (maxW ≤ 1.0+1e-8) && (diag_max ≤ 1e-6) && 
                  (sym_err ≤ 1e-6) && !has_invalid
        
        if band_ok
            upper_tri_vals = [W[i,j] for i in 1:n_ch for j in (i+1):n_ch]
            println("  ✔ $bname: OK (media wPLI = $(round(mean(upper_tri_vals), digits=4)))")
        else
            println("  ⚠ $bname: FALLÓ verificación")
            (minW < -1e-8 || maxW > 1.0+1e-8) && println("    - valores fuera de [0,1]")
            (diag_max > 1e-6) && println("    - diagonal ≠ 0")
            (sym_err > 1e-6) && println("    - no simétrica")
            has_invalid && println("    - contiene NaN/Inf")
            all_ok = false
        end
    end
    
    return all_ok
end

println("=" ^ 80)
println("🔍 VERIFICACIÓN DE RESULTADOS")
println("=" ^ 80)
println()

verificacion_ok = verificar_resultados_wpli()

println()
println(verificacion_ok ? "✅ Todas las verificaciones pasaron correctamente" : 
                          "⚠️  Algunas verificaciones fallaron - revisar implementación")
println()

# Registrar verificación en el log
open(log_path, "a") do io
    logmsg(io, "Verificación: " * (verificacion_ok ? "OK" : "FALLOS DETECTADOS"))
end

println("=" ^ 80)
println("✅ PROCESAMIENTO wPLI COMPLETADO")
println("=" ^ 80)

