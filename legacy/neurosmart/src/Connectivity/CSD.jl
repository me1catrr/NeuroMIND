#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# src/Connectivity/CSD.jl
#
# CONECTIVIDAD / REDUCCIÓN DE CONDUCCIÓN DE VOLUMEN (CSD)
# ========================================================
# Esta rutina realiza la reducción de conducción de volumen (Current Source Density, CSD)
# de los datos EEG mediante el método de Perrin (spherical spline Laplacian).
#
# DESCRIPCIÓN:
# El CSD es una técnica de transformación espacial que reduce los efectos de la conducción
# de volumen en las señales EEG. A diferencia de la referencia común o promedio, el CSD
# proporciona una estimación de la actividad eléctrica local en cada electrodo, eliminando
# la influencia de fuentes distantes y mejorando la resolución espacial.
#
# MÉTODO:
# Se utiliza el algoritmo de Perrin et al. (1989) basado en splines esféricos, que es
# equivalente al método implementado en BrainVision Analyzer. Este método:
#   - Normaliza las posiciones de los electrodos a una esfera unidad
#   - Construye matrices G y H usando polinomios de Legendre
#   - Aplica regularización (lambda) para estabilidad numérica
#   - Genera un operador Laplaciano que transforma los datos EEG
#
# PROCESO:
#   1. Carga datos EEG segmentados tras la segunda corrección de baseline
#   2. Carga información de posiciones de electrodos desde archivo TSV
#   3. Verifica consistencia entre canales EEG y posiciones de electrodos
#   4. Construye el operador CSD usando coordenadas esféricas
#   5. Aplica el operador a cada segmento temporal
#   6. Guarda resultados y metadatos
#
# PARÁMETROS CSD (valores por defecto estilo BrainVision):
#   - m: orden del spline (default: 4)
#   - degree: grado máximo de polinomios de Legendre (default: 10)
#   - lambda: parámetro de regularización (default: 1e-5)

# ------------------------------------------------------------------------------------
# IMPORTACIÓN DE LIBRERÍAS
# ------------------------------------------------------------------------------------
using CSV              # Lectura de archivos TSV con posiciones de electrodos
using DataFrames        # Manipulación de datos tabulares
using Dates             # Formateo de fechas para logs
using Serialization     # Serialización de datos binarios
using Statistics        # Funciones estadísticas básicas
using Printf            # Formateo de salida
using LinearAlgebra     # Operaciones de álgebra lineal (matrices, identidad)
# Importación selectiva de CairoMakie para evitar ambigüedades con otros módulos de visualización
# (Plots también exporta scatter!, por lo que importamos selectivamente solo lo necesario)
using CairoMakie: Figure, Axis, scatter!, lines!, text!, Colorbar, barplot!,
                  heatmap!, contour!, axislegend, save, hidespines!, hidedecorations!

# Si se ejecuta este script directamente (fuera del módulo EEG_Julia),
# cargamos utilidades de rutas para disponer de `stage_dir` y `electrodes_dir`.
if !@isdefined(stage_dir) || !@isdefined(electrodes_dir)
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end

# ------------------------------------------------------------------------------------
# 1. CONFIGURACIÓN DE RUTAS
# ------------------------------------------------------------------------------------
# Definición de rutas relativas desde el directorio del script
# Nota: @__DIR__ apunta a src/Connectivity/, por lo que subimos dos niveles (.., ..)
#       para llegar a la raíz del proyecto

# Rutas de entrada: datos después del 2º baseline correction
dir_baseline = stage_dir(:baseline)
path_dict_2nd_baseline = joinpath(dir_baseline, "dict_2nd_baseline_correction.bin")

# Rutas de entrada: información de posiciones de electrodos
dir_electrodes = electrodes_dir()
electrode_tsv  = joinpath(dir_electrodes, "sub-M05_ses-T2_electrodes.tsv")

# Rutas de salida: datos procesados con CSD
dir_csd      = stage_dir(:CSD)
out_eeg_path = joinpath(dir_csd, "eeg_csd.bin")      # Datos EEG transformados
out_dict_path = joinpath(dir_csd, "dict_csd.bin")    # Metadatos y diccionario completo

# Rutas de salida: logs y resultados
dir_logs    = stage_dir(:CSD; kind = :logs)      # Archivos de log
dir_tables  = stage_dir(:CSD; kind = :tables)    # Tablas de resultados
dir_figures = stage_dir(:CSD; kind = :figures)   # Figuras y visualizaciones

# ------------------------------------------------------------------------------------
# 2. LIMPIEZA DE RESULTADOS PREVIOS
# ------------------------------------------------------------------------------------
# Se eliminan todos los archivos previos antes de iniciar el procesamiento
# para evitar confusiones con resultados antiguos y asegurar reproducibilidad.
println("✔ Eliminando resultados previos")

# Función auxiliar para limpiar archivos de una carpeta por extensión
# Parámetros:
#   - dir: directorio a limpiar
#   - extension: extensión de archivos a eliminar (ej: ".png", ".bin")
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
limpiar_carpeta(dir_csd, ".bin", "datos previos")
limpiar_carpeta(dir_tables, ".csv", "tablas previas")
limpiar_carpeta(dir_logs, ".log", "logs previos")

println("✔ Limpieza completada")
println()

# Crear carpetas si no existen (necesario antes de escribir logs y resultados)
mkpath(dir_csd)
mkpath(dir_logs)
mkpath(dir_tables)
mkpath(dir_figures)

# Generar identificador único para esta ejecución (timestamp)
run_id   = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
log_path = joinpath(dir_logs, "CSD_" * run_id * ".log")

# ------------------------------------------------------------------------------------
# 3. SISTEMA DE LOGGING
# ------------------------------------------------------------------------------------
# Función auxiliar para escribir mensajes con timestamp en el archivo de log
function logmsg(io, s)
    println(io, "[$(Dates.format(now(), "HH:MM:SS"))] $s")
end

# Inicializar archivo de log con información de inicio
open(log_path, "w") do io
    logmsg(io, "Inicio CSD.jl")
    logmsg(io, "path_dict_2nd_baseline = $path_dict_2nd_baseline")
    logmsg(io, "electrode_tsv = $electrode_tsv")
    logmsg(io, "out_eeg_path = $out_eeg_path")
end

# ------------------------------------------------------------------------------------
# 4. CARGA DE DATOS EEG DESPUÉS DEL 2ND BASELINE CORRECTION
# ------------------------------------------------------------------------------------
# Se cargan los datos EEG segmentados tras la segunda corrección de baseline.
# Formato esperado: hipermatriz 3D (canales × muestras × segmentos)
println("=" ^ 80)
println("📊 CARGA DE DATOS DESPUÉS DEL 2ND BASELINE CORRECTION")
println("=" ^ 80)
println()

# Deserializar diccionario con datos y metadatos
dict_2nd_baseline_info = Serialization.deserialize(path_dict_2nd_baseline)

# Extraer datos e información del diccionario
eeg_2nd_baseline_corrected = dict_2nd_baseline_info["eeg_2nd_baseline_corrected"]  # Array 3D: (canales, muestras, segmentos)
channels = dict_2nd_baseline_info["channels"]                                      # Vector{String}: nombres de canales
fs = dict_2nd_baseline_info["fs"]                                                  # Float64: frecuencia de muestreo (Hz)
segment_length_samples = dict_2nd_baseline_info["segment_length_samples"]            # Int: muestras por segmento
n_segments = dict_2nd_baseline_info["n_segments"]                                  # Int: número de segmentos
n_channels = dict_2nd_baseline_info["n_channels"]                                   # Int: número de canales

# Mostrar información de los datos cargados
println("✔ Datos después del 2nd baseline correction cargados desde: $(basename(path_dict_2nd_baseline))")
println("  → Dimensiones: $(size(eeg_2nd_baseline_corrected))")
println("  → Formato: (canales × muestras × segmentos)")
println("  → Número de canales: $n_channels")
println("  → Muestras por segmento: $segment_length_samples")
println("  → Número de segmentos: $n_segments")
println("  → Frecuencia de muestreo: $(fs) Hz")
println("  → Duración de cada segmento: $(round(segment_length_samples/fs, digits=3)) s")
println()

# Para compatibilidad con el resto del código, usar nombres estándar
eeg = eeg_2nd_baseline_corrected
ch_names = channels

# ------------------------------------------------------------------------------------
# 5. CARGA DE INFORMACIÓN DE POSICIONES DE ELECTRODOS
# ------------------------------------------------------------------------------------
# Se carga el archivo TSV con las coordenadas 3D de los electrodos.
# Este archivo contiene información espacial necesaria para construir el operador CSD.
elec = CSV.read(electrode_tsv, DataFrame; delim='\t')

# Filtrar: usar solo los electrodos marcados como tipo "EEG"
# (el archivo TSV puede contener otros tipos como referencias, EOG, etc.)
elec_eeg = elec[elec.type .== "EEG", :]
elec_names = Vector{String}(elec_eeg.name)

# ------------------------------------------------------------------------------------
# 6. VERIFICACIÓN DE CONSISTENCIA (CRÍTICO PARA CSD CORRECTO)
# ------------------------------------------------------------------------------------
# Es fundamental verificar que todos los canales del EEG tengan posiciones definidas
# y que el orden coincida exactamente. Un error aquí resultaría en un CSD incorrecto.

# 6.1) Verificar que todos los canales del EEG tienen posición en el TSV
missing_in_tsv = setdiff(ch_names, elec_names)  # Canales en EEG pero no en TSV
extra_in_tsv   = setdiff(elec_names, ch_names) # Canales en TSV pero no en EEG

if !isempty(missing_in_tsv)
    open(log_path, "a") do io
        logmsg(io, "ERROR: Canales del EEG sin posición en TSV: $(join(missing_in_tsv, ", "))")
    end
    error("Hay canales del EEG que no están en el TSV (faltan posiciones). Revisa nombres/montaje.")
end

# 6.2) Aviso si hay electrodos en TSV que no están en el EEG (no es fatal)
# Esto puede ocurrir si el TSV contiene electrodos adicionales no usados
if !isempty(extra_in_tsv)
    open(log_path, "a") do io
        logmsg(io, "Aviso: Electrodos en TSV no presentes en EEG (se ignoran): $(join(extra_in_tsv, ", "))")
    end
end

# 6.3) Reordenar electrodos para que coincidan EXACTAMENTE con el orden de ch_names
# Esto es crítico: el operador CSD se construye asumiendo que el orden de coordenadas
# coincide con el orden de canales en los datos EEG.
idx = Dict(name => i for (i, name) in enumerate(elec_eeg.name))
elec_ordered = elec_eeg[[idx[n] for n in ch_names], :]

# 6.4) Extraer coordenadas 3D en el orden correcto del EEG
X = collect(elec_ordered.x)  # Coordenada X (normalmente en metros o unidades normalizadas)
Y = collect(elec_ordered.y)  # Coordenada Y
Z = collect(elec_ordered.z)  # Coordenada Z

# Registrar información en el log
open(log_path, "a") do io
    logmsg(io, "fs = $fs Hz")
    logmsg(io, "Canales EEG = $n_channels")
    logmsg(io, "Muestras por segmento = $segment_length_samples")
    logmsg(io, "Número de segmentos = $n_segments")
    logmsg(io, "Electrodos TSV(EEG) = $(nrow(elec_eeg))")
end

# ------------------------------------------------------------------------------------
# 7. REFERENCIA ANTES DE CSD (OPCIONAL)
# ------------------------------------------------------------------------------------
# Nota: El CSD es inherentemente referencia-libre (reference-free), ya que el operador
# Laplaciano elimina la componente constante. Sin embargo, algunos protocolos recomiendan
# aplicar una referencia promedio antes del CSD para estabilidad numérica.
#
# Si los datos ya vienen referenciados y quieres mantenerlo: no hagas nada.
# Si quieres re-referenciar a promedio ANTES de CSD, hazlo aquí y registra en el log.
# TODO: implementar opción de average reference si es necesario.

# ------------------------------------------------------------------------------------
# 8. ALGORITMO CSD: FUNCIONES AUXILIARES
# ------------------------------------------------------------------------------------
# Implementación del método de Perrin (1989) para CSD usando splines esféricos.
# Este método es equivalente al implementado en BrainVision Analyzer.

# ------------------------------------------------------------------------------------
# 8.1. Polinomios de Legendre
# ------------------------------------------------------------------------------------
# Calcula el polinomio de Legendre P_n(x) usando la relación de recurrencia de Bonnet.
# Esta implementación es numéricamente estable y no requiere dependencias externas.
#
# Parámetros:
#   - n: grado del polinomio (n >= 0)
#   - x: valor donde evaluar (debe estar en [-1, 1])
# Retorna:
#   - Float64: valor de P_n(x)
function legendreP(n::Int, x::Float64)
    n == 0 && return 1.0
    n == 1 && return x
    pnm2 = 1.0  # P_0(x) = 1
    pnm1 = x    # P_1(x) = x
    pn = 0.0
    # Relación de recurrencia: P_n(x) = ((2n-1)*x*P_{n-1}(x) - (n-1)*P_{n-2}(x)) / n
    @inbounds for k in 2:n
        pn = ((2k - 1) * x * pnm1 - (k - 1) * pnm2) / k
        pnm2, pnm1 = pnm1, pn
    end
    return pn
end

# ------------------------------------------------------------------------------------
# 8.2. Construcción del Operador CSD
# ------------------------------------------------------------------------------------
# Construye el operador Laplaciano esférico L tal que: EEG_CSD = L * EEG
# El operador se aplica en la dimensión de canales, transformando cada muestra temporal.
#
# Algoritmo (Perrin et al., 1989):
#   1. Normaliza coordenadas de electrodos a esfera unidad
#   2. Construye matrices G y H usando expansión en polinomios de Legendre
#   3. Aplica regularización (lambda) para estabilidad numérica
#   4. Aplica constraint de referencia-libre (elimina componente constante)
#   5. Calcula operador final: L = H * C
#
# Parámetros:
#   - X, Y, Z: vectores con coordenadas 3D de electrodos (mismo orden que canales EEG)
#   - m: orden del spline esférico (default: 4, típico en BrainVision)
#   - degree: grado máximo de polinomios de Legendre (default: 10)
#   - lambda: parámetro de regularización (default: 1e-5)
# Retorna:
#   - Matrix{Float64}: operador CSD de tamaño (n_channels × n_channels)
function csd_operator_perrin(X::AbstractVector, Y::AbstractVector, Z::AbstractVector;
                             m::Int=4, degree::Int=10, lambda::Float64=1e-5)

    n = length(X)
    @assert length(Y) == n && length(Z) == n "Coordenadas X,Y,Z deben tener misma longitud"

    # Normalizar coordenadas a esfera unidad
    # Esto asegura que las coordenadas estén en una esfera de radio 1, necesario para
    # el método de splines esféricos. BrainVision Analyzer suele proporcionar coordenadas
    # ya normalizadas, pero es buena práctica normalizarlas aquí también.
    R = zeros(Float64, n, 3)
    @inbounds for i in 1:n
        r = sqrt(X[i]^2 + Y[i]^2 + Z[i]^2)
        r == 0 && error("Electrodo con radio 0: no es válido para CSD. Revisa coordenadas.")
        R[i,1] = X[i]/r  # Coordenada X normalizada
        R[i,2] = Y[i]/r  # Coordenada Y normalizada
        R[i,3] = Z[i]/r  # Coordenada Z normalizada
    end

    # Inicializar matrices G y H (Perrin)
    # G: matriz de interpolación de splines esféricos
    # H: matriz de derivadas segundas (Laplaciano)
    G = Matrix{Float64}(undef, n, n)
    H = Matrix{Float64}(undef, n, n)

    # Precalcular coeficientes de la expansión en polinomios de Legendre
    # Estos coeficientes dependen solo del orden m y del grado, no de las posiciones
    # específicas, por lo que se calculan una sola vez.
    cg = [ (2k + 1) / ((k^m) * ((k+1)^m)) for k in 1:degree ]
    ch = [ (2k + 1) / ((k^(m-1)) * ((k+1)^(m-1))) for k in 1:degree ]

    # Construir matrices G y H elemento por elemento
    @inbounds for i in 1:n
        # Elementos diagonales: cuando i=j, cos(γ) = 1, por lo que P_k(1) = 1 para todo k
        G[i,i] = sum(cg)
        H[i,i] = sum(ch)
        
        # Elementos fuera de la diagonal: calcular usando polinomios de Legendre
        for j in (i+1):n
            # Coseno del ángulo entre vectores de posición normalizados
            cosγ = R[i,1]*R[j,1] + R[i,2]*R[j,2] + R[i,3]*R[j,3]
            cosγ = clamp(cosγ, -1.0, 1.0)  # Asegurar rango válido para acos

            # Sumar serie de polinomios de Legendre
            sG = 0.0
            sH = 0.0
            for k in 1:degree
                pk = legendreP(k, cosγ)
                sG += cg[k] * pk
                sH += ch[k] * pk
            end
            # Matrices simétricas
            G[i,j] = sG; G[j,i] = sG
            H[i,j] = sH; H[j,i] = sH
        end
    end

    # Regularización: añadir lambda * I a la matriz G
    # Esto estabiliza numéricamente la inversión y es equivalente al parámetro
    # de regularización en BrainVision Analyzer.
    Gλ = G .+ lambda .* Matrix{Float64}(I, n, n)

    # Constraint de referencia-libre (reference-free constraint)
    # Elimina la componente constante del espacio de soluciones, haciendo el CSD
    # independiente de la referencia. Esto es equivalente a proyectar fuera del
    # espacio nulo de G.
    one = ones(Float64, n)  # Vector de unos (componente constante)
    Ginv = inv(Gλ)
    C = Ginv - (Ginv*one*(one'*Ginv)) / (one' * Ginv * one)

    # Operador final Laplaciano esférico
    # L transforma los datos EEG referenciados a CSD (reference-free)
    L = H * C
    return L
end

# ------------------------------------------------------------------------------------
# 9. APLICACIÓN DEL OPERADOR CSD
# ------------------------------------------------------------------------------------
println("=" ^ 80)
println("🧠 APLICANDO CSD (Perrin spherical splines) — estilo BrainVision")
println("=" ^ 80)
println()

# Parámetros del algoritmo CSD (valores por defecto estilo BrainVision Analyzer)
m_csd      = 4      # Orden del spline esférico
deg_csd    = 10     # Grado máximo de polinomios de Legendre
lambda_csd = 1e-5   # Parámetro de regularización

# Construir operador CSD usando las coordenadas de los electrodos
L = csd_operator_perrin(X, Y, Z; m=m_csd, degree=deg_csd, lambda=lambda_csd)

# Aplicar el operador L a cada segmento temporal
# El operador se aplica en la dimensión de canales, manteniendo la forma original
# (canales × muestras × segmentos)
n_ch, n_samp, n_seg = size(eeg)
eeg_csd = similar(eeg)  # Preasignar matriz de salida con misma forma

# Aplicar transformación: eeg_csd[:,:,s] = L * eeg[:,:,s] para cada segmento s
@inbounds for s in 1:n_seg
    @views eeg_csd[:,:,s] .= L * eeg[:,:,s]
end

println("✔ CSD aplicado")
println("  → Dimensiones salida: $(size(eeg_csd))")
println()

# Registrar aplicación de CSD en el log
open(log_path, "a") do io
    logmsg(io, "CSD aplicado (Perrin)")
    logmsg(io, "Parámetros: m=$m_csd, degree=$deg_csd, lambda=$lambda_csd")
end

# ------------------------------------------------------------------------------------
# 10. GUARDADO DE RESULTADOS Y METADATOS
# ------------------------------------------------------------------------------------
# Se guardan tanto los datos transformados como un diccionario completo con metadatos
# que incluye información sobre el método CSD utilizado y sus parámetros.

# Construir diccionario con datos y metadatos
dict_csd = Dict{String, Any}()
dict_csd["eeg_csd"] = eeg_csd                                    # Datos EEG transformados
dict_csd["channels"] = ch_names                                   # Nombres de canales
dict_csd["fs"] = fs                                              # Frecuencia de muestreo
dict_csd["segment_length_samples"] = segment_length_samples      # Muestras por segmento
dict_csd["n_segments"] = n_segments                              # Número de segmentos
dict_csd["n_channels"] = n_channels                             # Número de canales
dict_csd["space"] = "CSD"                                        # Indicador de espacio transformado
dict_csd["csd_method"] = "Perrin spherical spline Laplacian"    # Método utilizado
dict_csd["csd_params"] = Dict("m" => m_csd, "degree" => deg_csd, "lambda" => lambda_csd)  # Parámetros

# Asegurar que el directorio existe
mkpath(dir_csd)

# Serializar y guardar resultados
Serialization.serialize(out_dict_path, dict_csd)  # Diccionario completo con metadatos
Serialization.serialize(out_eeg_path, eeg_csd)   # Datos EEG transformados (solo datos)

println("💾 Guardado:")
println("  → $(basename(out_dict_path))")
println("  → $(basename(out_eeg_path))")

# Registrar guardado en el log
open(log_path, "a") do io
    logmsg(io, "Guardado dict_csd en: $out_dict_path")
    logmsg(io, "Guardado eeg_csd  en: $out_eeg_path")
    logmsg(io, "Fin OK")
end

# ------------------------------------------------------------------------------------
# 11. CONTROL DE CALIDAD (QC): ESTADÍSTICAS Y VISUALIZACIONES
# ------------------------------------------------------------------------------------
# Se generan estadísticas por canal (antes vs después CSD) y figuras de control
# de calidad estilo BrainVision Analyzer para verificar la transformación CSD.

# ------------------------------------------------------------------------------------
# 11.1. Función: Estadísticas por Canal
# ------------------------------------------------------------------------------------
# Calcula estadísticas (media, RMS, varianza) por canal a partir de un array 3D.
# Parámetros:
#   - A: Array 3D de forma (canales × muestras × segmentos)
# Retorna:
#   - μ: Vector{Float64} con media por canal
#   - rms: Vector{Float64} con RMS por canal
#   - v: Vector{Float64} con varianza por canal
function channel_stats(A::AbstractArray{<:Real,3})
    n_ch = size(A, 1)
    μ  = Vector{Float64}(undef, n_ch)
    rms = Vector{Float64}(undef, n_ch)
    v  = Vector{Float64}(undef, n_ch)

    @inbounds for c in 1:n_ch
        x = vec(@view A[c, :, :])
        μ[c]  = mean(x)
        rms[c] = sqrt(mean(abs2, x))
        v[c]  = var(x)
    end
    return μ, rms, v
end

# ------------------------------------------------------------------------------------
# 11.2. Cálculo de Estadísticas QC
# ------------------------------------------------------------------------------------
# Calcular estadísticas antes y después de CSD para comparación
mean_before, rms_before, var_before = channel_stats(eeg)
mean_after,  rms_after,  var_after  = channel_stats(eeg_csd)

# Construir tabla de estadísticas QC
qc_csv = joinpath(dir_tables, "CSD_QC_channel_stats.csv")
qc = DataFrame(
    channel = ch_names,
    mean_before = mean_before,
    mean_after  = mean_after,
    rms_before  = rms_before,
    rms_after   = rms_after,
    var_before  = var_before,
    var_after   = var_after,
    rms_ratio   = rms_after ./ rms_before,
    var_ratio   = var_after ./ var_before
)

CSV.write(qc_csv, qc)
println("📄 QC guardado: $(basename(qc_csv))")

open(log_path, "a") do io
    logmsg(io, "QC CSV: $qc_csv")
end

# ------------------------------------------------------------------------------------
# 11.3. Funciones de Visualización QC
# ------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------
# 11.3.1. Función: Topografía de Dispersión (Scatter)
# ------------------------------------------------------------------------------------
# Genera una figura de topografía simple usando scatter plot con valores coloreados.
# Parámetros:
#   - path_png: ruta donde guardar la figura
#   - title_str: título de la figura
#   - X, Y: coordenadas 2D de los electrodos
#   - values: valores a colorear (p.ej. RMS por canal)
#   - labels: etiquetas de los canales
function topo_scatter_png(path_png::String, title_str::String, X::Vector, Y::Vector, values::Vector, labels::Vector{String})
    f = Figure(size=(600, 600))
    ax = Axis(f[1,1], title=title_str)
    hidespines!(ax); hidedecorations!(ax)

    # Cabeza simple (círculo)
    θ = range(0, 2π, length=300)
    lines!(ax, cos.(θ), sin.(θ))

    # Normalizar X,Y a radio ~1 para visualización (como hace Analyzer)
    r = maximum(sqrt.(X.^2 .+ Y.^2))
    x = X ./ r
    y = Y ./ r

    sc = scatter!(ax, x, y; markersize=18, color=values)
    Colorbar(f[1,2], sc)

    # Etiquetas (opcional, si no quieres ruido visual comenta este bloque)
    for i in eachindex(labels)
        text!(ax, labels[i], position=(x[i], y[i]), align=(:center, :bottom), fontsize=10)
    end

    save(path_png, f)
end

# ------------------------------------------------------------------------------------
# 11.4. Generación de Figuras QC
# ------------------------------------------------------------------------------------

# Definir rutas de salida para figuras QC
fig_rms_bar = joinpath(dir_figures, "CSD_QC_RMS_bar.png")
fig_topo_before = joinpath(dir_figures, "CSD_QC_topo_RMS_before.png")
fig_topo_after  = joinpath(dir_figures, "CSD_QC_topo_RMS_after.png")

# Figura 1: Gráfico de barras RMS (antes vs después)
let
    f = Figure(size=(1400, 450))
    ax = Axis(f[1,1], title="RMS por canal (antes vs después CSD)",
              xlabel="Canal", ylabel="RMS")

    x = 1:length(ch_names)
    barplot!(ax, x .- 0.2, rms_before; width=0.35, label="Antes")
    barplot!(ax, x .+ 0.2, rms_after;  width=0.35, label="Después")

    ax.xticks = (x, ch_names)
    ax.xticklabelrotation = π/2
    axislegend(ax, position=:rt)

    save(fig_rms_bar, f)
end
println("🖼️ Figura guardada: $(basename(fig_rms_bar))")

# Figuras 2-3: Topografías de dispersión RMS (antes y después)
topo_scatter_png(fig_topo_before, "CSD QC — RMS (antes)",  X, Y, rms_before, ch_names)
topo_scatter_png(fig_topo_after,  "CSD QC — RMS (después)", X, Y, rms_after,  ch_names)

println("🖼️ Figuras guardadas: $(basename(fig_topo_before)), $(basename(fig_topo_after))")

open(log_path, "a") do io
    logmsg(io, "QC figuras: $fig_rms_bar")
    logmsg(io, "QC figuras: $fig_topo_before")
    logmsg(io, "QC figuras: $fig_topo_after")
end

# ------------------------------------------------------------------------------------
# 12. MAPAS CSD: TOPOGRAFÍA INTERPOLADA (ESTILO BRAINVISION)
# ------------------------------------------------------------------------------------
# Genera mapas topográficos interpolados del CSD usando interpolación IDW (Inverse
# Distance Weighting). Estos mapas muestran la distribución espacial del CSD en la
# superficie del cuero cabelludo, similar a los mapas generados por BrainVision Analyzer.

# ------------------------------------------------------------------------------------
# 12.1. Función: Interpolación IDW (Inverse Distance Weighting)
# ------------------------------------------------------------------------------------
# Interpola valores en una grilla 2D usando el método de ponderación por distancia inversa.
# Parámetros:
#   - xe, ye: coordenadas 2D de los puntos de datos (electrodos)
#   - ve: valores en los puntos de datos
#   - gridN: tamaño de la grilla (default: 200)
#   - power: exponente de la distancia para ponderación (default: 2.0)
#   - eps: umbral mínimo de distancia para evitar división por cero (default: 1e-12)
# Retorna:
#   - xs, ys: rangos de coordenadas de la grilla
#   - Z: matriz interpolada (valores NaN fuera del círculo de la cabeza)
function idw_grid(xe::Vector{Float64}, ye::Vector{Float64}, ve::Vector{Float64};
                  gridN::Int=200, power::Float64=2.0, eps::Float64=1e-12)

    xs = range(-1.0, 1.0, length=gridN)
    ys = range(-1.0, 1.0, length=gridN)

    Z = Matrix{Float64}(undef, gridN, gridN)

    @inbounds for (iy, y) in enumerate(ys)
        for (ix, x) in enumerate(xs)
            # Máscara: solo dentro del círculo de la cabeza
            if x^2 + y^2 > 1.0
                Z[iy, ix] = NaN
                continue
            end

            num = 0.0
            den = 0.0
            # IDW
            for k in eachindex(ve)
                dx = x - xe[k]
                dy = y - ye[k]
                d2 = dx*dx + dy*dy
                if d2 < eps
                    num = ve[k]
                    den = 1.0
                    break
                end
                w = 1.0 / (d2^(power/2))
                num += w * ve[k]
                den += w
            end
            Z[iy, ix] = num / den
        end
    end

    return xs, ys, Z
end

# ------------------------------------------------------------------------------------
# 12.2. Función: Dibujo de Cabeza Simple
# ------------------------------------------------------------------------------------
# Dibuja un contorno simple de cabeza (círculo + nariz) en el eje proporcionado.
# Parámetros:
#   - ax: Axis de CairoMakie donde dibujar
function draw_head!(ax)
    θ = range(0, 2π, length=400)
    lines!(ax, cos.(θ), sin.(θ))
    # Nariz (triángulo simple arriba)
    lines!(ax, [0.0, 0.08, -0.08, 0.0], [1.0, 1.08, 1.08, 1.0])
end

# ------------------------------------------------------------------------------------
# 12.3. Función: Generación de Mapa CSD
# ------------------------------------------------------------------------------------
# Genera un mapa topográfico interpolado del CSD y lo guarda como PNG.
# El mapa muestra la distribución espacial del CSD usando interpolación IDW y contornos.
# Parámetros:
#   - path_png: ruta donde guardar la figura PNG
#   - values: Vector{Float64} con el valor por canal (p.ej. CSD en un sample o promedio)
#   - X, Y: coordenadas 2D de los electrodos (del TSV) en el orden de ch_names
#   - ch_names: nombres de los canales
#   - title: título de la figura (default: "CSD Map")
#   - gridN: tamaño de la grilla de interpolación (default: 220)
function csd_map_png(path_png::String;
                     values::Vector{Float64},
                     X::Vector, Y::Vector,
                     ch_names::Vector{String},
                     title::String="CSD Map",
                     gridN::Int=220)

    # Proyectar a 2D y normalizar a radio 1 para visualización
    r = maximum(sqrt.(X.^2 .+ Y.^2))
    xe = (X ./ r) |> collect
    ye = (Y ./ r) |> collect
    ve = values |> collect

    xs, ys, Z = idw_grid(xe, ye, ve; gridN=gridN, power=2.0)

    f = Figure(size=(700, 600))
    ax = Axis(f[1,1], title=title)
    hidespines!(ax); hidedecorations!(ax)
    draw_head!(ax)

    hm = heatmap!(ax, xs, ys, Z)          # mapa continuo
    contour!(ax, xs, ys, Z; levels=10)    # líneas de contorno (look Analyzer)

    scatter!(ax, xe, ye; markersize=10)   # electrodos
    # Etiquetas opcionales (si lo prefieres sin etiquetas, comenta este bucle)
    for i in eachindex(ch_names)
        text!(ax, ch_names[i], position=(xe[i], ye[i]), align=(:center, :bottom), fontsize=9)
    end

    Colorbar(f[1,2], hm)
    save(path_png, f)
end

# ------------------------------------------------------------------------------------
# 12.4. Generación de Mapas CSD
# ------------------------------------------------------------------------------------
# Se generan dos tipos de mapas CSD:
#   A) Mapa instantáneo: CSD en un momento temporal específico (una muestra)
#   B) Mapa promedio: CSD promedio temporal en un segmento completo

# A) CSD map instantáneo: un segmento y una muestra concreta
# (p.ej. sample 250 ~ mitad del segundo)
seg_idx = 1
sample_idx = 250
values_instant = collect(vec(@view eeg_csd[:, sample_idx, seg_idx]))  # Convertir a Vector{Float64}

path_map_instant = joinpath(dir_figures, "CSD_map_instant_seg$(seg_idx)_samp$(sample_idx).png")
csd_map_png(path_map_instant;
    values=values_instant, X=X, Y=Y, ch_names=ch_names,
    title="CSD map (instant) — seg=$seg_idx, sample=$sample_idx"
)
println("🖼️ CSD map instantáneo guardado: $(basename(path_map_instant))")

# B) CSD map promedio temporal: promedio en el segmento (1s)
values_mean = collect(vec(mean(@view(eeg_csd[:, :, seg_idx]); dims=2)))  # Convertir a Vector{Float64}
path_map_mean = joinpath(dir_figures, "CSD_map_mean_seg$(seg_idx).png")
csd_map_png(path_map_mean;
    values=values_mean, X=X, Y=Y, ch_names=ch_names,
    title="CSD map (mean 1s) — seg=$seg_idx"
)
println("🖼️ CSD map promedio guardado: $(basename(path_map_mean))")

open(log_path, "a") do io
    logmsg(io, "CSD map instant: $path_map_instant")
    logmsg(io, "CSD map mean:    $path_map_mean")
end
