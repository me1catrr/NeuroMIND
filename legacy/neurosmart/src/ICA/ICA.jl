#!/usr/bin/env julia
# -*- coding: utf-8 -*-
# src/ICA/ICA.jl
#
# DESCOMPOSICIÓN ICA (Independent Component Analysis) usando FastICA
# ===================================================================
# Esta rutina implementa el algoritmo FastICA simétrico para descomponer
# señales EEG multicanal en componentes estadísticamente independientes.
#
# PROCESO:
#   1. Carga datos filtrados (post-notch, bandreject, highpass, lowpass)
#   2. Organiza datos en matriz X (canales × muestras)
#   3. Aplica FastICA simétrico:
#      - Centrado de datos (media cero por canal)
#      - Blanqueo mediante PCA (whitening)
#      - Extracción de componentes independientes mediante optimización
#        de no-gaussianidad (función tanh)
#   4. Obtiene:
#      - S: componentes independientes (n_comp × n_samples)
#      - W_total: matriz de desmezcla total (n_comp × n_channels)
#      - A: matriz de mezcla aproximada (n_channels × n_comp)
#   5. Guarda resultados en dict_EEG_ICA.bin
#
# NOTA: El número de componentes ICA es igual al número de canales.

using Serialization
using LinearAlgebra
using Statistics
using Random

# Si se ejecuta este script directamente (fuera del módulo EEG_Julia),
# cargamos utilidades de rutas para disponer de `stage_dir`.
if !@isdefined(stage_dir)
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end

# ------------------------------------------------------------------------------------
# 1. Carga de datos filtrados (tras notch y bandreject, hp y lp filtrados)
# ------------------------------------------------------------------------------------
# Se cargan los datos EEG que ya han sido filtrados en pasos anteriores:
# - Filtrado notch (eliminación de ruido de línea eléctrica)
# - Filtrado bandreject (eliminación de bandas específicas)
# - Filtrado highpass (eliminación de deriva de baja frecuencia)
# - Filtrado lowpass (eliminación de frecuencias altas)
#
# Estos datos filtrados son la entrada para el análisis ICA.

println("========================================")
println("Paso 1: Carga de datos filtrados    ")
println("========================================")

# Directorio base para datos ICA (donde se guardarán los resultados)
dir_ica = stage_dir(:ICA)
path_dict_ica = joinpath(dir_ica, "dict_EEG_ICA.bin")

# Directorio de datos filtrados (entrada para ICA)
dir_filtering     = stage_dir(:filtering)
path_dict_lowpass = joinpath(dir_filtering, "dict_EEG_Lowpass.bin")

# Verificamos que el resultado del filtrado esté disponible
if !isfile(path_dict_lowpass)
    error("No se encontró $(abspath(path_dict_lowpass)). Ejecuta antes src/Preprocessing/filtering.jl para generar dict_EEG_Lowpass.bin (y src/Setup/IO.jl si parte desde cero).")
end

# Cargar diccionario con señales por canal
dict_EEG = Serialization.deserialize(path_dict_lowpass)

println("✓ Archivo: $(basename(path_dict_lowpass))")
println()
display(dict_EEG)
println()

# Orden estable y reproducible de canales (orden alfabético)
# Esto asegura que la matriz X tenga un orden consistente
channels   = sort(collect(keys(dict_EEG)))
n_channels = length(channels)
first_ch   = first(channels)
n_samples  = length(dict_EEG[first_ch])

println("✓ Dimensiones: $n_channels canales × $n_samples muestras")
println("✓ Canales (ordenados): $(join(channels, ", "))")
println()

# Construir matriz X: canales × muestras
# Cada fila corresponde a un canal, cada columna a una muestra temporal
X = Array{Float64}(undef, n_channels, n_samples)

# Rellenar la matriz X con los datos de cada canal
# Se verifica que todos los canales tengan la misma longitud
for (i, ch) in enumerate(channels)
    sig = dict_EEG[ch]
    @assert length(sig) == n_samples "Canal $ch con longitud distinta"
    X[i, :] = sig
end

# Número de componentes ICA a extraer
# En ICA, típicamente se extraen tantos componentes como canales hay
k = n_channels  

println("ICA: $n_channels canales, $n_samples muestras, $k componentes")
println()

# ------------------------------------------------------------------------------------
# 2. FUNCIONES AUXILIARES: BLANQUEO Y ORTONORMALIZACIÓN SIMÉTRICA
# ------------------------------------------------------------------------------------
# Estas funciones son necesarias para el preprocesamiento de datos antes de aplicar
# FastICA. El blanqueo reduce la dimensionalidad y elimina correlaciones de segundo
# orden, mientras que la ortonormalización simétrica mantiene las filas de W
# ortogonales durante la optimización.

"""
    whiten_pca(X, k)

Blanqueo (PCA whitening / probabilistic sphering).
Transforma los datos X de manera que la matriz de covarianza de los datos
blanqueados sea aproximadamente la identidad (cov(Z) ≈ I).

PROCESO:
  1. Calcula la matriz de covarianza de X
  2. Realiza descomposición en autovalores/autovectores (PCA)
  3. Selecciona las k primeras componentes principales
  4. Aplica transformación de blanqueo: Z = Λ^(-1/2) * E' * X
     donde Λ contiene los autovalores y E los autovectores

Entrada:
  X :: Matrix{Float64}    # (n_channels × n_samples)
  k :: Int                # nº de componentes a mantener (≤ n_channels)

Salida:
  Z       :: Matrix{Float64}  # datos blanqueados (k × n_samples), cov(Z) ≈ I
  V_whit  :: Matrix{Float64}  # matriz de blanqueo (k × n_channels)
"""
function whiten_pca(X::AbstractMatrix{<:Real}, k::Int)
    n_channels, n_samples = size(X)

    @assert 1 ≤ k ≤ n_channels "k debe estar entre 1 y n_channels"

    # Matriz de covarianza: (n_channels × n_channels)
    # cov = (1/N) * X * X'
    # Mide las correlaciones entre canales
    CovX = (X * X') / n_samples

    # Descomposición en autovalores/autovectores (PCA)
    # Los autovectores (E) son las direcciones de máxima varianza
    # Los autovalores (λ) son las varianzas en esas direcciones
    F = eigen(Symmetric(CovX))  # aseguramos simetría numérica
    eigvals  = F.values         # vector (n_channels)
    eigvecs  = F.vectors        # matriz (n_channels × n_channels), columnas = autovectores

    # Ordenar autovalores y autovectores de mayor a menor
    # Las componentes principales más importantes primero
    idx = sortperm(eigvals; rev = true)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Nos quedamos con las k primeras componentes principales
    # Esto reduce la dimensionalidad de n_channels a k
    λ_k = eigvals[1:k]
    E_k = eigvecs[:, 1:k]   # (n_channels × k)

    # Matriz diagonal de λ^(-1/2)
    # Normaliza las varianzas a 1 (blanqueo)
    Λ_inv_sqrt = Diagonal(1 ./ sqrt.(λ_k))   # (k × k)

    # Matriz de blanqueo: V = Λ^(-1/2) * E'
    # Tamaño: (k × n_channels)
    # Esta matriz transforma X en datos blanqueados Z
    V_whit = Λ_inv_sqrt * E_k'

    # Datos blanqueados: Z = V * X  (k × n_samples)
    # Z tiene covarianza aproximadamente igual a la identidad
    Z = V_whit * X

    return Z, V_whit
end

"""
    sym_decorrelation(W)

Ortonormalización simétrica de W (symmetric decorrelation),
tal y como se usa en FastICA.

Esta función asegura que las filas de W sean ortonormales (W*W' ≈ I),
lo cual es necesario para que FastICA extraiga componentes independientes
y no correlacionadas.

PROCESO:
  1. Calcula C = W * W' (matriz de productos internos entre filas)
  2. Encuentra C^(-1/2) mediante descomposición en autovalores
  3. Aplica transformación: W_ortho = C^(-1/2) * W

Entrada:
  W :: Matrix{Float64}  # (k × k)

Salida:
  W_ortho :: Matrix{Float64}  # filas ortonormales, W*W' ≈ I
"""
function sym_decorrelation(W::AbstractMatrix{<:Real})
    # C = W * W'
    # Matriz de productos internos entre las filas de W
    C = W * W'

    # Descomposición propia de C
    # Necesaria para calcular C^(-1/2)
    F = eigen(Symmetric(C))
    E = F.vectors          # autovectores (k × k)
    d = F.values           # autovalores  (k)

    # C^(-1/2) = E * diag(1/sqrt(d)) * E'
    # Raíz cuadrada inversa de C mediante descomposición espectral
    C_inv_sqrt = E * Diagonal(1 ./ sqrt.(d)) * E'

    # W_ortho = C^(-1/2) * W
    # Transformación que hace que las filas de W sean ortonormales
    W_ortho = C_inv_sqrt * W
    return W_ortho
end

# ------------------------------------------------------------------------------------
# 3. IMPLEMENTACIÓN DE FASTICA SIMÉTRICO (SUPER-GAUSIANO, tanh)
# ------------------------------------------------------------------------------------
# FastICA es un algoritmo iterativo que encuentra componentes estadísticamente
# independientes maximizando la no-gaussianidad mediante el método de Newton.
#
# PRINCIPIO: El Teorema del Límite Central establece que la suma de variables
# aleatorias independientes tiende a ser gaussiana. Por tanto, para encontrar
# componentes independientes, buscamos señales que sean lo más no-gaussianas posible.
#
# ALGORITMO:
#   1. Centrado: eliminar la media de cada canal
#   2. Blanqueo: reducir dimensionalidad y eliminar correlaciones (PCA whitening)
#   3. Inicialización: matriz W aleatoria y ortonormalizada
#   4. Iteración (hasta convergencia o max_iter):
#      a) Proyectar: Y = W * Z
#      b) Aplicar función no lineal: g(Y) = tanh(a*Y)
#      c) Actualizar W según gradiente de no-gaussianidad
#      d) Ortonormalizar W (symmetric decorrelation)
#      e) Verificar convergencia
#   5. Construir salidas: S, W_total, A

"""
    fastica(X; n_comp, max_iter, tol, a, seed)

Aplica FastICA simétrico a la matriz X para extraer componentes independientes.

El algoritmo busca maximizar la no-gaussianidad de las componentes extraídas
usando la función tanh como aproximación de la negentropía (medida de no-gaussianidad).

Entrada:
  X         :: Matrix{Float64}     # (n_channels × n_samples), ya filtrado
  n_comp    :: Int                 # nº de componentes ICA a extraer
  max_iter  :: Int                 # máximo nº de iteraciones
  tol       :: Float64             # tolerancia de convergencia
  a         :: Float64             # parámetro de tanh (g(u) = tanh(a*u))
  seed      :: Int                 # semilla para reproducibilidad

Salida:
  S        :: Matrix{Float64}  # (n_comp × n_samples), señales independientes
  W_total  :: Matrix{Float64}  # (n_comp × n_channels), matriz de desmezcla total
  A        :: Matrix{Float64}  # (n_channels × n_comp), matriz de mezcla
"""
function fastica(
    X::AbstractMatrix{<:Real};
    n_comp::Int,
    max_iter::Int = 512,
    tol::Float64 = 1e-7,
    a::Float64 = 1.0,
    seed::Int = 1234,
)

    # --------------------------
    # 0) Preparativos
    # --------------------------
    n_channels, n_samples = size(X)
    @assert 1 ≤ n_comp ≤ n_channels "n_comp debe estar entre 1 y n_channels"

    # --------------------------
    # 1) Centrado por canal
    #    (media cero en cada fila)
    # --------------------------
    # ICA requiere que los datos tengan media cero. Se centra cada canal
    # restando su media, lo que elimina el componente DC (corriente continua).
    Xc = copy(Matrix{Float64}(X))  # copiamos a matriz densa de Float64
    for i in 1:n_channels
        μ = mean(Xc[i, :])
        Xc[i, :] .-= μ
    end

    # --------------------------
    # 2) Blanqueo (PCA whitening)
    #    Obtendremos:
    #    - Z: datos blanqueados (n_comp × n_samples)
    #    - V_whit: matriz de blanqueo (n_comp × n_channels)
    # --------------------------
    # El blanqueo reduce la dimensionalidad y elimina correlaciones de segundo orden.
    # Esto simplifica el problema: en lugar de buscar componentes independientes
    # en el espacio original, las buscamos en el espacio blanqueado donde solo
    # necesitamos una matriz ortogonal W (más fácil de optimizar).
    Z, V_whit = whiten_pca(Xc, n_comp)

    # --------------------------
    # 3) Inicialización de W (k × k)
    #    W trabaja en el espacio blanqueado
    # --------------------------
    # W es la matriz de desmezcla en el espacio blanqueado.
    # Se inicializa aleatoriamente y luego se ortonormaliza para que sus filas
    # sean ortogonales (esto asegura que las componentes extraídas sean ortogonales).
    Random.seed!(seed)    # para reproducibilidad
    W = randn(n_comp, n_comp)  # inicialización aleatoria

    # Ortonormalizar W inicialmente
    # Esto asegura que las filas de W sean ortogonales desde el inicio
    W = sym_decorrelation(W)

    # --------------------------
    # 4) Bucle principal de FastICA
    # --------------------------
    # En cada iteración, se actualiza W para maximizar la no-gaussianidad
    # de las componentes Y = W * Z. La función tanh aproxima la negentropía,
    # que es una medida de no-gaussianidad.
    for iter in 1:max_iter
        # Proyección actual: Y = W * Z
        # Y: (n_comp × n_samples)
        # Y contiene las componentes actuales en el espacio blanqueado
        Y = W * Z

        # No linealidad supergaussiana:
        #   g(u)  = tanh(a*u)  - función de contraste (aproxima negentropía)
        #   g'(u) = a * (1 - tanh(a*u)^2)  - derivada de g
        # La función tanh es adecuada para señales supergaussianas (con picos),
        # que es típico en artefactos como parpadeos o actividad muscular.
        GY = tanh.(a .* Y)                      # g(Y)
        Gp = a .* (1 .- tanh.(a .* Y).^2)       # g'(Y)

        # Medias de g'(Y) sobre las muestras (por fila)
        # Estas medias se usan en la actualización de W
        mean_Gp = mean(Gp, dims = 2)           # (n_comp × 1)
        D = Diagonal(vec(mean_Gp))             # diag de medias

        # Actualización tipo FastICA simétrico:
        #   W_new = (1/N) * GY * Z' - D * W
        # Esta fórmula proviene del método de Newton aplicado a la maximización
        # de la negentropía. El primer término es el gradiente, el segundo
        # es un término de normalización.
        W_new = (GY * Z') / n_samples - D * W

        # Ortonormalización simétrica
        # Asegura que las filas de W_new sean ortogonales, lo cual garantiza
        # que las componentes extraídas sean ortogonales (y por tanto independientes
        # en el espacio blanqueado).
        W_new = sym_decorrelation(W_new)

        # ----------------------
        # 4.1) Criterio de convergencia
        #      max_i |1 - |w_i_new ⋅ w_i_old|| < tol
        # ----------------------
        # Se verifica si las filas de W han cambiado significativamente.
        # Si el producto interno entre w_i_new y w_i_old es cercano a 1 (o -1),
        # significa que la dirección no ha cambiado mucho → convergencia.
        M = W_new * W'                         # productos internos entre filas
        diagM = diag(M)                         # productos internos de cada fila consigo misma
        max_diff = maximum(abs.(1 .- abs.(diagM)))  # diferencia máxima

        println("Iteración $iter  ->  max_diff = $max_diff")

        # Si convergió, actualizamos W y salimos
        if max_diff < tol
            println("✔ Convergencia alcanzada en $iter iteraciones (tol = $tol)")
            W = W_new
            break
        end

        # Preparar siguiente iteración
        W = W_new

        # Si llegamos al final sin romper, informamos
        if iter == max_iter
            println("⚠ Se alcanzó max_iter = $max_iter sin converger (max_diff = $max_diff)")
        end
    end

    # --------------------------
    # 5) Construir salidas
    # --------------------------
    # Una vez convergido el algoritmo, construimos las matrices finales:
    # - S: componentes independientes (señales fuente)
    # - W_total: matriz de desmezcla total (del espacio original al de componentes)
    # - A: matriz de mezcla (del espacio de componentes al original)

    # Componentes independientes en el espacio blanqueado:
    # S = W * Z   (n_comp × n_samples)
    # Cada fila de S es una componente independiente (IC)
    S = W * Z

    # Matriz de desmezcla total:
    #   S = W * Z = W * V_whit * Xc = W_total * Xc
    #   => W_total = W * V_whit
    # W_total transforma directamente los datos centrados Xc a componentes S
    # Tamaño: (n_comp × n_channels)
    W_total = W * V_whit

    # Matriz de mezcla aproximada:
    #   Xc ≈ A * S
    # A es la inversa de W_total y representa cómo se mezclan las componentes
    # para reconstruir los datos originales. Las columnas de A son los mapas
    # topográficos de cada componente (pesos espaciales).
    # Tamaño: (n_channels × n_comp)
    A = inv(W_total)

    return S, W_total, A
end

# ------------------------------------------------------------------------------------
# 4. EJECUTAR ICA SOBRE X
# ------------------------------------------------------------------------------------
# Se ejecuta el algoritmo FastICA con los parámetros especificados.
# Los parámetros pueden ajustarse según las características de los datos:
# - max_iter: número máximo de iteraciones (512 suele ser suficiente)
# - tol: tolerancia de convergencia (1e-7 es estricto, asegura buena convergencia)
# - a: parámetro de la función tanh (1.0 es el valor estándar)
# - seed: semilla para reproducibilidad

println("🚀 Ejecutando FastICA...")

n_comp    = k          # nº de componentes ICA = nº canales
max_iter  = 512        # máximo número de iteraciones permitidas
tol       = 1e-7       # criterio de convergencia (cuando max_diff < tol, se detiene)
a         = 1.0        # parámetro de tanh(a*u), valor típico en FastICA
seed      = 1234       # semilla para reproducibilidad (inicialización aleatoria)
println("n_comp: ", n_comp)
println("max_iter: ", max_iter)
println("tol: ", tol)
println("a: ", a)
println("seed: ", seed)

# Ejecutar FastICA
# Esta función realizará todo el proceso: centrado, blanqueo, optimización iterativa
S, W_total, A = fastica(
    X;
    n_comp   = n_comp,
    max_iter = max_iter,
    tol      = tol,
    a        = a,
    seed     = seed,
)

println()
println("✅ ICA completado")
println("  - S (componentes): size = ", size(S))
println("  - W_total (desmezcla): size = ", size(W_total))
println("  - A (mezcla): size = ", size(A))
println()

# ------------------------------------------------------------------------------------
# 5. GUARDAR RESULTADOS EN UN DICCIONARIO
# ------------------------------------------------------------------------------------
# Se guardan todos los resultados del análisis ICA en un diccionario serializado.
# Este diccionario será usado por la rutina de limpieza (ICA_cleaning.jl) para:
# - Visualizar las componentes (mapas topográficos, time courses, espectros)
# - Evaluar y etiquetar componentes artefactuales
# - Reconstruir los datos limpios eliminando componentes artefactuales

# Creamos el diccionario de resultados ICA
dict_EEG_ICA = Dict(
    "S"         => S,          # componentes independientes (n_comp × n_samples)
                                # Cada fila es una componente IC (señal fuente)
    "W_total"   => W_total,    # desmezcla total (n_comp × n_channels)
                                # Transforma datos centrados Xc a componentes S
    "A"         => A,          # mezcla aproximada (n_channels × n_comp)
                                # Columnas de A son mapas topográficos de cada IC
    "channels"  => channels,   # orden de canales (necesario para reconstrucción)
    "max_iter"  => max_iter,   # parámetros usados (para referencia)
    "tol"       => tol,         # parámetros usados (para referencia)
)

# Nos aseguramos de que el directorio exista
isdir(dir_ica) || mkpath(dir_ica)

# Serializar y guardar el diccionario
# El formato binario es eficiente para matrices grandes
Serialization.serialize(path_dict_ica, dict_EEG_ICA)

println("💾 Resultados ICA guardados en: $path_dict_ica")
