### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 2ff7ee54-7bec-4e21-b804-e9fcc2f86a11
begin
	using PlutoUI
	PlutoUI.TableOfContents(title = "Contenido")
end

# ╔═╡ 49ca7fca-11dc-4f81-b74f-b54e8b252e10
begin
	include(joinpath(@__DIR__, "..", "_template_base.jl"))
	using .PlutoTemplateBase
	using CSV, DataFrames, Serialization, Statistics, StatsBase, DSP, Plots, Dates, InlineStrings, LinearAlgebra, Random, FFTW
end

# ╔═╡ 789269a7-1c3a-41a2-811c-c5b8e08ee1f3
begin
# Si se ejecuta este script directamente (fuera del módulo EEG_Julia),
# cargamos utilidades de rutas para disponer de `stage_dir`.
if !@isdefined(stage_dir)
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end
end

# ╔═╡ 7a5d6ddb-33a3-422a-88a4-f76fc8f0c620
md"""
**PAQUETES CARGADOS**
"""

# ╔═╡ 8a0cde3b-88ab-4b56-bef6-8b9a394b1c2e
notebook_intro("Independent Component Analysis (ICA)")

# ╔═╡ 177911fa-6acc-4cd5-b862-81fa0fef1daa
begin
# -----------------------------------------------------------------------------------
# 1. Carga de datos filtrados (tras notch y bandreject, hp y lp filtrados)
# -----------------------------------------------------------------------------------
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
end

# ╔═╡ c5f19470-bda8-4420-b7f4-c3a1c94b4753
begin
# Orden estable y reproducible de canales (orden alfabético)
# Esto asegura que la matriz X tenga un orden consistente
channels   = sort(collect(keys(dict_EEG)))
n_channels = length(channels)
first_ch   = first(channels)
n_samples  = length(dict_EEG[first_ch])

println("✓ Dimensiones: $n_channels canales × $n_samples muestras")
println("✓ Canales (ordenados): $(join(channels, ", "))")
println()
end

# ╔═╡ 01b14950-89b8-4df1-981f-44a13f6026cb
md"""
# Algoritmo ICA

Después del preprocesamiento (**Fase 1**), el **EEG filtrado** se descompone en **componentes estadísticamente independientes** para separar la actividad neuronal de artefactos típicos (oculares, musculares y ruido de línea).

**ICA**, Análisis de Componentes Independientes, asume que la señal multicanal observada es una **mezcla lineal de fuentes latentes** que son mutuamente independientes y no gaussianas. Bajo este modelo, las observaciones pueden escribirse como ``X = A S`` donde ``X \in \mathbb{R}^{C \times N}`` es la matriz de C **canales** y N **muestras temporales**,``A \in \mathbb{R}^{C \times C}`` es la **matriz de mezcla** (pesos espaciales), y ``S \in \mathbb{R}^{C \times N}`` contiene las **señales fuente independientes**.

El objetivo es estimar una **matriz de desmezcla** ``W`` tal que ``S = W X`` recupere las fuentes originales (salvo permutación y signo).

Debido a que la suma de variables aleatorias independientes tiende a ser gaussiana (*teorema central del límite*), **ICA** identifica las fuentes **maximizando la no gaussianidad** de las señales proyectadas.
"""

# ╔═╡ 09a79ac1-4430-4729-bca5-5425a2ba330e
begin
# Construir matriz X: canales × muestras
# Cada fila corresponde a un canal, cada columna a una muestra temporal
X = Array{Float64}(undef, n_channels, n_samples)
end

# ╔═╡ 7523a05e-bb9e-48bc-9720-9a6178941273
begin
# Rellenar la matriz X con los datos de cada canal
# Se verifica que todos los canales tengan la misma longitud
for (i, ch) in enumerate(channels)
    sig = dict_EEG[ch]
    @assert length(sig) == n_samples "Canal $ch con longitud distinta"
    X[i, :] = sig
end
end

# ╔═╡ 580b258a-cde8-4c8a-bdd4-5c71f411cc6b
begin
# Número de componentes ICA a extraer
# En ICA, típicamente se extraen tantos componentes como canales hay
k = n_channels  

println("ICA: $n_channels canales, $n_samples muestras, $k componentes")
println()
end

# ╔═╡ d9936d83-6fae-47af-9890-4c0594e2c782
md"
---

## Algoritmo Fast-ICA simétrico

La implementación utiliza el algoritmo **FastICA simétrico**. Primero, los datos se **centran** (media cero por canal) y posteriormente se **blanquean** mediante PCA para que los datos blanqueados ``Z`` tengan **covarianza identidad**.

El número de componentes se fija igual al número de canales ``k = C``. Esto conserva todas las dimensiones espaciales y permite que el paso posterior de limpieza asigne varianza a componentes neuronales o artefactuales.

La proporción de varianza explicada por cada componente puede calcularse a partir de los autovalores de **PCA** y utilizarse para control de calidad (por ejemplo, detectar sujetos con varianza anormalmente baja en los primeros componentes).

En el espacio blanqueado, la matriz de desmezcla es **ortogonal**. **FastICA** actualiza iterativamente una matriz ortogonal ``W`` maximizando una función de contraste que aproxima la **negentropía**. En esta implementación se utiliza la no linealidad ``g(u) = \tanh(a u)`` con ``a = 1`` adecuada para fuentes supergaussianas (por ejemplo, artefactos oculares o musculares).

En cada iteración se calcula ``Y = W Z`` se aplica la función \(g\) y su derivada, y se actualiza \(W\) mediante un paso tipo Newton seguido de **decorrelación simétrica**, de modo que

```math
W W^{\top} \approx I
```

La convergencia se verifica mediante

```math
\max_i \left|1 - \left| w_i^{new} \cdot w_i^{old} \right|\right| < \tau
```

por ejemplo con ``\tau = 10^{-7}``.

Finalmente, las salidas del algoritmo son:

- los **componentes independientes**

```math
S
```

- la **matriz de desmezcla total**

```math
W_{\mathrm{total}} = W V_{\mathrm{whit}}
```

- la **matriz de mezcla**

```math
A = W_{\mathrm{total}}^{-1}
```

cuyas columnas representan los **mapas topográficos** de cada componente.

El **Algoritmo FastICA** resume el procedimiento **FastICA** simétrico utilizado en `src/ICA.jl`.
"

# ╔═╡ 9fd4518f-a809-4c1e-8a6f-ca26a557094e
md"
---

## Funciones auxiliares 

Estas funciones son necesarias para el preprocesamiento de datos antes de aplicar
**Fast-ICA**:

- [whiten_pca] El **blanqueo** reduce la dimensionalidad y elimina correlaciones de segundo orden
- [sym_decorrelation] La **ortonormalización simétrica** mantiene las filas de W ortogonales durante la optimización
"

# ╔═╡ 98207097-b070-4b07-8ff5-1d29e64e5142
begin
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
end

# ╔═╡ ad20b61f-27c5-46bd-b7c4-2d5125a523b7
begin
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
end

# ╔═╡ 3b9148cd-9027-42e5-9446-b1b1d1c39fdf
md"
---
## Pseudo-código Fast-ICA simétrico

**Fast-ICA** (súper gaussiano, tanh) es un algoritmo iterativo que encuentra componentes estadísticamente
independientes maximizando la no-gaussianidad mediante el método de Newton.

PRINCIPIO: El **Teorema del Límite Central** establece que la suma de variables
aleatorias independientes tiende a ser gaussiana. Por tanto, para encontrar
componentes independientes, buscamos señales que sean lo más no-gaussianas posible.

ALGORITMO:
   1. **Centrado**: eliminar la media de cada canal
   2. **Blanqueo**: reducir dimensionalidad y eliminar correlaciones (PCA whitening)
   3. **Inicialización**: matriz W aleatoria y ortonormalizada
   4. **Iteración** (hasta convergencia o max_iter):
      - Proyectar: Y = W * Z
      - Aplicar función no lineal: g(Y) = tanh(a*Y)
      - Actualizar W según gradiente de no-gaussianidad
      - Ortonormalizar W (symmetric decorrelation)
      - Verificar convergencia
   5. **Construir salidas**: S, W_total, A
"

# ╔═╡ 04fd2eda-02c4-49b5-ac82-c747a300bdf4
md"""
!!! info "Pipeline Fase 2a: Descomposición ICA"

    **Objetivo.**
    - Descomponer el **EEG filtrado** en componentes estadísticamente independientes utilizando **FastICA simétrico**.
    - Generar las matrices ``S``, ``W_{\text{total}}`` y ``A`` para su posterior evaluación y limpieza de artefactos.

    **Entrada.**
    - **EEG filtrado** procedente de la Fase 1: `data/filtering/dict_EEG_Lowpass.bin`
      (diccionario indexado por canal, ``C \times N`` muestras).
    - **Frecuencia de muestreo** ``f_s = 500`` Hz.
    - **Canales** ``C = 32``.
    - **Número de componentes** ``k = C``.

    **Salida:** (`data/ICA/dict_EEG_ICA.bin`)
    - ``S`` (ICs × muestras),
    - ``W_{\text{total}}``,
    - ``A``, lista de canales,
    - `max_iter`, `tol`.

    **Pasos de procesamiento (alineados con `src/ICA.jl`).**

    1. **Carga del EEG filtrado.**
       - Cargar el **EEG** desde `data/filtering/`.
       - Construir la matriz ``X \in \mathbb{R}^{C \times N}`` (canales × muestras) en un orden fijo de canales (por ejemplo, orden alfabético).

    2. **Centrado de la señal.**
       Restar la media de cada canal en ``X``.

    3. **Blanqueamiento (*whitening*).**
       - Calcular la matriz de covarianza de ``X``.
       - Realizar descomposición en autovalores mediante **PCA** y construir la señal blanqueada ``Z`` con ``k = C`` componentes.

    4. **Algoritmo FastICA.**
       - Inicializar ``W``.
       - Aplicar iteraciones hasta convergencia:
           - Proyectar ``Y = WZ``.
           - Aplicar la función de contraste `tanh`.
           - Actualizar ``W`` con decorrelación simétrica.

       El proceso se repite hasta que

       ```math
       \max_i \left|1 - |M_{ii}|\right| < \text{tol} = 10^{-7}
       ```

       o hasta alcanzar el número máximo de iteraciones ``\mathrm{max_iter} = 512``.

    5. **Resultados.**
       Calcular:

       ```math
       S = WZ
       ```

       ```math
       W_{\text{total}} = W \cdot V_{\text{whit}}
       ```

       ```math
       A = W_{\text{total}}^{-1}
       ```

       Serializar los resultados en `data/ICA/dict_EEG_ICA.bin`.
"""

# ╔═╡ 62db513e-efef-4f65-898e-4d1283e3c55a
begin
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
end

# ╔═╡ 818a6675-189e-44fe-abb5-94a1dab065ea
md"
## Ejecutar ICA sobre X

Se ejecuta el algoritmo **Fast-ICA** con los parámetros especificados. Los parámetros pueden ajustarse según las características de los datos:
- **max_iter**: número máximo de iteraciones (512 suele ser suficiente)
- **tol**: tolerancia de convergencia (1e-7 es estricto, asegura buena convergencia)
- **a**: parámetro de la función tanh (1.0 es el valor estándar)
- **seed**: semilla para reproducibilidad
"

# ╔═╡ 7e8b0401-0b43-43e4-a3f2-6f51c471da46
begin
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
end

# ╔═╡ c289b2e6-294e-48ea-a8c7-75f1d5448d13
begin
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
end

# ╔═╡ a036b633-2812-444e-9f72-7117faaac04e
md"
## Diccionario de datos ICA

Se guardan todos los resultados del análisis **ICA** en un diccionario serializado.
Este diccionario será usado por la rutina de limpieza (**ICA_cleaning.jl**) para:
- Visualizar las componentes (mapas topográficos, time courses, espectros)
- Evaluar y etiquetar componentes artefactuales
- Reconstruir los datos limpios eliminando componentes artefactuales
"

# ╔═╡ c5ddf477-9edd-4ce7-aefa-d9a17b0c0fe0
begin
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
end

# ╔═╡ eed453df-7646-41a7-af88-c5beb48902cd
md"
---
# ICA cleaning

Esta rutina evalúa, identifica y elimina componentes artefactuales de los datos
descompuestos mediante **ICA**. Utiliza un sistema automático de evaluación basado
en features espaciales, espectrales y estadísticas para clasificar componentes
como artefactos o señales neuronales.

**PROCESO COMPLETO**:

   1. Preliminares: cálculo de métricas de referencia (datos pre-ICA)
   2. Carga de resultados ICA (matrices S, A, W_total)
   3. Carga de posiciones de electrodos (para mapas topográficos)
   4. Generación automática de visualizaciones para todos los ICs
   5. Evaluación automática mediante sistema de features y scores
   6. Selección y eliminación de componentes artefactuales
   7. Reconstrucción de datos limpios
   8. Visualización comparativa antes/después
"

# ╔═╡ 76fb19ba-33ed-4214-992e-a3d1f0a1983f
md"

La etapa de limpieza (`src/ICA_cleaning.jl`) toma los resultados de la anterior y evalúa cada componente independiente utilizando características espaciales, espectrales y estadísticas.

Las **características espaciales** incluyen `frontal_ratio` y `temporal_ratio`, calculadas a partir de las columnas de la matriz de mezcla **A**.

Las **características espectrales** incluyen:

- `blink_ratio` (0.5–4 Hz frente a 4–40 Hz)
- `emg_ratio` (30–80 Hz frente a 1–30 Hz)
- `line_ratio` (48–52 Hz frente al espectro total)

Las **características estadísticas** incluyen:

- **kurtosis**
- `extreme_frac`, definido como la fracción de muestras que superan ±5σ.

Opcionalmente, la característica `corrEOG` utiliza canales EOG cuando están disponibles.

Todas las características se normalizan mediante **z-score** y se combinan para obtener las puntuaciones:

- `ocular_score`
- `muscle_score`
- `line_score`
- `jump_score`

La **puntuación global de artefacto**, `artifact_score`, se define como el máximo de estas puntuaciones.

Los componentes con

```math
artifact\_score > 1.5
```
"

# ╔═╡ f87ef5c2-1dfc-4d9a-be14-6238bf723025
md"
## Métricas cuantitativas

Esta sección calcula métricas de referencia sobre los datos filtrados  
(**antes del ICA**) para poder comparar posteriormente con los datos limpios.

El objetivo es cuantificar propiedades estadísticas de la señal EEG que permitan
evaluar la efectividad del proceso de limpieza mediante ICA.

Se emplean dos métricas principales:

**Kurtosis**

- Mide el grado de **no-gaussianidad** de la señal.
- Valores altos indican:
  - Presencia de picos pronunciados
  - Posibles outliers o artefactos
- Es especialmente útil para detectar componentes no neuronales.

**Energía**

- Definida como la **varianza de la señal**.
- Representa la **potencia** de la señal EEG.
- Permite comparar cambios globales en la amplitud tras el procesamiento.

Estas métricas sirven como referencia para:

- Comparar señales **antes vs después del ICA**
- Evaluar si el ICA ha:
  - Reducido artefactos (↓ kurtosis)
  - Modificado la potencia de la señal (energía)
- Validar cuantitativamente el proceso de limpieza

Trabajar con métricas cuantitativas complementa la inspección visual y permite
una evaluación más objetiva del preprocesamiento de señales EEG.
"

# ╔═╡ 1eb23ea8-5f04-4b28-b293-b8b376f259fd
begin
println("=========================================================================")
println("  Preliminares: Métricas cuantitativas (datos de dict_EEG_Lowpass.bin)")
println("========================================================================\n")

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
end

# ╔═╡ e2a46fba-a630-4be3-920b-1023ebbfb570
begin
# Obtener canales ordenados (orden alfabético para consistencia)
channels_lowpass = sort(collect(keys(dict_EEG_Lowpass)))
end

# ╔═╡ eef5e8c8-5415-4407-8ad8-5b15dafcfb96
begin
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
end

# ╔═╡ d290c6c6-eeca-4a72-87c6-6d5c4058736d
md"
## Carga resultados de ICA

Se cargan los resultados del análisis *ICA* realizado previamente:

   - **S**: componentes independientes (señales fuente)
   - **W_total**: matriz de desmezcla total
   - **A**: matriz de mezcla (columnas = mapas topográficos de cada IC)
   - **channels**: orden de canales usado en el ICA

También se calculan métricas básicas (**kurtosis** y **energía**) para cada componente, que servirán como referencia para la evaluación posterior.
"

# ╔═╡ d7d45a88-a256-462d-ada2-50622508b62c
begin
dict_ICA = Serialization.deserialize(path_dict_ica)
end

# ╔═╡ 86cfb6a8-c4d0-4fe7-9626-65f82eede2b6
begin
println("✔ Carga de datos de ICA:")
println("  S (n_comp × n_samples):        ", size(S))
println("  W_total (n_comp × n_channels):  ", size(W_total))
println("  A (n_channels × n_comp):        ", size(A))
println("  channels: ", channels, "\n")
end

# ╔═╡ f604db17-6269-4b09-b461-0f8f699d8485
md"
Calcular límite de color para mapas topográficos. Se usa el máximo absoluto de todos los pesos para normalizar la escala de colores
"

# ╔═╡ 5d86e143-86a3-4e46-98f2-e65266e9ea1a
begin
global topo_clim = maximum(abs, vec(A))
println("Límite color: ", topo_clim)
end

# ╔═╡ 73fbcedf-3755-4f28-a308-94f4965feaca


# ╔═╡ 0356c672-a4b1-4624-abcc-b97cf501cb8c
begin
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
end

# ╔═╡ 70d2ad9e-2a67-48b0-a5a3-a64b9fc974b0
md"
**CARGA FUNCIONES**
"

# ╔═╡ 282c6662-55e5-44a6-9af1-030b73b6545e
begin
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
end

# ╔═╡ 9f250578-81b4-486f-8d30-b54914158bf4
xs2d_all, ys2d_all, haspos = project_elec_xy(elec_pos, channels)

# ╔═╡ cda10a74-3790-499a-bc8a-2415bd3bb6cc
md"
Función para dibujar topografía de un componente **ICA** (A)  
"

# ╔═╡ 753841c5-5b71-4fcf-a482-5a30317ecb29
begin
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
end

# ╔═╡ 1ef7f377-0d74-4e14-879e-30c13dd54df9
begin
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
end

# ╔═╡ 2312ebe1-9c4c-4e26-856e-73c414d7441f
begin
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
end

# ╔═╡ 7077db9b-1233-4085-a6f3-e2e18d063c3c
begin
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
end

# ╔═╡ 77f56832-2314-4385-a847-571fa65eb65f
begin
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
end

# ╔═╡ 0f1d4728-d74e-4871-9845-bd470e9c990b
begin
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
end

# ╔═╡ c8fb5e5f-ca7f-41f6-ad13-4d7211da73fc
md"
---
## Visualizaciones componentes ICA

Generar y guardar `plot_ic_summary` para todas las componentes ICA

Esta sección genera automáticamente visualizaciones completas (tipo EEGLAB) para todas las componentes ICA. Cada visualización incluye:

   - Mapa topográfico (distribución espacial)
   - Serie temporal (actividad temporal)
   - Espectro de potencia (distribución espectral)

Estas visualizaciones se guardan como archivos PNG y permiten inspección manual de todas las componentes para verificar la evaluación automática.
"

# ╔═╡ e9ca5fa3-21fb-42cb-8114-8538a69ebc60
begin
println("=" ^ 60)
println("  Generando resúmenes visuales de componentes ICA")
println("=" ^ 60 * "\n")

# Crear directorio de salida si no existe
# Las figuras se guardan en results/figures/ICA_cleaning/
dir_output = stage_dir(:ICA_cleaning; kind = :figures)
isdir(dir_output) || mkpath(dir_output)

# Obtener número de componentes ICA
fs = 500.0  # frecuencia de muestreo (Hz)

println("Generando resúmenes para $n_comp componentes ICA...")
println("→ Directorio de salida: $dir_output")
println()
end

# ╔═╡ 95643a9e-59d2-4873-a8d6-ca3c07f8fb14
begin
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
end

# ╔═╡ 4615eb2e-9602-425e-899a-2a590e71e1cc
md"
---
## Scores para las componentes 

Estas funciones combinan features normalizadas (z-scores) para crear scores que indican la probabilidad de que un componente sea un tipo específico de artefacto.

Los scores se calculan como combinaciones lineales ponderadas de las features relevantes.

- Score_1. Ratio espacial-temporal
- Score_2. Potencias de banda PSD
- Score_3. Blink Ratio
- Score_4. EMG Ratio
- Score_5. Line Ratio
- Score_6. Kurtosis
- Score_7. Fracción por muestras
- Score_8. Correlación con EOG
"

# ╔═╡ 758bdd3b-9681-496a-b5ed-0844552b851b
md"
**Funciones Auxiliares Scores**
"

# ╔═╡ 1720bdbb-e56b-406f-a524-7ecef994c8f7
begin
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
end

# ╔═╡ 319b5a80-5f4c-4973-ae94-75d09f735876
begin
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
end

# ╔═╡ 48d49f69-58de-4fff-89f8-7f71eeb5adc9
begin
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
end

# ╔═╡ 58424b09-e5c9-431f-a971-1faa4a1d6104
begin
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
end

# ╔═╡ 6c8a21b0-a2d8-4d78-90d5-0acf83909a84
begin
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

end

# ╔═╡ 9bd532f8-dd85-460e-9d8b-c43f7251c6e3
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

# ╔═╡ 6e0aef9c-f860-45c7-abf7-28a7e3691308
"""
    line_score(line_z)

Score para identificar componentes de ruido de línea eléctrica (50/60 Hz).

Usa directamente el z-score del line_ratio (peso 1.0).
Un score alto indica que el componente contiene principalmente ruido de línea.
"""
function line_score(line_z)
    return line_z      # aquí el peso es 1 directamente
end

# ╔═╡ 452b7200-4636-42c8-91d1-660ff4497ae7
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

# ╔═╡ 0f66f445-fe02-4c27-b7a1-2dd5a15ac405
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

# ╔═╡ 53b80a6f-b5f4-434b-bec6-6ab699b4268c
md"
## Índices canales frontales y temporales

Esta sección identifica automáticamente los canales frontales y temporales basándose en la nomenclatura estándar 10-20. Estos canales se usan para calcular los ratios espaciales que ayudan a identificar tipos de artefactos:

   - Canales frontales: típicamente afectados por parpadeos y movimientos oculares
   - Canales temporales: típicamente afectados por actividad muscular (EMG)

La identificación se hace por nombre de canal usando prefijos estándar. 
"

# ╔═╡ 81c89737-2b24-48ab-8c9e-4fea3ea5fa78
md"
Identificar **canales frontales (F, FC, FT, FP)**. Estos canales están en la parte frontal del cuero cabelludo
"

# ╔═╡ aba0ac49-d995-45fb-a650-dcba81b8379f
begin
frontal_idxs = Int[]
for (i, ch) in enumerate(channels)
    ch_upper = uppercase(ch)  # convertir a mayúsculas para comparación
    # Buscar canales que empiecen con F, FC, FT o FP
    if startswith(ch_upper, "F") || startswith(ch_upper, "FC") || 
       startswith(ch_upper, "FT") || startswith(ch_upper, "FP")
        push!(frontal_idxs, i)
    end
end
println("📍 Canales frontales detectados: ", [channels[i] for i in frontal_idxs])
println()
end

# ╔═╡ bf4f70fc-c2b3-4932-9320-20d88cde92a6
md"
Identificar **canales temporales (T, TP)**. Estos canales están en las regiones temporales (lados de la cabeza)
"

# ╔═╡ d20e61ec-8992-4629-81f2-67390dad967e
begin
temporal_idxs = Int[]
for (i, ch) in enumerate(channels)
    ch_upper = uppercase(ch)
    # Buscar canales que empiecen con T o TP
    if startswith(ch_upper, "T") || startswith(ch_upper, "TP")
        push!(temporal_idxs, i)
    end
end
println("📍 Canales temporales detectados: ", [channels[i] for i in temporal_idxs])
println()
end

# ╔═╡ bfad913c-0ef8-46fc-a168-20d0f3eff740
md"
**Canales EOG** (no disponibles en estos datos). Si hubiera canales EOG, se usarían para calcular correlaciones que ayudarían a identificar componentes oculares
"

# ╔═╡ 2018768c-862a-44be-8f4e-28fa1ca54c76
eog = nothing

# ╔═╡ c500f42c-34ad-4256-8b5e-d04990f41fd7
md"
Evaluar todos los ICs usando el sistema automático. Esta función calcula features, scores y etiquetas para todas las componentes
"

# ╔═╡ 3ee89855-3d96-4107-bf43-c57c4ad96663
begin
df_all = evaluate_ics(A, S, fs;
                      frontal_idxs = frontal_idxs,
                      temporal_idxs = temporal_idxs,
                      eog = eog)
display(df_all)
end

# ╔═╡ 904430e1-10d6-40b0-8a08-a2325bc11ee1
md"Guardar resultados de evaluación en tabla
"

# ╔═╡ 0dc439c9-05ba-48ed-8150-3b50ef9dfcc9
begin
dir_tables = stage_dir(:ICA_cleaning; kind = :tables)
isdir(dir_tables) || mkpath(dir_tables)
path_table = joinpath(dir_tables, "ICA_components_evaluation.csv")
CSV.write(path_table, df_all)
println("\n💾 Tabla de evaluación de componentes ICA guardada en: $path_table\n")
end

# ╔═╡ c5e3850c-c33d-45da-847b-7727630392d9
md"
---
# Selección componentes a eliminar

Esta sección identifica automáticamente los componentes artefactuales basándose en los scores calculados y luego reconstruye los datos EEG eliminando esos componentes.

PROCESO:
   1. Identificar componentes con artifact_score > umbral
   2. Agrupar por tipo de artefacto (parpadeo, músculo, línea, salto)
   3. Anular las componentes artefactuales en S_clean (poner a cero)
   4. Reconstruir datos limpios: X_clean = A * S_clean
   5. Guardar datos limpios y metadatos
"

# ╔═╡ c8d481a7-6ffd-4ae1-b5a4-69fd5cac68f0
begin
println("====================================")
println("  Selección de Componentes Artefactuales")
println("====================================\n")
artifact_thresh = 1.5
bad_ics = Int[]  # Lista de componentes a eliminar
artifacts_by_type = Dict{String, Vector{Int}}(
    "parpadeo" => Int[],
    "músculo" => Int[],
    "línea" => Int[],
    "salto" => Int[],
    "otros" => Int[]
)
end

# ╔═╡ 50c0d9b2-71f1-401a-a471-2c81f93b9b31
begin
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
end

# ╔═╡ 32ee530f-e66a-44fb-ba43-96f59dde3616
begin
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
end

# ╔═╡ 308a57b0-5e60-41b9-8ed6-a2545ef434f8
md"
# Reconstrucción EEG limpio
"

# ╔═╡ dfc805d8-b213-4aa8-b4ce-87c5a7df53ca
begin
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
end

# ╔═╡ a2687d65-dd55-4ef1-a4ad-c6857c4a736c
begin
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
end

# ╔═╡ 3a4ad856-e5ae-431e-a010-f69f7b199f6c
md"
Guardar versión completa con info ICA + cleaning
"

# ╔═╡ 4ce98dc3-4644-4a70-8762-f94df502c485
begin
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
end

# ╔═╡ 9e09d533-942b-4a1e-a741-51d4eb82ae63
md"
# Comparativa (antes/después ICA)

Esta sección prepara los datos para comparar visualmente el EEG antes y después
de la limpieza ICA. Se cargan los datos originales (filtrados) y los datos limpios, y se convierten a formato de matriz para facilitar la visualización.
"

# ╔═╡ fd41bc1b-fea1-4bfd-bacd-3d953fda896a
begin
# Cargar datos después del ICA (limpio)
# Estos son los datos después de eliminar componentes artefactuales
dict_EEG_ICA_clean = Serialization.deserialize(path_dict_clean)
end

# ╔═╡ 04b9b52c-d9e6-4232-996a-bcf5f65cd2af
begin
# Crear matrices para comparación
# Cada fila es un canal, cada columna es una muestra temporal
eeg_orig = zeros(Float64, n_channels, n_samples)   # Antes del ICA (Lowpass)
eeg_clean = zeros(Float64, n_channels, n_samples)  # Después del ICA (limpio)
end

# ╔═╡ e8b98bb9-8bb1-4793-944d-384b4153cb6b
md" 
Rellenar matrices con datos de los diccionarios. Se mantiene el mismo orden de canales para comparación directa
"

# ╔═╡ f716dd11-c0a7-4a9a-8113-5d449255fad1
begin
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
end

# ╔═╡ bbe9e255-f695-4888-9163-ce8ae6726865
md"
## Funciones auxiliares
"

# ╔═╡ 8d50e6a0-0a2c-45bf-bfa8-88733e52ebaf
begin
# Alias para compute_psd usando la función existente
compute_psd(signal::AbstractVector, fs::Real) = calculate_PSD_signal(signal, fs)
end

# ╔═╡ a75dbfaf-da36-497b-926a-990872f74364
begin
# Función para calcular kurtosis de un canal
kurtosis_channel(signal::AbstractVector) = kurtosis(signal)
end

# ╔═╡ beec7fc7-1620-43c6-afb6-3fd7c098250a
begin
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
end

# ╔═╡ 51b05771-8410-40ca-b032-95792bf34935
begin
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
end

# ╔═╡ 856e9b3b-3261-4236-8b9f-3d8a9b60c2a4
begin
# Ejemplo: visualizar segmentos de EEG de los canales 1, 10 y 20
# desde el segundo 10 hasta el segundo 14
println("📊 Visualizando segmentos de EEG (canales 1, 10, 20)...")
plot_eeg_segment(eeg_orig, eeg_clean, fs;
                 channels = [1, 10, 20],
                 t_start = 10.0,
                 t_len = 4.0,
                 chan_labels = channels)   
end

# ╔═╡ 9eb77221-5d78-48e4-b7f5-0317c3e1b492
begin
# Ejemplo: visualizar segmentos de EEG de los canales 1, 10 y 20
# desde el segundo 10 hasta el segundo 14  
plot_channel_psd(eeg_orig, eeg_clean, fs;
                 chan = 1,
                 fmax = 80.0,
                 chan_labels = nothing)
end

# ╔═╡ 7d96aef0-87a7-4148-b095-ea33d2b5af78
begin
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
end

# ╔═╡ abf2aa4b-b279-458a-afd0-fba54065d893
begin
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
end

# ╔═╡ c0d7e4d3-4bef-47bc-bd35-88d6147a7613
begin
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
end

# ╔═╡ ff20dba2-7a6c-4853-a7b9-b47a0f60696d
begin
println("\n📊 Visualizando Kurtosis y Energía por canal...")
plot_kurtosis_energy(eeg_orig, eeg_clean; chan_labels = channels)
end

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
# ╠═7a5d6ddb-33a3-422a-88a4-f76fc8f0c620
# ╠═2ff7ee54-7bec-4e21-b804-e9fcc2f86a11
# ╠═8a0cde3b-88ab-4b56-bef6-8b9a394b1c2e
# ╠═49ca7fca-11dc-4f81-b74f-b54e8b252e10
# ╠═789269a7-1c3a-41a2-811c-c5b8e08ee1f3
# ╠═177911fa-6acc-4cd5-b862-81fa0fef1daa
# ╠═c5f19470-bda8-4420-b7f4-c3a1c94b4753
# ╠═01b14950-89b8-4df1-981f-44a13f6026cb
# ╠═09a79ac1-4430-4729-bca5-5425a2ba330e
# ╠═7523a05e-bb9e-48bc-9720-9a6178941273
# ╠═580b258a-cde8-4c8a-bdd4-5c71f411cc6b
# ╠═d9936d83-6fae-47af-9890-4c0594e2c782
# ╠═9fd4518f-a809-4c1e-8a6f-ca26a557094e
# ╠═98207097-b070-4b07-8ff5-1d29e64e5142
# ╠═ad20b61f-27c5-46bd-b7c4-2d5125a523b7
# ╠═3b9148cd-9027-42e5-9446-b1b1d1c39fdf
# ╠═04fd2eda-02c4-49b5-ac82-c747a300bdf4
# ╠═62db513e-efef-4f65-898e-4d1283e3c55a
# ╠═818a6675-189e-44fe-abb5-94a1dab065ea
# ╠═7e8b0401-0b43-43e4-a3f2-6f51c471da46
# ╠═c289b2e6-294e-48ea-a8c7-75f1d5448d13
# ╠═a036b633-2812-444e-9f72-7117faaac04e
# ╠═c5ddf477-9edd-4ce7-aefa-d9a17b0c0fe0
# ╠═eed453df-7646-41a7-af88-c5beb48902cd
# ╠═76fb19ba-33ed-4214-992e-a3d1f0a1983f
# ╠═f87ef5c2-1dfc-4d9a-be14-6238bf723025
# ╠═1eb23ea8-5f04-4b28-b293-b8b376f259fd
# ╠═e2a46fba-a630-4be3-920b-1023ebbfb570
# ╠═eef5e8c8-5415-4407-8ad8-5b15dafcfb96
# ╠═d290c6c6-eeca-4a72-87c6-6d5c4058736d
# ╠═d7d45a88-a256-462d-ada2-50622508b62c
# ╠═86cfb6a8-c4d0-4fe7-9626-65f82eede2b6
# ╠═f604db17-6269-4b09-b461-0f8f699d8485
# ╠═5d86e143-86a3-4e46-98f2-e65266e9ea1a
# ╠═73fbcedf-3755-4f28-a308-94f4965feaca
# ╠═0356c672-a4b1-4624-abcc-b97cf501cb8c
# ╠═70d2ad9e-2a67-48b0-a5a3-a64b9fc974b0
# ╠═282c6662-55e5-44a6-9af1-030b73b6545e
# ╠═9f250578-81b4-486f-8d30-b54914158bf4
# ╠═cda10a74-3790-499a-bc8a-2415bd3bb6cc
# ╠═753841c5-5b71-4fcf-a482-5a30317ecb29
# ╠═1ef7f377-0d74-4e14-879e-30c13dd54df9
# ╠═2312ebe1-9c4c-4e26-856e-73c414d7441f
# ╠═7077db9b-1233-4085-a6f3-e2e18d063c3c
# ╠═77f56832-2314-4385-a847-571fa65eb65f
# ╠═0f1d4728-d74e-4871-9845-bd470e9c990b
# ╠═c8fb5e5f-ca7f-41f6-ad13-4d7211da73fc
# ╠═e9ca5fa3-21fb-42cb-8114-8538a69ebc60
# ╠═95643a9e-59d2-4873-a8d6-ca3c07f8fb14
# ╠═4615eb2e-9602-425e-899a-2a590e71e1cc
# ╠═758bdd3b-9681-496a-b5ed-0844552b851b
# ╠═1720bdbb-e56b-406f-a524-7ecef994c8f7
# ╠═319b5a80-5f4c-4973-ae94-75d09f735876
# ╠═48d49f69-58de-4fff-89f8-7f71eeb5adc9
# ╠═58424b09-e5c9-431f-a971-1faa4a1d6104
# ╠═6c8a21b0-a2d8-4d78-90d5-0acf83909a84
# ╠═9bd532f8-dd85-460e-9d8b-c43f7251c6e3
# ╠═6e0aef9c-f860-45c7-abf7-28a7e3691308
# ╠═452b7200-4636-42c8-91d1-660ff4497ae7
# ╠═0f66f445-fe02-4c27-b7a1-2dd5a15ac405
# ╠═53b80a6f-b5f4-434b-bec6-6ab699b4268c
# ╠═81c89737-2b24-48ab-8c9e-4fea3ea5fa78
# ╠═aba0ac49-d995-45fb-a650-dcba81b8379f
# ╠═bf4f70fc-c2b3-4932-9320-20d88cde92a6
# ╠═d20e61ec-8992-4629-81f2-67390dad967e
# ╠═bfad913c-0ef8-46fc-a168-20d0f3eff740
# ╠═2018768c-862a-44be-8f4e-28fa1ca54c76
# ╠═c500f42c-34ad-4256-8b5e-d04990f41fd7
# ╠═3ee89855-3d96-4107-bf43-c57c4ad96663
# ╠═904430e1-10d6-40b0-8a08-a2325bc11ee1
# ╠═0dc439c9-05ba-48ed-8150-3b50ef9dfcc9
# ╠═c5e3850c-c33d-45da-847b-7727630392d9
# ╠═c8d481a7-6ffd-4ae1-b5a4-69fd5cac68f0
# ╠═50c0d9b2-71f1-401a-a471-2c81f93b9b31
# ╠═32ee530f-e66a-44fb-ba43-96f59dde3616
# ╠═308a57b0-5e60-41b9-8ed6-a2545ef434f8
# ╠═dfc805d8-b213-4aa8-b4ce-87c5a7df53ca
# ╠═a2687d65-dd55-4ef1-a4ad-c6857c4a736c
# ╠═3a4ad856-e5ae-431e-a010-f69f7b199f6c
# ╠═4ce98dc3-4644-4a70-8762-f94df502c485
# ╠═9e09d533-942b-4a1e-a741-51d4eb82ae63
# ╠═fd41bc1b-fea1-4bfd-bacd-3d953fda896a
# ╠═04b9b52c-d9e6-4232-996a-bcf5f65cd2af
# ╠═e8b98bb9-8bb1-4793-944d-384b4153cb6b
# ╠═f716dd11-c0a7-4a9a-8113-5d449255fad1
# ╠═bbe9e255-f695-4888-9163-ce8ae6726865
# ╠═8d50e6a0-0a2c-45bf-bfa8-88733e52ebaf
# ╠═a75dbfaf-da36-497b-926a-990872f74364
# ╠═beec7fc7-1620-43c6-afb6-3fd7c098250a
# ╠═51b05771-8410-40ca-b032-95792bf34935
# ╠═856e9b3b-3261-4236-8b9f-3d8a9b60c2a4
# ╠═9eb77221-5d78-48e4-b7f5-0317c3e1b492
# ╠═7d96aef0-87a7-4148-b095-ea33d2b5af78
# ╠═abf2aa4b-b279-458a-afd0-fba54065d893
# ╠═c0d7e4d3-4bef-47bc-bd35-88d6147a7613
# ╠═ff20dba2-7a6c-4853-a7b9-b47a0f60696d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
