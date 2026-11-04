### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 4d8ca8b0-149f-4f6f-93f4-4e2d7f1a0b01
begin
	using PlutoUI
	PlutoUI.TableOfContents(title = "Contenido")
end

# ╔═╡ d249d3f0-cf95-44e1-aec6-8ca4ad0fe102
begin
	include(joinpath(@__DIR__, "..", "_template_base.jl"))
	using .PlutoTemplateBase
	using CSV, DataFrames, Serialization, Statistics, StatsBase, DSP, Plots, Dates, InlineStrings
end

# ╔═╡ a4097250-6b5a-4415-8469-e2a8dab637ec
# Si se ejecuta este script directamente (fuera del módulo EEG_Julia),
# cargamos utilidades de rutas para disponer de `stage_dir`.
if !@isdefined(stage_dir)
    include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
end

# ╔═╡ fc2c1f3e-448d-4bd0-8208-4d531d1e5800
md"""
**PAQUETES CARGADOS**
"""

# ╔═╡ e857112f-df27-4f7f-bcc7-34421e0c3103
notebook_intro("PREPROCESSING")

# ╔═╡ e456bde5-04a7-4db8-a5d4-cf2068dfe29b
md"""

El **preprocesamiento** de las señales **EEG crudas** tiene como objetivo eliminar:

- ruido de línea eléctrica,
- interferencias relacionadas con el hardware,
- deriva lenta del baseline,

y limitar el ancho de banda antes del análisis de conectividad específico por banda.

La implementación se proporciona en `src/Preprocessing/filtering.jl`, que:

1. carga los datos serializados desde `data/Preprocessing/IO/dict_EEG.bin`,
2. aplica una **cascada de filtros digitales**,
3. guarda en cada paso el resultado en `data/Preprocessing/filtering/`,
4. genera gráficos de la respuesta del filtro (magnitud y fase),
5. compara la **PSD media** antes y después del filtrado.

Todos los filtros son del tipo **Butterworth**. Las etapas **pasa-altos** y **pasa-bajos** utilizan `filtfilt` (filtrado de fase cero), de modo que:

- el **orden efectivo es el doble** del orden de diseño
- se preserva la **alineación temporal** de la señal.

La cascada de filtrado se aplica en el siguiente orden:

- Filtro **Notch** (50 Hz): elimina interferencia de la red eléctrica.
- Filtro **Bandreject** (100 Hz): elimina posible interferencia de hardware.
- Filtro **Highpass** (0.5 Hz): elimina deriva (drift) preservando actividad lenta.
- Filtro **Lowpass** (150 Hz): limita alias y ruido de alta frecuencia.

---
"""

# ╔═╡ 204da7a3-9f60-4918-97f4-06a2d621dc15
md"
## Carga de datos (dict_EEG.bin)

Se carga el diccionario **EEG** desde el paso de IO (`dict_EEG.bin`) y se configuran las rutas de salida para los datos filtrados.
"

# ╔═╡ 6337d051-3f9e-4dd9-9925-30b188eb6ac0
md"""
!!! info "Pipeline Fase 1: Preprocesamiento"

    **Objetivo:**
    - Eliminar interferencias de red eléctrica y del hardware, deriva lenta y ruido fuera de banda del **EEG crudo** generado en la **Fase 0**.
    - Almacenar las señales filtradas para análisis posteriores específicos por banda.

    **Entrada.**
    - **EEG** serializado procedente de la **Fase** ``0`` (por ejemplo, `data/IO/dict_EEG.bin`).
    - Frecuencia de muestreo ``f_s = 500`` Hz.

    **Salida.**
    - Señales filtradas almacenadas en `data/filtering/` tras cada etapa (Notch, Bandreject, Highpass, Lowpass).
    - Gráficas de respuesta del **filtro** (magnitud y fase) y comparaciones de **PSD promedio** (antes/después de cada etapa) para control de calidad.

    **Pasos de procesamiento (alineados con `src/filtering.jl`).**

    1. **Carga de datos EEG.**  
       Cargar el diccionario EEG desde `data/IO/`; definir ``f_s`` y la lista de canales.

    2. **Filtro Notch 50 Hz.**  
       Aplicar un filtro Butterworth *bandstop* (orden 4, ancho 1 Hz).  
       Guardar `dict_EEG_Notch.bin`.  
       Representar la respuesta del filtro y comparar la PSD promedio (señal original vs señal con notch).

    3. **Rechazo de banda 100 Hz.**  
       Aplicar un filtro Butterworth *bandstop* (orden 4, ancho de banda 1 Hz) sobre la salida del notch.  
       Guardar `dict_EEG_Bandreject.bin`.  
       Representar la respuesta del filtro y comparar PSD.

    4. **Filtro paso alto 0.5 Hz.**  
       Aplicar un filtro Butterworth *high-pass* (orden 4, `filtfilt`).  
       Guardar `dict_EEG_Highpass.bin`.  
       Representar la respuesta del filtro y comparar PSD.

    5. **Filtro paso bajo 150 Hz.**  
       Aplicar un filtro Butterworth *low-pass* (orden 4, `filtfilt`).  
       Guardar `dict_EEG_Lowpass.bin`.  
       Representar la respuesta del filtro y comparar PSD.

    6. **Salida del preprocesamiento.**  
       Devolver el diccionario completamente filtrado (o la ruta a `dict_EEG_Lowpass.bin`) para la **Fase 2a (ICA)**.

    Opcionalmente, el re-referenciado, la eliminación de artefactos y la segmentación en épocas pueden realizarse en un paso separado antes del cálculo de conectividad.
"""

# ╔═╡ b1a100f5-2f4f-4d30-b4c4-1dca20c50e04
begin
# -----------------------------------------------------------------------------------
# 1. CARGA DE DATOS Y CONFIGURACIÓN
# -----------------------------------------------------------------------------------
# Se carga el diccionario EEG desde el paso de IO (dict_EEG.bin) y se configuran
# las rutas de salida para los datos filtrados.

dir_io = stage_dir(:IO)
path_dict = joinpath(dir_io, "dict_EEG.bin")
if !isfile(path_dict)
    error("No se encontró $(abspath(path_dict)). Ejecuta antes src/Setup/IO.jl para generar dict_EEG.bin.")
end
dict_EEG = Serialization.deserialize(path_dict)

# Directorio base para datos filtrados
dir_filtering = stage_dir(:filtering)
end

# ╔═╡ aa74685e-d70c-4b92-99c6-9c8f051d387e
begin
# -----------------------------------------------------------------------------------
# 2. INFORMACIÓN TEMPORAL DEL REGISTRO
# -----------------------------------------------------------------------------------
fs = 500                                        		# Frecuencia de muestreo (Hz)
n_muestras = length(dict_EEG[first(keys(dict_EEG))])  	# Muestras por canal
tiempo_seg = collect(0:(n_muestras-1)) ./ fs    		# En segundos
duracion_total = tiempo_seg[end]
end

# ╔═╡ 2a1d5d68-6618-467d-9701-1167dd2d54cc
begin
println("=" ^ 40)
println("⏱️  INFORMACIÓN TEMPORAL")
println("=" ^ 40)
println("  Frecuencia de muestreo: $(fs) Hz")
println("  Número de muestras: $(n_muestras)")
println("  Duración total: $(round(duracion_total, digits=2)) segundos ($(round(duracion_total/60, digits=2)) minutos)")
println()
end

# ╔═╡ 557acfe8-89d9-4060-b851-b0fef0412631
begin
# Mostramos los datos de los canales
println("📊 Datos de los canales")
println("-" ^ 80)
display(dict_EEG)
println()
end

# ╔═╡ 583ea341-04cf-4c52-ba27-5a1e9617d397
begin
	# Filtrado de las señales
	println("=" ^ 50)
	println("📊 Pasos a seguir en el filtrado de las señales:")
	println("=" ^ 50)
	println("STEP Nº1. Notch Filter")
	println("STEP Nº2. Bandreject filter")
	println("STEP Nº3. Highpass filter")
	println("STEP Nº4. Lowpass filter")
	println()
end

# ╔═╡ d620dfe1-1628-4da7-8e1b-b456d659dbe8
md"
## Definición de funciones

- Funciones para filtros Butterworth (notch, bandreject, highpass, lowpass)
- Visualización de respuesta en frecuencia, PSD y guardado de señales filtradas
"

# ╔═╡ 39690c4c-c5c8-4235-b832-4e78d8f72b3e
begin
"""Aplica un filtro Butterworth notch_filter a un vector x"""
function notch_filter(x::Vector{<:Real}, freq::Real, fs::Real, order::Int = 5, width::Real = 1.0)  
    # Calculamos las frecuencias inferior y superior del notch
    freq_low = freq - width/2
    freq_high = freq + width/2
    # Normalizamos las frecuencias (0–1)
    wn_low = freq_low / (fs/2)
    wn_high = freq_high / (fs/2)
    responsetype = Bandstop(wn_low, wn_high)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    y = filt(flt, x)                # señal filtrada
    return y, flt                   # devolvemos también el filtro para analizarlo
end
end

# ╔═╡ 722a26b7-1f85-489a-a9bb-66d8104f93ca
begin
"""Aplica un filtro Butterworth bandreject (bandstop) a un vector x"""
function bandreject_filter(x::Vector{<:Real}, freq::Real, fs::Real, order::Int = 4, bandwidth::Real = 1.0)  
    # Calculamos las frecuencias inferior y superior del bandreject
    freq_low = freq - bandwidth/2
    freq_high = freq + bandwidth/2
    # Normalizamos las frecuencias (0–1)
    wn_low = freq_low / (fs/2)
    wn_high = freq_high / (fs/2)
    responsetype = Bandstop(wn_low, wn_high)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    y = filt(flt, x)                # señal filtrada
    return y, flt                   # devolvemos también el filtro para analizarlo
end
end

# ╔═╡ 3249aa95-e924-4e96-bf72-30c9a02c8607
begin
"""Aplica un filtro Butterworth pasa-bajos a un vector x con zero phase shift (filtfilt)
   
   Nota: filtfilt aplica el filtro dos veces (adelante y atrás), duplicando el orden efectivo.
   Para obtener orden N efectivo, se debe diseñar un filtro de orden N/2.
   
   Parámetros:
   - x: señal de entrada
   - cutoff: frecuencia de corte en Hz
   - fs: frecuencia de muestreo en Hz
   - order: orden del diseño del filtro (el orden efectivo será order × 2 con filtfilt)
   - zero_phase: si es true (por defecto), usa filtfilt para zero phase shift
"""
function lowpass_filter(x::Vector{<:Real}, cutoff::Real, fs::Real, order::Int = 5; zero_phase::Bool = true)
    wn = cutoff / (fs/2)            # Frecuencia de corte normalizada (0–1)
    responsetype = Lowpass(wn)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    if zero_phase
        y = filtfilt(flt, x)        # señal filtrada con zero phase (orden efectivo = order × 2)
    else
        y = filt(flt, x)            # señal filtrada con desplazamiento de fase (orden efectivo = order)
    end
    return y, flt                   # devolvemos también el filtro para analizarlo
end
end

# ╔═╡ d359220a-c31c-4ce4-91e6-b58a7a9e7d08
begin
"""Aplica un filtro Butterworth pasa-altos a un vector x con zero phase shift (filtfilt)
   
   Nota: filtfilt aplica el filtro dos veces (adelante y atrás), duplicando el orden efectivo.
   Para obtener orden N efectivo, se debe diseñar un filtro de orden N/2.
   
   Parámetros:
   - x: señal de entrada
   - cutoff: frecuencia de corte en Hz
   - fs: frecuencia de muestreo en Hz
   - order: orden del diseño del filtro (el orden efectivo será order × 2 con filtfilt)
   - zero_phase: si es true (por defecto), usa filtfilt para zero phase shift
"""
function highpass_filter(x::Vector{<:Real}, cutoff::Real, fs::Real, order::Int = 5; zero_phase::Bool = true)
    wn = cutoff / (fs/2)            # Frecuencia de corte normalizada (0–1)
    responsetype = Highpass(wn)
    designmethod = Butterworth(order)
    flt = digitalfilter(responsetype, designmethod)
    if zero_phase
        y = filtfilt(flt, x)        # señal filtrada con zero phase (orden efectivo = order × 2)
    else
        y = filt(flt, x)            # señal filtrada con desplazamiento de fase (orden efectivo = order)
    end
    return y, flt                   # devolvemos también el filtro para analizarlo
end 
end

# ╔═╡ 8ade48ae-34e1-462b-9c0e-206547e7432b
begin
"""Plotea la respuesta en frecuencia de un filtro (magnitud y fase)"""
function plot_filter_response(flt, fs::Real; title::String = "Respuesta del Filtro", n_points::Int = 1024)
    # Frecuencias en rad/muestra donde evaluar (0..π)
    ω = range(0, stop = π, length = n_points)
    # Respuesta compleja del filtro en esas frecuencias
    H = freqresp(flt, ω)
    # Pasamos a Hz reales
    f = ω .* fs ./ (2π)
    # Magnitud en dB y fase (convertida a grados sexagesimales)
    mag_db   = 20 .* log10.(abs.(H))
    fase_rad = angle.(H)
    fase_deg = fase_rad .* (180 / π)
    # Aseguramos que la fase esté en el rango [-180°, 180°]
    fase_deg = mod.(fase_deg .+ 180, 360) .- 180
    # Plot de magnitud
    p_mag = Plots.plot(
        f, mag_db,
        xlabel = "Frecuencia [Hz]",
        ylabel = "Magnitud [dB]",
        title  = "$(title) - Magnitud",
        grid   = true,
        linewidth = 2,
        legend = false,
    )
    # Plot de fase
    p_phase = Plots.plot(
        f, fase_deg,
        xlabel = "Frecuencia [Hz]",
        ylabel = "Fase [degrees]",
        title  = "$(title) - Fase",
        grid   = true,
        linewidth = 2,
        ylim = (-180, 180),
        yticks = -180:90:180,
        legend = false,
    )
    # Combinamos ambos plots
    p_combined = Plots.plot(p_mag, p_phase, layout = (2, 1), size = (800, 600))
    
    return p_combined, p_mag, p_phase
end
end

# ╔═╡ dca7f5db-a29d-4028-aaef-8d51f4ebfb9c
begin
"""Calcula la PSD promedio de un diccionario de señales EEG (sin plotear)"""
function calculate_PSD_average(dict_EEG_data::Dict, fs::Real)
    channels = collect(keys(dict_EEG_data))
    n_muestras = length(dict_EEG_data[first(channels)])
    
    # Calcular PSD para cada canal
    PSD = Dict(channel => begin
        p = welch_pgram(dict_EEG_data[channel]; fs=fs, window=hamming, nfft=n_muestras)
        (; freq = DSP.freq(p), power = DSP.power(p))
    end for channel in channels)
    
    # Calcular PSD promedio
    freqs_avg = PSD[first(channels)].freq
    avg_power = mapreduce(ch -> PSD[ch].power, +, channels) ./ length(channels)
    
    return freqs_avg, avg_power, PSD
end
end

# ╔═╡ 81c91017-c887-4aa0-b5f0-12436d83b0b3
begin
"""Calcula y visualiza la PSD (Power Spectral Density) para un diccionario de señales EEG"""
function calculate_and_plot_PSD(dict_EEG_data::Dict, fs::Real; title_prefix::String = "PSD")
    println("📈 CÁLCULO DE PSD POR CANAL")
    println("-" ^ 80)
    println()
    
    # Calcular PSD promedio
    freqs_avg, avg_power, PSD = calculate_PSD_average(dict_EEG_data, fs)
    
    println("  → PSD promedio de todos los canales")
    
    power_min, power_max = extrema(avg_power)
    y_min_log = floor(log10(power_min))
    y_max_log = ceil(log10(power_max))
    ylim_log = (10.0^y_min_log, 10.0^y_max_log)
    
    p_psd_avg = plot(
        xlabel = "Frecuencia (Hz)",
        ylabel = "Potencia (µV²/Hz)",
        title = "$(title_prefix) - Promedio de todos los canales",
        legend = false,
        xlim = (0, fs / 2),
        yscale = :log10,
        ylim = ylim_log,
    )
    
    plot!(p_psd_avg, freqs_avg, avg_power; label = "Promedio", lw = 2, color = :black)
    
    display(p_psd_avg)
    println()
    
    return PSD, p_psd_avg
end
end

# ╔═╡ 6e0cdb1b-d1cb-48ac-9167-c4d0c94ef398
begin
"""Compara y grafica superpuestas las PSD promedio de dos diccionarios de señales EEG"""
function compare_PSD_averages(dict_EEG_1::Dict, dict_EEG_2::Dict, fs::Real;
                              label_1::String = "Original",
                              label_2::String = "Filtrado",
                              title::String = "Comparación PSD Promedio")
    println("📈 COMPARACIÓN DE PSD PROMEDIO")
    println()
    
    # Calcular PSD promedio de ambos diccionarios
    freqs_1, avg_power_1, _ = calculate_PSD_average(dict_EEG_1, fs)
    freqs_2, avg_power_2, _ = calculate_PSD_average(dict_EEG_2, fs)
    
    # Calcular límites Y para ambos conjuntos de datos
    all_powers = vcat(avg_power_1, avg_power_2)
    power_min, power_max = extrema(all_powers)
    y_min_log = floor(log10(power_min))
    y_max_log = ceil(log10(power_max))
    ylim_log = (10.0^y_min_log, 10.0^y_max_log)
    
    # Crear plot con ambos promedios superpuestos
    p_comparison = Plots.plot(
        xlabel = "Frecuencia (Hz)",
        ylabel = "Potencia (µV²/Hz)",
        title = title,
        legend = :topright,
        xlim = (0, fs / 2),
        yscale = :log10,
        ylim = ylim_log,
    )
    
    Plots.plot!(p_comparison, freqs_1, avg_power_1; label = label_1, lw = 2, color = :blue)
    Plots.plot!(p_comparison, freqs_2, avg_power_2; label = label_2, lw = 2, color = :red)
    
    display(p_comparison)
    println()
    
    return p_comparison
end
end

# ╔═╡ 19a25445-09f4-4f9e-8281-dd16357f965b
begin
"""Guarda un diccionario de señales filtradas en un archivo binario"""
function save_filtered_signals(dict_EEG_filtered::Dict{String, Vector{Float64}}, filename::String, base_dir::String = dir_filtering)
    # Crear el directorio si no existe
    if !isdir(base_dir)
        mkpath(base_dir)
        println("📁 Directorio creado: $(abspath(base_dir))")
    end
    # Construir el path completo
    save_path = joinpath(base_dir, filename)
    
    # Guardar el archivo
    Serialization.serialize(save_path, dict_EEG_filtered)
    println("✓ Diccionario EEG filtrado guardado en: $(abspath(save_path))")
    println()
end
end

# ╔═╡ 21ec0155-70ba-4535-943c-7f2fc48dbd1f
md"
## Filtro Notch 

Filtro **band-stop** ($50\,\mathrm{Hz}$) para eliminar la interferencia de la red eléctrica (Europa).

**Diseño:**

- orden: $4$
- ancho de banda: $1\,\mathrm{Hz}$ alrededor de $50\,\mathrm{Hz}$
---
"

# ╔═╡ 2eb37338-f2ee-430e-8872-611374277c1b
begin
println("=" ^ 30)
println("📊 Filtro Notch (50 Hz)")
println("=" ^ 30)

Notch_cutoff = 50
Notch_order = 4
Notch_width = 1.0

println()
println("Frecuencia de corte: $(Notch_cutoff) Hz")
println("Orden del filtro: $(Notch_order)")
println("Ancho de banda: $(Notch_width) Hz")
println()

# Inicializamos el diccionario para las señales filtradas
dict_EEG_Notch = Dict{String, Vector{Float64}}()

# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_notch = notch_filter(dict_EEG[first(keys(dict_EEG))], Notch_cutoff, fs, Notch_order, Notch_width)

# Aplicamos el filtro Notch a todos los canales
for (channel, signal) in dict_EEG
    y_notch, _ = notch_filter(signal, Notch_cutoff, fs, Notch_order, Notch_width)
    dict_EEG_Notch[channel] = y_notch
end
end

# ╔═╡ fa145f7e-4175-4aa9-ac93-8c1793b72b1b
begin
# Ploteamos la respuesta del filtro Notch
p_notch, _, _ = plot_filter_response(flt_notch, fs; title = "Filtro Notch (50 Hz)")
p_notch
end

# ╔═╡ ca639a99-34d1-4fec-b800-8c7a03c45f7c
save_filtered_signals(dict_EEG_Notch, "dict_EEG_Notch.bin")

# ╔═╡ e3d66280-fcca-4299-a59b-5187c4b05c57
begin
# Comparar PSD promedio del original vs Notch filtrado (superpuestos)
p_psd_comparison = compare_PSD_averages(
    dict_EEG, 
    dict_EEG_Notch, 
    fs;
    label_1 = "Original",
    label_2 = "Notch (50 Hz)",
    title = "Comparación PSD Promedio: Original vs Notch"
)
end

# ╔═╡ 3cfdf522-7ab0-4ca4-be68-a5cfe7a421a9
md"

## Rechazo de banda

Filtro **band-stop estrecho** ($100\,\mathrm{Hz}$) para eliminar posibles interferencias relacionadas con el hardware alrededor de $100\,\mathrm{Hz}$ (y múltiplos de $50\,\mathrm{Hz}$). Su ancho es muy estrecho (1 Hz) para evitar una eliminación generalizada de componentes de alta frecuencia.

**Diseño:**

- orden: $4$
- ancho de banda: $1\,\mathrm{Hz}$
"

# ╔═╡ a78f7513-f36d-4854-80d9-c0f6d3afd3e4
begin
println("=" ^ 40)
println("📊 Filtro Bandreject (100 Hz)")
println("=" ^ 40)
Bandreject_freq = 100
Bandreject_order = 4
Bandreject_bandwidth = 1.0
println()
println("Frecuencia central: $(Bandreject_freq) Hz")
println("Orden del filtro: $(Bandreject_order)")
println("Ancho de banda: $(Bandreject_bandwidth) Hz")
println()
end

# ╔═╡ 8ffaf690-abeb-4228-84b8-f31aa8b0fade
begin
# Inicializamos el diccionario para las señales filtradas
dict_EEG_Bandreject = Dict{String, Vector{Float64}}()
end

# ╔═╡ 1f79a883-9fc4-49d6-8119-d6f1344a0836
begin
# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_bandreject = bandreject_filter(dict_EEG_Notch[first(keys(dict_EEG_Notch))], Bandreject_freq, fs, Bandreject_order, Bandreject_bandwidth)

# Aplicamos el filtro Bandreject a todos los canales (usando dict_EEG_Notch como entrada)
println("Aplicando filtro Bandreject a todos los canales...")
for (channel, signal) in dict_EEG_Notch
    y_bandreject, _ = bandreject_filter(signal, Bandreject_freq, fs, Bandreject_order, Bandreject_bandwidth)
    dict_EEG_Bandreject[channel] = y_bandreject
    println("  ✓ Canal $(channel) filtrado")
end
println("✓ Filtro Bandreject aplicado a $(length(dict_EEG_Bandreject)) canales")
println()
end

# ╔═╡ ded8a18d-8a0f-4cbe-8eac-4e0eae3f1b31
begin
# Ploteamos la respuesta del filtro Bandreject
p_bandreject, _, _ = plot_filter_response(flt_bandreject, fs; title = "Filtro Bandreject (100 Hz)")
p_bandreject
end

# ╔═╡ 78673555-9272-4ac4-9236-9724e0386ad8
begin
# Comparar PSD promedio del Notch vs Bandreject filtrado (superpuestos)
p_psd_comparison_br = compare_PSD_averages(
    dict_EEG_Notch, 
    dict_EEG_Bandreject, 
    fs;
    label_1 = "Notch (50 Hz)",
    label_2 = "Notch + Bandreject (100 Hz)",
    title = "Comparación PSD Promedio: Notch vs Notch+Bandreject"
)
end

# ╔═╡ 59c6cbf2-9a13-45c5-9d65-c60ba88a8707
save_filtered_signals(dict_EEG_Bandreject, "dict_EEG_Bandreject.bin")

# ╔═╡ 8de2f84f-2f33-44b7-b9cb-4a88fd6f8f09
md"""
## Filtro Pasa-Alto

Elimina la deriva lenta preservando actividad de baja frecuencia ($0.5\,\mathrm{Hz}$)

**Diseño:**

- orden: $4$
- orden efectivo: $8$ con `filtfilt`

El punto de corte en $0.5\,\mathrm{Hz}$ se seleccionó como un compromiso entre reducir la deriva del baseline y conservar la señal fisiológica lenta. Esta elección permite mantener actividad **infra-lenta** potencialmente relevante en **EM (Esclerosis Múltiple)**, asociada a procesos como inflamación o reorganización cortical.

Aunque se ha reportado que ICA puede beneficiarse de filtros pasa-alto ≥1 Hz, se priorizó un umbral más bajo para evitar la pérdida de información clínicamente significativa en bajas frecuencias.
"""

# ╔═╡ 8eb74e92-e770-4bb6-80ab-38a12c59e0c2
begin
println("=" ^ 40)
println("📊 Filtro Highpass (0.5 Hz)")
println("=" ^ 40)
	
Highpass_cutoff = 0.5
Highpass_order_design = 4          				# Orden del diseño del filtro
Highpass_order_effective = 8       				# Orden efectivo (4 × 2 con filtfilt)
	
println("Frecuencia de corte: $(Highpass_cutoff) Hz")
println("Orden del diseño del filtro: $(Highpass_order_design)")
println("Orden efectivo (con filtfilt): $(Highpass_order_effective)")
println()
end

# ╔═╡ e4dba693-d8bf-4501-94e0-2e258026f006
begin
# Inicializamos el diccionario para las señales filtradas
dict_EEG_Highpass = Dict{String, Vector{Float64}}()
end

# ╔═╡ d9b01b57-e43f-4f08-865f-0a68f99c5c69
begin
# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_highpass = highpass_filter(dict_EEG_Bandreject[first(keys(dict_EEG_Bandreject))], Highpass_cutoff, fs, Highpass_order_design)

# Aplicamos el filtro Highpass a todos los canales (usando dict_EEG_Bandreject como entrada)
println("Aplicando filtro Highpass a todos los canales...")
for (channel, signal) in dict_EEG_Bandreject
    # zero_phase = true por defecto, así que no necesitamos especificarlo explícitamente
    y_highpass, _ = highpass_filter(signal, Highpass_cutoff, fs, Highpass_order_design)
    dict_EEG_Highpass[channel] = y_highpass
    println("  ✓ Canal $(channel) filtrado")
end
println("✓ Filtro Highpass aplicado a $(length(dict_EEG_Highpass)) canales")
println()
end

# ╔═╡ 071e25e0-d109-4eea-b8f7-6fd1b50e05b4
begin
# Ploteamos la respuesta del filtro Highpass
p_highpass, _, _ = plot_filter_response(flt_highpass, fs; title = "Filtro Highpass (0.5 Hz)")
p_highpass
end

# ╔═╡ 6ebbe85d-c23a-4c92-8dbf-e127644df2b0
begin
# Comparar PSD promedio del Bandreject vs Highpass filtrado (superpuestos)
p_psd_comparison_hp = compare_PSD_averages(
    dict_EEG_Bandreject, 
    dict_EEG_Highpass, 
    fs;
    label_1 = "Notch + Bandreject",
    label_2 = "Notch + Bandreject + Highpass (0.5 Hz)",
    title = "Comparación PSD Promedio: Bandreject vs Bandreject+Highpass"
)
end

# ╔═╡ 86d3608c-d508-4751-b95f-b1a9dda4b863
begin
save_filtered_signals(dict_EEG_Highpass, "dict_EEG_Highpass.bin")
end

# ╔═╡ 445104ec-c077-4e84-86cc-d4b21813bf07
md"
## Filtro Paso-Bajo 

Limita el ancho de banda superior ($150\,\mathrm{Hz}$).

**Diseño:**

- orden: $4$
- orden efectivo: $8$ con `filtfilt`

Se eligió un corte en $150\,\mathrm{Hz}$ para conservar un rango amplio de información, incluyendo la banda **gamma**, y permitir una inspección inicial sin restricciones excesivas.

Este enfoque evita perder posibles componentes de gamma alta, que en algunos estudios sobre **EM (Esclerosis Múltiple)** se han relacionado con procesos de plasticidad cortical. Además, situar el corte lejos de las bandas de interés reduce el riesgo de artefactos de *ringing*.

Si el análisis posterior lo requiere, la señal puede restringirse sin inconvenientes a frecuencias más bajas (por ejemplo, $<50\,\mathrm{Hz}$).
"

# ╔═╡ cb505061-f6c3-48da-a970-100d98c38c79
begin
# -----------------------------------------------------------------------------------
# 7. FILTRO LOWPASS (150 Hz)
# -----------------------------------------------------------------------------------
# Limita el ancho de banda superior; preserva gamma para análisis posteriores.
println("=" ^ 40)
println("📊 Filtro Lowpass (150 Hz)")
println("=" ^ 40)
	
Lowpass_cutoff = 150
Lowpass_order_design = 4          				# Orden del diseño del filtro
Lowpass_order_effective = 8       				# Orden efectivo (4 × 2 con filtfilt)
	
println("Frecuencia de corte: $(Lowpass_cutoff) Hz")
println("Orden del diseño del filtro: $(Lowpass_order_design)")
println("Orden efectivo (con filtfilt): $(Lowpass_order_effective)")
println()
end

# ╔═╡ b7125ab8-d7cd-4ae9-b096-5a8f82c965cd
begin
# Inicializamos el diccionario para las señales filtradas
dict_EEG_Lowpass = Dict{String, Vector{Float64}}()
end

# ╔═╡ 7b84b771-3ef5-4213-ae36-43408b1c4f0c
begin
# Obtenemos el filtro una vez (es el mismo para todos los canales)
_, flt_lowpass = lowpass_filter(dict_EEG_Highpass[first(keys(dict_EEG_Highpass))], Lowpass_cutoff, fs, Lowpass_order_design)

# Aplicamos el filtro Lowpass a todos los canales (usando dict_EEG_Highpass como entrada)
println("Aplicando filtro Lowpass a todos los canales...")
for (channel, signal) in dict_EEG_Highpass
    # zero_phase = true por defecto, así que no necesitamos especificarlo explícitamente
    y_lowpass, _ = lowpass_filter(signal, Lowpass_cutoff, fs, Lowpass_order_design)
    dict_EEG_Lowpass[channel] = y_lowpass
    println("  ✓ Canal $(channel) filtrado")
end
println("✓ Filtro Lowpass aplicado a $(length(dict_EEG_Lowpass)) canales")
println()
end

# ╔═╡ 0b6f3dbf-57fc-4e62-b931-1f6ca745f4fb
begin
# Ploteamos la respuesta del filtro Lowpass
p_lowpass, _, _ = plot_filter_response(flt_lowpass, fs; title = "Filtro Lowpass (150 Hz)")
p_lowpass
end

# ╔═╡ a09a6827-9d97-4304-92dc-676322cf1028
begin
# Comparar PSD promedio del Highpass vs Lowpass filtrado (superpuestos)
p_psd_comparison_lp = compare_PSD_averages(
    dict_EEG_Highpass, 
    dict_EEG_Lowpass, 
    fs;
    label_1 = "Notch + Bandreject + Highpass",
    label_2 = "Notch + Bandreject + Highpass + Lowpass (150 Hz)",
    title = "Comparación PSD Promedio: Highpass vs Highpass+Lowpass"
)
end

# ╔═╡ b6f85f7f-400e-4f4a-b079-6228fc279483
begin
save_filtered_signals(dict_EEG_Lowpass, "dict_EEG_Lowpass.bin")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
InlineStrings = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
CSV = "~0.10.16"
DSP = "~0.8.4"
DataFrames = "~1.8.1"
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
project_hash = "ffe39d6236f9a44feaecf63064aafa557da8e1b2"

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
# ╟─fc2c1f3e-448d-4bd0-8208-4d531d1e5800
# ╠═4d8ca8b0-149f-4f6f-93f4-4e2d7f1a0b01
# ╠═d249d3f0-cf95-44e1-aec6-8ca4ad0fe102
# ╠═e857112f-df27-4f7f-bcc7-34421e0c3103
# ╠═e456bde5-04a7-4db8-a5d4-cf2068dfe29b
# ╠═204da7a3-9f60-4918-97f4-06a2d621dc15
# ╟─6337d051-3f9e-4dd9-9925-30b188eb6ac0
# ╠═a4097250-6b5a-4415-8469-e2a8dab637ec
# ╠═b1a100f5-2f4f-4d30-b4c4-1dca20c50e04
# ╠═aa74685e-d70c-4b92-99c6-9c8f051d387e
# ╠═2a1d5d68-6618-467d-9701-1167dd2d54cc
# ╠═557acfe8-89d9-4060-b851-b0fef0412631
# ╠═583ea341-04cf-4c52-ba27-5a1e9617d397
# ╠═d620dfe1-1628-4da7-8e1b-b456d659dbe8
# ╠═39690c4c-c5c8-4235-b832-4e78d8f72b3e
# ╠═722a26b7-1f85-489a-a9bb-66d8104f93ca
# ╠═3249aa95-e924-4e96-bf72-30c9a02c8607
# ╠═d359220a-c31c-4ce4-91e6-b58a7a9e7d08
# ╠═8ade48ae-34e1-462b-9c0e-206547e7432b
# ╠═dca7f5db-a29d-4028-aaef-8d51f4ebfb9c
# ╠═81c91017-c887-4aa0-b5f0-12436d83b0b3
# ╠═6e0cdb1b-d1cb-48ac-9167-c4d0c94ef398
# ╠═19a25445-09f4-4f9e-8281-dd16357f965b
# ╠═21ec0155-70ba-4535-943c-7f2fc48dbd1f
# ╠═2eb37338-f2ee-430e-8872-611374277c1b
# ╠═fa145f7e-4175-4aa9-ac93-8c1793b72b1b
# ╠═ca639a99-34d1-4fec-b800-8c7a03c45f7c
# ╠═e3d66280-fcca-4299-a59b-5187c4b05c57
# ╠═3cfdf522-7ab0-4ca4-be68-a5cfe7a421a9
# ╠═a78f7513-f36d-4854-80d9-c0f6d3afd3e4
# ╠═8ffaf690-abeb-4228-84b8-f31aa8b0fade
# ╠═1f79a883-9fc4-49d6-8119-d6f1344a0836
# ╠═ded8a18d-8a0f-4cbe-8eac-4e0eae3f1b31
# ╠═78673555-9272-4ac4-9236-9724e0386ad8
# ╠═59c6cbf2-9a13-45c5-9d65-c60ba88a8707
# ╠═8de2f84f-2f33-44b7-b9cb-4a88fd6f8f09
# ╠═8eb74e92-e770-4bb6-80ab-38a12c59e0c2
# ╠═e4dba693-d8bf-4501-94e0-2e258026f006
# ╠═d9b01b57-e43f-4f08-865f-0a68f99c5c69
# ╠═071e25e0-d109-4eea-b8f7-6fd1b50e05b4
# ╠═6ebbe85d-c23a-4c92-8dbf-e127644df2b0
# ╠═86d3608c-d508-4751-b95f-b1a9dda4b863
# ╠═445104ec-c077-4e84-86cc-d4b21813bf07
# ╠═cb505061-f6c3-48da-a970-100d98c38c79
# ╠═b7125ab8-d7cd-4ae9-b096-5a8f82c965cd
# ╠═7b84b771-3ef5-4213-ae36-43408b1c4f0c
# ╠═0b6f3dbf-57fc-4e62-b931-1f6ca745f4fb
# ╠═a09a6827-9d97-4304-92dc-676322cf1028
# ╠═b6f85f7f-400e-4f4a-b079-6228fc279483
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
