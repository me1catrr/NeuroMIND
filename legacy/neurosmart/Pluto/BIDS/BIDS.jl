### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 8d84f7cc-d13b-4b49-a2e0-cbf70d2d75a5
begin
	using PlutoUI
	PlutoUI.TableOfContents(title = "Contenido")
end

# ╔═╡ 62346384-1eef-4e5a-b7a0-e1b8af341cbe
begin
	include(joinpath(@__DIR__, "..", "_template_base.jl"))
	# Si este archivo se ejecuta directamente en Main, recarga utilidades de rutas
	# para evitar definiciones antiguas en la sesión REPL.
	if (@__MODULE__) == Main
		include(joinpath(@__DIR__, "..", "modules", "paths.jl"))
	end
	using .PlutoTemplateBase
	using CSV, DataFrames, Serialization, Statistics, StatsBase, DSP, Plots, Dates, StatsPlots, FFTW
end

# ╔═╡ 6a3f59f9-df9d-44a6-9a0e-bbc5a87d9057
notebook_intro("EEG-BIDS-RAW DATA")

# ╔═╡ 0f14b248-b15f-4fa4-b2f0-ff5f6609d66a
md"""
# Electroencefalografia (EEG)

La **electroencefalografía (EEG)** es una técnica no invasiva que registra la actividad eléctrica del cerebro mediante electrodos colocados sobre el cuero cabelludo. Estos electrodos se distribuyen siguiendo sistemas estandarizados, como el **sistema internacional 10–20**, que garantiza un muestreo espacial reproducible de la actividad cortical.

El conjunto de datos comprende registros de **EEG** en **estado de reposo** procedentes de dos grupos: **controles sanos** y **pacientes diagnosticados con esclerosis múltiple (EM)**. La pertenencia a cada grupo se definió de acuerdo con los criterios clínicos proporcionados en los metadatos asociados.

Las variables demográficas (**edad y sexo**) se registraron para todos los participantes y se incorporaron a la estructura de metadatos a nivel de sujeto, con el fin de permitir análisis estadísticos posteriores.

En este sistema, los electrodos se posicionan a distancias equivalentes al **10 % o 20 % de determinadas medidas craneales**, lo que permite una cobertura sistemática de la superficie del cráneo. Las etiquetas de los electrodos reflejan aproximadamente la región cortical subyacente:

- **F**: frontal  
- **C**: central  
- **T**: temporal  
- **P**: parietal  
- **O**: occipital  

Además, los **números impares** indican posiciones en el **hemisferio izquierdo**, los **números pares** posiciones en el **hemisferio derecho**, y la letra **z** denota posiciones en la **línea media**.

## Origen fisiológico de la señal EEG

Las señales registradas mediante EEG reflejan principalmente la **suma de potenciales postsinápticos** generados por grandes poblaciones de neuronas corticales que se activan de forma relativamente sincronizada. En particular, las contribuciones más relevantes provienen de **neuronas piramidales** orientadas perpendicularmente a la superficie cortical.

Estas neuronas generan **campos eléctricos dipolares** que se propagan a través del tejido cerebral, el líquido cefalorraquídeo, el cráneo y el cuero cabelludo mediante un proceso conocido como **conducción de volumen** (*volume conduction*).

Como consecuencia, el **EEG de superficie** no mide potenciales de acción individuales, sino una **actividad sináptica integrada espacial y temporalmente** que proviene de áreas corticales relativamente amplias. Esto implica:

- **Resolución espacial limitada**
- **Posible mezcla de señales entre canales**

## EEG en estado de reposo

El **EEG en estado de reposo** (*resting-state EEG*, **RS-EEG**) es un paradigma en el que la actividad cerebral se registra mientras el sujeto permanece relajado y no realiza ninguna tarea estructurada. Este enfoque captura las dinámicas oscilatorias basales y las fluctuaciones espontáneas de la actividad neuronal, evitando los factores de confusión asociados al diseño de tareas y a la carga motora o cognitiva.

En neurofisiología clínica, los registros de **RS-EEG** se utilizan ampliamente para:

- evaluar **ritmos de fondo**
- analizar la **reactividad cerebral**
- detectar **asimetrías**
- identificar **anomalías focales o difusas**

## Relevancia en esclerosis múltiple

En el contexto de la **esclerosis múltiple (EM)**, el RS-EEG resulta particularmente interesante porque la enfermedad afecta a **redes cerebrales distribuidas**. Estas alteraciones pueden manifestarse como cambios en los **patrones oscilatorios** y en la **conectividad funcional**, incluso en ausencia de tareas cognitivas explícitas.

Además, el EEG constituye una herramienta **relativamente simple, portátil y de bajo coste**, lo que lo convierte en un complemento útil a técnicas de neuroimagen como la **resonancia magnética (RM)** en estudios longitudinales o comparativos.

## Condiciones experimentales

- **Ojos abiertos (EO)**  
  mantiene un nivel moderado de activación cortical y reduce la dominancia del ritmo alfa.

- **Ojos cerrados (EC)**  
  potencia típicamente el **ritmo alfa**, siendo sensible a cambios dependientes del estado cerebral y potencialmente relacionados con la enfermedad.

La comparación entre **EO** y **EC** permite evaluar tanto efectos asociados al **nivel de alerta** como posibles **alteraciones neurofisiológicas** entre controles sanos y pacientes con EM.
"""

# ╔═╡ c5a6ef95-7fb6-4368-9d4d-6a44d95dd95e
md"""
## Parametros de adquisicion EEG

Las señales de **EEG** se registraron a partir de $C = 32$ electrodos de cuero cabelludo colocados según el **sistema internacional (SI) 10–20**. La disposición incluye **31 canales de EEG**, además del **electrodo de referencia REF (FCz)** y el **electrodo de tierra GRND (Fpz)**.

Las señales se muestrearon a una frecuencia de $f_s = 500\,\text{Hz}$, lo que produce

$N = f_s \cdot T = 50\,000$

muestras por canal y por condición.

Los datos se representan como matrices

$X_{s,k} \in \mathbb{R}^{C \times N},$

donde $s$ indexa a los sujetos y $k \in \{\text{EO}, \text{EC}\}$ denota la condición de registro. Cada fila corresponde a la serie temporal de un electrodo, expresada en microvoltios ($\mu V$).

## Informacion tecnica de adquisicion BrainVision

### Sistema de grabacion EEG

- **Amplificador**: actiCHamp Base Unit (5001) + modulo 32 CH
- **Software**: BrainVision Recorder Professional v. 1.21.0303
- **Formato**: Brain Vision Data Exchange Header File v1.0

### Parametros de grabacion

- **Canales**: 31 electrodos (sistema 10-20 internacional)
- **Frecuencia de muestreo**: 500 Hz
- **Resolucion**: 0.0488281 uV por unidad digital
- **Duracion**: 3-5 minutos por condicion (ojos abiertos/cerrados)
- **Filtros hardware**: DC-140 Hz
- **Filtros software**: 0.63-70 Hz + notch 50 Hz
"""

# ╔═╡ 97e1f81e-e8fc-4232-a7f4-0f145f2eae11
md"""
## Software BrainVision

**BrainVision Analyzer** es un software comercial ampliamente utilizado para el procesamiento y análisis de señales de **EEG**, desarrollado por **Brain Products**.

En este dataset, la adquisicion se realizo con **actiCHamp Base Unit (5001)** y modulo de **32 canales**, mediante **BrainVision Recorder Professional v. 1.21.0303**, utilizando el formato **Brain Vision Data Exchange Header File v1.0**.

### Estructura típica de archivos BrainVision

Está diseñado para trabajar de forma nativa con el **formato BrainVision**, el cual separa la información del registro en tres archivos complementarios:

```text
recording.vhdr
recording.eeg
recording.vmrk
```

- **`recording.eeg`**  
  Contiene los datos binarios, las señales EEG crudas, almacenadas como datos binarios.  
  Cada canal representa la evolución temporal del potencial eléctrico medido en el cuero cabelludo.

- **`recording.vhdr`**  
  Archivo de texto estructurado en secciones que incluye:
  - información del sistema de adquisición
  - número de canales
  - frecuencia de muestreo
  - configuración de referencia
  - formato del archivo binario
  - **referencias a los archivos `.eeg` y `.vmrk`**

- **`recording.vmrk`**  
  Archivo de texto con la lista de eventos experimentales (por ejemplo estímulos o inicio de ensayo).  
  Cada marcador suele incluir:
  - tipo de evento
  - posición temporal
  - duración
  - canal asociado (si aplica)


| Archivo | Tipo | Función |
|-------|------|--------|
| `.eeg` | Binario | Contiene la **serie temporal cruda del EEG** (voltaje continuo en $\mu V$). |
| `.vhdr` | Texto | Archivo de **cabecera** con información sobre el registro: número de canales, frecuencia de muestreo, referencias, formato de datos y enlaces a los otros archivos. |
| `.vmrk` | Texto | Archivo de **marcadores de eventos** que describe los tiempos y tipos de eventos experimentales. |

Esta estructura modular permite una gestión clara de la información y facilita su integración en flujos de análisis reproducibles.

Uno de los aspectos clave de este formato es la representación de los datos en el archivo `.eeg`, donde las señales pueden almacenarse en distintos tipos numéricos. Entre ellos, el formato **`IEEE_FLOAT_32`** (coma flotante de 32 bits) es especialmente relevante, ya que permite una alta precisión en la representación de la señal **EEG**, evitando pérdidas de información durante el registro y el procesamiento. Este nivel de precisión es importante para análisis sensibles, como la detección de componentes, el filtrado o el modelado de fuentes.

El archivo `.vhdr` especifica, entre otros parámetros, el tipo de datos (`BinaryFormat=IEEE_FLOAT_32`), la frecuencia de muestreo, el número de canales y las referencias a los archivos asociados. Por su parte, el archivo `.vmrk` contiene la información temporal de los eventos experimentales, esencial para segmentar y analizar la señal en función del diseño del experimento.

En conjunto, BrainVision Analyzer no solo proporciona herramientas avanzadas de procesamiento (filtrado, ICA, análisis espectral, entre otros), sino que también se apoya en un formato de datos transparente y bien documentado. Esto facilita su compatibilidad con estándares como **EEG-BIDS** y su integración con otras herramientas de análisis dentro del ecosistema de neurociencia computacional.

Estos tres archivos forman un **conjunto inseparable**. El archivo `.vhdr` actúa como punto de entrada para la mayoría del software de análisis, ya que describe cómo interpretar el archivo binario `.eeg` y dónde encontrar los eventos en `.vmrk`.

---

### Enlaces internos entre archivos

Dentro de los archivos `.vhdr` y `.vmrk` existen **campos que apuntan a los otros archivos del conjunto**.

Ejemplo típico dentro de un archivo `.vhdr`:

```text
DataFile=recording.eeg
MarkerFile=recording.vmrk
```

En el archivo `.vmrk` también aparece una referencia al archivo de datos:

```text
DataFile=recording.eeg
```

Estos enlaces permiten que los programas de análisis localicen correctamente los datos y los eventos.
"""

# ╔═╡ 3c6fd072-198c-410a-a771-85bc8e8ac5a3
md"""
## BIDS

El **Brain Imaging Data Structure (BIDS)** es un estándar desarrollado por la comunidad de neuroimagen para organizar, describir y compartir datos experimentales de forma consistente entre laboratorios y herramientas de análisis. 

Su objetivo principal es reducir la heterogeneidad en los formatos y estructuras de datos, facilitando su **reutilización, interoperabilidad y reproducibilidad** en estudios científicos.

En el caso del **EEG**, la extensión **EEG-BIDS** adapta este estándar a las particularidades de los registros electroencefalográficos, definiendo convenciones claras para nombres de archivos, organización por sujeto/sesión y el uso de archivos de metadatos (`.json`, `.tsv`). Además, promueve el uso de formatos ampliamente soportados como **BrainVision** (`.vhdr`, `.eeg`, `.vmrk`), lo que permite integrar datos provenientes de distintos sistemas en un flujo de trabajo común.

---

### Validación de datasets con BIDS Validator

El **BIDS Validator** es una herramienta diseñada para comprobar automáticamente si un dataset cumple con la especificación **BIDS**. Su función principal es detectar errores en la estructura de carpetas, nombres de archivos y metadatos, así como identificar información faltante o inconsistente, lo que ayuda a garantizar que los datos sean interoperables y reutilizables.

El validador puede ejecutarse tanto desde la línea de comandos como directamente en el navegador, lo que facilita su uso en diferentes entornos de trabajo. Al utilizarlo durante el proceso de organización de datos, es posible corregir problemas de forma temprana y asegurar la compatibilidad con herramientas de análisis como **BrainVision Analyzer**, **EEGLAB** o **MNE-Python**.

Puedes acceder al validador en el siguiente enlace:  
[BIDS Validator](https://bids-standard.github.io/bids-validator/)

---

### Introducción a BIDS y su uso con BrainVision Analyzer

Al adaptar un dataset al estándar **BIDS (Brain Imaging Data Structure)**, no basta con renombrar los archivos siguiendo la convención establecida. En el caso del formato **BrainVision** (`.vhdr`, `.eeg`, `.vmrk`), es necesario también actualizar las **referencias internas** contenidas en los archivos de texto, ya que estos incluyen enlaces explícitos entre sí. Si los nombres cambian pero estas referencias no se actualizan, el software de análisis no podrá localizar correctamente los datos asociados, generando errores en la carga o interpretación del dataset.

La integración de **BIDS** con herramientas como **BrainVision Analyzer** es especialmente natural, dado que este software utiliza el formato BrainVision, uno de los formatos recomendados dentro de EEG-BIDS. No obstante, convertir un dataset a BIDS implica más que reorganizar archivos: supone estructurar la información de forma coherente, estandarizar metadatos y asegurar la consistencia entre todos los componentes del registro.

De este modo, **BIDS** no solo facilita el uso de los datos en BrainVision Analyzer, sino que también permite su interoperabilidad con otras herramientas ampliamente utilizadas como EEGLAB, FieldTrip o MNE. En conjunto, BIDS actúa como una capa de estandarización que favorece la reproducibilidad, mientras que BrainVision Analyzer se integra como una herramienta de procesamiento dentro de un ecosistema más amplio y compatible.

---

### Problema al convertir a BIDS

En el estándar **BIDS-EEG**, los archivos deben renombrarse siguiendo una convención estructurada, por ejemplo:

```text
sub-01_ses-01_task-rest_run-01_eeg.vhdr
sub-01_ses-01_task-rest_run-01_eeg.eeg
sub-01_ses-01_task-rest_run-01_eeg.vmrk
```

Sin embargo, cuando se cambian los nombres de archivo:

- los **enlaces internos dentro de `.vhdr` y `.vmrk` siguen apuntando al nombre antiguo**
- esto provoca que el software no encuentre los archivos correspondientes.

---

### Actualización necesaria de enlaces

Después de renombrar los archivos, se deben actualizar las referencias internas.

En el archivo `.vhdr`:

```text
DataFile=sub-01_ses-01_task-rest_run-01_eeg.eeg
MarkerFile=sub-01_ses-01_task-rest_run-01_eeg.vmrk
```

En el archivo `.vmrk`:

```text
DataFile=sub-01_ses-01_task-rest_run-01_eeg.eeg
```

---

### Importancia para el análisis de EEG

Actualizar estos enlaces garantiza que:

- el software compatible con **BrainVision** (por ejemplo EEGLAB, FieldTrip o MNE)
  pueda cargar correctamente los datos
- los **marcadores de eventos** se alineen correctamente con la señal EEG
- los **validadores BIDS** puedan interpretar el dataset sin errores

En resumen, la actualización de los enlaces internos asegura la **consistencia estructural del dataset** y permite que las herramientas de análisis y validación localicen correctamente tanto la señal EEG como los eventos experimentales asociados.

Para cada sujeto y sesión, **BIDS** proporciona un archivo `electrodes.tsv` (por ejemplo, `sub-M05_ses-T2_electrodes.tsv`) que contiene las columnas:

(`name`, `x`, `y`, `z`, `type`)

Las coordenadas $(x, y, z)$ están **normalizadas** (por ejemplo, sobre una **esfera unitaria** o proporcionalmente a las dimensiones de la cabeza) y se utilizan para:

- generar **mapas topográficos**
- realizar **modelado de fuentes**

La **Tabla** enumera los nombres de los electrodos, sus coordenadas y el tipo (`EEG`, `REF` o `GND`) tal como aparecen en el conjunto de datos.

En esta fase de Input-Ouput **Fase (IO)**, donde las señales *EEG* registradas mediante el equipo comercial (en este caso **BrainVision**) se han organizado como un diccionario o matriz indexada por canales, junto con metadatos temporales y de canal, y estadísticas opcionales de calidad obtenidas durante la inspección inicial.
"""

# ╔═╡ cc99f17e-dd5f-4911-ac73-5264b3235940
begin
	# Rutas por etapa BIDS sin hardcodeo
	cfg_paths = @isdefined(cfg) ? cfg : nothing
	dir_raw_bids = raw_dir(cfg_paths)
	dir_build_eeg_bids = stage_dir(:build_eeg_bids; cfg = cfg_paths)
	dir_build_participants = stage_dir(:build_participants; cfg = cfg_paths)
	dir_validate_bids = stage_dir(:validate_bids; cfg = cfg_paths)

	(
		dir_raw = dir_raw_bids,
		dir_build_eeg_bids = dir_build_eeg_bids,
		dir_build_participants = dir_build_participants,
		dir_validate_bids = dir_validate_bids,
	)
end

# ╔═╡ 601ed72d-3245-4922-b887-e9919e1fdbe7
md"### Tabla de electrodos (name, x, y, z, type)"

# ╔═╡ d0ec074d-59d7-4f8a-9146-7fe9a4ca870e
begin
	cfg_electrodes = @isdefined(cfg) ? cfg : nothing
	electrodes_path = joinpath(electrodes_dir(cfg_electrodes), "sub-M05_ses-T2_electrodes.tsv")
	electrodes = if isfile(electrodes_path)
		CSV.read(electrodes_path, DataFrame; delim = '\t')
	else
		@warn "No se encontro electrodes.tsv; usando tabla vacia (modo publico/export)." electrodes_path
		DataFrame(name = String[], x = Float64[], y = Float64[], z = Float64[], type = String[])
	end
	electrodes
end

# ╔═╡ e13eb5ef-c053-4dca-9a96-8509068937df
md"
Tipos de archivos **BIDS**: ubicación, nombre y contenido principal

| Localización | Archivo | Descripción |
|----------|------|-------------|
| Root | `dataset_description.json` | Metadatos generales del dataset: nombre, autores, versión de BIDS y licencia. |
| Root | `participants.tsv` | Demografía tabular: `participant_id` (p. ej., `sub-M05`), `age`, `sex`, `group` (p. ej., `control`, `MS`). |
| Root | `participants.json` | Descripción y unidades para cada columna en `participants.tsv`. |
| Root | `README.md` | Documentación del dataset legible por humanos. |
| `sub-*/ses-*/eeg/` | `sub-*_ses-*_task-*_run-*_eeg.eeg` | Serie temporal EEG binaria cruda (voltaje continuo, $\mu V$). |
| `sub-*/ses-*/eeg/` | `*_eeg.vhdr` | Cabecera BrainVision: referencias a `.eeg` y `.vmrk`, información de canales y frecuencia de muestreo. Tras renombrar en BIDS, `DataFile=` y `MarkerFile=` deben referenciar los nuevos nombres base. |
| `sub-*/ses-*/eeg/` | `*_eeg.vmrk` | Archivo de marcadores BrainVision: tiempos de eventos y tipos. `DataFile=` debe apuntar al nuevo nombre del archivo `.eeg`. |
| `sub-*/ses-*/eeg/` | `*_eeg.json` | *Sidecar* BIDS-EEG con metadatos: `SamplingFrequency` (p. ej., 500), `PowerLineFrequency` (50), `EEGReference` (p. ej., `TP9-TP10` linked), `EEGGround` (p. ej., `Fpz`), `EEGChannelCount` (31), `RecordingType` (continuous), entre otros campos requeridos por BIDS-EEG. |
| `sub-*/ses-*/eeg/` | `*_events.tsv` | Eventos del *run*: columnas `onset`, `duration`, `trial_type` (p. ej., `rest`); tiempos en segundos (decimal). |
| `sub-*/ses-*/eeg/` | `*_events.json` | Descripción y unidades de las columnas presentes en `events.tsv`. |
"

# ╔═╡ 8c6ad877-f4e2-422a-a86f-5ae47f76f2ba
begin
	# Archivos base del root BIDS
	cfg_bids = @isdefined(cfg) ? cfg : nothing
	bids_base = bids_root(cfg_bids)
	dataset_description = joinpath(bids_base, "dataset_description.json")
	participants_tsv = joinpath(bids_base, "participants.tsv")
	participants_json = joinpath(bids_base, "participants.json")

	(
		dataset_description = dataset_description,
		dataset_description_exists = isfile(dataset_description),
		participants_tsv = participants_tsv,
		participants_tsv_exists = isfile(participants_tsv),
		participants_json = participants_json,
		participants_json_exists = isfile(participants_json),
	)
end

# ╔═╡ 0381f849-d27f-4208-a63f-08b7ed99e1db
begin
	participants_preview = if isfile(participants_tsv)
		CSV.read(participants_tsv, DataFrame; delim = '\t')
	else
		DataFrame(
			participant_id = String[],
			age = Union{Missing, Float64}[],
			sex = String[],
			group = String[],
		)
	end
	first(participants_preview, min(nrow(participants_preview), 10))
end

# ╔═╡ 11f6f202-6ca0-4f51-a6d3-e59aeb4e8130
md"
# Raw data (datos crudos)

Antes de cualquier etapa de preprocesamiento es recomendable realizar una **inspección básica del EEG crudo**.  
Esto permite detectar problemas de adquisición, artefactos evidentes o canales defectuosos que podrían afectar al análisis posterior.

Las siguientes representaciones gráficas y estadísticas son útiles para realizar un **control de calidad inicial (QC)** de los datos.
"

# ╔═╡ 5fddf73d-9b85-47f7-b8eb-c4788de33484
md"""
## Pipeline fase 0 (Raw data)

!!! info "Pipeline Fase 0: Ingesta de datos y estandarizacion de metadatos"

    **Objetivo:**
    - Cargar datos EEG crudos desde archivos `.tsv`.
    - Validar metadatos y consistencia dimensional.
    - Construir representacion estructurada por sujeto.
    - Generar visualizaciones y estadisticas de calidad por canal.
"""

# ╔═╡ 77df8d14-c178-44e8-982f-71b38ffaf77b
md"
## Carga datos raw`.tsv`

Se carga el archivo `.tsv` que contiene los datos **EEG raw**, en este caso solamente para un sujeto denominado **M05**, en condiciones de Ojos Cerrados **(EC)**: 

`sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv`

El formato esperado es:

- Primera columna: nombres de canales **(Channel)**
- Columnas siguientes: muestras temporales (una columna por punto de tiempo)
- Los datos están en unidades de microvoltios **(µV)**

Este archivo es la entrada principal del pipeline y contiene los datos sin ningún procesamiento previo **(raw data)**. El `.tsv` es el resultado de transformar los archivos de **BrainVision** al formato **BIDS**.
"

# ╔═╡ f9f9a9a3-a2c3-4a93-a3de-80de403f0d5f
begin
# -----------------------------------------------------------------------------------
# 1. CARGA DE DATOS RAW DESDE ARCHIVO TSV
# -----------------------------------------------------------------------------------
# Se carga el archivo TSV que contiene los datos EEG raw. El formato esperado es:
# - Primera columna: nombres de canales (Channel)
# - Columnas siguientes: muestras temporales (una columna por punto de tiempo)
# - Los datos están en unidades de microvoltios (µV)
#
# Este archivo es la entrada principal del pipeline y contiene los datos sin
# ningún procesamiento previo.

println("=" ^ 80)
println("📊 CARGA DE DATOS EEG")
println("=" ^ 80)

# Ruta al archivo de datos raw
raw_filename = "sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv"
cfg_local = @isdefined(cfg) ? cfg : nothing
dir_raw = joinpath(bids_root(cfg_local), "raw", raw_filename)

if !isfile(dir_raw)
    error("No se encontró archivo EEG raw: $(abspath(dir_raw)). Revisa cfg.data_dir y la estructura BIDS (raw_dir).")
end

# Leer el archivo TSV como DataFrame
# CSV.read carga el archivo y lo convierte en una estructura tabular
data_raw = CSV.read(dir_raw, DataFrame)   

println("✓ Archivo cargado: $(basename(dir_raw))")
println("✓ Dimensiones: $(size(data_raw, 1)) canales × $(size(data_raw, 2) - 1) muestras")
println("✓ Tabla de datos:")
display(data_raw)
println()
end

# ╔═╡ a9f6c35a-7e7b-4383-80ff-1df3d49f2bf6
md"
**Información temporal del registro EEG**

En este caso se extrae información sobre las características temporales del registro **EEG**, como son:

- Frecuencia de muestreo **(fs)**: número de muestras por segundo
- Número de muestras: longitud temporal del registro
- Duración total **(t seg)**: tiempo total de registro en segundos y minutos
"

# ╔═╡ 6dcf7d5d-fc68-4f21-8f40-406f173122f2
begin
# -----------------------------------------------------------------------------------
# 2. INFORMACIÓN TEMPORAL DEL REGISTRO
# -----------------------------------------------------------------------------------
# Se extrae información sobre las características temporales del registro:
# - Frecuencia de muestreo (fs): número de muestras por segundo
# - Número de muestras: longitud temporal del registro
# - Duración total: tiempo total de registro en segundos y minutos
#
# Esta información es esencial para:
# - Cálculo de frecuencias en análisis espectral
# - Aplicación de filtros (requieren fs)
# - Segmentación temporal
# - Visualización de señales en el dominio temporal

println("⏱️  INFORMACIÓN TEMPORAL")
println("-" ^ 80)

# Frecuencia de muestreo (Hz)
# Define la resolución temporal: cuántas muestras se toman por segundo
fs = (cfg_local !== nothing && hasproperty(cfg_local, :filter) && hasproperty(cfg_local.filter, :fs)) ?
    cfg_local.filter.fs : 500

# Número de muestras temporales
# Se resta 1 porque la primera columna es "Channel" (nombres de canales)
n_muestras = size(data_raw, 2) - 1

# Vector de tiempo en segundos
# Se genera desde 0 hasta (n_muestras-1) y se divide por fs para obtener segundos
tiempo_seg = collect(0:(n_muestras-1)) ./ fs

# Duración total del registro
duracion_total = tiempo_seg[end]

println("  Frecuencia de muestreo: $(fs) Hz")
println("  Número de muestras: $(n_muestras)")
println("  Duración total: $(round(duracion_total, digits=2)) segundos ($(round(duracion_total/60, digits=2)) minutos)")
println()
end

# ╔═╡ 1ccb18f5-b6f8-431a-a670-a94efea7e988
md"
**Extraccion y organizacion de canales**

Se extraen los nombres de los canales y se organizan los datos en estructuras
convenientes para el procesamiento, de forma siguiente:

- **channels**: vector con nombres de canales
- **EEG_matrix**: matriz (canales × muestras) con todos los datos
- **dict_EEG**: diccionario que mapea nombre de canal → señal (vector)

El diccionario es especialmente útil porque permite acceso por nombre de canal
y facilita el procesamiento posterior en el pipeline.
"

# ╔═╡ 8b676680-8e27-49f7-991f-c8dffef2f62d
begin
# -----------------------------------------------------------------------------------
# 3. EXTRACCIÓN Y ORGANIZACIÓN DE CANALES
# -----------------------------------------------------------------------------------
# Se extraen los nombres de canales y se organizan los datos en estructuras
# convenientes para el procesamiento:
# - channels: vector con nombres de canales
# - EEG_matrix: matriz (canales × muestras) con todos los datos
# - dict_EEG: diccionario que mapea nombre de canal → señal (vector)
#
# El diccionario es especialmente útil porque permite acceso por nombre de canal
# y facilita el procesamiento posterior en el pipeline.

println("📡 INFORMACIÓN DE CANALES")
println("-" ^ 80)

# Nombres de canales (primera columna del DataFrame)
channels = data_raw.Channel

# Matriz de datos EEG (canales × muestras)
# Se excluye la columna "Channel" usando Not(:Channel)
# Cada fila corresponde a un canal, cada columna a una muestra temporal
EEG_matrix = Matrix(data_raw[:, Not(:Channel)])

println("  Número de canales: $(length(channels))")
println("  Canales: $(join(channels, ", "))")
println()

# Diccionario: canal → señal (vector de muestras)
# Esta estructura facilita el acceso por nombre de canal y es compatible
# con el resto del pipeline que espera datos organizados de esta manera
dict_EEG = Dict(data_raw.Channel[i] => EEG_matrix[i, :]
                for i in 1:size(data_raw, 1))
end

# ╔═╡ 70094435-8e89-4d22-a9de-fd8f0b2ac76d
md"
**Guardado de datos organizados**

Se guarda el diccionario de canales en formato binario nativo de *Julia* (`dict_EEG.bin`). Este formato es eficiente para matrices grandes y permite serialización rápida de estructuras complejas. El archivo guardado será la entrada para el siguiente paso del pipeline (filtrado, corrección de línea base, etc.).
"

# ╔═╡ 66fa1086-568d-4eb8-a8ea-eb96f19dc6ef
begin
# -----------------------------------------------------------------------------------
# 4. GUARDADO DE DATOS ORGANIZADOS
# -----------------------------------------------------------------------------------
# Se guarda el diccionario de canales en formato binario nativo de Julia.
# Este formato es eficiente para matrices grandes y permite serialización
# rápida de estructuras complejas.
#
# El archivo guardado será la entrada para el siguiente paso del pipeline
# (filtrado, corrección de línea base, etc.)

# Directorio de salida para datos IO
dir_io = stage_dir(:IO; cfg = cfg_local)
path_dict = joinpath(dir_io, "dict_EEG.bin")

# Asegurar que el directorio existe
isdir(dir_io) || mkpath(dir_io)

# Serializar y guardar el diccionario
Serialization.serialize(path_dict, dict_EEG)

println("💾 Diccionario EEG raw guardado en: $(abspath(path_dict))")
end

# ╔═╡ fef42c56-e8e2-4ec2-9ad0-6bc9d1f37772
md"
## Control de calidad senales EEG

Es necesario generar visualizaciones exploratorias en el dominio temporal/frecuencial. Estos gráficos se generan a partir de los **datos raw**. Estas visualizaciones ayudan a identificar:

- Artefactos obvios (parpadeos, movimiento, ruido de línea).
- Diferencias entre canales.
- Calidad general del registro.

Algunas notas prácticas para control de calidad de las señales **EEG**:

- Un **canal plano** suele indicar desconexión del electrodo o fallo del amplificador.
- **Picos estrechos** en **50/60 Hz** indican interferencia de red eléctrica.
- **Exceso de potencia <1 Hz** puede deberse a mal contacto electrodo–piel o movimiento.
- En registros con **ojos cerrados (EC)**, suele observarse aumento de potencia en **banda alfa (8–12 Hz)** en regiones occipitales.
- Es recomendable marcar o excluir **canales sospechosos** antes de realizar análisis espectral o conectividad.
"

# ╔═╡ b4de8852-c51f-423f-be51-74f67f8404d2
md"
Tabla x: Graficas y estadisticas para la inspeccion y control de calidad de **EEG crudo**

| Representación o estadístico | Propósito y uso |
|---|---|
| **Trazas crudas** (EEG vs. tiempo, vista multicanal apilada o desplazable) | Visualizar el comportamiento temporal de los canales, identificar artefactos o canales planos, comprobar la continuidad de la señal y localizar eventos; primera verificación de la integridad y sincronía de la grabación. |
| **Mapa de impedancias** | Refleja la calidad del contacto electrodo–piel. Impedancias altas (>10–20 kΩ) suelen anticipar ruido o inestabilidad. Un mapa topográfico ayuda a detectar áreas problemáticas (p. ej., electrodos frontales secos). |
| **Histograma de amplitud por canal** | Permite comprobar la amplitud típica por canal y detectar saturación, *clipping* o ruido excesivo. Una distribución simétrica centrada en cero y sin colas largas suele indicar buena calidad. |
| **Topografía de varianza** | Resume la actividad media en el dominio espacial. Los canales con varianza (o desviación estándar) mucho mayor que el resto suelen estar contaminados por ruido o desconectados. |
| **PSD media cruda por canal** | El espectro de potencia sin filtrar revela interferencias eléctricas (picos estrechos a 50/100 Hz) o mal acoplamiento (exceso de potencia en frecuencias muy bajas). Útil para planificar el filtrado. |
| **Amplitud media, varianza o desviación estándar por canal** | Caracteriza el nivel de actividad basal y su variabilidad; permite comparar condiciones (p. ej., ojos abiertos vs. ojos cerrados) y detectar canales ruidosos o planos. |
| **Potencia en bandas absoluta o relativa** | Cuantifica la potencia en bandas canónicas (δ, θ, α, β, γ); caracteriza el perfil espectral y facilita la comparación entre condiciones y la detección de artefactos. |
| **Implementación** | Los pasos anteriores (carga de TSV, extracción temporal y de canales, gráficos temporales crudos y apilados, PSD y potencia por bandas, estadísticas por canal y marcado de canales sospechosos) están implementados en la rutina de Julia `src/IO.jl`. |
"

# ╔═╡ f3cfdba1-8641-4c5e-b9e4-f3f1138859ea
md"
**Grafico de un canal individual** (ejemplo: canal Cz). El canal Cz (central) suele ser representativo de la actividad EEG general
"

# ╔═╡ 286ace1d-be59-4a0b-8b73-2ba751129850
begin
p_cz = plot(tiempo_seg, dict_EEG["Cz"], 
xlabel = "Tiempo (s)", ylabel= "Amplitud (µV)", title="EEG - Canal Cz", legend = false)
p_cz
end

# ╔═╡ 2a92ae15-c2c4-4b2c-8523-7f63475d5057
md"
En la representación de todos los **canales apilados**, cada canal se desplaza verticalmente para facilitar la comparación visual. Esta visualización permite detectar artefactos que afectan a múltiples canales. Una separación vertical entre canales (en µV) mayor aumenta la separación visual pero puede hacer que los detalles
de amplitud sean menos visibles.
"

# ╔═╡ ab98be9e-c06f-4a9f-929f-0ab2b8ce6edf
begin
# Separación vertical entre canales (en µV)
# Un valor mayor aumenta la separación visual pero puede hacer que los detalles
# de amplitud sean menos visibles
sep = 80.0

# Calcular offsets (posiciones Y) para cada canal
# Los canales se apilan de arriba a abajo, con el primer canal en la parte superior
offsets = [(length(channels) - i) * sep for i in 1:length(channels)]

# Gráfico apilado
# Se transpone EEG_matrix para que cada columna sea un canal
# Se suman los offsets para desplazar verticalmente cada señal
p = plot(
    tiempo_seg,
    EEG_matrix' .+ offsets',  # transpuesta + offsets para apilar
    xlabel = "Tiempo (s)",
    ylabel = "",
    yticks = (offsets, channels),  # etiquetas Y con nombres de canales
    legend = false,
    grid = false,
    size = (1000, 600)
)
p
end

# ╔═╡ 6668e05c-ec9f-44cb-a66f-ccb41b2f44dd
md"
En el **análisis estadístico** de control de calidad de las señales de los canales **EEG** se calculan estadísticas descriptivas para cada canal que permiten evaluar la calidad de los datos y detectar canales problemáticos:

- **Mean**: valor medio (indica offset DC)
- **Min/Max**: rango de amplitudes
- **Std**: desviación estándar (variabilidad, ruido)
- **RMS**: raíz cuadrada de la media de cuadrados (potencia promedio)
- **Kurtosis**: medida de picos (valores altos indican outliers/artefactos)
- **Skewness**: asimetría de la distribución

De manera que los canales con **kurtosis** muy alta (>6) o desviación estándar muy baja (<5 µV) pueden indicar problemas (artefactos, canales muertos, etc.).
"

# ╔═╡ 2e4c5804-bdc9-4f53-bbc4-b3df6f82c145
begin
# -----------------------------------------------------------------------------------
# 8. ANÁLISIS ESTADÍSTICO DE CALIDAD DE CANALES
# -----------------------------------------------------------------------------------
# Se calculan estadísticas descriptivas para cada canal que permiten evaluar
# la calidad de los datos y detectar canales problemáticos:
# - Mean: valor medio (indica offset DC)
# - Min/Max: rango de amplitudes
# - Std: desviación estándar (variabilidad, ruido)
# - RMS: raíz cuadrada de la media de cuadrados (potencia promedio)
# - Kurtosis: medida de "picos" (valores altos indican outliers/artefactos)
# - Skewness: asimetría de la distribución
#
# Canales con kurtosis muy alta (>6) o desviación estándar muy baja (<5 µV)
# pueden indicar problemas (artefactos, canales muertos, etc.)


# Calcular estadísticas para cada canal
# eachrow itera sobre cada fila de la matriz (cada canal)
stats_channels = DataFrame(
    :Channel                 => channels,
    Symbol("Mean (µV)")     => mean.(eachrow(EEG_matrix)),      # Media
    Symbol("Min (µV)")      => minimum.(eachrow(EEG_matrix)),   # Mínimo
    Symbol("Max (µV)")      => maximum.(eachrow(EEG_matrix)),  # Máximo
    Symbol("Std (µV)")      => std.(eachrow(EEG_matrix)),       # Desviación estándar
    Symbol("RMS (µV)")      => rms.(eachrow(EEG_matrix)),       # RMS (root mean square)
    :Kurtosis                => kurtosis.(eachrow(EEG_matrix)), # Curtosis (medida de picos)
    :Skewness                => skewness.(eachrow(EEG_matrix)), # Asimetría
)
stats_channels
end

# ╔═╡ 8b984de8-8ec4-46ac-adf0-0a2443774f24
md"
En un **histograma de la desviación estándar** **Std(μV)** por canales, se puede interpretar la variabilidad de los mismos. Valores muy bajos **(<5 µV)** pueden indicar canales muertos o con muy poca señal.
"

# ╔═╡ 3083f7e8-b7f8-4f91-a35d-b6720fc03fdc
begin
# Histograma de desviación estándar
# La desviación estándar mide la variabilidad. Valores muy bajos (<5 µV)
# pueden indicar canales muertos o con muy poca señal
histogram(
    stats_channels[!, Symbol("Std (µV)")],
    bins   = 20,
    xlabel = "Std (µV)",
    ylabel = "Número de canales",
    legend = false,
    title  = "Distribución de desviación típica por canal"
) 
end

# ╔═╡ 530abe8f-3dc1-424f-a04c-b6a57e723d2f
md"
En un **histograma de kurtosis** por canales, se puede ver qué tan **picuda** es la distribución. Como se ha dicho antes, valores altos **(>6)** pueden indicar presencia de artefactos (picos, spikes) o canales problemáticos.
"

# ╔═╡ 6a5025f6-c6a7-43cb-9cdd-b6ea8cd3935e
begin
# Histograma de kurtosis
# La kurtosis mide qué tan "picuda" es la distribución. Valores altos (>6)
# pueden indicar presencia de artefactos (picos, spikes) o canales problemáticos
histogram(
    stats_channels.Kurtosis,
    bins = 20,
    xlabel = "Kurtosis",
    ylabel = "Número de canales",
    legend = false,
    title  = "Distribución de kurtosis por canal"
)
end

# ╔═╡ ae0949bc-70b0-4cf2-9ab9-b5dc33e9b5ad
md"
**Detección de canales sospechosos** Se usa un gráfico de dispersión (Std vs Kurtosis) para identificar canales que se desvían de la distribución normal, lo cual puede indicar problemas.
"

# ╔═╡ 6ea65e96-4001-4795-b1f0-052db14cec52
begin
# -----------------------------------------------------------------------------------
# 8.2. DETECCIÓN DE CANALES SOSPECHOSOS
# -----------------------------------------------------------------------------------
# Se usa un gráfico de dispersión (Std vs Kurtosis) para identificar canales
# que se desvían de la distribución normal, lo cual puede indicar problemas.

# Extraer variables para el scatter plot
stds   = stats_channels[!, Symbol("Std (µV)")]
kurts  = stats_channels.Kurtosis
names_ = stats_channels.Channel

# Criterio de canales sospechosos:
# - Kurtosis > 6: distribución muy picuda (posibles artefactos)
# - Std < 5 µV: variabilidad muy baja (posible canal muerto)
idx_raros = (kurts .> 6) .| (stds .< 5)

# Colores: rojo para sospechosos, azul para normales
colors = ifelse.(idx_raros, :red, :blue)

# Crear scatter plot
p_scatter = scatter(
    stds, kurts;
    color  = colors,
    xlabel = "Std (µV)",
    ylabel = "Kurtosis",
    title  = "Calidad de canales EEG (Std vs Kurtosis)",
    legend = false,
)

# Añadir etiquetas SOLO a los canales sospechosos
# Se desplazan ligeramente para mejor legibilidad
offsetx = 1.5   # desplazamiento en X
offsety = 0.2   # desplazamiento en Y

for i in eachindex(names_)
    if idx_raros[i]
        annotate!(
            p,
            stds[i] + offsetx,
            kurts[i] + offsety,
            Plots.text(string(names_[i]), 6, :red)
        )
    end
end
p_scatter
end

# ╔═╡ 13177dee-e72f-4997-a204-e110ef7ea60c
md"
## Análisis Espectral

El **análisis espectral** constituye una de las herramientas fundamentales en el estudio de señales **EEG**, ya que permite descomponer la actividad cerebral en sus distintas componentes de frecuencia. A diferencia del análisis en el dominio temporal, el enfoque espectral facilita la identificación de patrones rítmicos que están directamente relacionados con **estados cognitivos**, **procesos fisiológicos** y posibles **alteraciones patológicas**.

En este contexto, la estimación de la **densidad de potencia espectral (PSD)** proporciona una medida cuantitativa de cómo se distribuye la **energía de la señal** a lo largo de las **frecuencias**. Esto resulta clave para caracterizar la dinámica cerebral, comparar condiciones experimentales y evaluar cambios en la actividad neuronal de forma objetiva.

El **análisis espectral** es ampliamente utilizado en múltiples aplicaciones, como el estudio de estados de vigilia y sueño, la detección de anomalías (por ejemplo, actividad epiléptica), el análisis de ritmos asociados a funciones cognitivas (atención, memoria) y el desarrollo de interfaces cerebro-computador (BCI). Además, el **cálculo de potencia en bandas específicas** permite resumir la información de la señal en métricas interpretables y comparables entre sujetos o condiciones.

En conjunto, estas herramientas convierten al **análisis espectral** en un componente esencial dentro de cualquier pipeline de procesamiento de **EEG**, proporcionando una base sólida para la interpretación funcional de la actividad cerebral.
"

# ╔═╡ 23b2f8f8-0133-46e9-a5d2-0c00aab6d2a4
md"
### Bandas espectrales

Los análisis espectrales y de conectividad se realizaron dentro de **siete bandas de frecuencia predefinidas**. Estas bandas se seleccionaron en función de las convenciones neurofisiológicas estándar y de su relevancia para la dinámica del **EEG en estado de reposo**:

| Banda | Rango de frecuencia | Interpretación neurofisiológica |
|------|---------------------|--------------------------------|
| $\delta$ | $0.5$–$4\ \text{Hz}$ | Actividad cortical lenta asociada al sueño profundo; incrementos anómalos en vigilia pueden indicar disfunción cortical difusa. |
| $\theta$ | $4$–$8\ \text{Hz}$ | Relacionada con somnolencia, estados de baja vigilancia y procesos de memoria. |
| $\alpha$ | $8$–$12\ \text{Hz}$ | Ritmo dominante en regiones occipitales durante reposo con **ojos cerrados**; marcador clásico del estado basal del EEG. |
| $\beta_{\text{low}}$ | $12$–$15\ \text{Hz}$ | Asociada a actividad sensoriomotora y a estados de alerta moderada. |
| $\beta_{\text{mid}}$ | $15$–$18\ \text{Hz}$ | Relacionada con procesos atencionales y actividad cortical sostenida. |
| $\beta_{\text{high}}$ | $18$–$30\ \text{Hz}$ | Vinculada a procesamiento cognitivo activo y actividad motora; puede aumentar con tensión muscular. |
| $\gamma$ | $30$–$50\ \text{Hz}$ | Oscilaciones rápidas asociadas al procesamiento cortical de orden superior e integración neuronal; susceptible a artefactos musculares. |

La **banda δ** refleja actividad cortical lenta y puede indicar disfunción difusa cuando aparece anormalmente elevada.  
El rango **θ** se asocia con somnolencia y con procesos relacionados con la memoria.  
El **ritmo α** domina en regiones posteriores durante el reposo con **ojos cerrados** y constituye un marcador clásico del estado de reposo.

Las **subbandas β** capturan procesos sensoriomotores y atencionales a frecuencias progresivamente más altas.  
Por último, la **banda γ** refleja actividad oscilatoria rápida vinculada al procesamiento cortical de orden superior, aunque es más susceptible a artefactos musculares.
"

# ╔═╡ 98c1de45-f9cf-43ee-b903-5f6ca9897269
begin
# Nomenclatura de bandas alineada con Pluto/BIDS/BIDS.jl
Δ = (0.5, 4.0)
δ = (4.0, 8.0)
α = (8.0, 12.0)
β_low = (12.0, 15.0)
β_medium = (15.0, 18.0)
β_high = (18.0, 30.0)
γ = (30.0, 50.0)

band_limits = Dict(:Δ => Δ, :δ => δ, :α => α, :β_low => β_low, :β_medium => β_medium, :β_high => β_high, :γ => γ)
ordered_bands = [:Δ, :δ, :α, :β_low, :β_medium, :β_high, :γ]
band_labels = Dict(
    :Δ => "DELTA (0.5-4)", :δ => "THETA (4-8)", :α => "ALPHA (8-12)",
    :β_low => "BETA LOW (12-15)", :β_medium => "BETA MID (15-18)",
    :β_high => "BETA HIGH (18-30)", :γ => "GAMMA (30-50)"
)

# Colores para visualización de bandas en gráficos
band_colors = [:lightcyan, :lavender, :lightgoldenrod, :lightgreen, :lightsalmon, :lightpink, :lightgray]
end

# ╔═╡ e2643d97-07f8-4fd4-b4ad-9e36b340f6e6
begin
PSD = Dict(channel => begin
p = welch_pgram(dict_EEG[channel]; fs=fs, window=hamming, nfft=n_muestras)
# Extraer frecuencias y potencias del periodograma
(; freq = DSP.freq(p), power = DSP.power(p))
end for channel in channels)
PSD
end

# ╔═╡ 8f52ce15-6075-4c01-8361-4c18b4aae41a
md"
### Densidad de Potencia Espectral (PSD)

Se calcula la **densidad de potencia espectral (PSD)** usando el **método de Welch**, que divide la señal en segmentos superpuestos, aplica ventanas (**Hamming**) y promedia los periodogramas.

La **PSD** muestra cómo se distribuye la potencia de la señal en el dominio frecuencial, lo cual es esencial para:

- Identificar bandas de frecuencia dominantes (α, β, etc.)
- Detectar ruido de línea (50/60 Hz)
- Evaluar la calidad espectral de cada canal
- Comparar características espectrales entre canales

Los parámetros fundamentales para calcular la **PSD** son los siguientes:

- **fs**: frecuencia de muestreo (necesaria para escalar frecuencias)
- **window**: tipo de ventana (**Hamming** reduce leakage espectral)
- **nfft**: tamaño de **FFT** (usar n_muestras para máxima resolución)

---

### Visualización de la PSD

Se superponen todas las **PSD**, canales superpuestos, para comparar características espectrales entre canales. La escala logarítmica facilita visualizar rangos amplios de potencia.
"

# ╔═╡ 50bb2eb0-d367-4f48-a8db-c5d4db35e7f4
begin
# Calcular límites Y amigables para escala logarítmica
# Se buscan potencias de 10 que enmarquen todos los datos
all_powers = vcat([PSD[ch].power for ch in channels]...)
all_powers_pos = filter(>(0), all_powers)
if isempty(all_powers_pos)
    error("No hay potencias positivas para graficar en escala log10.")
end
power_min_all, power_max_all = extrema(all_powers_pos)
if power_max_all <= power_min_all
    power_max_all = power_min_all * 10
end
y_min_log_all = floor(log10(power_min_all))
y_max_log_all = ceil(log10(power_max_all))
ylim_log_all = (10.0^y_min_log_all, 10.0^y_max_log_all)

# Gráfico de PSD superpuestas
p_psd = plot(
    xlabel = "Frecuencia (Hz)",
    ylabel = "Potencia (µV²/Hz)",
    title = "PSD por canal (Welch)",
    legend = :right,
    xlim = (0, fs / 2),  # Frecuencia de Nyquist = fs/2
    yscale = :log10,     # Escala logarítmica para mejor visualización
    ylim = ylim_log_all,
)

# Añadir cada canal al gráfico
for channel in channels
    plot!(p_psd, PSD[channel].freq, PSD[channel].power; label = channel, lw = 1)
end
p_psd
end

# ╔═╡ ad3ebb32-88af-4f4d-8aac-e35f4ee76aee
md"
### Cálculo de potencia por banda

El **cálculo de la potencia integrada por banda** se realiza a partir del cálculo de la potencia total en cada banda de frecuencia para cada canal.
La potencia se integra sumando las densidades espectrales dentro de cada banda y multiplicando por la **resolución espectral (df)**.

Esto permite cuantificar:

- Qué bandas dominan en cada canal
- Distribución relativa de potencia entre bandas
- Comparaciones entre canales
"

# ╔═╡ c4d34a8d-ad3d-44f6-934d-cf5c14670440
begin
# Calcular potencia integrada por banda para cada canal
# La integración se hace sumando las densidades de potencia en cada banda
# y multiplicando por la resolución espectral (df = diferencia entre frecuencias)
band_powers = Dict(channel => begin
    freqs = PSD[channel].freq
    powers = PSD[channel].power
    
    # Resolución espectral: diferencia entre frecuencias adyacentes
    # Esto es necesario para convertir densidad (µV²/Hz) en potencia total (µV²)
    df = freqs[2] - freqs[1]
    
    # Para cada banda, sumar las potencias dentro del rango de frecuencias
    Dict(
        band => begin
            # Índices de frecuencias dentro de la banda [f_min, f_max)
            idx = (freqs .>= limits[1]) .& (freqs .< limits[2])
            # Potencia integrada = suma de densidades × resolución
            sum(powers[idx]) * df
        end
        for (band, limits) in band_limits
    )
end for channel in channels)
end

# ╔═╡ 4d28c8d4-c1fc-4c7d-bcc8-09e8f568b0f6
begin
# Crear DataFrame con potencias por banda y canal
band_powers_df = DataFrame(Band = [band_labels[band] for band in ordered_bands])
channel_totals = Dict{String, Float64}()

# Añadir columnas para cada canal (potencia absoluta y porcentaje)
for channel in channels
    values = [band_powers[channel][band] for band in ordered_bands]
    total_power = sum(values)
    channel_totals[channel] = total_power
    
    # Potencia absoluta en µV²
    band_powers_df[!, Symbol(channel)] = values
    
    # Potencia relativa (porcentaje del total)
    # Evita división por cero si total_power es muy pequeño
    band_powers_df[!, Symbol(channel * "_pct")] = total_power ≈ 0 ? zeros(length(values)) : (values ./ total_power)
end

# Añadir fila de totales
total_row_pairs = Any[:Band => "TOTAL POWER"]
for channel in channels
    push!(total_row_pairs, Symbol(channel) => channel_totals[channel])
    push!(total_row_pairs, Symbol(channel * "_pct") => 1.0)
end
total_row = (; total_row_pairs...)

# Combinar fila de totales con el resto del DataFrame
band_powers_df = vcat(DataFrame([total_row]), band_powers_df; cols = :union)

band_powers_df
end

# ╔═╡ 2962c98d-6b35-42cf-9dfa-30ecebe4f38b
md"""
### PSD promedio de todos los canales

Visualización equivalente a la rutina `IO.jl`, con sombreado por bandas EEG para facilitar la interpretación fisiológica del espectro.
"""

# ╔═╡ 61ec2843-f706-41f3-9eb2-8bad3e1ea3ca
begin
# -----------------------------------------------------------------------------------
# 6.3. VISUALIZACIÓN DE PSD PROMEDIO CON BANDAS DE FRECUENCIA
# -----------------------------------------------------------------------------------
# Se calcula la PSD promedio de todos los canales y se visualiza con las bandas
# de frecuencia marcadas. Esto proporciona una visión general del contenido
# espectral del registro y ayuda a identificar qué bandas dominan.

# Crear gráfico base
p_psd_avg = plot(
    xlabel = "Frecuencia (Hz)",
    ylabel = "Potencia (µV²/Hz)",
    title = "PSD promedio de todos los canales",
    legend = :outerright,
    xlim = (0, fs / 2),
    yscale = :log10,
)

# Calcular PSD promedio
# Se asume que todos los canales tienen las mismas frecuencias (mismo nfft)
freqs_avg = PSD[first(channels)].freq

# Promediar potencias de todos los canales
# mapreduce suma todas las potencias y luego divide por el número de canales
avg_power = mapreduce(ch -> PSD[ch].power, +, channels) ./ length(channels)

# Calcular límites de potencia para escala Y
avg_power_pos = filter(>(0), avg_power)
if isempty(avg_power_pos)
    error("La PSD promedio no contiene valores positivos para escala log10.")
end
power_min, power_max = extrema(avg_power_pos)
if power_max <= power_min
    power_max = power_min * 10
end

# Calcular límites Y amigables para escala logarítmica (potencias de 10)
y_min_log = floor(log10(power_min))
y_max_log = ceil(log10(power_max))
ylim_log = (10.0^y_min_log, 10.0^y_max_log)

# Actualizar el gráfico con límites Y explícitos
plot!(p_psd_avg, ylim = ylim_log)

# Añadir regiones sombreadas para cada banda de frecuencia
# fillrange crea un área sombreada entre límites positivos en cada banda
for (idx, band) in enumerate(ordered_bands)
    limits = band_limits[band]
    plot!(
        p_psd_avg,
        [limits[1], limits[2]],      # límites de frecuencia en X
        [power_max, power_max];      # altura superior del sombreado
        fillrange = power_min,       # altura inferior del sombreado
        fillcolor = band_colors[idx], # color de la banda
        fillalpha = 0.25,            # transparencia
        linecolor = :transparent,    # sin borde
        label = band_labels[band],   # etiqueta para leyenda
    )
end

# Añadir línea de PSD promedio
plot!(p_psd_avg, freqs_avg, avg_power; label = "Promedio", lw = 2, color = :black)

p_psd_avg
end

# ╔═╡ 8ad9dcae-9ea1-4eb1-bf41-12db11170ad1
md"""
### Comparación con datos Javier Espuny

Carga opcional del archivo de referencia espectral (`M5T2cerrados.txt`) para contrastar resultados del notebook frente a una salida externa.
"""

# ╔═╡ 8bc13219-3baf-41d7-bf8d-3bc0e09be6bb
begin
# -----------------------------------------------------------------------------------
# 7. COMPARACIÓN CON RESULTADOS DE REFERENCIA (BrainVision Analyzer)
# -----------------------------------------------------------------------------------
# Se cargan resultados espectrales obtenidos con BrainVision Analyzer (software
# comercial de análisis EEG) para comparación y validación.
#
# NOTA: Los valores de referencia pueden representar:
# - Mean Activity: densidad media de potencia en la banda (µV²/Hz)
# - Area: potencia integrada en la banda (µV²)
#
# La comparación permite verificar que el procesamiento en Julia produce
# resultados consistentes con herramientas estándar del campo.
# "real" (µV²/Hz) o un simple "Area" (µV²) para la potencia total.

# Ruta al archivo de resultados de referencia
path_candidates = [
    joinpath(project_root(), "Javier_results", "M5T2cerrados.txt"),      # layout actual (raíz del proyecto)
    joinpath(@__DIR__, "..", "Javier_results", "M5T2cerrados.txt"),      # compatibilidad legacy
]
existing_paths = filter(isfile, path_candidates)
path_javier = isempty(existing_paths) ? "" : first(existing_paths)
if isempty(path_javier)
    error("No se encontró M5T2cerrados.txt en: $(join(path_candidates, " | "))")
end

# Leer líneas del archivo y filtrar líneas vacías
lines_javier = filter(x -> !isempty(strip(x)), readlines(path_javier))

# Parsear encabezado (primera línea)
# Se separa por espacios múltiples (2 o más)
header = split(strip(lines_javier[1]), r"\s{2,}")
channels_javier = header[2:end]  # Primera columna es "Band", resto son canales

# Parsear filas de datos
rows_javier = [split(strip(line), r"\s{2,}") for line in lines_javier[2:end]]
bands_javier = first.(rows_javier)  # Primera columna de cada fila es el nombre de banda

# Parsear números (convertir comas a puntos para formato decimal)
numbers = [
    parse(Float64, replace(rows_javier[i][j], ',' => '.'))
    for i in eachindex(rows_javier), j in 2:length(header)
]

# Construir DataFrame con resultados de referencia
Javier_results = DataFrame(Band = bands_javier)
for (col_idx, channel) in enumerate(channels_javier)
    Javier_results[!, Symbol(channel)] = numbers[:, col_idx]
end

Javier_results
end

# ╔═╡ f7193f9f-0d48-4e8e-bf11-0aaf67b9f04d
md"""
## Siguientes scripts de referencia

- `src/BIDS/build_eeg_bids.jl`
- `src/BIDS/build_participants.jl`
- `src/BIDS/validate_bids.jl`
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
CSV = "~0.10.16"
DSP = "~0.8.4"
DataFrames = "~1.8.1"
FFTW = "~1.10.0"
Plots = "~1.41.6"
PlutoUI = "~0.7.79"
StatsBase = "~0.34.10"
StatsPlots = "~0.15.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.3"
manifest_format = "2.0"
project_hash = "d530016d602be709656807f6199f79c1ac9f1ac6"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "35ea197a51ce46fcd01c4a44befce0578a1aaeca"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.5.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "libblastrampoline_jll"]
git-tree-sha1 = "7f54761502ff149a9d492e4acefe9805898e29b3"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.2+0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

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

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

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
weakdeps = ["OffsetArrays"]

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

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

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

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

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

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

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

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

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

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

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "4260cfc991b8885bf747801fb60dd4503250e478"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.11"

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

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariateStats]]
deps = ["Arpack", "Distributions", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "816620e3aac93e5b5359e4fdaf23ca4525b00ddf"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NearestNeighbors]]
deps = ["AbstractTrees", "Distances", "StaticArrays"]
git-tree-sha1 = "e2c3bba08dd6dedfe17a17889131b885b8c082f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.27"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

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

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

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

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

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

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

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

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

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

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

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
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

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

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "88cf3587711d9ad0a55722d339a013c4c56c5bbc"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.8"

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

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.8.3+2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

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

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "e9aeb174f95385de31e70bd15fa066a505ea82b9"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.7"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

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
# ╠═8d84f7cc-d13b-4b49-a2e0-cbf70d2d75a5
# ╠═62346384-1eef-4e5a-b7a0-e1b8af341cbe
# ╠═6a3f59f9-df9d-44a6-9a0e-bbc5a87d9057
# ╠═0f14b248-b15f-4fa4-b2f0-ff5f6609d66a
# ╟─c5a6ef95-7fb6-4368-9d4d-6a44d95dd95e
# ╟─97e1f81e-e8fc-4232-a7f4-0f145f2eae11
# ╟─3c6fd072-198c-410a-a771-85bc8e8ac5a3
# ╟─cc99f17e-dd5f-4911-ac73-5264b3235940
# ╟─601ed72d-3245-4922-b887-e9919e1fdbe7
# ╟─d0ec074d-59d7-4f8a-9146-7fe9a4ca870e
# ╟─e13eb5ef-c053-4dca-9a96-8509068937df
# ╟─8c6ad877-f4e2-422a-a86f-5ae47f76f2ba
# ╟─0381f849-d27f-4208-a63f-08b7ed99e1db
# ╠═11f6f202-6ca0-4f51-a6d3-e59aeb4e8130
# ╠═5fddf73d-9b85-47f7-b8eb-c4788de33484
# ╠═77df8d14-c178-44e8-982f-71b38ffaf77b
# ╠═f9f9a9a3-a2c3-4a93-a3de-80de403f0d5f
# ╠═a9f6c35a-7e7b-4383-80ff-1df3d49f2bf6
# ╠═6dcf7d5d-fc68-4f21-8f40-406f173122f2
# ╠═1ccb18f5-b6f8-431a-a670-a94efea7e988
# ╠═8b676680-8e27-49f7-991f-c8dffef2f62d
# ╠═70094435-8e89-4d22-a9de-fd8f0b2ac76d
# ╠═66fa1086-568d-4eb8-a8ea-eb96f19dc6ef
# ╠═fef42c56-e8e2-4ec2-9ad0-6bc9d1f37772
# ╠═b4de8852-c51f-423f-be51-74f67f8404d2
# ╠═f3cfdba1-8641-4c5e-b9e4-f3f1138859ea
# ╠═286ace1d-be59-4a0b-8b73-2ba751129850
# ╟─2a92ae15-c2c4-4b2c-8523-7f63475d5057
# ╠═ab98be9e-c06f-4a9f-929f-0ab2b8ce6edf
# ╠═6668e05c-ec9f-44cb-a66f-ccb41b2f44dd
# ╠═2e4c5804-bdc9-4f53-bbc4-b3df6f82c145
# ╟─8b984de8-8ec4-46ac-adf0-0a2443774f24
# ╠═3083f7e8-b7f8-4f91-a35d-b6720fc03fdc
# ╠═530abe8f-3dc1-424f-a04c-b6a57e723d2f
# ╠═6a5025f6-c6a7-43cb-9cdd-b6ea8cd3935e
# ╠═ae0949bc-70b0-4cf2-9ab9-b5dc33e9b5ad
# ╠═6ea65e96-4001-4795-b1f0-052db14cec52
# ╟─13177dee-e72f-4997-a204-e110ef7ea60c
# ╠═23b2f8f8-0133-46e9-a5d2-0c00aab6d2a4
# ╠═98c1de45-f9cf-43ee-b903-5f6ca9897269
# ╠═e2643d97-07f8-4fd4-b4ad-9e36b340f6e6
# ╠═8f52ce15-6075-4c01-8361-4c18b4aae41a
# ╠═50bb2eb0-d367-4f48-a8db-c5d4db35e7f4
# ╟─ad3ebb32-88af-4f4d-8aac-e35f4ee76aee
# ╠═c4d34a8d-ad3d-44f6-934d-cf5c14670440
# ╠═4d28c8d4-c1fc-4c7d-bcc8-09e8f568b0f6
# ╠═2962c98d-6b35-42cf-9dfa-30ecebe4f38b
# ╠═61ec2843-f706-41f3-9eb2-8bad3e1ea3ca
# ╠═8ad9dcae-9ea1-4eb1-bf41-12db11170ad1
# ╠═8bc13219-3baf-41d7-bf8d-3bc0e09be6bb
# ╟─f7193f9f-0d48-4e8e-bf11-0aaf67b9f04d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
