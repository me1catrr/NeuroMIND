# NeuroMIND

**Framework de análisis EEG para conectividad funcional en Esclerosis Múltiple**

NeuroMIND toma señales EEG de reposo en formato BrainVision, las procesa de principio a fin y genera matrices de conectividad wPLI / dwPLI con inferencia estadística opcional. El análisis cubre el dataset MINDEM-IMIBIC (41 pacientes EM + 37 controles sanos, ~206 grabaciones válidas) y produce resultados listos para comparar grupos y sesiones.

---

## Tabla de contenidos

1. [¿Qué hace el pipeline?](#1-qué-hace-el-pipeline)
2. [Requisitos previos](#2-requisitos-previos)
3. [Primeros pasos: preparar el entorno](#3-primeros-pasos-preparar-el-entorno)
4. [Flujo de trabajo completo](#4-flujo-de-trabajo-completo)
   - 4.1 [Ejecutar el pipeline en lote (todos los sujetos)](#41-ejecutar-el-pipeline-en-lote)
   - 4.2 [Revisar la tabla de QC](#42-revisar-la-tabla-de-qc)
   - 4.3 [Análisis transversal (EM vs controles)](#43-análisis-transversal)
   - 4.4 [Análisis longitudinal (T1 → T2)](#44-análisis-longitudinal)
   - 4.5 [Lanzar el dashboard](#45-lanzar-el-dashboard)
   - 4.6 [Compilar el informe PDF](#46-compilar-el-informe-pdf)
5. [Los 8 pasos del pipeline explicados](#5-los-8-pasos-del-pipeline-explicados)
6. [Política de calidad de señal (QC)](#6-política-de-calidad-de-señal-qc)
7. [Interpretar los resultados](#7-interpretar-los-resultados)
8. [Configuración avanzada](#8-configuración-avanzada)
9. [Estructura del proyecto](#9-estructura-del-proyecto)
10. [Referencia de dependencias](#10-referencia-de-dependencias)

---

## 1. ¿Qué hace el pipeline?

```
Señal EEG cruda (.vhdr + .eeg)
        │
        ▼
[1] Carga y validación BIDS
[2] Control de calidad (QC) de canales
[3] Filtrado (highpass + lowpass + notch + bandreject)
[4] ICA — eliminación de artefactos oculares y musculares
[5] Segmentación en épocas + rechazo de artefactos ±70 µV
[6] Espectro de potencia (PSD) por banda de frecuencia
[7] Conectividad funcional wPLI/dwPLI  [Hilbert | FourierCSD | Multitaper]
[8] Surrogates de desplazamiento circular + corrección FDR al 5 %
        │
        ▼
Tablas CSV  ·  Figuras PNG  ·  Dashboard interactivo  ·  Informe PDF
```

**¿Por qué wPLI?**  
El *weighted Phase Lag Index* mide la sincronización de fase entre pares de canales ignorando las contribuciones de campo de volumen y el ruido de amplitud. Es el estimador de conectividad más robusto para señales EEG de reposo. La variante **dwPLI** (debiased wPLI, Vinck 2011) elimina el sesgo por número variable de épocas entre sujetos y es la opción recomendada para comparaciones grupales.

**Bandas de frecuencia analizadas:**

| Banda | Rango | Relevancia clínica |
|-------|-------|-------------------|
| δ (Delta) | 0.5 – 4 Hz | Sueño, estados de baja vigilancia |
| θ (Theta) | 4 – 8 Hz | Memoria de trabajo, cognición |
| α (Alpha) | 7.8 – 11.7 Hz | Estado de reposo, inhibición cortical |
| β_low | 12 – 15 Hz | Control motor, atención sostenida |
| β_mid | 15 – 18 Hz | Actividad sensoriomotora |
| β_high | 18 – 30 Hz | Procesos cognitivos de alto nivel |
| γ (Gamma) | 30 – 50 Hz | Procesamiento sensorial integrado |

> **Nota sobre Delta:** con épocas de 1 s, la banda Delta (0.5 Hz) acumula solo 0.5 ciclos por época, por debajo del mínimo recomendado de 4. Si se trabaja con el perfil `eeg_julia` (épocas de 1 s), considerar activar `exclude_unreliable_bands = true` o usar el perfil `default` con `segment_length_seconds = 8.0`.

---

## 2. Requisitos previos

- **Julia 1.9 o superior** — [descargar en julialang.org](https://julialang.org/downloads/)
- **LaTeX** (TeX Live o MacTeX) — solo si quieres compilar el informe PDF
- **~10 GB de espacio libre** en disco para los resultados de todo el dataset

Verifica tu versión de Julia:

```bash
julia --version
# Debe mostrar julia version 1.9.x o superior
```

---

## 3. Primeros pasos: preparar el entorno

Ejecuta este comando **una sola vez** al clonar el repositorio (o al cambiar de ordenador). Descargará e instalará todas las dependencias de Julia:

```bash
cd NeuroMIND/
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Para verificar que todo está correctamente instalado:

```bash
julia --project=. -e 'include("src/NeuroMIND.jl"); println("✓ NeuroMIND cargado correctamente")'
```

Si ves `✓ NeuroMIND cargado correctamente`, el entorno está listo.

---

## 4. Flujo de trabajo completo

Sigue estos pasos en orden. Cada uno depende del anterior.

### 4.1 Ejecutar el pipeline en lote

Este comando procesa todas las grabaciones válidas del dataset y genera los resultados en `results/`.

```bash
julia --project=. scripts/run_batch_pipeline.jl
```

**Para procesar solo un subconjunto de sujetos** (útil para pruebas):

```bash
# Procesar 3 sujetos específicos
julia --project=. scripts/run_batch_pipeline.jl --subjects M43 M44 MC01

# Procesar solo sujetos EM
julia --project=. scripts/run_batch_pipeline.jl --group MS

# Procesar solo controles
julia --project=. scripts/run_batch_pipeline.jl --group Control
```

**¿Cuánto tarda?**

| Sujetos | Surrogates | Tiempo estimado |
|---------|------------|-----------------|
| 1 sujeto, 1 condición | OFF | ~2–5 min |
| 1 sujeto, 1 condición | 200 permutaciones | ~10–15 min |
| Dataset completo (206 grabaciones) | OFF | ~6–12 h |
| Dataset completo | 200 permutaciones | ~23–34 h |

> **Nota:** El pipeline es incremental — si se interrumpe, los sujetos ya procesados se saltan automáticamente en la siguiente ejecución.

**¿Qué verás en pantalla durante la ejecución?**

```
══════════════════════════════════════════════════════════════
NeuroMIND · Pipeline individual · 2026-05-26 10:32:15
══════════════════════════════════════════════════════════════
  Dataset : MINDEM-IMIBIC
  ▶  sub-M43  ·  ses-T1  ·  eyesclosed

  [1/8] Cargando EEG...       ✓  31 ch · 47460 muestras · 500.0 Hz · 94.9 s
  [2/8] QC canales...         ✓  σ̄=8.3 µV · ✓ filtro OK · bad: ninguno
  [3/8] Filtrado...           ✓  HP 0.5 Hz · LP 150.0 Hz · Notch 50.0 Hz
  [4/8] ICA...                ✓  30 comp · rechazados: IC01, IC05 [CACHÉ]
  [5/8] Segmentación + AR...  ✓  82/94 epochs válidos (87.2%) · ±70 µV · 30 ch
  [6/8] PSD...                ✓  ALPHA=1.23 · THETA=0.87 · DELTA=2.11 ... µV²
  [7/8] wPLI (FourierCSD)...  ✓  30×30 · 7 bandas · 435 aristas
  [8/8] Guardando resultados... ✓  tablas · figuras · QC table

──────────────────────────────────────────────────────────────
✅  sub-M43 / ses-T1 / eyesclosed  →  include
     Epochs  : 82/94 válidos (87.2%) · rechazados: 12
     Montaje : 30 canales (Fp2 excluido) · bad no-Fp2: ninguno
     PSD     : ALPHA=1.23 · THETA=0.87 · DELTA=2.11 ... µV²
──────────────────────────────────────────────────────────────
```

### 4.2 Revisar la tabla de QC

Después del pipeline, revisa este archivo para ver el estado de cada grabación:

```
results/qc/qc_decision_table.csv
```

Cada fila representa una grabación. Los valores posibles de `final_decision` son:

| Decisión | Significado | ¿Incluir en análisis grupal? |
|----------|------------|------------------------------|
| `include` | Grabación limpia sin alertas | ✅ Sí |
| `include_with_warning` | Válida, con alguna alerta menor | ✅ Sí (con cautela) |
| `manual_review` | Pocas épocas válidas o señal ruidosa | ⚠️ Revisar manualmente |
| `exclude` | 0 épocas válidas tras rechazo AR | ❌ No |

Para abrir la tabla en Julia:

```julia
using CSV, DataFrames
qc = CSV.read("results/qc/qc_decision_table.csv", DataFrame)
# Ver resumen
combine(groupby(qc, :final_decision), nrow => :n_grabaciones)
```

### 4.3 Análisis transversal

Compara la conectividad media entre el grupo EM y el grupo de controles:

```bash
julia --project=. scripts/run_transversal_analysis.jl
```

**Prerequisito:** que exista `data/bids/groups.csv` con columnas `subject_id, group, session_id`.

Los resultados se guardan en:
```
results/group/transversal/EC/    ← condición ojos cerrados
results/group/transversal/EO/    ← condición ojos abiertos
```

Archivos generados por banda (ej. para alpha):
- `group_connectivity_ctrl_ALPHA.csv` — matriz wPLI media del grupo control
- `group_connectivity_ms_ALPHA.csv` — matriz wPLI media del grupo EM
- `group_difference_ALPHA.csv` — diferencia EM − control
- `group_statistics_ALPHA.csv` — test estadístico por par de canales
- `significant_edges_ALPHA.csv` — pares significativos tras FDR

### 4.4 Análisis longitudinal

Compara la conectividad entre la primera visita (T1) y la segunda (T2) para pacientes EM:

```bash
julia --project=. scripts/run_longitudinal_analysis.jl
```

El script lee `data/bids/longitudinal_pairs.csv` (generado por `scripts/audit_full_dataset.jl`) y detecta automáticamente los pares T1/T2 disponibles. Los resultados se guardan en:
```
results/group/longitudinal/EC/
results/group/longitudinal/EO/
```

### 4.5 Lanzar el dashboard

El dashboard permite explorar interactivamente todos los resultados paso a paso:

```bash
julia --project=. scripts/launch_dashboard.jl
```

Abre `http://localhost:8080` en el navegador.

**¿Qué puedes hacer en el dashboard?**

| Panel | Qué muestra |
|-------|-------------|
| 0 – Proyecto | Estado general del dataset |
| 1–3 | Metadatos BIDS, señal cruda, QC de canales |
| 4 | Respuesta del filtro diseñado |
| 5 | Componentes ICA, topomaps, antes/después |
| 6–7 | Épocas segmentadas, artefactos rechazados |
| 8 | Espectro de potencia y topomapas por banda |
| 9 | Heatmaps y red de conectividad wPLI |
| 10 | Resultados de surrogates y FDR |
| 11 | Resumen global de calidad de la grabación |
| 12 | Inventario de archivos exportados |
| 13 | Comparación transversal EM vs controles |
| 14 | Comparación longitudinal T1 → T2 |
| 15 | Validación MNE-Python (pipeline mne_brain) |

### 4.6 Compilar el informe PDF

El informe oficial en LaTeX se encuentra en `report/`. Para compilarlo:

```bash
cd report
make build-es
```

El PDF se genera en:
```
report/build/pdf/main_es.pdf
```

Para validar errores de compilación:

```bash
rg -n "LaTeX Error|File .* not found|Undefined control sequence" report/build/pdf/main_es.log
```

---

## 5. Los 8 pasos del pipeline explicados

### [1/8] Carga EEG
Lee los archivos BrainVision (`.vhdr` + `.eeg`) y crea un objeto `EEGRecording` con la señal en formato canales × muestras.

- Frecuencia de muestreo: 500 Hz
- Canales: 31 (10-20 internacional, referencia FCz)
- Duración típica: 60–120 s de señal en reposo

### [2/8] Control de calidad de canales
Calcula la desviación estándar de cada canal. Los canales cuya amplitud sea estadísticamente atípica (z-score ≥ 3.0) se marcan como malos y se excluyen del análisis de conectividad.

También calcula la **amplitud media global** (σ̄):
- σ̄ < 20 µV → grabación con filtro online activo (normal)
- σ̄ > 20 µV → posible grabación sin filtro online activo → se emite `amplitude_warning`

### [3/8] Filtrado
Aplica un banco de filtros Butterworth de orden 4 con `filtfilt` (cero retardo de fase), en el orden del protocolo `eeg_julia`:

1. **Notch 50 Hz** — elimina interferencia de red eléctrica
2. **Bandreject 99.5–100.5 Hz** — elimina subarmónico de red
3. **Highpass 0.5 Hz** — elimina deriva lenta y artefactos DC
4. **Lowpass 150 Hz** — elimina ruido de alta frecuencia

### [4/8] ICA (Análisis de Componentes Independientes)
Aplica FastICA simétrico a la señal **continua** filtrada (no segmentada). Esto es crítico: aplicar ICA sobre la señal continua maximiza la información disponible para separar fuentes.

Cada componente se clasifica automáticamente como artefacto ocular, muscular, cardíaco o señal cerebral. Los componentes de artefacto se sustraen de la señal antes de continuar.

Los resultados ICA se guardan en caché (`results/{suj}/{ses}/cache/ica_result.jls`). Si el pipeline se relanza con la misma configuración, ICA se recupera del caché en segundos.

### [5/8] Segmentación y rechazo de artefactos
Corta la señal limpiada en épocas y aplica rechazo por amplitud. El perfil de segmentación se configura en `config/batch_pipeline.toml`:

| Perfil | Longitud de época | Ventaja |
|--------|-------------------|---------|
| `eeg_julia` | 1.0 s fijo (compatibilidad EEG_Julia) | Reproducibilidad con pipeline original |
| `default` | Configurable (recomendado: 2.0 s) | Más ciclos por época, mejor estimación wPLI |

Pasos internos:
1. **Corrección de baseline** — sustrae la media de los primeros 100 ms
2. **Rechazo por amplitud ±70 µV** — elimina épocas con excursiones extremas
3. Segunda pasada de corrección de baseline (perfil `eeg_julia`)

> **Umbral ±70 µV:** estándar del protocolo RS-MIND. Las grabaciones con `amplitude_warning` pueden tener pocas épocas válidas — esto se registra en la tabla QC.

**Montaje del análisis:** se excluye Fp2 (canal sistemáticamente problemático en el 46% del dataset) y los canales malos detectados en [2/8]. El análisis trabaja con **30 canales**.

### [6/8] Espectro de potencia (PSD)
Calcula la densidad espectral de potencia de cada época usando ventana Hanning (FFT de 1024 puntos). La potencia se integra en cada banda de frecuencia.

### [7/8] Conectividad wPLI / dwPLI
Calcula el *weighted Phase Lag Index* entre todos los pares de canales en cada banda de frecuencia. Con 30 canales, hay **435 pares de electrodos** por banda.

**Tres métodos de estimación configurables:**

| Método | Descripción | Cuándo usarlo |
|--------|-------------|---------------|
| `hilbert` | Filtrado Butterworth + señal analítica Hilbert | Análisis exploratorio rápido |
| `fourier_csd` | Espectro cruzado FFT con ventana Hanning/Hamming | **Recomendado para publicación** |
| `multitaper` | DPSS multitaper (equivalente a MNE `spectral_connectivity`) | Cuando se requiere bajo sesgo espectral |

El estimador se selecciona en `config/batch_pipeline.toml`:

```toml
[connectivity]
wpli_method = "fourier_csd"  # o "hilbert", "multitaper"
use_dwpli   = true           # true → dwPLI no sesgado (recomendado para grupos)
```

> La variante **dwPLI** (debiased wPLI) es matemáticamente insesgada respecto al número de épocas y produce rango [-1, 1]. Se recomienda para comparaciones grupales donde el número de épocas varía entre sujetos.

### [8/8] Surrogates y FDR
Para cada par de canales y cada banda, genera surrogates mediante **desplazamiento circular independiente por canal y época**. Esta técnica:
- Destruye la sincronía de fase inter-canal (lo que mide wPLI)
- Preserva el espectro de potencia y la autocorrelación dentro de cada canal
- Usa el **mismo estimador** configurado en [7/8] → distribución nula metodológicamente coherente

Los p-valores se calculan con corrección Monte Carlo (+1) y se corrigen por comparaciones múltiples con FDR de Benjamini-Hochberg al 5 %.

```toml
[surrogates]
enabled      = true
n_surrogates = 200    # ≥ 200 recomendado para publicación; 20 para exploración rápida
alpha        = 0.05
seed         = 42
```

---

## 6. Política de calidad de señal (QC)

### Montaje de 30 canales

El análisis principal trabaja con **30 canales** (Fp2 excluido en todos los sujetos). Fp2 aparece como canal problemático en el 46% del dataset (artefactos oculares y de contacto). Excluirlo garantiza matrices de conectividad homogéneas entre grabaciones.

Esto produce matrices wPLI de **30 × 30** con **435 aristas por banda**.

Para un análisis de sensibilidad con Fp2 activo, editar `config/batch_pipeline.toml`:

```toml
[montage]
exclude_fp2 = false
```

### Alerta de amplitud (`amplitude_warning`)

Si la amplitud media de la señal cruda supera 20 µV, el pipeline emite esta alerta. Indica probable ausencia de filtro online en la grabación. **No provoca exclusión automática** — la decisión final se toma tras ver el porcentaje de épocas válidas tras el rechazo AR ±70 µV.

### Criterios de inclusión

| Condición | Decisión |
|-----------|----------|
| Sin alertas, ≥10 épocas válidas | `include` |
| `amplitude_warning` o canales malos adicionales | `include_with_warning` |
| `amplitude_warning` + ≥2 canales malos **o** <50% épocas válidas | `manual_review` |
| 0 épocas válidas | `exclude` |

---

## 7. Interpretar los resultados

### Estructura de resultados por sujeto

```
results/subjects/sub-{ID}/ses-{SES}/{task}/
├── overview.csv              ← resumen general (canales, epochs, duración)
├── qc_summary.csv            ← QC de la grabación
├── channel_statistics.csv    ← estadísticas por canal
├── band_power_summary.csv    ← potencia por banda
├── wpli_{BANDA}.csv          ← matriz wPLI 30×30 para cada banda
├── wpli_pvalues_{BANDA}.csv  ← p-valor por par de canales (si surrogates ON)
├── wpli_qvalues_{BANDA}.csv  ← q-valor corregido FDR
├── wpli_significant_{BANDA}.csv  ← máscara de significancia
├── significant_connections.csv   ← top conexiones significativas
├── surrogate_summary.json    ← n_sig, % significativas, umbral FDR
├── ica_summary.json          ← componentes ICA, rechazados, varianza
├── figures/                  ← figuras PNG (señal, PSD, wPLI, ICA)
└── pipeline_log.txt          ← log completo de la ejecución
```

### ¿Cómo leer los resultados de conectividad?

El archivo `significant_connections.csv` contiene las conexiones estadísticamente significativas:

```
ch1, ch2, band,  wpli_observed, p_value, q_value
F3,  F4,  ALPHA, 0.234,         0.002,   0.038
...
```

- `wpli_observed` > 0 → sincronización en fase (leading-lagging relationship)
- `q_value` < 0.05 → significativo tras corrección por comparaciones múltiples

### Señales de alerta en los resultados

| Señal | Qué hacer |
|-------|-----------|
| `n_epochs_valid = 0` | Revisar manualmente; posible grabación corrupta |
| `valid_epochs_pct < 30%` | La grabación es marginal; usar con cautela en análisis grupal |
| `amplitude_warning = true` | Normal si el filtro online estaba desactivado; verificar filtrado offline |
| `n_sig_total = 0` en surrogates | Puede ser normal para sujeto individual; las diferencias aparecen a nivel grupal |
| `@warn [DELTA] solo N ciclos` | Con épocas cortas, Delta tiene pocos ciclos; considerar `exclude_unreliable_bands = true` |

---

## 8. Configuración avanzada

La configuración del pipeline en lote vive en `config/batch_pipeline.toml`. Los parámetros más relevantes:

### Segmentación

```toml
[segmentation]
profile                = "default"         # "eeg_julia" = 1s fijo; "default" = configurable
segment_length_seconds = 2.0               # épocas de 2 s → mejor estimación wPLI en Delta
overlap_seconds        = 0.0
min_segments           = 10
```

### Método de conectividad

```toml
[connectivity]
wpli_method              = "fourier_csd"   # "hilbert" | "fourier_csd" | "multitaper"
use_dwpli                = true            # debiased wPLI — recomendado para grupos
min_cycles_for_wpli      = 4.0             # ciclos mínimos por época para estimación fiable
exclude_unreliable_bands = true            # omitir bandas con < min_cycles (no solo advertir)

[connectivity.fourier_csd]
window = "hann"   # "hann" | "hamming" | "rect"
nfft   = 0        # 0 = longitud de la época completa

[connectivity.multitaper]
nw       = 4.0    # time-bandwidth product (mayor = más suavizado)
n_tapers = 0      # 0 = automático: floor(2×nw)−1
low_bias = true   # descartar tapers con concentración espectral λ < 0.9
```

### Surrogates

```toml
[surrogates]
enabled      = false    # true para inferencia estadística; false para análisis rápido
n_surrogates = 200      # ≥ 200 para publicación; 20 para exploración
method       = "circular_shift"
alpha        = 0.05
fdr_method   = "bh"
seed         = 42
```

### Montaje y rechazo de artefactos

```toml
[montage]
exclude_channels    = ["Fp2"]
exclude_fp2         = true      # false para análisis de sensibilidad
n_channels_analysis = 30

[artifact_rejection]
min_amplitude_uv = -70.0   # protocolo RS-MIND estándar
max_amplitude_uv =  70.0
```

### Ejecutar un sujeto individual (para pruebas)

```bash
julia --project=. scripts/run_single_subject.jl
```

Configuración en `config/single_subject.toml`:

```toml
[subject]
subject_id = "M05"
session_id = "T2"
task       = "eyesclosed"
```

### Tests unitarios

```bash
julia --project=. tests/runtests.jl
```

---

## 9. Estructura del proyecto

```
NeuroMIND/
├── config/
│   ├── single_subject.toml     ← Parámetros para análisis de un sujeto
│   └── batch_pipeline.toml     ← Parámetros para el análisis en lote
│
├── data/
│   ├── bids/                   ← Estructura BIDS ligera (metadata JSON, sin señales)
│   │   ├── raw/                   sub-{ID}_ses-{SES}_task-*_eeg_metadata.json
│   │   ├── groups.csv             sujeto → grupo (MS / Control)
│   │   └── longitudinal_pairs.csv pares T1/T2 detectados por audit_full_dataset.jl
│   └── full_data/
│       └── inventory.csv       ← Inventario de las ~212 grabaciones del dataset
│
├── results/                    ← Generado por el pipeline (NO en Git)
│   ├── subjects/               ← Resultados individuales por sujeto/sesión
│   ├── group/                  ← Resultados de análisis grupal
│   │   ├── transversal/        ← EM vs controles
│   │   └── longitudinal/       ← T1 vs T2
│   └── qc/
│       └── qc_decision_table.csv   ← Estado QC de cada grabación
│
├── src/                        ← Código fuente Julia
│   ├── NeuroMIND.jl            ← Entry point del módulo
│   ├── types.jl                ← Tipos de datos (EEGRecording, EpochSet, etc.)
│   ├── SingleSubjectPipeline.jl ← Pipeline de 8 pasos
│   ├── io/                     ← Carga de datos (BIDS, BrainVision)
│   ├── preprocessing/          ← Filtrado
│   ├── ica/                    ← FastICA, clasificación, inspección
│   ├── segmentation/           ← Épocas, baseline, rechazo AR
│   ├── spectral/               ← PSD
│   ├── connectivity/
│   │   ├── wPLI.jl             ← Estimadores Hilbert / FourierCSD / Multitaper
│   │   ├── GraphMetrics.jl     ← Métricas de grafo (strength, clustering, path length)
│   │   └── CSD.jl              ← Current Source Density (opcional)
│   ├── statistics/             ← Surrogates (circular_shift), FDR
│   ├── visualization/          ← Topomaps, heatmaps, espectros
│   └── webapp/                 ← Dashboard web (Genie.jl)
│
├── scripts/
│   ├── run_batch_pipeline.jl        ← Procesar todos los sujetos
│   ├── run_single_subject.jl        ← Procesar un sujeto
│   ├── run_transversal_analysis.jl  ← Comparación grupal EM vs controles
│   ├── run_longitudinal_analysis.jl ← Comparación temporal T1 → T2
│   ├── audit_full_dataset.jl        ← Generar longitudinal_pairs.csv
│   └── launch_dashboard.jl          ← Dashboard interactivo
│
├── mne_brain/                  ← Pipeline de validación MNE-Python (Panel 15)
│   └── ...
│
├── web/
│   └── views/dashboard.html    ← Interfaz del dashboard (SPA, 16 paneles)
│
├── report/                     ← Informe científico LaTeX
│   ├── main_es.tex
│   └── build/pdf/main_es.pdf   ← PDF compilado
│
├── tests/
│   └── runtests.jl
│
└── CLAUDE.md                   ← Contexto unificado para asistentes IA (Claude, Codex, Cursor)
```

---

## 10. Referencia de dependencias

| Paquete Julia | Uso en el pipeline |
|---------------|--------------------|
| `DSP` | Filtros Butterworth, `filtfilt`, DPSS tapers (multitaper wPLI) |
| `FFTW` | FFT, espectro cruzado (FourierCSD), transformada de Hilbert |
| `LinearAlgebra`, `Random` | FastICA (implementación propia, sin paquetes externos) |
| `Statistics`, `StatsBase` | Estadísticas de canales, z-scores |
| `CairoMakie` | Figuras PNG (señal, PSD, heatmaps, topomaps) |
| `DataFrames`, `CSV` | Tablas de resultados |
| `TOML` | Lectura de configuración |
| `Serialization` | Caché de resultados ICA intermedios |
| `Genie` | Servidor web del dashboard |
