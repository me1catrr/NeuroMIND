# NeuroMIND

**Framework de análisis EEG para conectividad funcional en Esclerosis Múltiple**  
Señal continua → ICA → Segmentación → PSD → wPLI · Dashboard web interactivo

---

## Estructura del proyecto

```
NeuroMIND/
├── config/
│   ├── single_subject.toml    ← Parámetros para análisis de un sujeto
│   ├── pipeline.toml          ← Parámetros del pipeline grupal
│   └── subjects.toml          ← Registro de sujetos y sesiones
│
├── data/
│   └── BIDS/
│       ├── raw/               ← Señales EEG (.tsv) + metadata (.json)
│       └── electrodes/        ← Posiciones de electrodos (.tsv)
│
├── results/                   ← Generado por el pipeline
│   ├── {subj}/{sess}/
│   │   ├── figures/           ← Figuras PNG (señal, PSD, wPLI)
│   │   ├── tables/            ← Tablas CSV (QC, PSD, band power)
│   │   ├── cache/             ← Resultados intermedios serializados
│   │   └── logs/
│   └── subjects/
│       └── sub-{id}/ses-{sess}/{task}/
│           ├── overview.csv
│           ├── qc_summary.csv
│           ├── channel_statistics.csv
│           ├── psd_by_channel.csv
│           ├── band_power_summary.csv
│           ├── connectivity_edges.csv
│           ├── wpli_{band}.csv
│           ├── ica_components.csv
│           ├── ica_summary.json
│           ├── ica_activations.csv
│           ├── ica_signal_before.csv
│           ├── ica_signal_after.csv
│           ├── pipeline_log.txt
│           └── config_snapshot.toml
│
├── scripts/
│   ├── run_single_subject.jl  ← Pipeline para un sujeto
│   ├── run_pipeline.jl        ← Pipeline grupal
│   └── launch_dashboard.jl    ← Lanzar dashboard web
│
├── src/
│   ├── NeuroMIND.jl           ← Módulo principal (entry point)
│   ├── types.jl               ← Tipos centrales del framework
│   ├── Pipeline.jl            ← Orquestador grupal con cache
│   ├── SingleSubjectPipeline.jl ← Pipeline sujeto individual (8 pasos)
│   │
│   ├── io/
│   │   ├── BIDSLoader.jl      ← Carga datasets BIDS
│   │   ├── Config.jl          ← Rutas y acceso a PipelineConfig
│   │   └── Serializer.jl      ← Cache versionado
│   │
│   ├── qc/
│   │   └── QualityControl.jl  ← Estadísticas y flags de canales malos
│   │
│   ├── preprocessing/
│   │   └── Filtering.jl       ← Highpass, lowpass, notch, bandreject
│   │
│   ├── ica/
│   │   ├── ICACore.jl         ← FastICA simétrico (PCA whitening + tanh)
│   │   └── ICAInspection.jl   ← Carga de labels y rechazo de componentes
│   │
│   ├── segmentation/
│   │   └── Epochs.jl          ← Segmentación + baseline + rechazo AR
│   │
│   ├── spectral/
│   │   └── PowerSpectrum.jl   ← PSD por ventana Hanning + potencia por banda
│   │
│   ├── connectivity/
│   │   ├── wPLI.jl            ← wPLI across-segments (Hilbert analítica)
│   │   ├── CSD.jl             ← CSD (opcional, use_csd = false por defecto)
│   │   └── GraphMetrics.jl    ← Métricas de teoría de grafos
│   │
│   ├── statistics/
│   │   ├── Surrogates.jl      ← Phase-shuffle surrogates
│   │   └── FDR.jl             ← Benjamini-Hochberg
│   │
│   ├── visualization/
│   │   ├── Topomaps.jl        ← Mapas topográficos IDW
│   │   ├── Heatmaps.jl        ← Heatmaps de conectividad
│   │   └── Spectra.jl         ← Espectros de potencia
│   │
│   └── webapp/
│       └── App.jl             ← Dashboard web (Genie.jl + HTML/JS)
│
├── web/
│   ├── views/
│   │   └── dashboard.html     ← UI completa del dashboard (single-page)
│   └── public/                ← Assets estáticos
│
├── tests/
│   └── runtests.jl            ← Suite de tests unitarios
│
├── Project.toml
└── README.md
```

---

## Inicio rápido

### 1. Preparar el entorno

```bash
cd NeuroMIND/
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### 2. Ejecutar el pipeline (sujeto individual)

```bash
julia --project=. scripts/run_single_subject.jl
```

El sujeto se detecta automáticamente desde `data/BIDS/raw/`  
(primer archivo `sub-*_ses-*_task-*_run-*_eeg_data.tsv` por orden alfabético).

Para especificar el sujeto, editar `config/single_subject.toml`:

```toml
[subject]
subject_id = "M05"      # "auto" para detección automática
session_id = "T2"
task       = "eyesclosed"
run        = 1
```

### 3. Lanzar el dashboard

```bash
julia --project=. scripts/launch_dashboard.jl
```

Abre `http://localhost:8080` en el navegador.

---

## Pipeline de 8 pasos

El pipeline (`run_single_subject_pipeline`) ejecuta los siguientes pasos en orden:

| Paso | Descripción | Salida |
|------|-------------|--------|
| **[1/8] Carga** | Lee señal EEG (TSV), metadata (JSON) y posiciones de electrodos | `EEGRecording` |
| **[2/8] QC** | Estadísticas por canal; detecta canales malos por z-score | `qc_channels.csv` |
| **[3/8] Filtrado** | Highpass + Lowpass + Notch + Bandreject (Butterworth, `filtfilt`) | `EEGRecording` filtrado |
| **[4/8] ICA** | FastICA simétrico sobre señal **continua** filtrada; guarda activaciones y señal antes/después | 5 ficheros ICA |
| **[5/8] Segmentación** | Corta en epochs sobre señal **limpiada** por ICA; baseline + rechazo AR | `EpochSet` |
| **[6/8] Espectral** | PSD por ventana Hanning; potencia por banda (DELTA…GAMMA) | `SpectralResult` |
| **[7/8] wPLI** | Conectividad funcional entre pares de canales por banda; CSD opcional | `ConnectivityMatrix` |
| **[8/8] Guardado** | CSV/PNG/JSON en `results/`; log del pipeline; snapshot de config | ver carpeta `results/` |

> **Orden ICA → Segmentación**: ICA se ejecuta sobre la señal continua filtrada (paso 4), antes de segmentar (paso 5). Si se han cargado labels de rechazo, la señal limpiada alimenta directamente la segmentación y el análisis posterior.

### ICA — detalle técnico

- Algoritmo: **FastICA simétrico** (implementación propia en Julia puro)
- Blanqueo: PCA whitening (primeras `n_components` componentes principales)
- No linealidad: `tanh` (G supergaussiana)
- Decorrelación: simétrica en cada iteración
- Dependencias: solo `LinearAlgebra` + `Random` (sin paquetes externos)
- Defaults: 30 componentes, 500 iteraciones, tol=1e-5, seed=42

Para marcar componentes como artefactos, crear un archivo de labels y volver a ejecutar el pipeline (ver `ICAInspection.jl`).

---

## Dashboard web — 13 paneles

El dashboard muestra el estado de cada etapa del pipeline y los resultados calculados. Los paneles con fondo implementados están marcados con ✓.

| Panel | Nombre | Estado |
|-------|--------|--------|
| **0** | Proyecto / Dataset | ✓ implementado |
| **1** | BIDS & Metadata | ✓ implementado |
| **2** | Señal cruda | ✓ implementado |
| **3** | QC inicial | ✓ implementado |
| **4** | Preprocessing / Filtrado | ✓ implementado |
| **5** | ICA | ✓ implementado |
| **6** | Segmentación | placeholder |
| **7** | Rechazo artefactos | placeholder |
| **8** | Análisis espectral | placeholder |
| **9** | Conectividad wPLI | placeholder |
| **10** | Surrogates / Inferencia | placeholder |
| **11** | Resultados finales | placeholder |
| **12** | Exportación / Reporte | placeholder |

### Panel 5 — ICA

Si no se ha ejecutado ICA, el panel muestra un placeholder elegante con botón de acción. Una vez ejecutado el pipeline, muestra:

- **Resumen**: número de componentes, rechazados, varianza retenida, duración
- **Tabla de componentes**: % varianza, tipo de artefacto, estado (aceptado/rechazado)
- **Donut**: distribución aceptados vs rechazados
- **Señal del componente seleccionado**: activación temporal (Canvas)
- **Espectro del componente**: PSD con bandas sombreadas (SVG)
- **Comparación antes/después**: butterfly plot multicanal (Canvas)
- **Métricas de calidad y resumen de artefactos**

---

## Configuración — `config/single_subject.toml`

```toml
[subject]
subject_id = "auto"          # "auto" o ID explícito (ej. "M05")
session_id = "T2"
task       = "eyesclosed"    # "eyesclosed" → condición "EC"
run        = 1

[recording]
sampling_rate = 500.0        # Hz (sobreescrito por metadata si existe)

[filtering]
highpass_hz   = 0.5          # Highpass Butterworth
lowpass_hz    = 48.0         # Lowpass Butterworth
notch_hz      = 50.0         # Notch (red eléctrica)
notch_bw_hz   = 2.0          # Ancho de banda del notch
bandreject_lo = 100.0        # Inicio del rechazo de banda alta
bandreject_hi = 120.0        # Fin del rechazo de banda alta
filter_order  = 4

[segmentation]
segment_length_seconds = 2.0
overlap_seconds        = 0.0
min_segments           = 10

[artifact_rejection]
amplitude_threshold_uv = 100.0
gradient_threshold_uv  = 50.0
enabled                = true

[spectral]
nfft       = 1024
window_pct = 10.0            # Cosine taper en % de cada extremo

[bands]
DELTA     = [0.5,  4.0]
THETA     = [4.0,  8.0]
ALPHA     = [7.8, 11.7]
BETA_LOW  = [12.0, 15.0]
BETA_MID  = [15.0, 18.0]
BETA_HIGH = [18.0, 30.0]
GAMMA     = [30.0, 50.0]

[connectivity]
filter_order = 8
use_csd      = false         # CSD opcional — false por defecto

[qc]
bad_channel_zscore_threshold = 3.0

[paths]
bids_root = "data/BIDS"
results   = "results"

[dashboard]
port         = 8080
open_browser = true
```

---

## Formato de datos BIDS

```
data/BIDS/
├── raw/
│   ├── sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv
│   └── sub-M05_ses-T2_task-eyesclosed_run-01_eeg_metadata.json
└── electrodes/
    └── sub-M05_ses-T2_electrodes.tsv
```

**TSV de señal**: filas = canales, primera columna = nombre del canal, resto = muestras  
**JSON de metadata**: campos `fs`, `channel_names`, etc.  
**TSV de electrodos**: columnas `name`, `x`, `y` (posiciones en plano 2D)

---

## Tipos centrales

| Tipo | Descripción |
|------|-------------|
| `EEGRecording` | Señal EEG continua (channels × samples) con metadata |
| `EpochSet` | Epochs segmentados (channels × samples × epochs) + máscara AR |
| `ICAResult` | Matrices de mezcla/desmezclado, activaciones, componentes rechazados, varianza |
| `SpectralResult` | PSD + potencia por banda por canal |
| `ConnectivityMatrix` | Matrices wPLI por banda de frecuencia |
| `PipelineConfig` | Configuración completa del pipeline (cargada desde TOML) |
| `RecordingMeta` | Metadatos: sujeto, sesión, fs, nombres de canales, posiciones |

---

## Dependencias

| Paquete | Uso |
|---------|-----|
| `Genie` | Servidor web del dashboard |
| `CairoMakie` | Figuras PNG (señal, PSD, heatmaps wPLI) |
| `DSP` | Filtros Butterworth, `filtfilt` |
| `FFTW` | FFT + transformada de Hilbert para wPLI |
| `DataFrames` + `CSV` | Tablas de resultados |
| `LinearAlgebra` + `Random` | FastICA (implementación propia) |
| `Statistics` + `StatsBase` | Estadísticas por canal |
| `TOML` | Lectura de configuración |
| `Serialization` | Cache de resultados intermedios |

---

## Tests

```bash
julia --project=. tests/runtests.jl
```

Tests incluidos: tipos, filtrado, segmentación, PSD, wPLI, FDR, configuración.

---

## Reproducibilidad

- Toda la configuración versionable en `config/single_subject.toml`
- `seed = 42` configurable para FastICA
- Log completo del pipeline en `results/.../pipeline_log.txt`
- Snapshot de la configuración copiado junto con los resultados
- Cache versionado: detecta incompatibilidades entre versiones del framework
