# NeuroMIND

**Framework de conectividad funcional EEG basado en wPLI**
Conectividad funcional en Esclerosis Múltiple
Rafael Castro Triguero
*Última modificación: 22 Julio 2026*

NeuroMIND toma señales EEG de reposo en formato BrainVision, las procesa de principio a fin con un pipeline reproducible de 8 pasos y genera matrices de conectividad **weighted Phase Lag Index (wPLI)** con inferencia estadística opcional. Está diseñado para el dataset **MINDEM-IMIBIC** (41 pacientes con EM + 37 controles sanos, sesiones T1/T2, condiciones ojos cerrados / abiertos) y produce resultados listos para comparar grupos y sesiones.

> **Configuración única.** Desde el 21 Julio 2026 todo el proyecto se controla desde un solo fichero, [`config/pipeline.toml`](config/pipeline.toml).

---

## Tabla de contenidos

**Parte I — Guía de uso**

1. [¿Qué hace el pipeline?](#1-qué-hace-el-pipeline)
2. [Contexto científico](#2-contexto-científico)
3. [Requisitos e instalación](#3-requisitos-e-instalación)
4. [Cadena de ejecución: los 7 scripts](#4-cadena-de-ejecución-los-7-scripts)
5. [Configuración: `config/pipeline.toml`](#5-configuración-configpipelinetoml)
6. [Los 8 pasos del pipeline](#6-los-8-pasos-del-pipeline)
7. [Salidas del pipeline](#7-salidas-del-pipeline)
8. [Estructura de `results/` y estado actual](#8-estructura-de-results-y-estado-actual)
9. [Política de calidad de señal (QC)](#9-política-de-calidad-de-señal-qc)
10. [Dashboard e informe](#10-dashboard-e-informe)
11. [Estructura del proyecto y dependencias](#11-estructura-del-proyecto-y-dependencias)
12. [Reglas de Git y changelog](#12-reglas-de-git-y-changelog)

**Parte II — Verificación detallada (caso de referencia sub-M05)**

Recorrido fase por fase con los valores numéricos reales de la ejecución del 2026-07-09, que reproduce exactamente la configuración vigente. Ver [Anexo](#anexo-caso-de-referencia-sub-m05).

---

## 1. ¿Qué hace el pipeline?

```
Fase A/B (preparación)   audit_full_dataset.jl  →  build_bids_full.jl
        │
        ▼
Señal EEG cruda (.vhdr + .eeg BrainVision, o TSV BIDS)
        │
        ▼
[1/8] Carga y validación BIDS
[2/8] Control de calidad (QC) de canales
[3/8] Filtrado (notch + bandreject + highpass + lowpass)
[4/8] ICA — eliminación de artefactos (señal continua, ANTES de segmentar)
[5/8] Segmentación (1 s) + baseline + rechazo de artefactos ±70 µV
[6/8] Espectro de potencia (PSD) por banda
[7/8] Conectividad wPLI  [Hilbert | FourierCSD | Multitaper]
        └─ opcional: surrogates + FDR (si [surrogates] enabled = true)
[8/8] Guardado de tablas, figuras, índices QC y snapshot de configuración
        │
        ▼
Tablas CSV · Figuras PNG · Dashboard · Análisis de grupo · Informe PDF
```

**Regla crítica:** el paso 4 (ICA) se ejecuta siempre sobre la **señal continua filtrada** y **antes** de segmentar (paso 5). El orden no es negociable.

---

## 2. Contexto científico

| Concepto | Descripción |
|----------|-------------|
| **EEG de reposo** | Señal continua, 31 canales 10-20, ~500 Hz, condiciones EC (ojos cerrados) / EO (ojos abiertos) |
| **wPLI** | *Weighted Phase Lag Index*: conectividad funcional entre pares de canales por banda; robusto frente a conducción de volumen y ruido de amplitud |
| **dwPLI** | Variante *debiased* (Vinck 2011): elimina el sesgo por número variable de épocas. **Disponible** (`use_dwpli = true`), pero la configuración vigente usa wPLI clásico — ver [§5](#5-configuración-configpipelinetoml) |
| **ICA** | *Independent Component Analysis*: elimina artefactos oculares / musculares / cardíacos |
| **CSD** | *Current Source Density*: reduce conducción de volumen; `use_csd = false` por defecto |
| **Hipótesis** | Los pacientes con EM muestran alteraciones de conectividad en bandas α y β frente a controles sanos |

**Bandas de frecuencia analizadas:**

| Banda | Rango | Relevancia |
|-------|-------|-----------|
| δ (Delta) | 0.5 – 4 Hz | Sueño, baja vigilancia |
| θ (Theta) | 4 – 8 Hz | Memoria de trabajo, cognición |
| α (Alpha) | 7.8 – 11.7 Hz | Reposo, inhibición cortical |
| β_low | 12 – 15 Hz | Control motor, atención sostenida |
| β_mid | 15 – 18 Hz | Actividad sensoriomotora |
| β_high | 18 – 30 Hz | Procesos cognitivos de alto nivel |
| γ (Gamma) | 30 – 50 Hz | Procesamiento sensorial integrado |

> **Nota sobre Delta.** Con épocas de 1 s (perfil `eeg_julia` vigente), δ (0.5 Hz) acumula solo 0.5 ciclos/época, por debajo del mínimo `min_cycles_for_wpli = 4.0`. Dispara un `@warn` en cada lanzamiento pero no se excluye salvo que se active `exclude_unreliable_bands = true`.

---

## 3. Requisitos e instalación

- **Julia** 1.9 – 1.12
- Datos BrainVision en `data/full_data/`
- Para el informe: LaTeX (`latexmk` + XeLaTeX)

```bash
# 1. Comprobar Julia
julia --version                # 1.9.x o superior

# 2. Instalar dependencias (una vez)
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 3. Comprobar que el módulo carga
julia --project=. -e 'include("src/NeuroMIND.jl"); println("OK")'
```

---

## 4. Cadena de ejecución: los 7 scripts

Todos los lanzadores activos leen la misma configuración: `config/pipeline.toml`.

```bash
# A. Preparación del dataset (una vez por dataset)
julia --project=. scripts/audit_full_dataset.jl
julia --project=. scripts/build_bids_full.jl

# B. Pipeline por sujeto (recomendado antes del lote)
julia --project=. scripts/run_single_subject.jl
julia --project=. scripts/run_single_subject.jl --config config/pipeline.toml --force

# C. Lote completo
julia --project=. scripts/run_batch_pipeline.jl
julia --project=. scripts/run_batch_pipeline.jl --condition EC --group MS --session T1
julia --project=. scripts/run_batch_pipeline.jl --subjects M11,M12 --skip-done --dry-run

# D. Análisis de grupo (post-hoc, tras el lote)
julia --project=. scripts/run_transversal_analysis.jl
julia --project=. scripts/run_longitudinal_analysis.jl

# E. Inspección interactiva
julia --project=. scripts/launch_dashboard.jl --port 8080
```

**Orden de dependencia:**

```
audit_full_dataset → build_bids_full → run_single_subject (validación)
                                     → run_batch_pipeline (dataset completo)
                                            ↓
                       run_transversal_analysis / run_longitudinal_analysis
                                            ↓
                                     launch_dashboard
```

### Estado de los scripts

| Script | Estado | Fase | Orquestador | Salida principal |
|--------|--------|------|-------------|------------------|
| `audit_full_dataset.jl` | ✅ Activo | A — inventario | autónomo | `inventory.csv`, `participants.tsv`, `groups.csv`, `longitudinal_pairs.csv` |
| `build_bids_full.jl` | ✅ Activo | B — metadata BIDS | autónomo | `*_eeg_metadata.json`, `electrodes.tsv`, `dataset_description.json` |
| `run_single_subject.jl` | ✅ Activo | pipeline individual | `SingleSubjectPipeline.jl` (8 pasos) | `results/subjects/sub-{ID}/ses-{SES}/{task}/` |
| `run_batch_pipeline.jl` | ✅ Activo | C — lote | idem, en bucle sobre `inventory.csv` | misma ruta BIDS + `results/logs/batch_run_*.csv` |
| `run_transversal_analysis.jl` | ✅ Activo | post-hoc grupal | autónomo | `results/transversal/{EC\|EO}/` |
| `run_longitudinal_analysis.jl` | ✅ Activo | post-hoc longitudinal | autónomo | `results/longitudinal/{EC\|EO}/` |
| `launch_dashboard.jl` | ✅ Activo | visualización | `webapp/App.jl` (Genie) | `http://localhost:8080` |

**Opciones de `run_batch_pipeline.jl`** (combinables): `--condition` EC\|EO\|ALL · `--group` MS\|HC\|ALL · `--session` T1\|T2\|ALL · `--subjects` lista · `--max-subjects` N · `--skip-done` · `--dry-run`.

---

## 5. Configuración: `config/pipeline.toml`

Fichero único que controla todo el pipeline. Cada grabación guarda además su propio `config_snapshot.toml`: **esa es la fuente fiable** de con qué parámetros se produjo un resultado, por encima de lo que diga el fichero de configuración en un momento dado.

En modo lote, `run_batch_pipeline.jl` copia este fichero y sobreescribe `[subject]` y `[paths]` por cada trabajo; el resto de secciones se propaga sin cambios.

### Catálogo de secciones

| Sección | Parámetros y decisiones vigentes |
|---------|----------------------------------|
| `[subject]` | `subject_id` / `session_id` / `task` / `run`. **Solo aplica en modo single.** En lote, `run_batch_pipeline.jl` los sobreescribe desde `inventory.csv`. `subject_id = "auto"` autodetecta el primer `*_eeg_data.tsv` de `{bids_root}/raw/`. |
| `[recording]` | `sampling_rate = 500.0` (fallback; el fs real se lee del `.vhdr`). `reference = "average"` (informativo — no re-referencia la señal). |
| `[qc]` | `bad_channel_zscore_threshold = 3.0` · `amplitude_warning_sigma_uv = 20.0`. |
| `[filtering]` | `profile = "eeg_julia"`: Notch 50 Hz (49.5–50.5) y Bandreject 99.5–100.5 Hz **causales** (`filt`); HP 0.5 Hz y LP 150 Hz **zero-phase** (`filtfilt`). Orden 4. |
| `[ica]` | `profile = "eeg_julia"` → `n_components = nº canales`, `max_iter = 512`, `tol = 1e-7`, `seed = 1234` (fijados en código; **no** se leen del TOML con este perfil). ⚠️ `artifact_threshold` se parsea pero está **inerte**: el umbral real está hardcodeado a 1.5 en `SingleSubjectPipeline.jl`. |
| `[segmentation]` | `profile = "eeg_julia"` → épocas de 1.0 s sin solape. ⚠️ Los nombres de clave del TOML (`segment_length_seconds`, `overlap_seconds`, `min_segments`) los **traduce** `load_ss_config` a `epoch_length_s`/`epoch_overlap`/`min_epochs`. `min_segments = 10` es criterio de **exclusión dura** en la tabla QC. |
| `[baseline]` | `method = "first_window_mean"` (media de 0–100 ms), `n_passes = 2` (antes y después del AR). Alternativas: `"mean"`, `"median"`. ⚠️ `baseline_start_s` no afecta al cálculo (siempre arranca en t=0). |
| `[artifact_rejection]` | `profile = "eeg_julia"`, ±70 µV, sin gradiente. **`n_channels_used = 31`** (todos los canales, coherente con `exclude_fp2 = false`). ⚠️ Es posicional. `before_event_ms`/`after_event_ms` se guardan pero no se aplican. |
| `[spectral]` | FFT con ventana Hamming-taper, `nfft = 1024` (~0.488 Hz/bin, 513 bins), `window_pct = 10.0`. |
| `[bands]` | Las 7 bandas. ⚠️ Intervalo semiabierto `[flo, fhi)` e independiente por banda → hay solape THETA/ALPHA (7.8–8.0 Hz) y hueco ALPHA/BETA_LOW (11.7–12.0 Hz). La suma de bandas no iguala la potencia total. |
| `[connectivity]` | **`wpli_method = "hilbert"`** · **`use_dwpli = false`** (wPLI clásico, rango 0–1) · `filter_order = 8` · `use_csd = false` · `min_cycles_for_wpli = 4.0` · `exclude_unreliable_bands = false`. Sub-tablas `[connectivity.fourier_csd]` y `[connectivity.multitaper]` para los otros métodos. |
| `[montage]` | **`exclude_fp2 = false`** → Fp2 SE CONSERVA: **31 canales, 465 aristas** por banda. ⚠️ Quien elimina canales es `exclude_channels`; `exclude_fp2` solo filtra Fp2 de esa lista. `n_channels_analysis` es informativo (se calcula, no se lee). |
| `[graph]` | `density = 0.1`, `threshold_method = "proportional"`. Solo afecta a métricas binarias (path_length, efficiency, degree); strength y clustering usan la matriz completa. |
| `[surrogates]` | `enabled = false`. Si se activa: `n_surrogates = 200`, `alpha = 0.05`, `fdr_method = "bh"`, `seed = 42`. ⚠️ `method` es inerte (siempre `circular_shift`). |
| `[output]` | `figure_format = "png"`, `figure_dpi = 150`. |
| `[paths]` | `bids_root = "data/bids"` (minúscula), `results = "results"`. |
| `[dashboard]` | `port = 8080`, `open_browser = true`. Solo lo lee `launch_dashboard.jl`; `--port` en CLI tiene prioridad. |

---

## 6. Los 8 pasos del pipeline

### [1/8] Carga EEG
`load_single_subject` lee el `.vhdr`/`.eeg` BrainVision (o TSV BIDS) → `EEGRecording` (canales × muestras). `validate_channels` contrasta los nombres contra `electrodes.tsv`. No genera salida propia.

### [2/8] Control de calidad de canales
`compute_channel_stats` + `flag_bad_channels` (z-score de RMS > 3.0). Calcula `amplitude_warning` si σ̄ de la señal cruda > 20 µV. Salida: `channel_statistics.csv` / `qc_summary.csv` (columna `is_bad`).

### [3/8] Filtrado
`filter_recording`, perfil `eeg_julia`: Notch 50 Hz causal → Bandreject 99.5–100.5 Hz causal → HP 0.5 Hz zero-phase → LP 150 Hz zero-phase (Butterworth, orden 4). Salida: `filtered_signal_preview.png`.

### [4/8] ICA
`run_ica` — FastICA simétrica con PCA whitening, en Julia puro (solo `LinearAlgebra` + `Random`). Con `profile = "eeg_julia"` usa todos los canales; matriz de mezcla `A = inv(W_total)` (invertible por ser cuadrada). `load_ica_labels` + `apply_ica_rejection` eliminan los componentes marcados. Salidas: `ica_summary.json`, `ica_components.csv`, `ica_component_features.csv`, matrices de mezcla/separación, activaciones, `raw_signal.csv`, `figures/ica_topomap_NNN.png`.

> El rechazo de componentes **no es automático por score**: `apply_ica_rejection` actúa sobre los índices de un CSV de inspección manual (`ica_labels.csv`). Sin ese CSV se conservan todos los componentes.

### [5/8] Segmentación + baseline + rechazo de artefactos
`segment_recording` (épocas de 1.0 s sin solape) → `apply_baseline` (`first_window_mean`, pre-AR) → `reject_artifacts` (±70 µV sobre los 31 canales) → segunda pasada de baseline (`n_passes = 2`). Salidas: `segmentation_summary.json`, `artifact_rejection_summary.json`, `segments_table.csv`, `channel_coverage.csv`, `rejected_segments.csv`, `channel_artifact_summary.csv`.

### [6/8] Espectro de potencia (PSD)
`compute_psd` — FFT con ventana Hamming-taper. Salidas: `spectral_summary.json`, `band_power_summary.csv`, `psd_by_channel.csv`, `regional_psd.csv`, `spectral_indices.csv`, `psd_all_channels.png`, `band_power_summary.png`.

### [7/8] Conectividad wPLI
`compute_wpli` con el estimador `hilbert` (Butterworth bandpass + señal analítica Hilbert). Con 31 canales → **465 aristas** por banda. `compute_graph_metrics` deriva métricas de red. Salidas: `connectivity_summary.json`, `network_metrics.csv`, `connectivity_edges.csv`, `wpli_{banda}.csv`, `wpli_{banda}.png`.

**Surrogates (opcional, `enabled = true`):** `surrogate_test` genera la distribución nula por *circular shift* con el mismo estimador que el observado, y `fdr_correction` aplica Benjamini-Hochberg. Salidas: `wpli_{observed,pvalues,qvalues,significant}_{banda}.csv`, `surrogate_null_stats_{banda}.csv`, `significant_connections.csv`, `surrogate_{quality,summary}` y `surrogate_null_{banda}.png`.

### [8/8] Guardado de resultados
`_save_all_results` + `_save_config_snapshot`. Actualiza dos índices globales: `results/subjects_index.csv` y `results/qc/qc_decision_table.csv`. Salidas por sujeto: `overview.csv`, `config_snapshot.toml`, `pipeline_log.txt`.

---

## 7. Salidas del pipeline

### Árbol de salida único (BIDS)

Desde el 2026-07-21 cada grabación se escribe **una sola vez**, en el árbol BIDS. El antiguo árbol heredado con sufijos (`results/{ID}/{SES}/tables/overview_EC.csv`) se eliminó: `_save_all_results` escribe ahora directamente en `export_dir` sin copias.

```
results/subjects/sub-{ID}/ses-{SES}/{task}/
├── overview.csv · channel_statistics.csv · qc_summary.csv
├── segmentation_summary.json · segments_table.csv · channel_coverage.csv · …
├── ica_summary.json · ica_components.csv · ica_*.csv · raw_signal.csv
├── spectral_summary.json · psd_by_channel.csv · band_power_summary.csv · …
├── connectivity_summary.json · network_metrics.csv · connectivity_edges.csv
├── wpli_{banda}.csv                  ← matriz por banda (nombre canónico)
├── (surrogates opcionales: wpli_{pvalues,qvalues,significant}_{banda}.csv, …)
├── config_snapshot.toml · pipeline_log.txt
├── cache/ica_result.jls              ← caché ICA (antes en el árbol heredado)
└── figures/                          ← todas las figuras (sin sufijo _EC/_EO)
```

El nombre BIDS ya distingue la condición por el nivel `{task}` (`eyesclosed`/`eyesopen`), así que los sufijos `_EC`/`_EO` eran redundantes. Las tablas de edges por banda (`wpli_edges_{banda}`) se descartaron: su información está en `connectivity_edges.csv` (todas las bandas, con `rank`).

### Análisis de grupo

Ambos scripts leen los `wpli_{banda}.csv` ya generados y escriben una carpeta por condición (39 ficheros por condición con 7 bandas):

| Análisis | Salida | Ficheros clave |
|----------|--------|----------------|
| **Transversal** (MS vs Control) | `results/transversal/{EC\|EO}/` | `group_{connectivity_ms,connectivity_control,difference,statistics}_{banda}.csv`, `significant_edges_{banda}.csv`, `band_statistics.csv`, `subject_inclusion.csv`, `transversal_summary.json` |
| **Longitudinal** (T1 → T2) | `results/longitudinal/{EC\|EO}/` | `longitudinal_{connectivity_t1,connectivity_t2,difference,statistics}_{banda}.csv`, `significant_longitudinal_edges_{banda}.csv`, `paired_subjects.csv`, `longitudinal_summary.json` |

El transversal requiere `data/bids/groups.csv`; el longitudinal usa `longitudinal_pairs.csv` o autodetecta los pares T1/T2 en `results/subjects/`.

---

## 8. Estructura de `results/` y estado actual

Reorganizada el 2026-07-21. Cada modo de ejecución tiene un área propia.

```
results/
├── README.md
├── subjects/                run_single_subject.jl · run_batch_pipeline.jl
│   └── sub-{ID}/ses-{T1|T2}/{eyesclosed|eyesopen}/
├── transversal/             run_transversal_analysis.jl
│   └── {EC|EO}/
├── longitudinal/            run_longitudinal_analysis.jl
│   └── {EC|EO}/
├── qc/qc_decision_table.csv
├── logs/batch_run_{timestamp}.csv
└── subjects_index.csv

deprecated/                  fuera del árbol de producción
├── code/                    ← VERSIONADO en git (ficheros pequeños)
│   ├── config/  scripts/  src/    (run_pipeline.jl, Pipeline.jl, configs antiguas)
└── results/                 ← IGNORADO por git (GB de derivados EEG)
    └── 2026-05-26_pre-unificacion/
```

El pipeline individual y el lote escriben en la **misma** ruta: reprocesar un sujeto suelto sobrescribe su carpeta y nada más.

### Estado de los resultados

En producción solo está el caso de referencia `subjects/sub-M05/ses-T2/eyesclosed/` (reprocesado el 2026-07-09), cuyo snapshot coincide con la configuración vigente. Las 204 grabaciones restantes y los análisis de grupo de mayo están archivados porque cuatro parámetros cambian los valores numéricos:

| Parámetro | Corridas de mayo (archivadas) | Configuración vigente |
|-----------|-------------------------------|-----------------------|
| `connectivity.use_dwpli` | `true` — dwPLI, rango −1…1 | `false` — wPLI clásico, 0…1 |
| `baseline.method` | `mean` — media de la época | `first_window_mean` — 0–100 ms |
| `montage.exclude_fp2` | `true` — 30 canales, 435 aristas | `false` — 31 canales, 465 aristas |
| `artifact_rejection.n_channels_used` | `30` | `31` |

Las matrices wPLI tienen además dimensión distinta (30×30 vs 31×31): no son comparables. Orden previsto de regeneración: **M05 (validación) → lote completo → transversal y longitudinal**.

---

## 9. Política de calidad de señal (QC)

### Montaje de 31 canales
Con la configuración vigente (`exclude_fp2 = false`) el análisis trabaja con los **31 canales**, incluido Fp2, produciendo matrices wPLI de **31 × 31** con **465 aristas por banda**. Fp2 fue históricamente problemático (artefactos oculares/contacto en el 46% del dataset); se conserva porque la ICA (paso 4) corre sobre los 31 canales y le resta la activación del componente ocular dominante antes de la conectividad. `n_channels_used = 31` garantiza que Fp2 también pase la criba de artefactos ±70 µV. Para un análisis de sensibilidad sin Fp2: `exclude_fp2 = true` y `n_channels_used = 30`.

### Alerta de amplitud (`amplitude_warning`)
Si σ̄ de la señal cruda > 20 µV → probable grabación sin filtro online activo. **No excluye automáticamente**; la decisión final se toma tras el AR.

### Tabla de decisión QC (`results/qc/qc_decision_table.csv`)

| `final_decision` | Condición | ¿Incluir en grupo? |
|------------------|-----------|--------------------|
| `include` | Sin alertas, ≥10 épocas válidas | ✅ Sí |
| `include_with_warning` | `amplitude_warning` o ≥1 canal malo | ✅ Con cautela |
| `manual_review` | `amplitude_warning` + ≥2 canales malos, **o** <50% épocas válidas | ⚠️ Revisar |
| `exclude` | 0 épocas válidas, o menos de `min_segments` (10) | ❌ No |

---

## 10. Dashboard e informe

### Dashboard web (Genie.jl) — 16 paneles

```bash
julia --project=. scripts/launch_dashboard.jl --port 8080   # → http://localhost:8080
```

| Panel | Contenido | | Panel | Contenido |
|-------|-----------|--|-------|-----------|
| 0 | Proyecto / Dataset | | 8 | Análisis espectral |
| 1 | BIDS y metadata | | 9 | Conectividad wPLI |
| 2 | Señal cruda | | 10 | Surrogates / Inferencia |
| 3 | QC inicial | | 11 | Resultados finales |
| 4 | Preprocesado / Filtrado | | 12 | Exportación / Informe |
| 5 | ICA | | 13 | Evaluación transversal |
| 6 | Segmentación | | 14 | Evaluación longitudinal |
| 7 | Rechazo de artefactos | | 15 | Validación MNE-Python |

### Informe PDF
El informe científico LaTeX vive en `report/` (`main_es.tex` → `report/build/pdf/main_es.pdf`). Requiere `latexmk` + XeLaTeX.

### Validación cruzada (mne_brain/)
`mne_brain/` reimplementa el pipeline en MNE-Python para validar los resultados de forma independiente. Ver `mne_brain/README.md`.

---

## 11. Estructura del proyecto y dependencias

```
NeuroMIND/
├── config/
│   └── pipeline.toml           ← ⭐ CONFIGURACIÓN ÚNICA (la leen los 7 scripts)
├── data/
│   ├── bids/                   ← BIDS ligero (metadata JSON, sin señales)
│   └── full_data/              ← .vhdr crudos + inventory.csv
├── src/                        ← Código Julia
│   ├── NeuroMIND.jl            ← Entry point del módulo
│   ├── types.jl                ← EEGRecording, EpochSet, ICAResult, …
│   ├── SingleSubjectPipeline.jl ← Pipeline de 8 pasos (canónico) + load_ss_config
│   ├── io/                     ← BrainVisionLoader, BIDSLoader, Config
│   ├── preprocessing/          ← Filtering
│   ├── ica/                    ← ICACore (FastICA), ICAClassification, ICAInspection
│   ├── segmentation/           ← Epochs (segment, baseline, AR)
│   ├── spectral/               ← PowerSpectrum
│   ├── connectivity/           ← wPLI, CSD, GraphMetrics
│   ├── statistics/             ← Surrogates, FDR, GroupStats
│   ├── visualization/          ← Topomaps, Heatmaps, Spectra
│   ├── report/                 ← HTMLReport
│   └── webapp/                 ← App.jl (dashboard Genie)
├── scripts/                    ← 7 lanzadores activos (ver §4)
├── results/                    ← Generado por el pipeline (NO en git)
├── deprecated/                 ← Archivado (code/ en git, results/ ignorado)
├── mne_brain/                  ← Pipeline de validación MNE-Python
├── report/                     ← Informe científico LaTeX
├── tests/ · test/              ← runtests.jl
├── config/pipeline.toml · Project.toml · Manifest.toml
├── README.md · AGENTS.md · CLAUDE.md
└── NeuroMIND Claude Code/      ← Referencia de arquitectura (PDF)
```

**Orquestador en `src/`:**

| Módulo | Usado por | Estado |
|--------|-----------|--------|
| `SingleSubjectPipeline.jl` | `run_single_subject.jl`, `run_batch_pipeline.jl` | ✅ Canónico — 8 pasos, BrainVision, surrogates, export BIDS |

### Dependencias (`Project.toml`)

| Paquete | Uso |
|---------|-----|
| `DSP` | Filtros Butterworth, `filtfilt`, DPSS (multitaper) |
| `FFTW` | Transformadas para PSD, Hilbert y espectro cruzado |
| `CairoMakie` | Figuras PNG (topomaps, heatmaps, espectros) |
| `Genie` | Servidor web del dashboard |
| `CSV`, `DataFrames` | Lectura/escritura de tablas |
| `StatsBase`, `Statistics` | Estadística descriptiva y correlaciones |
| `TOML` | Parseo de `config/pipeline.toml` |
| `Serialization` | Caché de ICA (`.jls`) |

---

## 12. Reglas de Git y changelog

### Nunca versionar

```bash
# Antes de cualquier commit, verificar que no se cuelan datos:
git ls-files | grep -E '(^data/|^results/|^deprecated/results/|\.DS_Store$|^\.claude/|^\.vscode/)'
# → debe devolver vacío
```

`data/`, `results/`, `deprecated/results/`, `reports/`, `.claude/`, `.vscode/`, `.cursor/`, `.DS_Store`, `.env`, logs clínicos y secretos. `deprecated/code/` **sí** se versiona (ficheros pequeños, histórico útil).

**Identidad:** Rafael Castro Triguero · `me1catrr@uco.es`. Nunca commitear directamente en `main`; usar ramas `feat/<nombre>`.

### Changelog

| Fecha | Cambios |
|-------|---------|
| **2026-07-21** | Configuración unificada en `config/pipeline.toml`; `run_pipeline.jl` archivado; `results/` reorganizado (subjects/transversal/longitudinal); archivo movido a `deprecated/`; montaje a 31 canales; wPLI clásico |
| **2026-07-09** | Reprocesado del caso de referencia M05 con la configuración actual |
| **2026-05-26** | wPLI multi-método (Hilbert / FourierCSD / Multitaper) |
| **2026-05-25** | dwPLI, surrogates válidos, GraphMetrics, QC v2 |
| **2026-05-24** | Dataset completo, BrainVision loader, pipeline en lote |
| **2026-05-23** | Dashboard completo (paneles 0–14), corrección PSD |

---
## Anexo: Caso de referencia sub-M05

> Verificación rutina a rutina del pipeline individual: código, tablas y figuras para `sub-M05 / ses-T2 / eyesclosed` (EC). Asume que has leído la [Parte I](#1-qué-hace-el-pipeline), en especial las secciones [4](#4-cadena-de-ejecución-los-7-scripts) y [6](#6-los-8-pasos-del-pipeline).

Sujeto de trabajo para verificar rutina a rutina el código, las tablas y las figuras generadas por el pipeline individual.

| Campo | Valor |
|-------|-------|
| **Sujeto** | `sub-M05` |
| **Sesión** | `ses-T2` |
| **Tarea** | `eyesclosed` (ojos cerrados, condición **EC**) |
| **Último lanzamiento** | **2026-07-09** (inicio 10:49:50, duración 521.7 s) |
| **Comando** | `julia --project=. scripts/run_single_subject.jl --config <toml>` (ver `pipeline_log.txt` para el TOML exacto) |
| **Traza** | `results/subjects/sub-M05/ses-T2/eyesclosed/pipeline_log.txt` |
| **Config aplicada** | `results/subjects/sub-M05/ses-T2/eyesclosed/config_snapshot.toml` (copia fiel del TOML usado) |

> Esta fecha y esta config deben actualizarse cada vez que se relance el pipeline sobre M05.
> El `config/pipeline.toml` tiene `[surrogates] enabled = false` por defecto; la ejecución del 2026-07-09 activó surrogates para regenerar el ejemplo con inferencia estadística.

### Directorio de salida

Salida **única** en el árbol BIDS (desde 2026-07-21; ver [§7](#7-salidas-del-pipeline)). Todo lo de M05 vive en una sola carpeta:

```
results/subjects/sub-M05/ses-T2/eyesclosed/
  overview.csv · channel_statistics.csv · *.csv · *.json   ← tablas y resúmenes
  wpli_{banda}.csv · connectivity_edges.csv               ← conectividad
  config_snapshot.toml · pipeline_log.txt                 ← parámetros y traza
  figures/                                                ← todas las figuras (sin sufijo)
  cache/ica_result.jls                                    ← caché ICA (interno, no citar)
```

Los nombres siguen BIDS (`sub-{ID}/ses-{SES}/{task}/`) y no llevan sufijo `_EC`/`_EO`: el nivel `{task}` (`eyesclosed`/`eyesopen`) ya distingue la condición.

| Uso | Ruta |
|-----|------|
| Informe, revisión científica, citas en `Report_Pre/` | `results/subjects/sub-M05/ses-T2/eyesclosed/` (raíz o `figures/`) |
| Dashboard web (paneles 0–15) | misma ruta — el dashboard lee de `results/subjects/` |
| Caché ICA (no versionar, no citar) | `results/subjects/sub-M05/ses-T2/eyesclosed/cache/ica_result.jls` |

> **Nota histórica.** La corrida del 2026-07-09 se hizo con el código previo, que además escribía un árbol heredado `results/M05/T2/` con sufijos `_EC` y copiaba al BIDS. Ese doble árbol se eliminó; las corridas actuales producen solo la carpeta de arriba.

### Fase 0 — Preparación del dataset (verificado M05)

La Fase 0 es **previa al pipeline de 8 pasos**. No procesa señal EEG; solo inventaría los `.vhdr` crudos y crea la metadata BIDS ligera que alimenta el paso 1/8 (`load_single_subject`).

```
data/full_data/Pacientes MINDEM_IMIBIC_…/*.vhdr
        │
        ▼
[Fase A] scripts/audit_full_dataset.jl
        │  inventory.csv, participants.tsv, groups.csv, longitudinal_pairs.csv
        ▼
[Fase B] scripts/build_bids_full.jl
        │  *_eeg_metadata.json, electrodes.tsv, dataset_description.json
        ▼
[Paso 1/8] load_single_subject → lee metadata JSON (+ .vhdr o TSV si existe)
```

**Comandos (ejecutar una vez por dataset, o al añadir grabaciones nuevas):**

```bash
julia --project=. scripts/audit_full_dataset.jl
julia --project=. scripts/build_bids_full.jl
```

**Salidas globales (sin figuras):**

| Script | Fichero | Ubicación |
|--------|---------|-----------|
| Fase A | `inventory.csv` | `data/full_data/inventory.csv` |
| Fase A | `participants.tsv` | `data/bids/participants.tsv` |
| Fase A | `groups.csv` | `data/bids/groups.csv` |
| Fase A | `longitudinal_pairs.csv` | `data/bids/longitudinal_pairs.csv` |
| Fase B | `*_eeg_metadata.json` | `data/bids/raw/` (uno por grabación válida) |
| Fase B | `*_electrodes.tsv` | `data/bids/electrodes/` (uno por sujeto/sesión) |
| Fase B | `dataset_description.json` | `data/bids/dataset_description.json` |

> La Fase B **no copia** los binarios BrainVision (`.eeg`/`.vhdr`); el JSON apunta al `.vhdr` original vía `vhdr_path`.

**Comprobación M05 — `eyesclosed` / ses-T2 (2026-07-13):**

| Comprobación | Resultado |
|--------------|-----------|
| Filas en `inventory.csv` | 3: T1 excluida (paradigma ODDBALL), T2 EC ✅, T2 EO ✅ |
| `bids_id` | `M05` (ID crudo del estudio: `M5`) |
| `participants.tsv` | `sub-M05` · grupo MS · sexo F · edad 51 · `has_t1=false` · `has_t2=true` |
| `longitudinal_pairs.csv` | `include_longitudinal=false` (sin T1 resting EC+EO; solo T2 válido) |
| Metadata BIDS | `data/bids/raw/sub-M05_ses-T2_task-eyesclosed_run-01_eeg_metadata.json` |
| `vhdr` fuente | `…/M5_T2_MLLERGAS_011221_OJOS CERRADOS.vhdr` — **existe en disco** |
| `electrodes.tsv` | `data/bids/electrodes/sub-M05_ses-T2_electrodes.tsv` — 31 canales EEG + REF/GND |
| Parámetros en metadata | `fs=500 Hz`, `n_channels=31`, `task=eyesclosed`, `group=MS` |
| `dataset_description.json` | Presente en `data/bids/` |
| Figuras Fase 0 | Ninguna (esperado) |

**Hallazgos a tener en cuenta:**

1. **Carga real en paso 1/8:** existe un `sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv` (~16 MB, 31 ch × 50 180 muestras). `load_single_subject` prioriza el TSV si está presente, aunque el metadata diga `data_format=brainvision`. Para M05 EC, el pipeline del 2026-07-09 cargó por **TSV**, no por `.vhdr` directo.
2. **Ruta absoluta en metadata:** el campo `vhdr_path` del JSON guarda una ruta de otro usuario (`/Users/rafa/…`). El inventario tiene la ruta actual correcta; conviene regenerar metadata si se migra de máquina.
3. **T1 de M05:** la grabación T1 es paradigma cognitivo (ODDBALL), marcada `excluded=true` — no entra al pipeline de reposo.

**Estado Fase 0 M05:** ✅ preparación verificada — listo para el paso 1/8.

### Fase 1/8 — Carga de datos (verificado M05)

El paso 1/8 lee la señal cruda y valida el montaje contra `electrodes.tsv`. **No genera ficheros propios** en el momento de la carga; deja un objeto `EEGRecording` en memoria que alimenta todos los pasos siguientes.

**Script de entrada:** `run_single_subject.jl` → `run_single_subject_pipeline()` (`src/SingleSubjectPipeline.jl`)

```
config/pipeline.toml  →  load_ss_config()
        │
        ▼
[subject] subject_id / session_id / task / run
        │  ("auto" → detect_first_subject en data/bids/raw/)
        ▼
load_single_subject()  ──►  EEGRecording  (canales × muestras)
        │                    ├─ ruta TSV si existe *_eeg_data.tsv
        │                    └─ ruta BrainVision si no hay TSV (load_eeg_brainvision)
        ▼
validate_channels()  ──►  compara nombres de canal vs electrodes.tsv
```

**Rutinas y módulos:**

| Rutina | Archivo | Función |
|--------|---------|---------|
| `load_ss_config` | `src/SingleSubjectPipeline.jl` | Parsea el TOML → `PipelineConfig` |
| `detect_first_subject` | `src/SingleSubjectPipeline.jl` | Solo si `subject_id = "auto"`: primer `*_eeg_data.tsv` en `data/bids/raw/` |
| `load_single_subject` | `src/SingleSubjectPipeline.jl` | Orquesta carga TSV o BrainVision |
| `load_eeg_brainvision` | `src/io/BrainVisionLoader.jl` | Lee `.vhdr` + `.eeg` (IEEE_FLOAT_32, escala µV) |
| `validate_channels` | `src/SingleSubjectPipeline.jl` | Cruza nombres de canal con `electrodes.tsv` |

**Entradas M05 (desde Fase 0):**

| Fichero | Rol |
|---------|-----|
| `data/bids/raw/sub-M05_ses-T2_task-eyesclosed_run-01_eeg_metadata.json` | `fs`, `n_channels`, `data_format` |
| `data/bids/raw/sub-M05_ses-T2_task-eyesclosed_run-01_eeg_data.tsv` | **Señal usada** (único TSV del dataset; 205 metadata JSON, 1 TSV) |
| `data/bids/electrodes/sub-M05_ses-T2_electrodes.tsv` | Posiciones 10-20 para validación y topomaps |

**Salidas relacionadas con la carga** (escritas en pasos posteriores, no en 1/8):

| Fichero | Escrito en | Contenido |
|---------|------------|-----------|
| `overview.csv` | Paso 8/8 | Ficha resumen: canales, muestras, duración, validación electrodos |
| `raw_signal.csv` | Paso 4/8 | Serie temporal cruda completa (31 ch × 100.36 s) — copia de lo cargado en 1/8 |
| `pipeline_log.txt` | Toda la ejecución | Líneas `[1/8] Carga de datos` |

> El paso 1/8 **no produce figuras**.

**Comprobación M05 — ejecución 2026-07-09:**

| Comprobación | Esperado | Verificado |
|--------------|----------|------------|
| Canales | 31 | ✅ `pipeline_log.txt` y `overview.csv` |
| Muestras | 50 180 | ✅ |
| Frecuencia de muestreo | 500 Hz | ✅ |
| Duración | 100.36 s | ✅ (50 180 / 500) |
| Validación electrodos | OK | ✅ `electrode_validation = Validación OK` |
| Formato de carga | TSV (prioridad sobre BrainVision) | ✅ único `*_eeg_data.tsv` del dataset |
| Nombres de canal | Fz…Fp2 (orden 10-20) | ✅ coinciden en TSV, `raw_signal.csv` y `electrodes.tsv` |
| Primer valor Fz (t=0 s) | −7.1946 µV | ✅ TSV ≡ `raw_signal.csv` |
| `subject_id = "auto"` | Detecta M05 | ✅ es el único TSV en `data/bids/raw/` |

**Traza en `pipeline_log.txt`:**

```
[1/8] Carga de datos
  Canales: 31 | Muestras: 50180 | fs: 500.0 Hz
  Duración: 100.36 s
  Electrodos: Validación OK
```

**Hallazgos a tener en cuenta:**

1. **M05 es un caso especial de carga:** de 205 grabaciones con metadata JSON, solo M05 tiene además un `*_eeg_data.tsv` legacy. El resto del cohorte cargará por **BrainVision** (`BrainVisionLoader.jl`) cuando se ejecute el batch.
2. **Prioridad TSV > BrainVision:** si existe `*_eeg_data.tsv`, `load_single_subject` lo usa aunque el JSON diga `data_format=brainvision` (comportamiento en `SingleSubjectPipeline.jl` L217–227).
3. **`detect_first_subject` solo busca TSV:** con `subject_id="auto"` fallará si no hay ningún `*_eeg_data.tsv` en `data/bids/raw/`. Para otros sujetos, fijar `subject_id = "M05"` (o el ID concreto) en el TOML.
4. **`raw_signal.csv` no es salida del paso 1/8:** se escribe en el paso 4/8 (`_save_raw_signal`) como copia de auditoría de la señal cruda cargada al inicio.

**Estado Fase 1/8 M05:** ✅ carga y validación verificadas — listo para el paso 2/8 (QC).

### Fase 2/8 — QC de canales (verificado M05)

El paso 2/8 analiza la **señal cruda** (`EEGRecording` del paso 1/8) antes de filtrar. Calcula estadísticos por canal, marca canales sospechosos y evalúa si la amplitud global sugiere grabación sin filtro online.

```
EEGRecording (raw)
        │
        ▼
compute_channel_stats()   →  DataFrame por canal (media, RMS, σ, rango, z-score RMS)
        │
        ▼
flag_bad_channels()       →  canales con |z-score RMS| > umbral (defecto 3.0)
        │
        ▼
amplitude_warning         →  σ̄_raw > 20 µV  →  aviso, no exclusión automática
```

**Rutinas y módulos (código actual en repo):**

| Rutina | Archivo | Función |
|--------|---------|---------|
| `compute_channel_stats` | `src/qc/QualityControl.jl` | 6 métricas básicas por canal |
| `flag_bad_channels` | `src/qc/QualityControl.jl` | Detección por z-score de RMS |
| Paso 2/8 (orquestación) | `src/SingleSubjectPipeline.jl` L355–379 | Añade `is_bad`, calcula `amplitude_warning`, escribe log |
| Guardado tablas QC | `src/SingleSubjectPipeline.jl` L1324–1328 (`_save_all_results`, **paso 8/8**) | Copia `qc_stats` → `qc_summary.csv` + `channel_statistics.csv` |

**Parámetros M05** (`config_snapshot.toml` → `[qc]`):

| Parámetro | Valor | Efecto |
|-----------|-------|--------|
| `bad_channel_zscore_threshold` | 3.0 | Umbral para marcar canal malo |
| `amplitude_warning_sigma_uv` | 20.0 (defecto código) | Media de σ por canal; si supera → aviso |

**Comprobación M05 — ejecución 2026-07-09:**

| Comprobación | Resultado |
|--------------|-----------|
| σ̄_raw (media de σ por canal) | **14.4 µV** → `amplitude_warning = false` |
| Canales malos (z > 3.0) | **Fp2** (`rms_zscore = 4.63`, `std_uv = 54.6`) |
| Resto de canales | `is_bad = false` (ej. Cz: `std_uv = 7.7`, `rms_zscore = −0.77`) |
| Decisión QC global (paso 8/8) | `include` en `results/qc/qc_decision_table.csv` |

**Traza en `pipeline_log.txt`:**

```
[2/8] QC de canales
  amplitude_warning: false (σ̄_raw=14.4 µV)
  Canales sospechosos (z>3.0): Fp2
  Figuras QC: amplitude_histograms, psd_raw_average, raw_butterfly, variance_topomap  (18.0 s)
```

**Salidas en disco (M05):**

| Fichero | Tipo | Escrito en | Notas |
|---------|------|------------|-------|
| `channel_statistics.csv` | Tabla (31 filas) | Paso 8/8 | **Duplicado byte a byte** de `qc_summary.csv` |
| `qc_summary.csv` | Tabla (31 filas) | Paso 8/8 | Misma tabla |
| `overview.csv` | Tabla (1 fila) | Paso 8/8 | Incluye `bad_channels=Fp2`, `amplitude_warning=false` |
| `qc_psd_raw_mean.csv` | Tabla auxiliar | Paso 2/8 (ext.) | PSD media Welch de la señal raw |
| `qc_band_power_raw.csv` | Tabla auxiliar | Paso 2/8 (ext.) | Potencia por banda (raw) |
| `figures/qc_amplitude_histograms.png` | Figura | Paso 2/8 (ext.) | Rejilla 31 histogramas de amplitud |
| `figures/qc_amplitude_histograms_p01.png` | Figura | Paso 2/8 (ext.) | Página 1 de histogramas (16 ch) |
| `figures/qc_amplitude_histograms_p02.png` | Figura | Paso 2/8 (ext.) | Página 2 (15 ch) |
| `figures/qc_raw_butterfly.png` | Figura | Paso 2/8 (ext.) | 31 canales, 0–10 s |
| `figures/qc_psd_raw_average.png` | Figura | Paso 2/8 (ext.) | PSD media + bandas sombreadas |
| `figures/qc_variance_topomap.png` | Figura | Paso 2/8 (ext.) | Topografía de `std_uv` |

> **(ext.)** = salidas de la **versión extendida de QC** usada en la ejecución del 2026-07-09. El `pipeline_log.txt` confirma que se generaron (18 s). Ver hallazgo 1 abajo.

**Columnas en `channel_statistics.csv` (M05, 15 columnas):**

`channel`, `mean_uv`, `rms_uv`, `std_uv`, `range_uv`, `rms_zscore`, `min_uv`, `max_uv`, `skewness`, `kurtosis`, `hfnoise`, `snr_db`, `mean_abs_corr`, `max_corr`, `min_corr`, `is_bad`

El módulo actual `QualityControl.jl` solo calcula las **6 primeras** (+ `is_bad` en el pipeline). Las columnas 7–14 (mín/máx, asimetría, curtosis, HFNoise, SNR, correlaciones) provienen del bloque extendido de la ejecución verificada.

**Hallazgos a tener en cuenta:**

1. **Desfase código ↔ salidas M05:** el `src/qc/QualityControl.jl` actual solo implementa QC básico (6 métricas). Las figuras QC y las 15 columnas de `channel_statistics.csv` en disco fueron generadas por una versión del pipeline con QC extendido (referenciada en `config_snapshot.toml` → `[qc.figures]`, `welch_nfft`, etc.) que **ya no está cableada** en el `SingleSubjectPipeline.jl` del repo. Si se relanza hoy, las tablas saldrían con **7 columnas** (6 + `is_bad`) y **sin figuras QC**.
2. **Fp2 marcado malo pero no excluido aún:** el QC solo **etiqueta**; la exclusión de Fp2 del montaje de análisis ocurre más adelante (antes del paso 5/8, política `exclude_fp2` + `bad_ch`).
3. **`channel_statistics.csv` = `qc_summary.csv`:** duplicación en `_save_all_results` — hallazgo menor pendiente de limpiar.
4. **`channel_statistics_compare.csv` / `_filtered.csv`:** generados aparte (walkthrough, 2026-07-09 11:07) — **no** son salida estándar del pipeline; comparan raw vs filtrado.

**Estado Fase 2/8 M05:** ✅ lógica de QC verificada en log y tablas — ⚠️ figuras y columnas extendidas dependen de reintegrar QC extendido en código para reproducirlas en un relanzamiento.

### Fase 3/8 — Filtrado (verificado M05)

El paso 3/8 transforma la señal cruda en señal filtrada (`rec_filt`) que alimenta ICA (paso 4/8). **No escribe tablas CSV**; el efecto se materializa en memoria y en figuras generadas más tarde.

```
EEGRecording (raw)
        │
        ▼
filter_recording()          →  EEGRecording (filtrada), misma meta y tiempos
        │
        ├─ describe_filter_chain()  →  log de la cadena (solo traza)
        └─ rec_filt  ─────────────────►  entrada del paso 4/8 (ICA)
```

**Rutinas y módulos:**

| Rutina | Archivo | Función |
|--------|---------|---------|
| `filter_recording` | `src/preprocessing/Filtering.jl` | Cadena Butterworth según perfil |
| `describe_filter_chain` | `src/preprocessing/Filtering.jl` | Lista orden/frecuencias/método para el log |
| `_plot_signal_preview` | `src/SingleSubjectPipeline.jl` L1465 | Raw vs filtrado (figura, **paso 8/8**) |
| Kernel `_filt` | `src/preprocessing/Filtering.jl` | `filt` (causal) o `filtfilt` (fase cero) vía `DSP` |

**Perfil M05:** `eeg_julia` (reproducibilidad con `EEG_Julia/src/Preprocessing/filtering.jl`)

| Paso | Filtro | Banda | Orden | Método |
|------|--------|-------|-------|--------|
| 1 | Notch | 49.5–50.5 Hz | 4 | `filt` (causal) |
| 2 | Bandreject | 99.5–100.5 Hz | 4 | `filt` (causal) |
| 3 | High-pass | 0.5 Hz | 4 | `filtfilt` (fase cero) |
| 4 | Low-pass | 150.0 Hz | 4 | `filtfilt` (fase cero) |

**Parámetros** (`config_snapshot.toml` → `[filtering]`): `profile=eeg_julia`, `filter_order=4`, `notch_bw_hz=1.0`, `bandreject_lo/hi=99.5/100.5`.

**Comprobación M05 — ejecución 2026-07-09:**

| Comprobación | Resultado |
|--------------|-----------|
| Perfil activo | `eeg_julia` ✅ |
| Cadena en log | 4 pasos, frecuencias y métodos coinciden con `Filtering.jl` |
| Canales tras filtrar | 31 (misma geometría que raw) |
| Efecto en Fp2 (σ raw → filtrado) | 54.6 µV → **13.8 µV** (`channel_statistics_compare.csv`, walkthrough) |
| Efecto en Cz (σ raw → filtrado) | 7.7 µV → **3.0 µV** |
| Tiempo de ejecución | < 1 s (mismo timestamp en log: 10:50:11) |

**Traza en `pipeline_log.txt`:**

```
[3/8] Filtrado
  Perfil: eeg_julia
    [1] Notch 49.5–50.5 Hz  ord=4  método=filt
    [2] Bandreject 99.5–100.5 Hz  ord=4  método=filt
    [3] High-pass 0.5 Hz  ord=4  método=filtfilt
    [4] Low-pass 150.0 Hz  ord=4  método=filtfilt
```

**Salidas en disco (M05):**

| Fichero | Tipo | Escrito en | Contenido |
|---------|------|------------|-----------|
| `filtered_signal_preview.png` | Figura | Paso 8/8 | 5 primeros canales (Fz, F3, F7, FT9, FC5), 10 s, raw gris vs filtrado azul |
| `results/subjects/sub-M05/ses-T2/eyesclosed/figures/filtered_signal_preview.png` | Figura | Paso 8/8 | Preview señal cruda vs. filtrada |
| `config_snapshot.toml` | Config | Paso 8/8 | Parámetros `[filtering]` aplicados |
| `ica_signal_before.csv` | Serie temporal | Paso 4/8 | Primeros 10 s de **señal filtrada** pre-ICA (derivada de `rec_filt`) |

> El paso 3/8 **no genera CSV propio**. Los parámetros de la cascada solo quedan en `config_snapshot.toml` y en `pipeline_log.txt`.

**Figuras de walkthrough (no estándar del pipeline):**

| Fichero | Origen | Uso |
|---------|--------|-----|
| `walkthrough_M05_filt_compare_Cz.png` | Script manual / walkthrough | Comparación raw vs filtrado en Cz |
| `walkthrough_M05_filt_compare_F3.png` | idem | Comparación en F3 |
| `walkthrough_M05_psd_raw_vs_filt.png` | idem | PSD antes/después del filtrado |
| `channel_statistics_filtered.csv` | Walkthrough (2026-07-09 11:07) | Estadísticos post-filtro por canal |

**Hallazgos a tener en cuenta:**

1. **Orden del perfil importa:** con `eeg_julia`, Notch y Bandreject son **causales** (`filt`); HP/LP son **fase cero** (`filtfilt`). El perfil `default` usa `filtfilt` en todos — no son equivalentes.
2. **La figura canónica se guarda en el paso 8/8**, no en el 3/8: `_plot_signal_preview(rec, rec_filt)` corre al final, cuando ya existen `rec` y `rec_filt`.
3. **Fp2 sigue siendo el canal más variable tras filtrar** (σ = 13.8 µV vs ~3–8 µV en la mayoría), coherente con el flag QC del paso 2/8.
4. **No confundir con filtrado wPLI:** el `[connectivity] filter_order = 8` del TOML es solo para el estimador Hilbert del paso 7/8; es independiente de esta cadena de preprocesado.

**Estado Fase 3/8 M05:** ✅ filtrado verificado en log, config y figura — listo para el paso 4/8 (ICA).

### Fase 4/8 — ICA (verificado M05)

> **Regla crítica:** ICA se aplica sobre la **señal continua filtrada**, **antes** de segmentar (paso 5/8). No alterar este orden.

El paso 4/8 descompone la señal filtrada en componentes independientes, clasifica artefactos y reconstruye la señal limpia (`rec_ica`) que alimenta la segmentación.

```
rec_filt (continua, filtrada)
        │
        ├─ run_ica()  [o caché ica_result.jls]  →  ICAResult (31 comp.)
        │
        ├─ load_ica_labels()  →  índices a rechazar (manual o auto*)
        │
        ├─ apply_ica_rejection()  →  rec_ica (señal limpia)
        │
        └─ _save_ica_results()  →  tablas + topomaps + figuras ICA
```

\* En la ejecución M05 del 2026-07-09 el rechazo fue **automático** (`artifact_threshold = 1.5`). Ver hallazgo 1.

**Rutinas y módulos:**

| Rutina | Archivo | Función |
|--------|---------|---------|
| `run_ica` | `src/ica/ICACore.jl` | FastICA simétrico + PCA whitening (`tanh`, perfil `eeg_julia`) |
| `load_ica_labels` | `src/ica/ICAInspection.jl` | Lee `ica_labels.csv` manual (si existe) |
| `apply_ica_rejection` | `src/ica/ICAInspection.jl` | Reconstruye señal sin componentes rechazados |
| `compute_ica_features` | `src/ica/ICAClassification.jl` | 7 features por componente |
| `evaluate_ica_components` | `src/ica/ICAClassification.jl` | Scores + etiqueta (`brain`, `jump`, `line_noise`, …) |
| `_save_ica_results` | `src/SingleSubjectPipeline.jl` L1539 | Exporta tablas, topomaps y resumen JSON |
| Caché | `results/subjects/sub-M05/ses-T2/eyesclosed/cache/ica_result.jls` | Evita re-ejecutar FastICA si la config ICA no cambia |

**Parámetros M05** (`config_snapshot.toml` → `[ica]`):

| Parámetro | Valor | Notas |
|-----------|-------|-------|
| `profile` | `eeg_julia` | 31 componentes (= nº canales), `seed=1234`, `max_iter=512`, `tol=1e-7` |
| `artifact_threshold` | 1.5 | Umbral de score para rechazo automático (ejecución 2026-07-09) |
| Matriz de mezcla | `A = inv(W_total)` | Perfil `eeg_julia` (no `pinv`) |

**Comprobación M05 — ejecución 2026-07-09:**

| Comprobación | Resultado |
|--------------|-----------|
| Componentes estimados | 31 |
| ICA desde caché | ✅ (`ica_result.jls`, 0.9 s) |
| Componentes rechazados | **5, 9, 13, 14, 21** |
| Tipos de artefacto | 2 × `jump` (IC5, IC13) · 3 × `line_noise` (IC9, IC14, IC21) |
| Varianza retenida | **82.6 %** (17.4 % eliminada) |
| Señal before ≠ after | ✅ `ica_signal_before.csv` ≠ `ica_signal_after.csv` |
| IC de mayor varianza aceptada | IC7 (12.0 %) — `brain` |

**Componentes rechazados (detalle):**

| IC | `variance_pct` | `artifact_type` | Motivo |
|----|----------------|-----------------|--------|
| IC5 | 11.44 % | jump | Mayor contribución a varianza eliminada |
| IC9 | 3.13 % | line_noise | |
| IC13 | 1.42 % | jump | |
| IC14 | 0.46 % | line_noise | |
| IC21 | 0.97 % | line_noise | |

**Traza en `pipeline_log.txt`:**

```
[4/8] ICA
  ICA cargado desde caché
  Componentes: 31 | Método: FastICA + PCA whitening
  Componentes rechazados (auto, umbral=1.5): 5, 9, 13, 14, 21
  Señal limpiada · 5 componente(s) eliminado(s)
  Duración ICA: 0.9 s
  Topomaps: 31/31 guardados en figures/
  Detalle visual de 5 componente(s) rechazado(s): 10 figuras
  Headplot varianza antes/después: guardado
  Butterfly antes/después: guardado
```

**Salidas tablas (M05):**

| Fichero | Contenido |
|---------|-----------|
| `ica_summary.json` | `n_comp=31`, `n_rejected=5`, `variance_retained=82.6`, tipos de artefacto |
| `ica_components.csv` | 31 filas: `variance_pct`, `rejected`, `artifact_type` |
| `ica_component_features.csv` | 31 filas × 7 features + scores + clasificación |
| `ica_mixing_matrix.csv` | Matriz **A** (31 ch × 31 IC) |
| `ica_unmixing_matrix.csv` | Matriz **W** (31 IC × 31 ch) |
| `ica_activations.csv` | Traza temporal IC1–IC31 (primeros 10 s) |
| `ica_signal_before.csv` | Señal filtrada pre-limpieza (10 s, 31 ch) |
| `ica_signal_after.csv` | Señal ICA-limpiada (10 s, 31 ch) |
| `raw_signal.csv` | Señal cruda completa (~100 s) — guardada aquí por conveniencia |

**Salidas figuras (`figures/`, M05):**

| Patrón | Cantidad | Descripción |
|--------|----------|-------------|
| `ica_topomap_NNN.png` | 31 | Topografía de cada componente |
| `ica_rejected_IC*_timecourse.png` | 5 | Traza temporal de IC rechazados |
| `ica_rejected_IC*_spectrum.png` | 5 | Espectro de IC rechazados |
| `ica_before_after_butterfly.png` | 1 | 31 canales antes/después |
| `ica_headplot_before.png` / `_after.png` | 2 | Mapa de varianza por canal |
| `ica_artifact_scores.png` | 1 | Scores de clasificación por componente |

**Hallazgos a tener en cuenta:**

1. **Desfase código ↔ ejecución M05 (rechazo automático):** el `config_snapshot.toml` define `artifact_threshold = 1.5` y el log confirma rechazo automático, pero el `SingleSubjectPipeline.jl` **actual** solo llama a `load_ica_labels()` (CSV manual). La lógica de auto-rechazo que aplicó M05 el 2026-07-09 **no está cableada** en el código vigente — pendiente de reintegrar (referencia en config: `ICACleaning.jl`, archivo ausente).
2. **Figuras ICA extendidas:** butterfly, headplots y detalle de rechazados (`ica_rejected_*`) se generaron en la ejecución verificada pero **no están** en el `_save_ica_results` actual (solo topomaps + tablas). Mismas figuras en disco = versión anterior del pipeline.
3. **Caché ICA:** `results/subjects/sub-M05/ses-T2/eyesclosed/cache/ica_result.jls` + `ica_config.hash`. Si cambian parámetros `[ica]`, se invalida y se re-ejecuta FastICA (~varios segundos).
4. **`ica_components.csv` vs rechazo real:** la columna `rejected` refleja `ica.rejected_components` al guardar; la clasificación en `ica_component_features.csv` puede etiquetar más ICs como artefacto sin rechazarlos si el umbral no se aplica en código.

**Estado Fase 4/8 M05:** ✅ ICA verificada en log, tablas y figuras en disco — ⚠️ rechazo automático y figuras extendidas requieren reintegración en código para reproducir en un relanzamiento.

### Fase 5/8 — Segmentación, baseline y AR (verificado M05)

El paso 5/8 opera sobre `rec_ica` (señal ICA-limpiada). Parte la señal continua en épocas, corrige baseline, rechaza artefactos por amplitud y vuelve a corregir baseline en las épocas válidas.

```
rec_ica (continua, ICA-limpiada, 31 ch)
        │
        ├─ [montaje] exclude_fp2=false  →  31 canales activos (config M05)
        │
        ▼
segment_recording()     →  100 épocas × 1.0 s × 500 muestras
        │
        ▼
apply_baseline()  [pasada 1]   →  ventana 0–100 ms, método first_window_mean
        │
        ▼
reject_artifacts()        →  ±70 µV, 31 canales, sin gradiente
        │
        ▼
apply_baseline()  [pasada 2]   →  misma ventana, solo épocas válidas
        │
        ▼
epochs (99 válidas)  →  PSD (paso 6) y wPLI (paso 7)
```

**Rutinas y módulos:**

| Rutina | Archivo | Función |
|--------|---------|---------|
| `segment_recording` | `src/segmentation/Epochs.jl` | Corta señal continua en `EpochSet` |
| `apply_baseline` | `src/segmentation/Epochs.jl` | Corrección por canal/época |
| `reject_artifacts` | `src/segmentation/Epochs.jl` | Rechazo ±µV (perfil `eeg_julia`) |
| `compute_epoch_quality_report` | `src/segmentation/Epochs.jl` | Tabla de calidad por época (pre-AR) |
| `compute_channel_coverage` | `src/segmentation/Epochs.jl` | % de épocas válidas por canal |
| `_save_segmentation_results` | `src/SingleSubjectPipeline.jl` L639 | JSON + CSVs de segmentación y AR |

**Parámetros M05** (`config_snapshot.toml`):

| Bloque | Parámetro | Valor |
|--------|-----------|-------|
| `[segmentation]` | `profile` | `eeg_julia` → épocas de **1.0 s**, sin solape |
| `[baseline]` | `method` | `first_window_mean` (0–100 ms) |
| `[baseline]` | `n_passes` | **2** (antes y después del AR) |
| `[artifact_rejection]` | `profile` | `eeg_julia` → ±70 µV |
| `[artifact_rejection]` | `n_channels_used` | **31** (actualizado 2026-07-08; antes 30) |
| `[montage]` | `exclude_fp2` | **false** — Fp2 permanece en el montaje tras limpieza ICA |

**Comprobación M05 — ejecución 2026-07-09:**

| Comprobación | Resultado |
|--------------|-----------|
| Épocas totales | 100 (100.36 s ÷ 1 s, paso 1 s) |
| Épocas válidas | **99** (99.0 %) |
| Épocas rechazadas | **1** (época 12, t = 11–12 s) |
| Motivo rechazo | Amplitud en **C4**: max = 72.56 µV > +70 µV |
| Calidad media / mediana | 0.708 / 0.715 |
| Cobertura C4 | 99.0 % (único canal afectado) |
| Entrada a segmentación | `ICA-limpiada` (`segmentation_summary.json`) |
| Duración paso 5 | 0.1 s |

**Época rechazada (detalle):**

| Campo | Valor |
|-------|-------|
| `epoch` | 12 |
| `start_s` / `end_s` | 11.0 / 12.0 |
| `worst_channel` | C4 |
| `max_amp_uv` | 72.56 |
| `p2p_uv` | 104.18 |
| `quality` | 0.482 |

**Traza en `pipeline_log.txt`:**

```
[5/8] Segmentación
  Perfil segmentación: eeg_julia | Perfil AR: eeg_julia | Baseline passes: 2
  Segmentos totales: 100 | Válidos: 99 | Rechazados: 1 (99.0%)
  Duración segmentación: 0.1 s
  Guardado: figuras de segmentación (figures/)
```

**Salidas tablas (M05):**

| Fichero | Contenido |
|---------|-----------|
| `segmentation_summary.json` | Resumen global: perfiles, conteos, histograma calidad, cobertura |
| `segments_table.csv` | 100 filas — estado de cada época (calidad, motivo, p2p, peor canal) |
| `channel_coverage.csv` | 31 filas — `coverage_pct` por canal |
| `rejected_segments.csv` | 1 fila — época 12 / C4 |
| `channel_artifact_summary.csv` | 1 fila — C4: 1 época mala (1.0 %) |
| `artifact_rejection_summary.json` | Umbrales, histograma p2p, estadísticos (p2p_mean=72.4 µV) |

**Salidas figuras (`figures/`, M05):**

| Fichero | Descripción |
|---------|-------------|
| `epoch_quality_histogram.png` | Distribución de `quality` (pico en 0.6–0.7) |
| `epochs_overlay_Oz.png` | Superposición de épocas válidas en Oz |
| `epoch_012_rejected_stacked.png` | Época rechazada (stacked por canal) |
| `baseline_before_after_C4_ep089.png` | Efecto baseline en C4, época 89 |
| `p2p_histogram.png` | Histograma pico-a-pico (diagnóstico, no criterio de rechazo) |

> Las figuras de segmentación se generaron en la ejecución del 2026-07-09 (`pipeline_log` lo confirma) pero **no están** en el `_save_segmentation_results` actual del repo — mismo patrón que QC/ICA extendido (ver hallazgo 1).

**Hallazgos a tener en cuenta:**

1. **Desfase código ↔ figuras:** el `_save_segmentation_results` vigente guarda tablas JSON/CSV pero **no** llama a funciones de figura (`epoch_figures.jl` referenciado en el informe, ausente en `src/visualization/`). Las 5 figuras en disco proceden de la ejecución verificada.
2. **Fp2 incluido en montaje M05:** `exclude_fp2=false`, que es la decisión vigente (desde 2026-07-08, tras la limpieza ICA). Con `exclude_fp2=true` el montaje sería de 30 canales (análisis de sensibilidad).
3. **`n_channels_used=31` en AR:** actualizado para que Fp2 (canal 31) entre en el chequeo de amplitud cuando permanece en el montaje.
4. **Baseline doble pasada:** la 2.ª pasada recalcula la media 0–100 ms solo sobre épocas que sobrevivieron al AR — coherente con EEG_Julia.
5. **Pico a pico (p2p):** se registra en `artifact_rejection_summary.json` como estadístico diagnóstico; el criterio de rechazo es solo amplitud instantánea ±70 µV, no p2p.

**Estado Fase 5/8 M05:** ✅ segmentación, baseline y AR verificados en log y tablas — ⚠️ figuras de épocas requieren reintegración en código para reproducir en un relanzamiento.

### Fase 6/8 — Análisis espectral (PSD) (verificado M05)

El paso **6/8** calcula el PSD sobre las **99 épocas válidas** post-baseline y post-AR. El **guardado en disco** de tablas y figuras espectrales ocurre en el paso **8/8** (`_save_all_results` + `_save_spectral_extras`), no en el 6/8.

```
epochs (99 válidas × 1.0 s × 31 ch)
        │
        ▼
compute_psd()          →  FFT + ventana Hamming-taper (10 %)
        │                 rFFT nfft=1024 (zero-padding desde 500 muestras)
        │                 promedio sobre épocas → PSD (ch × 513 bins)
        ▼
_band_power_from_psd() →  potencia media por banda (7 bandas)
        │
        ▼
SpectralResult  →  wPLI (paso 7) + guardado (paso 8)
```

**Rutinas y módulos:**

| Rutina | Archivo | Función |
|--------|---------|---------|
| `compute_psd` | `src/spectral/PowerSpectrum.jl` | PSD por canal (µV²/Hz), promedio inter-época |
| `_band_power_from_psd` | `src/spectral/PowerSpectrum.jl` | Integración por banda (media de bins) |
| `plot_spectrum_grid` | `src/visualization/Spectra.jl` | Grid PSD todos los canales |
| `_plot_band_power_summary` | `src/SingleSubjectPipeline.jl` L1498 | Barplot potencia por banda |
| `_save_spectral_extras` | `src/SingleSubjectPipeline.jl` L916 | `spectral_summary.json`, `regional_psd.csv`, `spectral_indices.csv` |
| `_save_all_results` | `src/SingleSubjectPipeline.jl` L1314 | `psd_by_channel.csv`, `band_power_summary.csv`, figuras |

**Parámetros M05** (`config_snapshot.toml`):

| Bloque | Parámetro | Valor |
|--------|-----------|-------|
| `[spectral]` | `nfft` | **1024** (≥ 500 muestras/época → zero-padding) |
| `[spectral]` | `window_pct` | **10.0** % (Hamming-taper en extremos) |
| `[bands]` | DELTA … GAMMA | 0.5–4, 4–8, 7.8–11.7, 12–15, 15–18, 18–30, 30–50 Hz |
| Entrada | `n_epochs_used` | **99** (épocas válidas tras AR) |
| Entrada | `n_channels` | **31** (montaje con Fp2) |

**Comprobación M05 — ejecución 2026-07-09:**

| Comprobación | Resultado |
|--------------|-----------|
| Resolución espectral | **0.488 Hz** (fs/nfft = 500/1024) |
| Bins de frecuencia | **513** (0–250 Hz) |
| Potencia media por banda (31 ch) | DELTA **3.64** · THETA **1.83** · ALPHA **3.37** · β_low **0.86** · β_mid **0.86** · β_high **1.00** · GAMMA **0.19** µV² |
| Potencia relativa (media bandas) | DELTA **31.0 %** · ALPHA **28.7 %** · THETA **15.6 %** · β **23.1 %** · GAMMA **1.6 %** |
| Potencia total integrada | **57.56 µV²** (`spectral_summary.json`) |
| Canal pico ALPHA | **O2** (11.76 µV²); pico espectral α en O2: **10.25 Hz**, 21.19 µV²/Hz |
| Canal pico DELTA | **C4** (21.70 µV²) — coherente con época 12 rechazada en AR |
| Duración paso 6/8 | < 0.1 s (solo cálculo, sin I/O) |

**Traza en `pipeline_log.txt` (paso 6/8):**

```
[6/8] Espectral
  nfft=1024 | Bins: 513 | Resolución: 0.488 Hz
  Potencia     DELTA: 3.6429 μV² (media canales)
  Potencia     ALPHA: 3.3715 μV² (media canales)
  … (7 bandas)
```

**Salidas tablas (guardadas en paso 8/8, M05):**

| Fichero | Contenido |
|---------|-----------|
| `psd_by_channel.csv` | 15 903 filas — PSD largo: `channel × freq_hz × power_uv2` (31×513) |
| `band_power_summary.csv` | 31 filas × 7 bandas (µV² por canal) |
| `spectral_summary.json` | Metadatos FFT + estadísticos globales por banda + `total_power_uv2` |
| `regional_psd.csv` | 35 filas — media ± std por región (frontal/central/parietal/occipital) y banda |
| `spectral_indices.csv` | 31 filas — ratios α/θ, β/α, θ/β, γ/α + pico α (Hz, µV²/Hz) |

**Salidas figuras (paso 8/8, M05):**

| Fichero (árbol BIDS, actual) | Nombre en el log 2026-07-09 (histórico, con sufijo) | Descripción |
|------------------------------|------------------------------------------------------|-------------|
| `figures/psd_all_channels.png` | `psd_all_channels_EC.png` | Grid PSD 31 canales (0–50 Hz) |
| `figures/band_power_summary.png` | `band_power_summary_EC.png` | Barplot potencia media por banda |
| `figures/band_topomap_grid.png` | `band_topomap_grid_EC.png` | Topomapas por banda (7 paneles) |

> `band_topomap_grid_EC.png` aparece en el log del 2026-07-09 pero **no hay** llamada a topomapas espectrales en el `_save_all_results` actual — posible módulo extendido no presente en el repo (mismo patrón que figuras QC/segmentación).

**Índices espectrales destacados (M05):**

| Canal | α/θ | β/α | pico α (Hz) | Nota |
|-------|-----|-----|-------------|------|
| O2 | 2.87 | 0.16 | 10.25 | Máxima potencia α absoluta |
| Oz | 2.30 | 0.19 | 10.25 | Occipital, α dominante |
| C4 | 0.99 | 0.55 | 9.77 | Elevada δ (21.7 µV²); canal implicado en AR |

**Hallazgos a tener en cuenta:**

1. **Cálculo vs guardado:** el paso 6/8 solo ejecuta `compute_psd` y escribe en log; las tablas/figuras espectrales se persisten en el paso 8/8 junto con wPLI y QC.
2. **Normalización PSD:** `PowerSpectrum.jl` usa `Pseg ./= (n_samp * fs)` (corrección 2026-05-23; antes se dividía por `nfft²`).
3. **Ventana Hamming-taper:** estilo BrainVision Analyzer — 10 % de cada extremo de la época atenuado; corrección de varianza vía `mw2`.
4. **PSD crudo vs procesado:** `qc_psd_raw_mean.csv` y `qc_band_power_raw.csv` (Fase 2) miden la señal **sin filtrar**; los ficheros de esta fase reflejan la señal **ICA-limpiada, segmentada y con baseline**.
5. **C4 y DELTA:** C4 concentra la mayor potencia δ (21.7 µV²), coherente con el rechazo de la época 12 en el paso 5 (amplitud 72.6 µV en C4).
6. **Salida única:** las figuras se escriben una sola vez, sin sufijo, en `results/subjects/sub-M05/ses-T2/eyesclosed/figures/` (el doble árbol se eliminó el 2026-07-21).

**Estado Fase 6/8 M05:** ✅ PSD y potencia por banda verificados en log, `spectral_summary.json` y tablas — ⚠️ `band_topomap_grid` requiere verificar módulo de topomapas espectrales para reproducir en relanzamiento.

### Fase 7/8 — Conectividad wPLI y surrogates (verificado M05)

El paso **7/8** calcula matrices de conectividad funcional wPLI en espacio sensor (31×31, 465 aristas/banda) y, si `[surrogates] enabled = true`, ejecuta inferencia por permutaciones. En M05 los surrogates estuvieron **activados** (verificación puntual). El guardado masivo ocurre en el paso **8/8**.

```
epochs (99 válidas × 1.0 s × 31 ch)
        │
        ├─ use_csd=false  →  espacio sensor (sin CSD)
        │
        ▼
compute_wpli()           →  Hilbert + Butterworth ord=8, 7 bandas
        │                    wPLI (no dwPLI en M05)
        ▼
ConnectivityMatrix       →  7 matrices 31×31
        │
        ▼  [surrogates enabled]
surrogate_test() × 7     →  circular_shift, N=200, FDR-BH por banda
        │
        ▼
SurrogateResult[]        →  p/q-values, máscaras, distribución nula
```

**Rutinas y módulos:**

| Rutina | Archivo | Función |
|--------|---------|---------|
| `compute_wpli` | `src/connectivity/wPLI.jl` | wPLI/dwPLI multi-método (Hilbert / FourierCSD / Multitaper) |
| `_build_estimator` | `src/connectivity/wPLI.jl` | Construye estimador desde `[connectivity]` |
| `apply_csd` | `src/connectivity/CSD.jl` | CSD opcional (`use_csd=false` en M05) |
| `surrogate_test` | `src/statistics/Surrogates.jl` | Permutaciones + p-valores Monte Carlo (+1) |
| `_circular_shift_surrogate` | `src/statistics/Surrogates.jl` | Desplazamiento circular independiente por canal/época |
| `fdr_correction` | `src/statistics/FDR.jl` | Benjamini-Hochberg por banda |
| `_save_connectivity_extras` | `src/SingleSubjectPipeline.jl` L1071 | `connectivity_summary.json`, `network_metrics.csv` |
| `_save_surrogate_results` | `src/SingleSubjectPipeline.jl` L1161 | Tablas p/q, `significant_connections.csv`, `surrogate_summary.json` |
| `_save_all_results` | `src/SingleSubjectPipeline.jl` L1314 | `wpli_{band}.csv`, `connectivity_edges.csv`, heatmaps PNG |

**Parámetros M05** (`config_snapshot.toml`):

| Bloque | Parámetro | Valor |
|--------|-----------|-------|
| `[connectivity]` | `wpli_method` | **`hilbert`** (Butterworth ord **8** + señal analítica) |
| `[connectivity]` | `use_dwpli` | **false** (wPLI estándar, no debiased) |
| `[connectivity]` | `use_csd` | **false** |
| `[connectivity]` | `min_cycles_for_wpli` | **4.0** |
| `[connectivity]` | `exclude_unreliable_bands` | **false** (DELTA se calcula aunque < 4 ciclos) |
| `[surrogates]` | `enabled` | **true** (habilitado para regenerar ejemplo M05) |
| `[surrogates]` | `n_surrogates` | **200** |
| `[surrogates]` | `method` | `circular_shift` *(informativo — `surrogate_test` siempre usa circular_shift)* |
| `[surrogates]` | `alpha` | **0.05** |
| `[surrogates]` | `fdr_method` | **`bh`** |
| `[surrogates]` | `seed` | **42** |

**Comprobación M05 — wPLI (ejecución 2026-07-09):**

| Comprobación | Resultado |
|--------------|-----------|
| Espacio | **sensor** (sin CSD) |
| Matrices | **7 bandas** × **31×31** canales |
| Aristas por banda | **465** (31×30/2) |
| Aristas totales (`connectivity_edges.csv`) | **3 255** (7×465) |
| Duración wPLI | **~3 s** (10:50:24 → 10:50:27) |
| wPLI medio (por banda) | ALPHA **0.327** · THETA **0.133** · DELTA **0.121** · β_low **0.151** · β_mid **0.120** · β_high **0.083** · GAMMA **0.094** |
| wPLI máximo observado | ALPHA **0.858** (O1–CP6); GAMMA **0.264** (Fz–Fp2) |
| Top-5 ALPHA (wPLI) | O1–CP6, Oz–P4, Oz–CP6, O1–P4, O1–P8 (red occipital-parietal) |

**Comprobación M05 — surrogates:**

| Comprobación | Resultado |
|--------------|-----------|
| Método efectivo | **circular_shift** (independiente por canal/época) |
| Permutaciones | **200** × 7 bandas |
| Duración total surrogates | **~480 s** (10:50:27 → 10:58:27) — **92 %** del tiempo de pipeline |
| Pares significativos (FDR 5 %, por banda) | **ALPHA: 258** (55.5 %) · **GAMMA: 189** (40.7 %) · resto: **0** |
| Total conexiones significativas | **447** (`significant_connections.csv`) |
| Umbral FDR | ALPHA **0.0249** · GAMMA **0.0199** |
| Semilla | **42** (+ offset por índice de banda) |

**Traza en `pipeline_log.txt`:**

```
[7/8] Conectividad wPLI
  Espacio: sensor | Bandas: ALPHA, BETA_HIGH, BETA_LOW, BETA_MID, DELTA, GAMMA, THETA

[SUR] Inferencia por surrogates
  Método: circular_shift (fijo — Surrogates.jl::surrogate_test no lee cfg.surrogates["method"]) | N=200 | FDR: bh
  Banda     ALPHA: 258 pares significativos  FDR-thr=0.0249
  Banda     GAMMA: 189 pares significativos  FDR-thr=0.0199
  … (5 bandas con 0 significativos)
  Surrogates guardados: 447 conexiones significativas en 7 bandas
```

**Salidas tablas wPLI (paso 8/8, M05):**

| Fichero | Contenido |
|---------|-----------|
| `wpli_{band}.csv` | Matriz simétrica 31×31 por banda (7 ficheros) |
| `connectivity_edges.csv` | 3 255 filas — aristas ordenadas por wPLI descendente |
| `connectivity_summary.json` | Estadísticos por banda (media, std, densidad > umbral 0.1) |
| `network_metrics.csv` | 31 filas — `strength`, `degree`, `norm_strength` agregados sobre bandas |

**Salidas tablas surrogates (M05):**

| Fichero | Contenido |
|---------|-----------|
| `wpli_observed_{band}.csv` | Matriz observada (duplicado de `wpli_{band}.csv`) |
| `wpli_pvalues_{band}.csv` | p-valores Monte Carlo (+1) |
| `wpli_qvalues_{band}.csv` | q-valores FDR-BH |
| `wpli_significant_{band}.csv` | Máscara binaria 31×31 |
| `surrogate_null_stats_{band}.csv` | 465 filas/banda — obs, null_mean, null_std, p, q, z |
| `significant_connections.csv` | 447 filas — conexiones con q < 0.05 |
| `surrogate_quality.csv` | 7 filas — QC por banda (% sig, media nula, etc.) |
| `surrogate_summary.json` | Resumen global de inferencia |

**Salidas figuras (paso 8/8, M05):**

| Fichero | Descripción |
|---------|-------------|
| `wpli_{band}.png` | Heatmap matriz wPLI (7 bandas) |
| `surrogate_null_{band}.png` | Distribución nula agregada por banda (7 bandas) |

> Los heatmaps wPLI **sí** se generan en `_save_all_results` actual. Los PNG `surrogate_null_*.png` están en disco pero **no hay** llamada en el `_save_surrogate_results` vigente — posible extensión de la ejecución del 2026-07-09.

**Métricas de red destacadas (M05):**

| Canal | Strength | Degree | Nota |
|-------|----------|--------|------|
| C4 | 5.86 (max) | 19 | Canal con mayor δ espectral; implicado en AR |
| Fp2 | 5.27 | 20 | Incluido en montaje M05 (`exclude_fp2=false`) |
| FT10 | 5.04 | 20 | Alta conectividad en GAMMA (Fz–FT10: 0.255) |

**Hallazgos a tener en cuenta:**

1. **Tiempo de cómputo:** con `n_surrogates=200` y 7 bandas, los surrogates dominan el coste (~8 min de 8.7 min totales en M05). Para el batch completo conviene evaluar `enabled=false` o reducir N.
2. **DELTA poco fiable:** con épocas de 1.0 s, DELTA tiene **0.5 ciclos/época** (< `min_cycles_for_wpli=4.0`); se calcula igualmente porque `exclude_unreliable_bands=false`. Interpretar con cautela.
3. **FDR por banda:** la corrección BH se aplica **independientemente** en cada banda (465 tests/banda), no globalmente sobre las 3 255 aristas.
4. **`method` en config no se lee:** `surrogate_test` siempre usa `circular_shift`; la clave `[surrogates].method` es solo informativa (el log del 2026-07-09 lo documenta explícitamente).
5. **Alta proporción de significativos en α y γ:** 55 % y 41 % de pares en reposo EC es esperable en análisis within-subject con FDR permisivo; la comparación grupal a nivel de grupo se hace con análisis de segundo nivel; `use_dwpli=true` está disponible como alternativa (la config vigente usa wPLI clásico).
6. **CSD desactivado:** valores wPLI no son directamente comparables con EEG_Julia (que aplica CSD antes de wPLI).
7. **Fp2 en conectividad:** con la config vigente (`exclude_fp2=false`) participa en las 465 aristas; con `exclude_fp2=true` (análisis de sensibilidad) serían 435 aristas (30 ch).

**Estado Fase 7/8 M05:** ✅ wPLI y surrogates verificados en log, matrices, `significant_connections.csv` y `surrogate_summary.json` — ⚠️ figuras `surrogate_null_*.png` y mensaje de log extendido sobre `method` requieren verificar alineación código ↔ ejecución 2026-07-09.

### Fase 8/8 — Guardado global e índices (verificado M05)

El paso **8/8** consolida tablas y figuras en el árbol BIDS, copia la configuración aplicada, actualiza índices globales del dataset y cierra el log. **No calcula** nada nuevo: persiste resultados de pasos 2–7 y genera las figuras finales que faltaban.

```
resultados en memoria (rec, spectra, conn, surr_results…)
        │
        ├─ _save_all_results()  →  escritura directa al árbol BIDS
        │        results/subjects/sub-M05/ses-T2/eyesclosed/  (sin sufijo)
        │
        ├─ _save_config_snapshot()  →  config_snapshot.toml
        ├─ _update_subjects_index() →  results/subjects_index.csv
        └─ _update_qc_decision_table() →  results/qc/qc_decision_table.csv
        │
        ▼
pipeline_log.txt cerrado · resumen en consola
```

> **Guardado incremental:** ICA (paso 4), segmentación (paso 5) y surrogates (paso 7b) ya escriben en `export_dir` **antes** del 8/8. El paso final añade tablas espectrales/conectividad consolidadas, figuras PNG y metadatos de trazabilidad.

**Rutinas y módulos:**

| Rutina | Archivo | Función |
|--------|---------|---------|
| `_save_all_results` | `src/SingleSubjectPipeline.jl` | Tablas + figuras; escritura **directa** al árbol BIDS (desde 2026-07-21; antes copiaba desde un árbol heredado) |
| `_save_config_snapshot` | `src/SingleSubjectPipeline.jl` L1719 | Copia TOML usado |
| `_update_subjects_index` | `src/SingleSubjectPipeline.jl` L1728 | Catálogo `subjects_index.csv` |
| `_update_qc_decision_table` | `src/SingleSubjectPipeline.jl` L1806 | Decisión include/exclude por grabación |

**Ruta de salida (M05):**

`results/subjects/sub-M05/ses-T2/eyesclosed/` — árbol **BIDS único**, sin sufijos.

> **Nota histórica.** La ejecución del 2026-07-09 se hizo con el código previo, que escribía además un árbol heredado `results/M05/T2/` con sufijos `_EC` y luego copiaba al BIDS. Ese doble árbol se eliminó el 2026-07-21 (ver [§7](#7-salidas-del-pipeline)); las corridas actuales producen solo el árbol BIDS. Por eso la traza de log de abajo menciona nombres con sufijo (`qc_channels_EC.csv`) que hoy serían `qc_summary.csv`.

**Comprobación M05 — ejecución 2026-07-09:**

| Comprobación | Resultado |
|--------------|-----------|
| Duración paso 8/8 | **~5 s** (10:58:27 → 10:58:32) |
| Duración pipeline total | **521.7 s** (~8.7 min) |
| Ficheros en export BIDS (raíz) | **98** (69 CSV · 21 PNG · 6 JSON · 1 TOML · 1 TXT) |
| Figuras en `figures/` | **67 PNG** (ICA, QC, segmentación, wPLI, surrogates…) |
| Decisión QC global | **`include`** |
| `subjects_index.csv` | `M05,T2,eyesclosed,31 ch,99 válidas,1 rechazada` |
| `overview.csv` | `qc_schema_version=2`, Fp2 en `bad_channels`, sin `amplitude_warning` |

**Traza en `pipeline_log.txt` (paso 8/8):**

```
[8/8] Guardando resultados
  Guardado: qc_channels_EC.csv
  Guardado: psd_by_channel_EC.csv
  Guardado: band_power_EC.csv
  Guardado: spectral_summary.json
  Guardado: regional_psd.csv
  Guardado: spectral_indices.csv
  Guardado: matrices y edges wPLI por banda
  Guardado: connectivity_summary.json + network_metrics.csv
  Guardado: signal_preview_EC.png
  Guardado: psd_all_channels_EC.png
  Guardados: heatmaps wPLI por banda
  Guardado: band_power_summary_EC.png
  Guardado: band_topomap_grid_EC.png
  QC decision: include

✓ Pipeline completado en 521.7 s
```

**Salidas del paso 8/8 (generadas o copiadas en M05):**

| Categoría | Ficheros principales |
|-----------|---------------------|
| QC / overview | `qc_summary.csv`, `channel_statistics.csv`, `overview.csv` |
| Espectral | `psd_by_channel.csv`, `band_power_summary.csv`, `spectral_summary.json`, `regional_psd.csv`, `spectral_indices.csv` |
| Conectividad | `wpli_{band}.csv` (×7), `connectivity_edges.csv`, `connectivity_summary.json`, `network_metrics.csv` |
| Figuras finales | `filtered_signal_preview.png`, `psd_all_channels.png`, `wpli_{band}.png` (×7), `band_power_summary.png` |
| Trazabilidad | `config_snapshot.toml`, `pipeline_log.txt` |

**Índices globales (fuera de `export_dir`):**

| Fichero | Contenido M05 |
|---------|---------------|
| `results/subjects_index.csv` | Fila actualizada: 31 ch, 99 épocas válidas, 1 rechazada, `processed_at=2026-07-09T10:58:30` |
| `results/qc/qc_decision_table.csv` | `M05,T2,EC` → `include`, 99.0 % válidos, Fp2 en QC original pero `fp2_removed=false`, 31 ch análisis |

**Criterio de decisión QC (M05 → `include`):**

| Condición | M05 | Efecto |
|-----------|-----|--------|
| `n_epochs_valid ≥ min_epochs` | 99 ≥ 10 | ✅ |
| `amplitude_warning` | false (σ̄=14.4 µV) | ✅ |
| `bad_ch_non_fp2` | vacío (solo Fp2 malo en QC) | ✅ |
| `fp2_removed` | false (`exclude_fp2=false`) | Fp2 permanece en análisis |

**Inventario completo por fase (export BIDS canónico):**

| Fase | Artefactos clave ya en disco |
|------|------------------------------|
| 2 QC | `qc_psd_raw_mean.csv`, `qc_band_power_raw.csv`, figuras `qc_*` en `figures/` |
| 4 ICA | `ica_*.csv/json`, `raw_signal.csv`, topomaps y butterfly en `figures/` |
| 5 Seg | `segmentation_summary.json`, `segments_table.csv`, figuras épocas |
| 6–7 | PSD, wPLI, surrogates (tablas en raíz + PNG en `figures/`) |
| 8 | Consolidación, `config_snapshot.toml`, índices globales |

**Hallazgos a tener en cuenta:**

1. **Salida única (resuelto 2026-07-21):** el doble árbol se eliminó; `_save_all_results` escribe directo al árbol BIDS `results/subjects/…`, sin copias ni riesgo de desincronización.
2. **`band_topomap_grid_EC.png`:** aparece en el log del 8/8 pero no en el `_save_all_results` actual — misma discrepancia código ↔ ejecución 2026-07-09.
3. **Caché ICA** (`results/subjects/sub-M05/ses-T2/eyesclosed/cache/ica_result.jls`): dentro del árbol BIDS del sujeto; interna, no se cita en el informe.
4. **`channel_statistics_compare.csv` / `_filtered.csv`:** timestamps 11:07 (post-pipeline); no forman parte del paso 8/8 estándar — posible análisis manual posterior.
5. **Config original ausente:** el TOML `_scratch_surrogates_verification.toml` ya no está en el repo; `config_snapshot.toml` es la única referencia reproducible de la corrida 2026-07-09.
6. **Desfase código ↔ salidas extendidas:** figuras QC (paso 2), segmentación (paso 5), `surrogate_null_*.png` y topomapas espectrales requieren reintegración para reproducir el inventario completo en un relanzamiento.

**Estado Fase 8/8 M05:** ✅ guardado global, `config_snapshot.toml`, `subjects_index.csv` y `qc_decision_table.csv` verificados — decisión **`include`**.

#### Cierre revisión pipeline M05 (fases 0–8)

| Fase | Estado | Nota principal |
|------|--------|----------------|
| 0 Preparación | ✅ | Dataset auditado; M05 T2 EC+EO válidas |
| 1/8 Carga | ✅ | 31 ch × 100.36 s vía TSV BIDS |
| 2/8 QC | ✅⚠️ | Fp2 bad; figuras QC en disco, código simplificado |
| 3/8 Filtrado | ✅ | Perfil `eeg_julia`, 4 filtros |
| 4/8 ICA | ✅⚠️ | 5 IC rechazados; auto-rechazo no en código actual |
| 5/8 Segmentación | ✅⚠️ | 99/100 épocas; figuras épocas no en código actual |
| 6/8 Espectral | ✅⚠️ | PSD OK; `band_topomap_grid` ausente en código |
| 7/8 wPLI+SUR | ✅⚠️ | 447 conexiones sig.; surrogates ~92 % del tiempo |
| 8/8 Guardado | ✅⚠️ | Índices OK; doble ruta + figuras extendidas pendientes |

**Próximos pasos sugeridos:** (1) unificar rutas de salida en `SingleSubjectPipeline.jl`; (2) reintegrar módulos de figuras extendidas; (3) alinear config repo con snapshot M05 o documentar TOML de verificación; (4) ejecutar batch completo con `[surrogates] enabled` evaluado conscientemente.

