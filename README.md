# NeuroMIND

**Framework de análisis EEG para conectividad funcional en Esclerosis Múltiple**

NeuroMIND toma señales EEG de reposo en formato BrainVision, las procesa de principio a fin y genera matrices de conectividad wPLI / dwPLI con inferencia estadística opcional. El análisis cubre el dataset MINDEM-IMIBIC (41 pacientes EM + 37 controles sanos, ~206 grabaciones válidas) y produce resultados listos para comparar grupos y sesiones.

---

## Tabla de contenidos

### Parte I — Guía de uso

1. [¿Qué hace el pipeline?](#1-qué-hace-el-pipeline)
2. [Requisitos previos](#2-requisitos-previos)
3. [Primeros pasos: preparar el entorno](#3-primeros-pasos-preparar-el-entorno)
4. [Flujo de trabajo completo](#4-flujo-de-trabajo-completo)
   - [4.0 Scripts de entrada](#40-scripts-de-entrada--estado-actual)
   - [4.1 Pipeline en lote](#41-ejecutar-el-pipeline-en-lote)
   - [4.2 Tabla de QC](#42-revisar-la-tabla-de-qc)
   - [4.3 Análisis transversal](#43-análisis-transversal)
   - [4.4 Análisis longitudinal](#44-análisis-longitudinal)
   - [4.5 Dashboard](#45-lanzar-el-dashboard)
   - [4.6 Informe PDF](#46-compilar-el-informe-pdf)
5. [Los 8 pasos del pipeline](#5-los-8-pasos-del-pipeline-explicados)
6. [Política de calidad (QC)](#6-política-de-calidad-de-señal-qc)
7. [Interpretar los resultados](#7-interpretar-los-resultados)
8. [Configuración avanzada](#8-configuración-avanzada)
9. [Estructura del proyecto](#9-estructura-del-proyecto)
10. [Referencia de dependencias](#10-referencia-de-dependencias)

### Parte II — Verificación detallada (caso M05)

- [Caso de referencia: sub-M05](#anexo-caso-de-referencia-sub-m05)
  - [Directorios de salida](#directorios-de-salida)
  - [Fases 0–8 verificadas](#fase-0--preparación-del-dataset-verificado-m05)
  - [Cierre y pendientes](#cierre-revisión-pipeline-m05-fases-08)

> **Cómo leer este documento:** empieza por la **Parte I** si es tu primer contacto con NeuroMIND. La **Parte II** documenta la verificación rutina a rutina del sujeto M05 (código ↔ salidas en disco) y sirve como plantilla para auditar otras ejecuciones.

---

## 1. ¿Qué hace el pipeline?

```
Fase 0 (previa)     audit_full_dataset.jl  →  build_bids_full.jl
        │
        ▼
Señal EEG cruda (.vhdr + .eeg, o TSV BIDS)
        │
        ▼
[1/8] Carga y validación BIDS
[2/8] Control de calidad (QC) de canales
[3/8] Filtrado (highpass + lowpass + notch + bandreject)
[4/8] ICA — eliminación de artefactos (señal continua, antes de segmentar)
[5/8] Segmentación + baseline + rechazo de artefactos ±70 µV
[6/8] Espectro de potencia (PSD) por banda
[7/8] Conectividad wPLI/dwPLI  [Hilbert | FourierCSD | Multitaper]
        └─ opcional: surrogates + FDR (si [surrogates] enabled = true)
[8/8] Guardado de tablas, figuras, índices QC y config snapshot
        │
        ▼
Tablas CSV  ·  Figuras PNG  ·  Dashboard  ·  Análisis grupal  ·  Informe PDF
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

| Paso | Acción | Sección |
|------|--------|---------|
| A | Inventariar dataset | [4.0](#40-scripts-de-entrada--estado-actual) · detalle en [Anexo M05, Fase 0](#fase-0--preparación-del-dataset-verificado-m05) |
| B | Generar metadata BIDS | idem |
| C | Probar un sujeto (`run_single_subject.jl`) | [Anexo M05](#anexo-caso-de-referencia-sub-m05) (verificación completa) |
| D | Pipeline en lote | [4.1](#41-ejecutar-el-pipeline-en-lote) |
| E | Revisar QC global | [4.2](#42-revisar-la-tabla-de-qc) |
| F | Análisis grupal / longitudinal | [4.3](#43-análisis-transversal) · [4.4](#44-análisis-longitudinal) |
| G | Explorar resultados | [4.5](#45-lanzar-el-dashboard) · [4.6](#46-compilar-el-informe-pdf) |

> **Recomendación:** antes del batch completo, ejecuta y revisa el [caso M05](#anexo-caso-de-referencia-sub-m05) con `run_single_subject.jl`. La Parte II del README documenta qué esperar en cada paso del pipeline.

### 4.0 Scripts de entrada — estado actual

La carpeta `scripts/` contiene **8 lanzadores**. Solo **uno está obsoleto** para el trabajo con el dataset MINDEM-IMIBIC; el resto forma la cadena operativa actual.

> **Regla práctica:** para M05 y el dataset real, usar siempre la ruta
> `audit → build_bids → run_single_subject` (o `run_batch_pipeline` en lote).
> **No usar** `run_pipeline.jl`.

| Script | Estado | Fase | Orquestador interno | Config | Salida principal |
|--------|--------|------|---------------------|--------|------------------|
| `audit_full_dataset.jl` | ✅ **Activo** | A — inventario | script autónomo | — | `data/full_data/inventory.csv`, `data/bids/participants.tsv`, `groups.csv`, `longitudinal_pairs.csv` |
| `build_bids_full.jl` | ✅ **Activo** | B — metadata BIDS | script autónomo | — | `data/bids/raw/*_eeg_metadata.json`, `electrodes.tsv`, `dataset_description.json` |
| `run_single_subject.jl` | ✅ **Activo** | pipeline individual | `SingleSubjectPipeline.jl` (**8 pasos**) | `config/pipeline.toml` | `results/subjects/sub-{ID}/ses-{SES}/{task}/` |
| `run_batch_pipeline.jl` | ✅ **Activo** | C — lote | idem, en bucle sobre `inventory.csv` | `config/pipeline.toml` | misma ruta BIDS + `logs/batch_run_*.csv` |
| `launch_dashboard.jl` | ✅ **Activo** | visualización | `webapp/App.jl` (Genie) | `config/pipeline.toml` | servidor en `http://localhost:8080` |
| `run_transversal_analysis.jl` | ✅ **Activo** | post-hoc grupal | script autónomo | `config/pipeline.toml` | `results/group/transversal/` |
| `run_longitudinal_analysis.jl` | ✅ **Activo** | post-hoc longitudinal | script autónomo | `config/pipeline.toml` | `results/group/longitudinal/` |
| `run_pipeline.jl` | ⚠️ **Obsoleto — archivado** | legacy | `Pipeline.jl` (**7 pasos**, sin surrogates ni export BIDS) | `legacy/config/pipeline.toml` + `legacy/config/subjects.toml` | caché serializada en `results/{ID}/{SES}/` |

**Por qué `run_pipeline.jl` está obsoleto** (archivado en `legacy/` el 2026-07-21 junto con su config y su orquestador):

1. Llama a `run_pipeline!` (`legacy/src/Pipeline.jl`), un orquestador **anterior** al actual `SingleSubjectPipeline.jl`.
2. Solo implementa **7 pasos** (sin paso 8 de surrogates/FDR ni el guardado completo de tablas/figuras en `results/subjects/`).
3. Lee sujetos desde `legacy/config/subjects.toml`, que contiene **entradas sintéticas de plantilla** (`SYN_MS_001`, etc.), no el registro real MINDEM-IMIBIC.
4. Carga datos con `load_eeg_bids` (TSV obligatorio); no integra `BrainVisionLoader` ni el flujo de metadata ligera de Fase B.
5. El caso M05 verificado (2026-07-09) se ejecutó con `run_single_subject.jl`, no con este script.

**Cadena recomendada (orden de ejecución):**

```bash
# 1–2. Preparación (una vez por dataset)
julia --project=. scripts/audit_full_dataset.jl
julia --project=. scripts/build_bids_full.jl

# 3. Prueba individual (recomendado antes del batch)
julia --project=. scripts/run_single_subject.jl
# → revisar salidas en results/subjects/… (ver Anexo M05)

# 4. Batch completo
julia --project=. scripts/run_batch_pipeline.jl
```

```
audit_full_dataset.jl  →  build_bids_full.jl  →  run_single_subject.jl  (prueba M05)
                                              ↘  run_batch_pipeline.jl  (dataset completo)
                                                    ↓
                              run_transversal_analysis.jl  /  run_longitudinal_analysis.jl
                                                    ↓
                              launch_dashboard.jl  (inspección interactiva)
```

**Dos orquestadores en `src/` (no confundir):**

| Módulo | Usado por | Estado |
|--------|-----------|--------|
| `SingleSubjectPipeline.jl` | `run_single_subject.jl`, `run_batch_pipeline.jl` | ✅ **Canónico** — 8 pasos, BrainVision, surrogates, export BIDS |
| `legacy/src/Pipeline.jl` | `legacy/scripts/run_pipeline.jl` | ⚠️ **Archivado (2026-07-21)** — fuera de `src/`, incluido desde `NeuroMIND.jl` solo por compatibilidad |

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

### [7/8] Conectividad wPLI / dwPLI (+ surrogates opcional)

Calcula el *weighted Phase Lag Index* entre todos los pares de canales en cada banda de frecuencia. Con 30 canales (montaje estándar), hay **435 pares de electrodos** por banda.

**Tres métodos de estimación configurables:**

| Método | Descripción | Cuándo usarlo |
|--------|-------------|---------------|
| `hilbert` | Filtrado Butterworth + señal analítica Hilbert | Análisis exploratorio rápido |
| `fourier_csd` | Espectro cruzado FFT con ventana Hanning/Hamming | **Recomendado para publicación** |
| `multitaper` | DPSS multitaper (equivalente a MNE `spectral_connectivity`) | Cuando se requiere bajo sesgo espectral |

```toml
[connectivity]
wpli_method = "fourier_csd"  # o "hilbert", "multitaper"
use_dwpli   = true           # true → dwPLI no sesgado (recomendado para grupos)
```

> La variante **dwPLI** (debiased wPLI) es matemáticamente insesgada respecto al número de épocas y produce rango [-1, 1]. Se recomienda para comparaciones grupales.

**Inferencia por surrogates** (sub-paso opcional, activar en config):

Si `[surrogates] enabled = true`, tras calcular wPLI se generan permutaciones mediante **desplazamiento circular independiente por canal y época**. Esta técnica destruye la sincronía de fase inter-canal pero preserva el espectro de cada canal. Los p-valores usan corrección Monte Carlo (+1) y FDR Benjamini-Hochberg al 5 %.

```toml
[surrogates]
enabled      = false   # true para inferencia; false en batch exploratorio
n_surrogates = 200     # ≥ 200 para publicación
alpha        = 0.05
seed         = 42
```

> Con `enabled = true`, los surrogates pueden consumir >90 % del tiempo de ejecución (ver [Anexo M05, Fase 7](#fase-78--conectividad-wpli-y-surrogates-verificado-m05)).

### [8/8] Guardado de resultados

Consolida tablas CSV, figuras PNG, `config_snapshot.toml`, `pipeline_log.txt` y actualiza índices globales (`subjects_index.csv`, `qc_decision_table.csv`). Escribe en dos rutas: export BIDS (`results/subjects/…`) y dashboard (`results/{ID}/{SES}/`). Ver [sección 7](#7-interpretar-los-resultados) y [Directorios de salida](#directorios-de-salida) en el anexo M05.

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

**Ruta canónica (informe, análisis grupal, BIDS):**

```
results/subjects/sub-{ID}/ses-{SES}/{task}/
├── overview.csv              ← resumen general (canales, epochs, duración)
├── qc_summary.csv            ← QC de la grabación
├── channel_statistics.csv    ← estadísticas por canal
├── band_power_summary.csv    ← potencia por banda
├── wpli_{BANDA}.csv          ← matriz wPLI por banda
├── wpli_pvalues_{BANDA}.csv  ← p-valor por par (si surrogates ON)
├── wpli_qvalues_{BANDA}.csv  ← q-valor FDR
├── wpli_significant_{BANDA}.csv
├── significant_connections.csv
├── surrogate_summary.json
├── ica_summary.json
├── config_snapshot.toml      ← TOML usado en la ejecución
├── figures/                  ← figuras PNG
└── pipeline_log.txt
```

**Ruta dashboard** (paneles web, sufijo `_EC` / `_EO`):

```
results/{ID}/{SES}/
├── tables/    ← mismas tablas con sufijo de condición
├── figures/   ← figuras con sufijo _EC.png
└── cache/     ← caché ICA (no citar en informe)
```

> El pipeline escribe primero en `results/{ID}/{SES}/` y copia a `results/subjects/…`. Para el informe, usar siempre la ruta BIDS. Detalle y reglas en [Directorios de salida](#directorios-de-salida).

### ¿Cómo leer los resultados de conectividad?

El archivo `significant_connections.csv` contiene las conexiones estadísticamente significativas:

```
ch_a, ch_b, band, wpli_obs, p_value, q_value, z_score
Fz,   P4,   ALPHA, 0.388,  0.005,   0.0115, 4.179
...
```

- `wpli_obs` > 0 → sincronización en fase (relación leading-lagging)
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
│   └── pipeline.toml           ← ⭐ CONFIGURACIÓN ÚNICA del proyecto
│                                  (unificada 2026-07-21; la leen los 5 scripts activos)
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
│   ├── README.md               ← Estructura y qué comando escribe dónde
│   ├── subjects/               ← Nivel sujeto: sub-{ID}/ses-{S}/{task}/
│   ├── transversal/            ← Nivel grupo: EM vs controles ({EC|EO})
│   ├── longitudinal/           ← Nivel grupo: T1 vs T2 ({EC|EO})
│   ├── qc/
│   │   └── qc_decision_table.csv   ← Estado QC de cada grabación
│   └── logs/                   ← batch_run_{timestamp}.csv
│
├── deprecated/                 ← Resultados archivados (NO en Git)
│   └── results/
│       └── 2026-05-26_pre-unificacion/  ← 204 grabaciones + grupo de mayo,
│                                           config incompatible (ver su README)
│
├── src/                        ← Código fuente Julia
│   ├── NeuroMIND.jl            ← Entry point del módulo
│   ├── types.jl                ← Tipos de datos (EEGRecording, EpochSet, etc.)
│   ├── SingleSubjectPipeline.jl ← Pipeline de 8 pasos (canónico)
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
│   ├── audit_full_dataset.jl        ← Fase A: inventario del dataset
│   ├── build_bids_full.jl           ← Fase B: metadata BIDS ligera
│   ├── run_single_subject.jl        ← Pipeline individual (8 pasos) ← USAR
│   ├── run_batch_pipeline.jl        ← Pipeline en lote (Fase C)
│   ├── run_transversal_analysis.jl  ← Comparación grupal EM vs controles
│   ├── run_longitudinal_analysis.jl ← Comparación temporal T1 → T2
│   └── launch_dashboard.jl          ← Dashboard interactivo
│
├── legacy/                     ← ⚠️ Archivado 2026-07-21 — no usar con MINDEM
│   ├── config/
│   │   ├── single_subject.toml      ← Config individual previa (generó el caso M05)
│   │   ├── batch_pipeline.toml      ← Config de lote previa (generó las 205 grabaciones)
│   │   ├── pipeline.toml            ← Config del orquestador legacy (≠ config/pipeline.toml)
│   │   └── subjects.toml            ← Registro de sujetos (entradas sintéticas de plantilla)
│   ├── scripts/
│   │   └── run_pipeline.jl          ← Lanzador legacy (7 pasos, sin surrogates ni BIDS)
│   └── src/
│       └── Pipeline.jl              ← Orquestador legacy de 7 pasos
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

Para la verificación detallada paso a paso (código ↔ salidas en disco), continúa con el [Anexo M05](#anexo-caso-de-referencia-sub-m05).

---

## Anexo: Caso de referencia sub-M05

> Verificación rutina a rutina del pipeline individual: código, tablas y figuras para `sub-M05 / ses-T2 / eyesclosed` (EC). Asume que has leído la [Parte I](#1-qué-hace-el-pipeline), en especial las secciones [4](#4-flujo-de-trabajo-completo) y [5](#5-los-8-pasos-del-pipeline-explicados).

Sujeto de trabajo para verificar rutina a rutina el código, las tablas y las figuras generadas por el pipeline individual.

| Campo | Valor |
|-------|-------|
| **Sujeto** | `sub-M05` |
| **Sesión** | `ses-T2` |
| **Tarea** | `eyesclosed` (ojos cerrados, condición **EC**) |
| **Último lanzamiento** | **2026-07-09** (inicio 10:49:50, duración 521.7 s) |
| **Comando** | `julia --project=. scripts/run_single_subject.jl --config <toml>` (ver `pipeline_log.txt` para el TOML exacto) |
| **Traza** | `results/subjects/sub-M05/ses-T2/eyesclosed/pipeline_log.txt` |
| **Config aplicada** | `results/subjects/sub-M05/ses-T2/eyesclosed/config_snapshot.toml` (copia fiel del TOML usado; el archivo `config/_scratch_surrogates_verification.toml` ya no está en el repo) |

> Esta fecha y esta config deben actualizarse cada vez que se relance el pipeline sobre M05.
> El `config/single_subject.toml` del repo tiene `[surrogates] enabled = false` por defecto; la ejecución del 2026-07-09 activó surrogates para regenerar el ejemplo con inferencia estadística.

### Directorios de salida

> **Pendiente de corrección:** el pipeline escribe los mismos resultados en **dos rutas distintas**. Hay que unificar esto en el código para que solo exista una carpeta canónica por sujeto/sesión/tarea.

Hay dos rutas de salida. Ambas contienen los mismos datos:

```
results/M05/T2/
  tables/     ← tablas con sufijo _EC (versión dashboard)
  figures/    ← imágenes con sufijo _EC
  cache/      ← caché ICA (uso interno, no para el informe)
  logs/

results/subjects/sub-M05/ses-T2/eyesclosed/    ← export BIDS (usar en informe)
  figures/
  (todos los CSV y JSON en la raíz)
```

**Los ficheros en `results/subjects/` son los canónicos:** mismos datos, sin sufijo `_EC`, nombres alineados con BIDS (`sub-{ID}/ses-{SES}/{task}/`).

| Uso | Ruta |
|-----|------|
| Informe, revisión científica, citas en `Report_Pre/` | `results/subjects/sub-M05/ses-T2/eyesclosed/` |
| Dashboard web (paneles 0–15) | `results/M05/T2/` (`tables/`, `figures/`) |
| Caché ICA (no versionar, no citar) | `results/M05/T2/cache/ica_result.jls` |

**Regla práctica:** al incorporar figuras o tablas al informe, copiar siempre desde `results/subjects/…` (raíz o `figures/` sin sufijo `_EC`). Las copias en `results/M05/T2/figures/*_EC.png` pueden estar desfasadas respecto a la última ejecución.

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
config/single_subject.toml  →  load_ss_config()
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
| `results/M05/T2/figures/signal_preview_EC.png` | Figura | Paso 8/8 | Copia dashboard con sufijo `_EC` |
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
| Caché | `results/M05/T2/cache/ica_result.jls` | Evita re-ejecutar FastICA si la config ICA no cambia |

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
3. **Caché ICA:** `results/M05/T2/cache/ica_result.jls` + `ica_config.hash`. Si cambian parámetros `[ica]`, se invalida y se re-ejecuta FastICA (~varios segundos).
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
2. **Fp2 incluido en montaje M05:** `exclude_fp2=false` en la config de esa ejecución (decisión 2026-07-08 tras limpieza ICA). Con `exclude_fp2=true` (default en código si no hay bloque `[montage]`) el montaje sería de 30 canales.
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

| Fichero (export BIDS) | Fichero (dashboard `results/M05/T2/figures/`) | Descripción |
|-----------------------|-----------------------------------------------|-------------|
| `psd_all_channels.png` | `psd_all_channels_EC.png` | Grid PSD 31 canales (0–50 Hz) |
| `band_power_summary.png` | `band_power_summary_EC.png` | Barplot potencia media por banda |
| `band_topomap_grid.png` | `band_topomap_grid_EC.png` | Topomapas por banda (7 paneles) |

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
6. **Rutas duplicadas:** figuras con sufijo `_EC` en `results/M05/T2/figures/`; copias sin sufijo en `results/subjects/sub-M05/ses-T2/eyesclosed/`.

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
5. **Alta proporción de significativos en α y γ:** 55 % y 41 % de pares en reposo EC es esperable en análisis within-subject con FDR permisivo; la comparación grupal requiere `use_dwpli=true` y análisis de segundo nivel.
6. **CSD desactivado:** valores wPLI no son directamente comparables con EEG_Julia (que aplica CSD antes de wPLI).
7. **Fp2 en conectividad:** al estar en el montaje, participa en las 465 aristas; con `exclude_fp2=true` (default batch) serían 435 aristas (30 ch).

**Estado Fase 7/8 M05:** ✅ wPLI y surrogates verificados en log, matrices, `significant_connections.csv` y `surrogate_summary.json` — ⚠️ figuras `surrogate_null_*.png` y mensaje de log extendido sobre `method` requieren verificar alineación código ↔ ejecución 2026-07-09.

### Fase 8/8 — Guardado global e índices (verificado M05)

El paso **8/8** consolida tablas y figuras en las dos rutas de salida, copia la configuración aplicada, actualiza índices globales del dataset y cierra el log. **No calcula** nada nuevo: persiste resultados de pasos 2–7 y genera las figuras finales que faltaban.

```
resultados en memoria (rec, spectra, conn, surr_results…)
        │
        ├─ _save_all_results()
        │     ├─ tables/  →  qc, overview, PSD, band_power, wPLI (sufijo _EC)
        │     └─ cp()     →  export BIDS sin sufijo (results/subjects/…)
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
| `_save_all_results` | `src/SingleSubjectPipeline.jl` L1314 | Tablas + figuras; `cp()` dashboard → BIDS |
| `_save_config_snapshot` | `src/SingleSubjectPipeline.jl` L1719 | Copia TOML usado |
| `_update_subjects_index` | `src/SingleSubjectPipeline.jl` L1728 | Catálogo `subjects_index.csv` |
| `_update_qc_decision_table` | `src/SingleSubjectPipeline.jl` L1806 | Decisión include/exclude por grabación |
| `load_dashboard_data` | `src/SingleSubjectPipeline.jl` L1761 | Lectura desde `results/{ID}/{SES}/tables/` |

**Doble ruta de salida (M05):**

| Ruta | Rol | Convención nombres |
|------|-----|-------------------|
| `results/M05/T2/` | Dashboard (`App.jl`, paneles 0–15) | `tables/*_EC.csv`, `figures/*_EC.png` |
| `results/subjects/sub-M05/ses-T2/eyesclosed/` | **Canónica** (informe, BIDS) | Sin sufijo `_EC`; CSV/JSON en raíz |

El patrón es **escribir en dashboard → copiar a BIDS** (`cp(...; force=true)`). Ver [Directorios de salida](#directorios-de-salida) al inicio del caso M05.

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
| 6–7 | PSD, wPLI, surrogates (tablas + PNG en raíz y `figures/`) |
| 8 | Consolidación, `config_snapshot.toml`, índices globales |

**Hallazgos a tener en cuenta:**

1. **Doble escritura pendiente de unificar:** el mismo dato vive en `results/M05/T2/` y `results/subjects/…`; riesgo de desincronización si solo se lee una ruta.
2. **`band_topomap_grid_EC.png`:** aparece en el log del 8/8 pero no en el `_save_all_results` actual — misma discrepancia código ↔ ejecución 2026-07-09.
3. **Caché ICA** (`results/M05/T2/cache/ica_result.jls`): solo dashboard; no se copia a export BIDS (correcto para informe).
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

