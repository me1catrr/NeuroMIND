# mne_brain

**Pipeline de réplica en MNE-Python para los análisis EEG del proyecto NeuroMIND/BRAIN.**

Este subproyecto reproduce, paso a paso y con herramientas nativas de Python/MNE,
el pipeline Julia de [NeuroMIND](../). Su finalidad es doble:

1. **Validar científicamente** los resultados de NeuroMIND comparándolos con una
   implementación independiente basada en MNE-Python (preprocesado, ICA, epochs,
   PSD, wPLI, surrogates).
2. **Servir de banco de pruebas** para evaluar diferencias metodológicas
   (perfil `eeg_julia` vs. perfil `default`, ICA propia vs. ICA de MNE,
   filtrado causal vs. zero-phase, etc.) sobre el mismo dataset MINDEM-IMIBIC.

---

## 1. Identidad del proyecto

| Campo | Valor |
|-------|-------|
| **Investigador** | Rafael Castro Triguero (`me1catrr@uco.es`) |
| **Repo padre** | https://github.com/me1catrr/NeuroMIND |
| **Lenguaje** | Python 3.10+ (testado con 3.14) |
| **Librerías clave** | MNE 1.12, mne-connectivity, numpy, scipy, pandas, matplotlib |
| **Dataset** | MINDEM-IMIBIC — 41 pacientes EM + 37 controles sanos, sesiones T1/T2, condiciones ojos cerrados (EC) y abiertos (EO) |
| **Formato origen** | BrainVision binario (`.vhdr`/`.eeg`/`.vmrk`), 31 canales 10-20, 500 Hz |

---

## 2. Flujo general del proyecto

El proyecto se organiza en **tres etapas** secuenciales:

```
┌────────────────────────────────┐
│  A. Construcción del BIDS      │   build_bids_full.py
│     (metadata por grabación)   │   → data/BIDS/raw/*.json
└────────────────────────────────┘
              ↓
┌────────────────────────────────┐
│  B. Pipeline por sujeto        │   run_full_pipeline.py
│     (Fases 3 → 7)              │   → results/subjects/sub-{id}/...
└────────────────────────────────┘
              ↓
┌────────────────────────────────┐
│  C. Análisis grupal y          │   (transversal / longitudinal)
│     comparación con NeuroMIND  │   → validation/
└────────────────────────────────┘
```

Cada etapa puede ejecutarse de forma independiente y todas dejan log estructurado
en `results/` para trazabilidad.

---

## 3. Pipeline EEG — orden de las fases

`mne_brain` aplica las mismas 7 fases científicas que NeuroMIND/Julia, pero
implementadas con MNE-Python.

> **Regla crítica:** ICA (Fase 4) se ejecuta sobre la señal **continua filtrada**
> y siempre **antes** de segmentar (Fase 5). La señal limpiada por ICA es la
> que alimenta el resto del pipeline.

### Fase 3 — Carga, control de calidad y filtrado

- Lee el `.vhdr` original mediante `mne.io.read_raw_brainvision`.
- Genera un informe QC por canal (`z-score > 3.0` → canal marcado como `bad`).
- Aplica la cadena de filtros del perfil `eeg_julia`:
  1. **Notch 50 Hz**, causal (red eléctrica europea)
  2. **Bandreject 99.5–100.5 Hz**, causal (artefacto del actiCHamp)
  3. **Pasa-alta 0.5 Hz**, zero-phase (`filtfilt`)
  4. **Pasa-baja 150 Hz**, zero-phase (`filtfilt`)
- Guarda la señal filtrada en `cache/filtered_{COND}.npz`.

### Fase 4 — ICA (FastICA)

- Ajusta `mne.preprocessing.ICA(method="fastica", n_components=30, random_seed=42)`.
- Calcula 7 features por componente (frontal, temporal, blink, EMG, line, kurtosis…)
  y propone candidatos a rechazar con dos criterios:
  - **Custom NeuroMIND** (`evaluate_ica_components`, umbral 1.5)
  - **MNE built-in** (`find_bads_eog`, `find_bads_muscle`)
- Aplica la reconstrucción excluyendo los componentes sugeridos.
- Genera topomapas (PNG, 10 por figura).

### Fase 5 — Epochs + baseline + rechazo de artefactos

- Segmenta la señal limpiada en **épocas de 1.0 s sin solape**.
- Aplica baseline `first_window_mean` con doble pase (perfil `eeg_julia`).
- Rechazo de artefactos (`±70 µV` por canal, perfil `eeg_julia`).
- Guarda `cache/epochs_{COND}-epo.fif` (formato nativo MNE).

### Fase 6 — Análisis espectral (PSD)

- PSD por época y canal con **Welch + ventana Hamming + zero-padding** (`nfft=512`).
- Potencia integrada por las 7 bandas: δ, θ, α, β_low, β_mid, β_high, γ.
- Salida en µV² (acorde a NeuroMIND).

### Fase 7 — Conectividad wPLI

- Conectividad funcional **weighted Phase Lag Index** par a par entre los 31 canales.
- Agregación **across-segments** (idéntica a NeuroMIND).
- Una matriz `wpli_{BANDA}_{COND}.csv` por banda + heatmap PNG.
- Resumen estadístico (media, máximo, top-edges) en `connectivity_summary_{COND}.json`.

> Las **Fases 1 y 2** son el _scaffolding_ del proyecto (paquete instalable,
> loaders BIDS, configuración) y no aparecen en el pipeline de ejecución.

---

## 4. Configuración

Toda la configuración científica vive en
[`config/pipeline_config.yaml`](config/pipeline_config.yaml). Es un espejo del
`config/single_subject.toml` de NeuroMIND/Julia para garantizar parámetros
idénticos (el antiguo `config/pipeline.toml` se archivó el 2026-07-21 en
`NeuroMIND/deprecated/code/config/`; ya no es la fuente activa).

Cambia ahí (no en el código) cosas como:

```yaml
filtering:        # cadena de filtros y perfil (eeg_julia / default)
segmentation:     # longitud y solape de épocas
artifact_rejection:  # umbrales de AR
ica:              # nº componentes, método, semilla
spectral:         # nfft, ventana, unidades
connectivity:     # método (wpli), agregación, CSD on/off
bands:            # límites de las 7 bandas EEG
```

---

## 5. Instalación

### Requisitos

- Python ≥ 3.10
- Acceso a los datos BrainVision originales en
  `../data/full_data/Pacientes MINDEM_IMIBIC_27 03 25/`
- El inventario auditado `../data/full_data/inventory.csv` (generado por
  NeuroMIND/Julia, fase A)

### Pasos

```bash
cd NeuroMIND/mne_brain
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest                                # ~150 tests unitarios
```

### Comprobación rápida del loader

```python
from mne_brain import load_config
from mne_brain.bids import load_eeg_bids

cfg = load_config()
rec = load_eeg_bids(cfg, "M05", "T2", "EC")
print(rec.data.shape, rec.meta.fs, rec.meta.channel_names[:5])
# (31, 47350) 500.0 ['Fz', 'F3', 'F7', 'FT9', 'FC5']
```

---

## 6. Cómo lanzar el pipeline

### Etapa A — Construir los metadatos BIDS (una sola vez)

Genera un JSON por grabación en `data/BIDS/raw/` leyendo los `.vhdr` originales
con MNE. **No copia datos binarios**: solo guarda la ruta `vhdr_path`, los
parámetros del header (`fs`, `n_channels`, `channel_names`) y la información
clínica (`group`, `sex`, `age`) tomada de `participants.tsv`.

```bash
# Todas las grabaciones válidas del inventario (~205 JSON en ~2 s)
python scripts/build_bids_full.py

# Restringido a controles, T1, ojos cerrados
python scripts/build_bids_full.py --groups HC --sessions T1 --conditions EC

# Previsualizar sin escribir
python scripts/build_bids_full.py --limit 5 --dry-run

# Forzar reescritura tras un cambio en la plantilla hardware/software
python scripts/build_bids_full.py --overwrite
```

**Salidas:**
- `data/BIDS/raw/sub-{ID}_ses-{T1|T2}_task-{eyesclosed|eyesopen}_run-01_eeg_metadata.json`
- `data/BIDS/{participants.tsv, groups.csv, longitudinal_pairs.csv, dataset_description.json}` (copiados desde NeuroMIND)
- `results/bids_build_log.csv` — una fila por grabación (ok / error)

### Etapa B — Pipeline EEG

Un único script (`run_full_pipeline.py`) cubre **ambos modos**:

#### B.1 — Un único sujeto (modo verboso)

```bash
# Ambas condiciones (EC + EO)
python scripts/run_full_pipeline.py --subject M07 --session T1

# Solo una condición
python scripts/run_full_pipeline.py --subject M20 --session T2 --condition EC

# Reanudar desde una fase concreta (reaprovecha cache)
python scripts/run_full_pipeline.py --subject M30 --session T1 --from-phase 5
```

#### B.2 — Todos los sujetos seguidos (modo batch)

```bash
# Procesar todo el dataset (saltando lo ya terminado)
python scripts/run_full_pipeline.py --all --skip-done

# En background, con log de stdout y posibilidad de cerrar la terminal
mkdir -p results/logs
nohup python -u scripts/run_full_pipeline.py --all --skip-done \
    > results/logs/batch_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Restringido a un subgrupo
python scripts/run_full_pipeline.py --all --groups HC --sessions T1 --conditions EC

# Subset arbitrario de sujetos
python scripts/run_full_pipeline.py --subjects M07 M20 MC10

# Solo recalcular wPLI (Fase 7) con epochs ya cacheados
python scripts/run_full_pipeline.py --all --from-phase 7 --skip-done

# Reanudar a partir de un sujeto concreto
python scripts/run_full_pipeline.py --all --start-from M20 --skip-done

# Previsualizar sin ejecutar
python scripts/run_full_pipeline.py --all --limit 5 --dry-run
```

#### El script detecta el modo automáticamente

| Si pasas… | Modo |
|-----------|------|
| `--subject X --session Y` (y nada más de batch) | Single subject |
| `--all` | Batch sobre todo el dataset |
| `--subjects M07 M20 MC10` (≥2 ids) | Batch sobre el subset |
| Cualquier filtro (`--groups`, `--sessions`, `--conditions`, `--start-from`, `--limit`, `--dry-run`, `--skip-done`) | Batch |

#### Argumentos del pipeline

| Argumento | Tipo | Descripción |
|-----------|------|-------------|
| `--subject ID` | str | (single) BIDS id del sujeto (p.ej. `M07`, `MC10`) |
| `--session T1\|T2` | str | (single) Sesión |
| `--condition EC\|EO` | str | (single) Una condición o ambas si se omite |
| `--all` | flag | (batch) Procesar todo `data/BIDS/raw/` |
| `--subjects ID ID ...` | lista | (batch) Restringir a una lista de sujetos |
| `--groups MS HC` | lista | (batch) Restringir a uno o más grupos |
| `--sessions T1 T2` | lista | (batch) Restringir a una o más sesiones |
| `--conditions EC EO` | lista | (batch) Restringir a una o más condiciones |
| `--start-from ID` | str | (batch) Solo sujetos con id `>= ID` (lexicográfico) |
| `--limit N` | int | (batch) Solo las primeras N parejas (sujeto, sesión) |
| `--from-phase 3..7` | int | (común) Fase de inicio (reutiliza cache previa) |
| `--skip-done` | flag | (batch) Saltar grabaciones con `connectivity_summary_*.json` existente |
| `--dry-run` | flag | (batch) Listar el trabajo planeado sin ejecutarlo |
| `--log-path PATH` | str | (batch) Ruta del CSV de log (default `results/batch_pipeline_log.csv`) |

#### Comportamiento del modo batch

- **Una fila por (sujeto, sesión, condición)** en `results/batch_pipeline_log.csv` con
  `timestamp`, `subject`, `session`, `task`, `condition`, `group`, `from_phase`,
  `status` (`ok`/`skip`/`error`), `elapsed_s`, `error_type`, `error_message`.
- **Un fallo no detiene el batch:** el siguiente sujeto continúa.
- **Skip-done inteligente:** comprueba si existe `tables/connectivity_summary_{COND}.json`
  (artefacto de Fase 7); útil para reanudar tras interrupción.
- **Tiempo típico:** ~3-10 s por condición → 200 grabaciones en ~15-30 minutos
  (mucho más rápido que NeuroMIND/Julia gracias al backend C de MNE).

### Comandos útiles durante el batch

```bash
# Seguir el log de stdout en vivo
tail -f results/logs/batch_*.log

# Resumen rápido por estado
awk -F, 'NR>1{print $8}' results/batch_pipeline_log.csv | sort | uniq -c

# Ver solo los errores
awk -F, 'NR>1 && $8=="error"' results/batch_pipeline_log.csv
```

---

## 7. Estructura de resultados

Cada `(sujeto, sesión, condición)` genera una carpeta independiente:

```text
results/subjects/sub-{ID}/ses-{T1|T2}/{eyesclosed|eyesopen}/
├── cache/                              ← intermedios reusables entre fases
│   ├── filtered_{COND}.npz             ← Fase 3 (señal filtrada cruda)
│   ├── ica_{COND}-ica.fif              ← Fase 4 (objeto ICA serializado)
│   ├── ica_activations_{COND}.npz      ← Fase 4 (activaciones IC)
│   ├── cleaned_{COND}.npz              ← Fase 4 (señal limpiada por ICA)
│   └── epochs_{COND}-epo.fif           ← Fase 5 (epochs MNE)
├── tables/                             ← métricas tabulares
│   ├── qc_channels_{COND}.csv          ← Fase 3
│   ├── filter_chain_{COND}.json
│   ├── signal_preview_{COND}.json
│   ├── ica_summary_{COND}.json         ← Fase 4
│   ├── ica_component_features_{COND}.csv
│   ├── ica_component_maps_{COND}.csv
│   ├── ica_mixing_matrix_{COND}.csv
│   ├── ica_unmixing_matrix_{COND}.csv
│   ├── epoch_summary_{COND}.json       ← Fase 5
│   ├── psd_by_channel_{COND}.csv       ← Fase 6
│   ├── band_power_summary_{COND}.csv
│   ├── spectral_params_{COND}.json
│   ├── wpli_{BANDA}_{COND}.csv         ← Fase 7 (una por banda × 7)
│   ├── connectivity_edges_{COND}.csv
│   └── connectivity_summary_{COND}.json
└── figures/                            ← visualizaciones PNG
    ├── ica/
    │   └── ica_topomaps_{COND}_*.png   ← topomapas (10 IC por figura)
    ├── psd_spectrum_{COND}.png         ← espectro promedio
    ├── wpli_all_bands_{COND}.png       ← parrilla 7-bandas
    └── wpli_heatmap_{BANDA}_{COND}.png ← heatmap par-canal por banda (×7)
```

### Logs globales (en `results/`)

| Fichero | Contenido |
|---------|-----------|
| `bids_build_log.csv` | Una fila por metadata JSON construido (Etapa A) |
| `batch_pipeline_log.csv` | Una fila por grabación procesada (Etapa B) |
| `logs/batch_*.log` | Stdout completo de cada ejecución batch |

---

## 8. Estructura del repositorio

```text
mne_brain/
├── config/
│   └── pipeline_config.yaml        ← parámetros científicos (espejo de NeuroMIND)
├── data/
│   └── BIDS/                       ← metadatos BIDS ligeros (no datos binarios)
│       ├── raw/                    ← *_eeg_metadata.json por grabación
│       ├── electrodes/             ← posiciones 10-20 (TSV)
│       ├── participants.tsv
│       ├── groups.csv
│       ├── longitudinal_pairs.csv
│       └── dataset_description.json
├── scripts/
│   ├── build_bids_full.py          ← Etapa A: construye metadatos BIDS
│   ├── run_full_pipeline.py        ← Etapa B: pipeline single + batch unificado
│   ├── run_phase3_m05.py …         ← scripts auxiliares por fase (debug)
│   └── run_phase8_validation.py    ← comparación con NeuroMIND/Julia
├── src/mne_brain/
│   ├── bids/                       ← loaders BIDS y BrainVision
│   ├── common/                     ← tipos y config
│   ├── preprocessing/              ← filtrado y QC (Fase 3)
│   ├── ica/                        ← ICA + features + suggestions (Fase 4)
│   ├── processing/                 ← epochs + baseline + AR (Fase 5)
│   ├── spectral/                   ← PSD + band power (Fase 6)
│   ├── connectivity/               ← wPLI (Fase 7)
│   └── surrogate/                  ← surrogates y FDR (Fase 8, in progress)
├── tests/
│   └── unit/                       ← tests unitarios por módulo (pytest)
├── validation/
│   └── reports/                    ← informes diff Python vs Julia (PNG + JSON)
└── results/                        ← resultados por sujeto + logs (NO va a git)
```

---

## 9. Flujo recomendado para reproducir todo el dataset

```bash
# 0. Activar entorno
cd NeuroMIND/mne_brain
source .venv/bin/activate

# 1. (Una sola vez) Construir metadatos BIDS para los ~205 recordings
python scripts/build_bids_full.py

# 2. Lanzar pipeline completo en background con skip-done por seguridad
mkdir -p results/logs
nohup python -u scripts/run_full_pipeline.py --all --skip-done \
    > results/logs/batch_$(date +%Y%m%d_%H%M%S).log 2>&1 &
echo "PID=$!"

# 3. Monitorizar
tail -f results/logs/batch_*.log

# 4. (Opcional) Cuando termine, revisar sujetos con errores y reintentar
awk -F, 'NR>1 && $8=="error"' results/batch_pipeline_log.csv \
    | cut -d, -f2 | sort -u > /tmp/failed.txt
python scripts/run_full_pipeline.py --subjects $(cat /tmp/failed.txt)
```

---

## 10. Comparación con NeuroMIND/Julia

`mne_brain` está diseñado para producir resultados **numéricamente comparables**
con NeuroMIND/Julia. Los outputs `.csv` y `.json` tienen los mismos nombres y
columnas que sus equivalentes en Julia (`NeuroMIND/results/subjects/...`), lo
que permite hacer `diff` directos.

Diferencias intencionadas conocidas:

| Aspecto | NeuroMIND/Julia | mne_brain |
|---------|-----------------|-----------|
| Reader BrainVision | Lector binario propio | `mne.io.read_raw_brainvision` |
| ICA | FastICA puro Julia (`ICACore.jl`) | `mne.preprocessing.ICA` (sklearn FastICA) |
| Filtros | Butterworth `DSP.jl` | Butterworth `mne.filter` |
| Hilbert para wPLI | FFTW propio | `scipy.signal.hilbert` |
| Surrogates | A implementar | A implementar (módulo `surrogate/`) |

La carpeta `validation/` contiene scripts de diff par-a-par para detectar
regresiones.

---

## 11. Información complementaria

- Reglas Git, convenciones y contexto científico ampliado: ver
  [`../CLAUDE.md`](../CLAUDE.md) y [`../AGENTS.md`](../AGENTS.md) en la raíz
  de NeuroMIND.
- Documentación del informe hardware BrainVision (que alimenta la plantilla
  de `build_bids_full.py`):
  [`data/BIDS/BrainVision_hardware_software_metadata.md`](data/BIDS/BrainVision_hardware_software_metadata.md).
- Reglas de exclusión del dataset (artefactos Oddball, duplicados, etc.):
  ver columna `excluded` en `../data/full_data/inventory.csv`.
