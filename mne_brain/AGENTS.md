# AGENTS.md — mne_brain

> **Documento de onboarding para agentes IA** (Claude Code, Codex CLI, Cursor, etc.)
> que vayan a trabajar en este subproyecto. Lee primero esto y luego el `README.md`.
>
> Última actualización: **2026-05-26**.

---

## 1. ¿Qué es este proyecto?

`mne_brain` es la **réplica en MNE-Python** del pipeline EEG de
[NeuroMIND/Julia](../). Sirve para validar científicamente los resultados de
NeuroMIND comparándolos con una implementación independiente.

- **Investigador:** Rafael Castro Triguero (`me1catrr@uco.es`)
- **Repo:** https://github.com/me1catrr/NeuroMIND (subcarpeta `mne_brain/`)
- **Lenguaje:** Python 3.10+ (probado con 3.14)
- **Librerías clave:** MNE 1.12, mne-connectivity, numpy, scipy, pandas, matplotlib
- **Dataset:** MINDEM-IMIBIC — 41 EM + 37 controles, BrainVision binario, 31 ch @ 500 Hz

---

## 2. Reglas estrictas — NUNCA hacer

```
❌ NO subir data/, results/, .venv/ al repo (.gitignore se encarga)
❌ NO subir señales EEG reales ni datos derivados por sujeto
❌ NO commitear directamente en main (salvo cambios triviales de 1 línea)
❌ NO modificar config/pipeline_config.yaml sin comparar con NeuroMIND/config/single_subject.toml
   (fuente de verdad activa; NeuroMIND/deprecated/code/config/pipeline.toml es una copia archivada,
   ver nota 2026-07-21 en la sección 8)
❌ NO duplicar dependencias en requirements.txt (fue eliminado; pyproject.toml es la única fuente)
❌ NO cambiar lógica científica (filtros, ICA, wPLI, PSD) sin comparar contra NeuroMIND/Julia
❌ NO importar de scripts/ desde otros scripts (cada script es entrypoint autónomo)
```

### Identidad Git

```bash
git config --global user.name "Rafael Castro Triguero"
git config --global user.email "me1catrr@uco.es"
```

---

## 3. Pipeline EEG — 7 fases

```
Etapa A (una vez)        Etapa B (cada sujeto × condición)         Etapa C (grupal)
                            ┌─[3] QC + Filtrado ─────────┐
build_bids_full.py          │                            │
  → data/BIDS/raw/*.json    │  [4] ICA (FastICA, 30 IC)  │           transversal/
                            │                            │ ────────→ longitudinal
                            │  [5] Epochs + AR (±70 µV)  │           (pendiente)
                            │                            │
                            │  [6] PSD (Welch, Hamming)  │
                            │                            │
                            └─[7] wPLI (7 bandas) ───────┘
                              run_full_pipeline.py
```

**Regla crítica:** ICA (Fase 4) se ejecuta sobre la **señal continua filtrada**,
SIEMPRE antes de segmentar (Fase 5). El orden no es negociable.

| Fase | Implementación | Output canónico |
|------|----------------|-----------------|
| 3 | `mne.filter` Butterworth, perfil `eeg_julia` (notch causal + bandreject causal + HP/LP zero-phase) | `cache/filtered_{COND}.npz` |
| 4 | `mne.preprocessing.ICA(method="fastica")` + features custom (`compute_ica_features`) | `cache/ica_{COND}-ica.fif` + `tables/ica_summary_{COND}.json` |
| 5 | `mne.make_fixed_length_epochs` (1 s, sin solape) + baseline `first_window_mean` + AR ±70 µV | `cache/epochs_{COND}-epo.fif` |
| 6 | `mne.time_frequency.psd_array_welch` Hamming + zero-padding | `tables/psd_by_channel_{COND}.csv` + `tables/band_power_summary_{COND}.csv` |
| 7 | `mne_connectivity.spectral_connectivity_epochs(method="wpli")` across-segments | `tables/wpli_{BANDA}_{COND}.csv` (×7) + `tables/connectivity_summary_{COND}.json` |

---

## 4. Archivos clave

| Archivo | Rol |
|---------|-----|
| `config/pipeline_config.yaml` | Parámetros científicos (espejo de `NeuroMIND/config/single_subject.toml`) |
| `pyproject.toml` | Dependencias + setup del paquete instalable (`mne-brain 0.1.0`) |
| `scripts/build_bids_full.py` | **Etapa A** — Lee inventory.csv + .vhdr → genera 205 JSON BIDS en `data/BIDS/raw/` |
| `scripts/run_full_pipeline.py` | **Etapa B** — Entrypoint unificado (single subject + batch). Detecta modo automáticamente según flags |
| `scripts/run_phase3..7_m05.py` | Scripts auxiliares por fase (debug, no usar en producción) |
| `scripts/run_phase8_validation.py` | Comparación cross-pipeline NeuroMIND vs mne_brain → `validation/reports/` |
| `src/mne_brain/bids/loader.py` | `load_eeg_bids` con resolución de `vhdr_path` y fallback TSV |
| `src/mne_brain/preprocessing/` | `filter_recording` + `qc_report` (Fase 3) |
| `src/mne_brain/ica/` | `run_ica` + `compute_ica_features` + `apply_ica_rejection` (Fase 4) |
| `src/mne_brain/processing/epochs.py` | `make_epochs` + `apply_baseline` + `reject_artifacts` (Fase 5) |
| `src/mne_brain/spectral/psd.py` | `compute_psd` + `save_spectral_results` (Fase 6) |
| `src/mne_brain/connectivity/wpli.py` | `compute_wpli` + `save_connectivity_results` (Fase 7) |
| `src/mne_brain/surrogate/` | **Vacío** — Fase 8 (surrogates + FDR) por implementar |
| `tests/unit/` | 34 tests pytest (preprocessing, processing, spectral, connectivity, ica, config, bids) |

---

## 5. Comandos frecuentes

```bash
# Activar entorno (siempre primero)
cd NeuroMIND/mne_brain
source .venv/bin/activate

# Tests
pytest                                         # 34 tests
pytest --co -q                                 # solo discovery, sin ejecutar
pytest tests/unit/test_spectral.py -v          # un módulo

# Verificar sintaxis de un script sin ejecutarlo
python -c "import importlib.util; importlib.util.spec_from_file_location('x','scripts/run_full_pipeline.py').loader.exec_module(__import__('types').ModuleType('x'))"

# Etapa A — (re)construir metadatos BIDS
python scripts/build_bids_full.py              # default: todos los recordings
python scripts/build_bids_full.py --overwrite  # forzar reescritura

# Etapa B — un sujeto
python scripts/run_full_pipeline.py --subject M07 --session T1
python scripts/run_full_pipeline.py --subject M20 --session T2 --from-phase 5

# Etapa B — todos los sujetos en background
mkdir -p results/logs
nohup python -u scripts/run_full_pipeline.py --all --skip-done \
    > results/logs/batch_$(date +%Y%m%d_%H%M%S).log 2>&1 &

# Monitorizar batch
tail -f results/logs/batch_*.log
awk -F, 'NR>1{print $8}' results/batch_pipeline_log.csv | sort | uniq -c

# Comparar con NeuroMIND/Julia
python scripts/run_phase8_validation.py --subject M05 --session T2 --condition EC

# Verificar antes de commit (debe devolver vacío)
git ls-files | grep -E '(^data/|^results/|\.DS_Store$|^\.venv/|\.egg-info)'
```

---

## 6. Estructura del repositorio (post-limpieza 2026-05-26)

```text
mne_brain/
├── config/
│   └── pipeline_config.yaml        ← parámetros científicos
├── data/
│   └── BIDS/                       ← metadatos ligeros (no datos binarios)
│       ├── raw/                    ← 205 *_eeg_metadata.json
│       ├── electrodes/             ← posiciones 10-20 (TSV)
│       ├── participants.tsv
│       ├── groups.csv
│       ├── longitudinal_pairs.csv
│       ├── dataset_description.json
│       └── BrainVision_hardware_software_metadata.md
├── scripts/
│   ├── build_bids_full.py          ← Etapa A
│   ├── run_full_pipeline.py        ← Etapa B (single + batch unificado)
│   ├── run_phase3_m05.py …         ← auxiliares por fase (debug)
│   └── run_phase8_validation.py    ← comparación cross-pipeline
├── src/mne_brain/                  ← paquete instalable
│   ├── bids/        common/        preprocessing/
│   ├── ica/         processing/    spectral/
│   ├── connectivity/  surrogate/   (surrogate vacío, pendiente)
│   └── __init__.py
├── tests/
│   └── unit/                       ← 34 tests pytest
├── validation/
│   └── reports/                    ← informes diff Python vs Julia
├── results/                        ← NO va a git
│   ├── batch_pipeline_log.csv
│   ├── bids_build_log.csv
│   ├── logs/                       ← stdout de batches
│   └── subjects/                   ← outputs por sujeto
├── pyproject.toml
├── README.md
└── AGENTS.md                       ← este documento
```

**Carpetas eliminadas en la limpieza del 2026-05-26:**
`mne_brain/mne_brain/`, `validation/{julia_reference,python_outputs}/`,
`tests/{integration,regression}/`, `results/{Connectivity,Preprocessing,Processing,Spectral}/`.
Eran scaffolding vacío del 7 de abril, sin referencias en código.
También se eliminó `requirements.txt` (duplicaba `pyproject.toml`).

---

## 7. Outputs por grabación

```text
results/subjects/sub-{ID}/ses-{T1|T2}/{eyesclosed|eyesopen}/
├── cache/      filtered_*.npz · *-ica.fif · cleaned_*.npz · *-epo.fif
├── tables/     qc_*.csv · ica_*.csv · ica_summary_*.json ·
│               psd_*.csv · band_power_*.csv · spectral_params_*.json ·
│               wpli_{BANDA}_*.csv (×7) · connectivity_summary_*.json
└── figures/    psd_spectrum_*.png · wpli_heatmap_{BANDA}_*.png (×7) ·
                wpli_all_bands_*.png · ica/ica_topomaps_*.png
```

**Nombres CSV/JSON intencionalmente iguales** a los de NeuroMIND/Julia para
permitir `diff` directo entre las dos implementaciones.

---

## 8. Configuración (`pipeline_config.yaml`)

> **Nota (2026-07-21):** `NeuroMIND/config/pipeline.toml` se archivó en
> `NeuroMIND/deprecated/code/config/pipeline.toml` (solo lo usa el orquestador legacy
> `run_pipeline.jl`/`Pipeline.jl`, ya no la cadena activa). La fuente de
> verdad viva para estos parámetros es ahora `NeuroMIND/config/single_subject.toml`
> — mismos valores en la fecha de este archivado, pero es ese fichero el que
> hay que vigilar para futuros cambios, no `deprecated/code/config/pipeline.toml`.

Espejo de `NeuroMIND/config/single_subject.toml` (histórico: antes `pipeline.toml`). Secciones clave:

```yaml
filtering:
  profile: eeg_julia          # filtra exactamente como Julia
  highpass_hz: 0.5
  lowpass_hz: 150.0
  notch_hz: 50.0
  bandreject_lo: 99.5
  bandreject_hi: 100.5

ica:
  method: fastica
  n_components: 30
  random_seed: 42

artifact_rejection:
  profile: eeg_julia
  max_amplitude_uv: 70.0      # ±70 µV, agresivo en sujetos ruidosos

bands:
  DELTA: [0.5, 4.0]
  THETA: [4.0, 8.0]
  ALPHA: [7.8, 11.7]          # límites exactos de NeuroMIND, no 8-12 estándar
  BETA_LOW: [12.0, 15.0]
  BETA_MID: [15.0, 18.0]
  BETA_HIGH: [18.0, 30.0]
  GAMMA: [30.0, 50.0]
```

Cualquier cambio en estos valores debe sincronizarse con NeuroMIND/Julia.

---

## 9. Estado actual del proyecto (2026-05-26)

### Hecho ✓

- [x] Paquete instalable bajo `src/mne_brain` (8 submódulos)
- [x] Configuración YAML espejada de NeuroMIND
- [x] Loaders BIDS + BrainVision (`mne.io.read_raw_brainvision`)
- [x] **Pipeline Fases 3–7 completo y funcional**
- [x] **Etapa A** — `build_bids_full.py`: 205 metadata JSON generados sin errores en ~2 s, usando MNE Python para leer cabeceras
- [x] **Etapa B** — `run_full_pipeline.py`: entry point unificado single + batch
  - Detección automática de modo (single si `--subject X --session Y`, batch si `--all` o filtros)
  - Tolerancia a fallos (un sujeto roto no detiene el batch)
  - Skip-done basado en `connectivity_summary_{COND}.json`
  - Log CSV por fila (`results/batch_pipeline_log.csv`)
- [x] **Primer run completo del dataset:** 185 grabaciones OK, ~21 minutos totales
- [x] 34 tests unitarios pasando (`tests/unit/`)
- [x] Script de validación cruzada `run_phase8_validation.py` (compara con NeuroMIND/Julia para M05/T2/EC)
- [x] README.md en español con flujo, comandos, estructura completa
- [x] Limpieza de carpetas vacías residuales (ver §6)
- [x] Consolidación: 2 scripts (`run_full_pipeline_m05.py` + `run_batch_pipeline.py`) → 1 único (`run_full_pipeline.py`)
- [x] Eliminación de `requirements.txt` (`pyproject.toml` es la fuente única de deps)

### Resultado del último batch (referencia)

```
results/batch_pipeline_log.csv (último estado conocido):
  185 ok
   21 error  → todos RuntimeError "Epochs-object is empty"
              (sujetos donde AR ±70 µV descarta todas las épocas;
               mismo comportamiento esperable en NeuroMIND/Julia)
    2 user interrupted
  236 skip   → relanzamientos con --skip-done

Tiempo total primer run: 20.6 minutos para 205 grabaciones
```

Los errores afectan principalmente a M11 y otros sujetos con muchos artefactos
musculares/oculares. No son bugs del pipeline.

### Pendiente

- [ ] **Fase 8 — Surrogates + FDR** (módulo `src/mne_brain/surrogate/` está vacío).
  Equivalente al `[SUR]` de NeuroMIND/Julia: phase-shuffle, BH-FDR par-a-par por banda.
- [ ] **Análisis transversal grupal** (MS vs HC por banda)
- [ ] **Análisis longitudinal pareado** (T1 vs T2 dentro de pacientes MS)
- [ ] **Tests de integración** end-to-end del pipeline completo
- [ ] **CI con GitHub Actions** (syntax check + pytest)
- [ ] **Tests de regresión numérica** que comparen mne_brain vs NeuroMIND/Julia bit a bit (módulo `tests/regression/` aún por crear si se necesita)
- [ ] Investigar por qué la cadena de filtrado da `mean wPLI` ligeramente distinto al de NeuroMIND/Julia en GAMMA (esperado por `mne.filter` vs `DSP.jl`, pero documentar)

---

## 10. Diferencias intencionadas con NeuroMIND/Julia

| Aspecto | NeuroMIND/Julia | mne_brain |
|---------|-----------------|-----------|
| Reader BrainVision | Lector binario propio (`BrainVisionLoader.jl`) | `mne.io.read_raw_brainvision` |
| ICA | FastICA puro Julia (`ICACore.jl`) | `mne.preprocessing.ICA` (sklearn FastICA) |
| Filtros | Butterworth `DSP.jl` | Butterworth `mne.filter` |
| Hilbert (wPLI) | FFTW propio | `scipy.signal.hilbert` |
| PSD | Welch + Hanning custom | `mne.time_frequency.psd_array_welch` Hamming |
| CSD | `apply_csd` opcional, off por defecto | No implementado (se asume `use_csd=false`) |
| Surrogates | Implementado (`[SUR]`) | **Pendiente** |
| Pares ALPHA | 7.8–11.7 Hz (no 8–12) | 7.8–11.7 Hz (idem, fijado en config) |

**Política:** No equiparar implementaciones a nivel de algoritmo; mantener cada
una idiomática a su ecosistema. Comparar siempre **outputs** (`.csv`, `.json`)
con `validation/reports/`, nunca código línea a línea.

---

## 11. Workflow Git

```bash
# Verificar antes de cualquier commit
git status --short
git diff --stat
git ls-files | grep -E '(^data/|^results/|\.venv/|\.egg-info|\.DS_Store)'
# debe devolver NADA

# Rama de trabajo
git checkout -b feat/<descripcion>      # nueva feature
git checkout -b fix/<descripcion>       # bug fix
git checkout -b docs/<descripcion>      # docs/README

# Commit + push
git add scripts/... src/... README.md
git commit -m "feat(mne_brain): descripción corta"
git push -u origin feat/<descripcion>
gh pr create --title "..." --body "..."
```

**Lo que NUNCA va al repo:**

- `data/` — señales EEG reales (los JSONs de BIDS sí, son metadata sin payload)
- `results/` — derivados por sujeto, logs
- `.venv/`, `.pytest_cache/`, `__pycache__/`, `*.egg-info/`
- `.DS_Store`, `.env`

---

## 12. Dependencias (`pyproject.toml`)

```toml
mne>=1.7
mne-connectivity>=0.7
numpy>=1.24
scipy>=1.10
pandas>=2.0
matplotlib>=3.7
pyyaml>=6.0
h5py>=3.9

# dev (opcional)
ipykernel>=6.25
jupyterlab>=4.0
pytest>=8.0
ruff>=0.5
```

Instalar con `pip install -e ".[dev]"` desde la raíz de `mne_brain/`.
**No** mantener un `requirements.txt` paralelo: ya se eliminó por duplicar.

---

## 13. Notas para el agente IA

- El **README.md** (en español) es la referencia primaria para usuarios.
  Este AGENTS.md es la referencia para agentes IA.
- **Antes de tocar lógica científica** (filtros, ICA, wPLI, PSD), comparar el código
  equivalente en `../src/SingleSubjectPipeline.jl` (NeuroMIND/Julia) y el output
  canónico en `../results/subjects/sub-{ID}/...`.
- **El batch en background** (`nohup python ... &`) puede tardar 15-30 min en
  procesar los 205 recordings. No interrumpir sin necesidad.
- **`scripts/run_full_pipeline.py` es el único entry point** del pipeline.
  Los `scripts/run_phase{3..7}_m05.py` son auxiliares legacy para debug por
  fase aislada; preferir `run_full_pipeline.py --from-phase N` cuando sea posible.
- **Errores típicos al ejecutar batch:**
  - `RuntimeError: Epochs-object is empty` → sujeto ruidoso, AR descarta todo. Documentado.
  - `ModuleNotFoundError: mne` → falta activar `.venv`.
  - `PermissionError: ~/.mne/mne-python.json` → en sandbox; ejecutar con permisos
    `all` o `full_network` según el agente.
- **Sincronización con NeuroMIND/Julia:** si NeuroMIND cambia parámetros en
  `config/single_subject.toml` (fuente activa; `deprecated/code/config/pipeline.toml` es
  copia archivada), actualizar **inmediatamente** `mne_brain/config/pipeline_config.yaml`.
  Las dos implementaciones deben usar exactamente los mismos parámetros.

> Para contexto completo del proyecto padre, leer también `../CLAUDE.md` y
> `../AGENTS.md` en la raíz de NeuroMIND.
