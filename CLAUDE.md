# CLAUDE.md — NeuroMIND Project Brain

> Este archivo es leído automáticamente por Claude Code al inicio de cada sesión.
> Contiene todo el contexto necesario para trabajar en el proyecto desde cualquier ordenador.
> **Mantenerlo actualizado es prioritario.** Última actualización: 2026-05-22 (Phase 6 EEG_Julia completa).

---

## 0. Identidad del desarrollador

| Campo | Valor |
|-------|-------|
| **Nombre** | Rafael Castro Triguero |
| **Email** | me1catrr@uco.es |
| **GitHub** | https://github.com/me1catrr |
| **Repo** | https://github.com/me1catrr/NeuroMIND |

```bash
# Ejecutar en cualquier ordenador nuevo antes del primer commit:
git config --global user.name "Rafael Castro Triguero"
git config --global user.email "me1catrr@uco.es"
```

---

## 1. Identidad del proyecto

**NeuroMIND** es un framework de análisis EEG para estudio de conectividad funcional en
**Esclerosis Múltiple (EM)** usando **weighted Phase Lag Index (wPLI)**.

- Repositorio: https://github.com/me1catrr/NeuroMIND
- Lenguaje: Julia 1.9+
- Dashboard: web local via Genie.jl (no GitHub Pages)
- Datos: formato BIDS, señales EEG en reposo (ojos cerrados `EC` / ojos abiertos `EO`)
- Proyecto relacionado (original): `EEG_Julia` / `NeuroSmart-EEG` (mismo directorio padre)

**Rutas locales típicas** (en cualquier Mac del usuario):
```
.../EEG_Julia/NeuroMIND/          ← este repo
.../EEG_Julia/                    ← proyecto original (solo referencia, no subir)
```

---

## 2. Contexto científico

| Concepto | Descripción |
|----------|-------------|
| **EEG en reposo** | Señal continua, 31 canales, ~500 Hz, condición ojos cerrados |
| **wPLI** | Conectividad funcional entre pares de canales por banda; robusto a campo de volumen |
| **ICA** | Separación de fuentes para eliminar artefactos oculares/musculares/cardíacos |
| **CSD** | Current Source Density — reduce campo de volumen; **opcional** (`use_csd = false` por defecto) |
| **Bandas** | δ(0.5–4), θ(4–8), α(7.8–11.7), β_low(12–15), β_mid(15–18), β_high(18–30), γ(30–50) Hz |
| **Hipótesis** | Pacientes EM presentan alteraciones en conectividad α y β respecto a controles sanos |

**Diferencia clave EEG_Julia vs NeuroMIND:**
- EEG_Julia aplica CSD antes de wPLI → el wPLI difiere significativamente
- NeuroMIND: `use_csd = false`, filtrado diferente, ICA propio, rechazo AR diferente
- El espectral (PSD) conserva bien el patrón espacial pero cambia la escala absoluta
- **Tarea futura**: modo reproducibilidad `EEG_Julia/BrainVision` para alinear parámetros

---

## 3. Pipeline de 8 pasos (`run_single_subject_pipeline`)

```
[1/8] Carga          BIDSLoader → EEGRecording (channels × samples)
[2/8] QC             flag_bad_channels por z-score (umbral 3.0)
[3/8] Filtrado       HP 0.5 Hz + LP 48 Hz + Notch 50 Hz + Bandreject (Butterworth filtfilt)
[4/8] ICA            FastICA simétrico sobre señal CONTINUA filtrada → ICAResult
[5/8] Segmentación   segment_recording(rec_ICA) → EpochSet + baseline + rechazo AR
[6/8] Espectral      PSD Hanning + potencia por banda → SpectralResult
[7/8] wPLI           Hilbert analítica across-segments → ConnectivityMatrix (± CSD)
[8/8] Guardado       CSV + PNG + JSON + log en results/subjects/sub-{id}/ses-{sess}/{task}/
```

**Regla crítica:** ICA se ejecuta sobre señal continua (paso 4), ANTES de segmentar (paso 5).
La señal limpiada `rec_ica` alimenta directamente la segmentación y todo lo posterior.

### ICA — implementación propia (sin MultivariateStats)
- Archivo: `src/ica/ICACore.jl`
- PCA whitening → FastICA simétrico con no-linealidad `tanh`
- Solo `LinearAlgebra` + `Random` (stdlib Julia, sin paquetes externos)
- Defaults: 30 componentes, 500 iter, tol=1e-5, seed=42
- Mixing matrix: `A = pinv(W_total)` — necesario porque n_comp < n_channels
- Varianza explicada: `||A[:,i]||² × var(S[i,:]) / total_var`

---

## 4. Dashboard web (Genie.jl)

13 paneles (Phase 0–12). Implementados completamente: **0, 1, 2, 3, 4, 5**.
Paneles 6–12: placeholder vacío (próximos a implementar).

| Panel | Nombre | Estado |
|-------|--------|--------|
| 0 | Proyecto / Dataset | ✓ |
| 1 | BIDS & Metadata | ✓ |
| 2 | Señal cruda | ✓ |
| 3 | QC inicial | ✓ |
| 4 | Preprocessing / Filtrado | ✓ |
| 5 | ICA | ✓ |
| 6 | Segmentación | ✓ |
| 7 | Rechazo artefactos | ✓ |
| 8 | Análisis espectral | placeholder |
| 9 | Conectividad wPLI | placeholder |
| 10 | Surrogates / Inferencia | placeholder |
| 11 | Resultados finales | placeholder |
| 12 | Exportación / Reporte | placeholder |

**Panel 5 ICA** muestra: resumen stats, **grid de topomaps paginado** (6 por página),
**tabla de features de clasificación** (frontal/temporal/blink/emg/line ratio, kurtosis),
tabla de componentes, donut SVG, señal temporal del componente (Canvas), espectro (SVG DFT),
butterfly antes/después (Canvas), componentes rechazados, métricas de calidad, vista de artefactos.
Profile badge indica si el perfil es `eeg_julia` o `default`.

Si ICA no se ha ejecutado → placeholder elegante con botón de acción.

---

## 5. Archivos clave

| Archivo | Rol |
|---------|-----|
| `src/NeuroMIND.jl` | Entry point del módulo; lista de `include()` en orden |
| `src/types.jl` | Tipos centrales: `EEGRecording`, `EpochSet`, `ICAResult`, `SpectralResult`, `ConnectivityMatrix`, `PipelineConfig`, `RecordingMeta` |
| `src/SingleSubjectPipeline.jl` | Pipeline 8 pasos + `load_ss_config` + helpers de guardado |
| `src/webapp/App.jl` | Servidor Genie + todas las rutas API (`/api/phase*`, `/api/ica_*`, etc.) |
| `web/views/dashboard.html` | SPA completa del dashboard (~6300 líneas; HTML + CSS + JS inline) |
| `src/ica/ICACore.jl` | FastICA puro Julia (PCA whitening + tanh); perfiles `eeg_julia`/`default` |
| `src/ica/ICAClassification.jl` | `compute_ica_features` (7 features) + `evaluate_ica_components` (scores → labels) |
| `src/ica/ICAInspection.jl` | `load_ica_labels` (multi-path), `apply_ica_rejection` |
| `src/preprocessing/Filtering.jl` | `filter_recording` (HP/LP/Notch/Bandreject, Butterworth) |
| `src/segmentation/Epochs.jl` | `segment_recording`, `apply_baseline`, `reject_artifacts` |
| `src/spectral/PowerSpectrum.jl` | `compute_psd`, `plot_spectrum_grid` |
| `src/connectivity/wPLI.jl` | `compute_wpli` (Hilbert analítica, across-segments) |
| `src/connectivity/CSD.jl` | `apply_csd` (opcional, `use_csd = false` por defecto) |
| `config/single_subject.toml` | Parámetros para `run_single_subject.jl` |
| `scripts/run_single_subject.jl` | Entrada CLI del pipeline |
| `scripts/launch_dashboard.jl` | Lanza Genie en http://localhost:8080 |
| `tests/runtests.jl` | Suite de tests unitarios |
| `docs/GITHUB_WORKFLOW.md` | Flujo de ramas y commits |
| `docs/MIGRATION_GUIDE.md` | Guía de migración desde EEG_Julia |

---

## 6. Resultados generados por el pipeline

```
results/subjects/sub-{id}/ses-{sess}/{task}/
├── overview.csv
├── qc_summary.csv
├── channel_statistics.csv
├── psd_by_channel.csv
├── band_power_summary.csv
├── connectivity_edges.csv
├── wpli_{band}.csv              ← una por banda
├── ica_components.csv           ← índice, varianza, artifact_type, rechazado
├── ica_component_features.csv   ← 7 features + 4 scores + artifact_type por IC
├── ica_mixing_matrix.csv        ← A (n_ch × n_comp)
├── ica_unmixing_matrix.csv      ← W_total (n_comp × n_ch)
├── ica_summary.json             ← n_comp, n_rej, varianza retenida, has_features, n_topomaps
├── ica_activations.csv          ← primeros 10s; columnas: t_s, IC1…IC30
├── ica_signal_before.csv        ← primeros 10s señal filtrada (todos los canales)
├── ica_signal_after.csv         ← primeros 10s señal limpiada por ICA
├── figures/ica_topomap_001.png  ← topomaps per IC (si hay ch_positions)
├── pipeline_log.txt
└── config_snapshot.toml
```

**Estos archivos NUNCA van al repo.** Están en `.gitignore` vía `/results/`.

---

## 7. Reglas Git estrictas

```bash
# Antes de cualquier commit:
git status --short
git diff --stat

# Antes de push, verificar que no hay datos:
git ls-files | grep -E '(^data/|^results/|\.DS_Store$|^\.claude/|^\.vscode/)'
# → debe devolver NADA

# Ramas: nunca commitear directamente en main salvo cambios triviales de 1 línea
# Rama de trabajo: feat/<nombre>, fix/<nombre>, docs/<nombre>
# PR desde rama → main (o push directo si es solo docs/README)
```

**Lo que NUNCA va al repo:**
- `data/` — señales EEG, BIDS raw
- `results/` — derivados por sujeto
- `reports/`, `exports/`
- `.claude/`, `.vscode/`, `.cursor/`
- `.DS_Store`, `.env`, logs clínicos, secretos

---

## 8. Comandos frecuentes

```bash
# Instalar dependencias
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Verificar sintaxis
julia --project=. -e 'include("src/NeuroMIND.jl"); println("OK")'

# Tests
julia --project=. tests/runtests.jl

# Pipeline (sujeto detectado automáticamente)
julia --project=. scripts/run_single_subject.jl

# Dashboard local
julia --project=. scripts/launch_dashboard.jl
# → http://localhost:8080

# GitHub
gh auth status
gh repo view me1catrr/NeuroMIND
git push -u origin <rama>
```

---

## 9. Estado actual del proyecto (2026-05-22)

### Implementado y funcionando
- [x] Pipeline 8 pasos completo con ICA
- [x] FastICA puro Julia (sin MultivariateStats); perfiles `eeg_julia` / `default`
- [x] Dashboard paneles 0–7
- [x] API routes: `/api/phase5_ica_info`, `/api/ica_activation`, `/api/ica_signal`, `/api/ica_features`
- [x] API route: `/api/phase6_segmentation` → segmentation_summary.json (enriquecido), segments_table.csv, channel_coverage.csv
- [x] API route: `/api/phase7_ar` → artifact_rejection_summary.json, rejected_segments.csv, channel_artifact_summary.csv
- [x] `compute_ica_features` — 7 features por componente (portado de EEG_Julia)
- [x] `evaluate_ica_components` — scores ocular/muscle/line/jump → labels
- [x] `_save_ica_topomaps` — genera PNGs por IC (requiere ch_positions en BIDS)
- [x] Dashboard Phase 5: topomap grid paginado (6 columnas, 36/pág), features table, profile badge
- [x] 241 tests unitarios pasando
- [x] `compute_epoch_quality_report` con worst_channel + p2p_uv, respeta perfil AR
- [x] **Phase 6 EEG_Julia profile** (rama `feat/phase6-eeg-julia-segmentation`):
  - `segment_recording` perfil `"eeg_julia"` → fuerza 1 s, sin solapamiento
  - `apply_baseline` método `"first_window_mean"` → media [0, 0.10 s] por canal (EEG_Julia exacto)
  - `reject_artifacts` perfil `"eeg_julia"` → ±70 µV, primeros 30 canales, sin gradiente
  - Pipeline: doble baseline (`n_passes=2`) — pre-AR + post-AR
  - `segmentation_summary.json` enriquecido con 10 nuevos campos metodológicos
  - Dashboard Phase 6: badges perfil, nota metodológica EEG_Julia, criterio AR detallado
- [x] README.md completo; `.gitignore` estricto; CLAUDE.md actualizado

### Ramas activas
- `feat/phase5-ica-classification` — commits de Phase 5, pendiente push/PR
- `feat/phase6-eeg-julia-segmentation` — Phase 6 EEG_Julia, pendiente push/PR

### Pendiente / próximo
- [ ] Dashboard paneles 8–12 (espectral → exportación)
- [ ] Tests de integración pipeline completo
- [ ] GitHub Actions CI (syntax check + tests)
- [ ] Push ramas + PR a main

### Bugs conocidos resueltos
- `JSON3.read` → regex parsing en `App.jl` (JSON3 no era dependencia)
- `inv(W_total)` → `pinv(W_total)` para caso n_comp < n_channels
- `MultivariateStats` → implementación propia en `ICACore.jl`
- Varianza explicada: cálculo corregido con norma de columnas de A

---

## 10. Workflow Git recomendado para esta sesión

```bash
# 1. Arrancar siempre así:
cd <ruta-local-NeuroMIND>
git fetch origin
git status --short

# 2. Crear rama para la tarea:
git checkout -b feat/nombre-descriptivo

# 3. Trabajar, luego commitear:
git add src/... web/...
git commit -m "feat: descripción corta"

# 4. Push y PR:
git push -u origin feat/nombre-descriptivo
gh pr create --title "..." --body "..."

# 5. Merge a main cuando esté listo
```

---

## 11. Dependencias Julia (`Project.toml`)

```toml
Base64, CSV, CairoMakie, DSP, DataFrames, Dates, FFTW,
Genie, LinearAlgebra, Random, Serialization, Statistics,
StatsBase, TOML
```

**No hay** `MultivariateStats`, `Pluto`, `JSON3`.

---

## 12. Notas para Claude Code

- El dashboard es una **SPA de ~6300 líneas** (`web/views/dashboard.html`). Siempre leer la sección relevante antes de editar; no releer completo.
- `App.jl` tiene las rutas API y helpers de filtrado para el viewer interactivo de la Fase 4.
- `SingleSubjectPipeline.jl` contiene el pipeline + visualizaciones + helpers de guardado en un solo archivo.
- Al tocar lógica científica (filtrado, ICA, wPLI, PSD), **comparar siempre** con el original en `EEG_Julia/` antes de cambiar comportamiento.
- Julia syntax check rápido: `julia --project=. -e 'include("src/NeuroMIND.jl"); println("OK")'`
