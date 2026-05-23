# CLAUDE.md — NeuroMIND Project Brain

> Este archivo es leído automáticamente por Claude Code al inicio de cada sesión.
> Contiene todo el contexto necesario para trabajar en el proyecto desde cualquier ordenador.
> **Mantenerlo actualizado es prioritario.** Última actualización: 2026-05-23 (Paneles 11-12 completos — todos los paneles del dashboard implementados).

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

15 paneles (Phase 0–14). **Todos completamente implementados** ✅.

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
| 8 | Análisis espectral | ✓ |
| 9 | Conectividad wPLI | ✓ |
| 10 | Surrogates / Inferencia | ✓ |
| 11 | Resultados finales | ✓ |
| 12 | Exportación / Reporte | ✓ |
| 13 | Evaluación Transversal | ✓ |
| 14 | Evaluación Longitudinal | ✓ |

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
| `web/views/dashboard.html` | SPA completa del dashboard (~12 600 líneas; HTML + CSS + JS inline) |
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
├── wpli_{band}.csv                      ← una por banda (wPLI observado)
├── connectivity_summary.json            ← resumen wPLI por banda
├── wpli_observed_{band}.csv             ← matriz wPLI observada por banda (surrogates)
├── wpli_pvalues_{band}.csv              ← p-valores par×par por banda
├── wpli_qvalues_{band}.csv              ← q-valores BH par×par por banda
├── wpli_significant_{band}.csv          ← máscara binaria de significancia
├── surrogate_null_stats_{band}.csv      ← (null_mean, null_std) por par (aprox. Gaussiana)
├── significant_connections.csv          ← top conexiones sig. (ch1, ch2, wPLI, p, q)
├── surrogate_summary.json               ← resumen por banda: n_sig, pct_sig, mean_p, fdr_thr
├── surrogate_quality.csv                ← QC checks por banda
├── ica_components.csv                   ← índice, varianza, artifact_type, rechazado
├── ica_component_features.csv           ← 7 features + 4 scores + artifact_type por IC
├── ica_mixing_matrix.csv                ← A (n_ch × n_comp)
├── ica_unmixing_matrix.csv              ← W_total (n_comp × n_ch)
├── ica_summary.json                     ← n_comp, n_rej, varianza retenida, has_features, n_topomaps
├── ica_activations.csv                  ← primeros 10s; columnas: t_s, IC1…IC30
├── ica_signal_before.csv                ← primeros 10s señal filtrada (todos los canales)
├── ica_signal_after.csv                 ← primeros 10s señal limpiada por ICA
├── figures/ica_topomap_001.png          ← topomaps per IC (si hay ch_positions)
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
- [x] Pipeline 8 pasos completo con ICA + paso opcional [SUR] surrogates
- [x] FastICA puro Julia (sin MultivariateStats); perfiles `eeg_julia` / `default`
- [x] Dashboard paneles 0–10 completamente implementados
- [x] API routes: `/api/phase5_ica_info`, `/api/ica_activation`, `/api/ica_signal`, `/api/ica_features`
- [x] API route: `/api/phase6_segmentation`
- [x] API route: `/api/phase7_ar` + `/api/phase7_epoch_signal`
- [x] API route: `/api/phase8_spectral` — PSD, topomaps IDW en Canvas, bandpower
- [x] API route: `/api/phase9_wpli` — ConnectivityMatrix, heatmap, network, band tabs
- [x] API route: `/api/phase10_surrogates` — surr_run flag, matrix_obs, sig_connections, null_stats, QC, p-value hist, timing
- [x] `compute_ica_features`, `evaluate_ica_components`, `_save_ica_topomaps`
- [x] Dashboard Phase 8: espectro PSD, topomaps IDW Canvas, bandpower por banda, comparativa
- [x] Dashboard Phase 9: heatmap wPLI Canvas (reusa `_p9Colormap`), red de conectividad con posiciones EEG, band tabs
- [x] **Phase 11 Resultados Finales** (tema #059669 emerald):
  - `/api/phase11_summary`: agrega overview, ica_summary, band_power_summary, connectivity_summary, surrogate_summary, significant_connections; calcula `quality_score` ponderado (canales, epochs, varianza ICA, fases completas) y `statuses` por fase
  - Dashboard: 6 KPIs, checklist 0–10 con iconos (ok/warn/pend), gauge SVG de calidad, perfil espectral SVG, tabla top 8 conexiones sig., grid wPLI por banda con barras comparativas
- [x] **Phase 12 Exportación / Reporte** (tema #0369a1 sky blue):
  - `/api/phase12_files`: escanea el directorio del sujeto, agrupa archivos en 7 categorías (QC, ICA, Espectral, wPLI, Surrogates, wPLI-por-banda, Figuras), devuelve tamaños y estado
  - Dashboard: 4 botones de acción (copiar JSON, CSV inventario, ver resultados, actualizar), grid de categorías con barra de progreso y lista de archivos, panel de texto exportable con inventario completo
- [x] **Phase 13 Evaluación Transversal** (tema #0891b2, 13 secciones):
  - `scripts/run_transversal_analysis.jl`: lee `groups.csv` + wPLI individuales → genera archivos grupales
  - Archivos: `group_connectivity_{ctrl|ms}_{band}.csv`, `group_difference_{band}.csv`, `group_statistics_{band}.csv`, `significant_edges_{band}.csv`, `band_statistics.csv`, `subject_inclusion.csv`, `subject_band_means.csv`, `transversal_summary.json`
  - `/api/phase13_transversal` → matrices, edges, estadísticas, distribución, inclusión
  - Dashboard: 3 heatmaps (ctrl/ms/diff divergente), red significativa, top edges, distribución SVG, tabla bandas, inclusión QC, 6 acciones
- [x] **Phase 14 Evaluación Longitudinal** (tema #d97706 ámbar, 13 secciones):
  - `scripts/run_longitudinal_analysis.jl`: auto-detecta pares T1/T2 o lee `longitudinal_pairs.csv` → paired t-test + Cohen's dz + BH-FDR por banda
  - Archivos: `longitudinal_connectivity_t1_{band}.csv`, `longitudinal_connectivity_t2_{band}.csv`, `longitudinal_difference_{band}.csv`, `longitudinal_statistics_{band}.csv`, `significant_longitudinal_edges_{band}.csv`, `band_statistics_longitudinal.csv`, `subject_band_means.csv`, `paired_subjects.csv`, `longitudinal_summary.json`
  - `/api/phase14_longitudinal` → matrices t1/t2/diff, edges significativos top-50, estadísticas, distribución, sujetos pareados
  - Dashboard: timeline T1→T2, 3 heatmaps (T1/T2/diff ámbar-azul), red significativa, top edges, distribución SVG, tabla bandas, tabla sujetos pareados
- [x] Dashboard Phase 10 (tema #7c3aed, 13 secciones):
  - Config surrogates, métricas inferencia, estado fase, archivos generados
  - Heatmap observado, red significativa (solo edges sig., ancho/color por wPLI)
  - Histograma p-valores SVG (barras verde/gris, línea alpha, línea uniforme)
  - Tabla top conexiones, estadísticas por banda/globales
  - Distribución par (Gaussiana N(null_mean,null_std), área p-value, línea observado)
  - 8 checks QC, 6 botones exportación
- [x] `_bh_qvalues` helper en Pipeline — q-values individuales BH con monotonicidad
- [x] `_save_surrogate_results` — 8 archivos por banda + surrogate_summary.json + significant_connections.csv + surrogate_quality.csv
- [x] 241 tests unitarios pasando (Types, Filtering, Segmentation, Spectral, wPLI, FDR, Config, ICA, …)
- [x] **Phase 6 EEG_Julia profile**: `"first_window_mean"` baseline, ±70 µV AR, doble baseline pass
- [x] **Phase 7 AR EEG_Julia profile**: enriquecimiento completo de AR summary, señal real epoch
- [x] README.md completo; `.gitignore` estricto; CLAUDE.md actualizado

### Ramas activas (pendientes push/PR)
- `feat/phase5-ica-classification` — Phase 5 ICA
- `feat/phase6-eeg-julia-segmentation` — Phase 6 segmentation
- `feat/phase7-ar-eeg-julia` — Phase 7 AR
- `feat/phase8-spectral` — Phase 8 espectral
- `feat/phase9-wpli` — Phase 9 conectividad wPLI
- `feat/phase10-surrogates` — Phases 10–14 completos ← **rama actual**

### Pendiente / próximo
- [ ] Push ramas feat/phase8, feat/phase9, feat/phase10 + PR a main
- [ ] Tests de integración pipeline completo
- [ ] GitHub Actions CI (syntax check + tests)

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

- El dashboard es una **SPA de ~13 500 líneas** (`web/views/dashboard.html`). Siempre leer la sección relevante antes de editar; no releer completo.
- `App.jl` tiene las rutas API y helpers de filtrado para el viewer interactivo de la Fase 4.
- `SingleSubjectPipeline.jl` contiene el pipeline + visualizaciones + helpers de guardado en un solo archivo.
- Al tocar lógica científica (filtrado, ICA, wPLI, PSD), **comparar siempre** con el original en `EEG_Julia/` antes de cambiar comportamiento.
- Julia syntax check rápido: `julia --project=. -e 'include("src/NeuroMIND.jl"); println("OK")'`
