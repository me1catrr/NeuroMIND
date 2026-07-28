# AGENTS.md — NeuroMIND

> **Single source of truth for AI agents.**
> Read automatically by Claude Code, Cursor, OpenAI Codex CLI, and any other AI assistant
> on project open. **Keep this file up to date — it is the authoritative project context.**
> Last updated: 2026-07-28.
>
> **README.md is the authoritative source for the pipeline, config/pipeline.toml,
> outputs and results structure.** This file keeps agent-specific rules (Git, dashboard,
> changelog) and cross-references the README to avoid drift.

---

## 0. Developer identity

| Field | Value |
|-------|-------|
| **Name** | Rafael Castro Triguero |
| **Email** | me1catrr@uco.es |
| **GitHub** | https://github.com/me1catrr |
| **Repo** | https://github.com/me1catrr/NeuroMIND |

```bash
git config --global user.name "Rafael Castro Triguero"
git config --global user.email "me1catrr@uco.es"
```

---

## 1. Project identity

**NeuroMIND** is an EEG connectivity analysis framework for studying **Multiple Sclerosis (MS)**
using **weighted Phase Lag Index (wPLI)**.

- Language: Julia 1.9+
- Dashboard: local web server via Genie.jl (`http://localhost:8080`)
- Data: BIDS format, resting-state EEG (eyes-closed `EC` / eyes-open `EO`)
- Related project (original): `../EEG_Julia/` — **do not push to repo**

**Local paths:**
```
.../EEG_Julia/NeuroMIND/   ← this repo
.../EEG_Julia/             ← original project (reference only)
```

---

## 2. Scientific context

| Concept | Description |
|---------|-------------|
| **Resting EEG** | Continuous signal, 31 channels, ~500 Hz, EC/EO conditions |
| **wPLI** | Functional connectivity between channel pairs per frequency band; robust to volume conduction |
| **dwPLI** | Debiased wPLI² estimator (Vinck 2011); available (`use_dwpli = true`) but current config uses classic wPLI |
| **ICA** | Independent Component Analysis — removes ocular/muscle/cardiac artifacts |
| **CSD** | Current Source Density — reduces volume conduction; `use_csd = false` by default |
| **Bands** | δ(0.5–4), θ(4–8), α(7.8–11.7), β_low(12–15), β_mid(15–18), β_high(18–30), γ(30–50) Hz |
| **Hypothesis** | MS patients show connectivity alterations in α and β bands vs healthy controls |

**Key difference EEG_Julia vs NeuroMIND:**
- EEG_Julia applies CSD before wPLI → wPLI values differ significantly
- NeuroMIND: `use_csd = false`, different filtering, custom ICA, different AR rejection
- **Never change scientific logic (ICA, wPLI, PSD) without cross-checking against EEG_Julia/**

---

## 3. 8-step pipeline (`run_single_subject_pipeline`)

```
[1/8] Load          BIDSLoader / BrainVisionLoader → EEGRecording (channels × samples)
[2/8] QC            flag_bad_channels by z-score (threshold 3.0) + amplitude_warning
[3/8] Filter        HP 0.5 Hz + LP 150 Hz + Notch 50 Hz + Bandreject (Butterworth filtfilt)
[4/8] ICA           Symmetric FastICA on CONTINUOUS filtered signal → ICAResult
[5/8] Segment       segment_recording(rec_ICA) → EpochSet + baseline + AR rejection
[6/8] Spectral      Hanning PSD + band power → SpectralResult
[7/8] wPLI          3 configurable methods: Hilbert / FourierCSD / Multitaper → ConnectivityMatrix
[8/8] Save          CSV + PNG + JSON + log → results/subjects/sub-{id}/ses-{sess}/{task}/
```

**Critical rule:** ICA (step 4) MUST always precede segmentation (step 5).

### ICA — custom implementation (no MultivariateStats)
- File: `src/ica/ICACore.jl`
- PCA whitening → symmetric FastICA with `tanh` nonlinearity
- Only `LinearAlgebra` + `Random` (Julia stdlib, no external packages)
- Mixing matrix: `A = inv(W_total)` with profile `eeg_julia` (square, all channels); `pinv` only in reduced `default` profile

### wPLI — multi-method architecture (since 2026-05-26)
Configured in `config/pipeline.toml`:
```toml
[connectivity]
wpli_method = "hilbert"      # ACTIVE — Butterworth + Hilbert analytic signal
# wpli_method = "fourier_csd" # Cross-power spectrum FFT + Hann window
# wpli_method = "multitaper"  # DPSS multitaper (MNE standard)
use_dwpli   = false          # ACTIVE — classic wPLI (dwPLI available)
```
All three methods support `use_dwpli = true`.
Surrogates automatically use the **same estimator** as the observed computation.

---

## 4. Web dashboard (Genie.jl) — 16 panels (Phase 0–15)

| Panel | Name | Status |
|-------|------|--------|
| 0 | Project / Dataset | ✓ |
| 1 | BIDS & Metadata | ✓ |
| 2 | Raw signal | ✓ |
| 3 | Initial QC | ✓ |
| 4 | Preprocessing / Filtering | ✓ |
| 5 | ICA | ✓ |
| 6 | Segmentation | ✓ |
| 7 | Artifact rejection | ✓ |
| 8 | Spectral analysis | ✓ |
| 9 | wPLI connectivity | ✓ |
| 10 | Surrogates / Inference | ✓ |
| 11 | Final results | ✓ |
| 12 | Export / Report | ✓ |
| 13 | Transversal evaluation | ✓ |
| 14 | Longitudinal evaluation | ✓ |
| 15 | MNE-Python validation | ✓ |

---

## 5. Key files

| File | Role |
|------|------|
| `src/NeuroMIND.jl` | Module entry point; ordered `include()` list |
| `src/types.jl` | Core types: `EEGRecording`, `EpochSet`, `ICAResult`, `SpectralResult`, `ConnectivityMatrix` |
| `src/SingleSubjectPipeline.jl` | 8-step pipeline + `load_ss_config` + save helpers |
| `src/webapp/App.jl` | Genie server + all API routes (`/api/phase*`). **Lazy-loaded**: included only on first `launch_webapp` call, so the pipeline never compiles Genie |
| `web/views/dashboard.html` | Full dashboard SPA (~13 500 lines; inline HTML+CSS+JS) |
| `src/ica/ICACore.jl` | Pure-Julia FastICA (PCA whitening + tanh); profiles `eeg_julia`/`default` |
| `src/ica/ICAClassification.jl` | `compute_ica_features` (7 features) + `evaluate_ica_components` |
| `src/preprocessing/Filtering.jl` | `filter_recording` (HP/LP/Notch/Bandreject, Butterworth) |
| `src/segmentation/Epochs.jl` | `segment_recording`, `apply_baseline`, `reject_artifacts` |
| `src/spectral/PowerSpectrum.jl` | `compute_psd`, `plot_spectrum_grid` |
| `src/connectivity/wPLI.jl` | `compute_wpli` — 3 estimators: Hilbert, FourierCSD, Multitaper |
| `src/connectivity/CSD.jl` | `apply_csd` (optional, `use_csd = false` by default) |
| `src/io/BrainVisionLoader.jl` | Native BrainVision binary reader (.vhdr + .eeg, IEEE_FLOAT_32) |
| `config/pipeline.toml` | **Single config** for all active scripts (see README §5) |
| `scripts/run_batch_pipeline.jl` | **Phase C** — batch pipeline with CLI filters |
| `scripts/run_transversal_analysis.jl` | Thin launcher — `include`s `Transversal.jl` and calls `Transversal.run(config_path)` |
| `scripts/run_longitudinal_analysis.jl` | Thin launcher — `include`s `Longitudinal.jl` and calls `Longitudinal.run(config_path)` |
| `src/transversal/Transversal.jl` | **Single module**: stats (Mann–Whitney/Welch/FDR-BH) + all figures (exploratory, every run; forest/raincloud/matrices/network+power, `eyesclosed` only) in one pass. No separate CSV-only regeneration mode — figures always reflect data just computed |
| `src/interactive/plot_transversal.jl` | Interactive MS vs Control viewer → `:8781` (standalone CLI) |
| `src/longitudinal/Longitudinal.jl` | Single module symmetric to `Transversal.jl`. Selects T1–T2 pairs independently for EC/EO, uses Wilcoxon exacto condicional for N≤30, and captures per-subject `mean_strength` T1/T2 on the common montage during `network_global_statistics.csv` (`tables/mean_strength_scores.csv`) — no separate reconstruction pass. Also generates the transversal↔longitudinal summary figure into `results/summary/` when both analyses' data exist on disk (checked from either script, whichever runs second) |
| `src/interactive/plot_longitudinal.jl` | Interactive T1→T2 viewer → `:8780` (standalone CLI) |
| `src/interactive/` | Home for every interactive viewer (not batch PNG generation): the two above + `viewer_support.jl` (shared IO/JSON helpers, ex `GroupVizCommon.jl`) + `viewer_common.js` + `aux/sub-M05_ses-T2_eyesclosed/` (12 standalone single-subject viewers, ports `:8765`–`:8775`, hardcoded to sub-M05/ses-T2/EC — see README §"Visores auxiliares") |
| `mne_brain/` | MNE-Python cross-validation pipeline |
| `test/runtests.jl` | Canonical Julia test suite (36 unit, integration, numerical and regression `@testset`s; direct runner + `Pkg.test()`) |

---

## 6. Pipeline output structure

```
results/subjects/sub-{id}/ses-{sess}/{task}/
├── overview.csv
├── qc_summary.csv
├── channel_statistics.csv
├── band_power_summary.csv
├── connectivity_edges.csv
├── wpli_{band}.csv                    ← one per band (observed wPLI)
├── wpli_pvalues_{band}.csv
├── wpli_qvalues_{band}.csv
├── wpli_significant_{band}.csv
├── surrogate_null_stats_{band}.csv
├── significant_connections.csv
├── surrogate_summary.json
├── ica_summary.json
├── ica_signal_before.csv              ← first 10s filtered signal
├── ica_signal_after.csv               ← first 10s cleaned signal
├── raw_signal.csv                     ← full raw signal (100s)
├── figures/
├── pipeline_log.txt
└── config_snapshot.toml
```

**These files never go to the repo** (`.gitignore` via `/results/`). Group results go to
`results/transversal/{eyesclosed,eyesopen}/{tables,figures}/` and
`results/longitudinal/{eyesclosed,eyesopen}/{tables,figures}/` (siblings of `subjects/`,
same `eyesclosed`/`eyesopen` naming), plus `results/summary/` for the one
cross-analysis figure. See README §7–8
for CSV layout, figure names (`heatmap_triplet_*`, `explore_network_topN_*`), and interactive viewers.

> **Single output tree (since 2026-07-21).** `_save_all_results` writes directly to the
> BIDS `export_dir`; the old dual tree `results/{ID}/{SES}/` with `_EC`/`_EO` suffixes was
> removed, along with the dead readers (`load_dashboard_data`, App.jl "Legacy API" routes).
> ICA cache now lives in `export_dir/cache/`.

> **Cohort viewers (since 2026-07-25).** `plot_*.jl` are standalone CLIs (not included in
> `NeuroMIND.jl`). They serve `/static/group_viewer_common.js` and send full `edge_stats`
> for heatmap tooltips (not only the filtered Top-N/FDR subset). Diverging colormap is
> `Reverse(:RdBu)` (positive Δ = red) in both PNG and UI. Do not change ICA/wPLI science
> when editing viewers — only visualization.

> **Single-subject aux viewers (since 2026-07-25; moved to `src/interactive/aux/` on
> 2026-07-27).** `src/interactive/aux/sub-M05_ses-T2_eyesclosed/*.jl`
> (12 CLIs, ports `:8765`–`:8775`) are the single-subject counterpart of the cohort
> viewers — same pattern (bare `Sockets` HTTP server, no HTTP.jl/Genie), but each reads
> `results/subjects/sub-M05/ses-T2/eyesclosed/` directly and is hardcoded to that
> subject/session/task (no CLI args). Versioned outside `results/` precisely so a future
> clean-slate regeneration of that folder (`docs/TDD_ejecucion_rutinas.md`) can't delete
> them; living in `src/` (not `scripts/`) because each is a full standalone feature, not
> a thin launcher — same reasoning as the cohort viewers. Full table (port ↔ purpose ↔
> input) in README, Anexo §"Visores auxiliares".

---

## 7. Signal quality policy (QC v2 — 2026-05-24)

### Main montage: 31 channels (current config)
`exclude_fp2 = false` → Fp2 is KEPT as a matter of static montage policy: wPLI matrices of
**31 × 31** with **465 edges per band** *when no channel is flagged by dynamic QC*.
`n_channels_used = 31` so Fp2 also passes the ±70 µV artifact rejection.
For sensitivity analysis without Fp2: `exclude_fp2 = true` + `n_channels_used = 30`. See README §9.

**Important — dynamic QC exclusion is independent of `exclude_fp2`.** Per-subject
channels flagged by the z-score QC check (`bad_ch`, threshold 3.0, see below) are
**always** merged into the exclusion list, regardless of `exclude_fp2`
(`src/SingleSubjectPipeline.jl:975-990`, `all_excl = unique(vcat(montage_excl, bad_ch))`).
So `exclude_fp2 = false` only guarantees Fp2 is not excluded *by static config* — if Fp2
(or any other channel) has z > 3.0 for a given subject, it is still dropped from the
analysis montage for that subject specifically. The real channel count therefore varies
per subject (e.g. sub-M05/ses-T2/EC: Fp2 has z=4.63 → 30 channels / 435 edges per band
analyzed, not 31/465) and downstream group analyses that hard-intersect channels across
the cohort can end up with noticeably fewer than 31 (confirmed: 24/31 in the transversal
EC cohort, ~59% of possible edges — see `docs/audits/batch_and_group_analysis_audit_2026-07-25.md`
§2.3/§5.3). The config field `n_channels_analysis` in `config/pipeline.toml` is
documentation-only and not read by the pipeline; the real per-subject count is reported
in `qc_decision_table.csv`'s `n_channels_analysis` column.

### Amplitude warning (`amplitude_warning`)
If mean raw signal σ̄ > 20 µV → probable recording without online filter active.
**Does not trigger automatic exclusion.** Final decision made after AR ±70 µV.

### QC decision table (`results/qc/qc_decision_table.csv`)

| `final_decision` | Condition | Include in group? |
|------------------|-----------|-------------------|
| `include` | No alerts, ≥10 valid epochs | ✅ Yes |
| `include_with_warning` | amplitude_warning or ≥1 bad channel | ✅ With caution |
| `manual_review` | amplitude_warning + ≥2 bad channels **or** <50% valid epochs | ⚠ Review |
| `exclude` | 0 valid epochs | ❌ No |

### wPLI reliability per band (minimum cycles per epoch)
With epochs of 1.0 s (`profile = eeg_julia`) and `min_cycles_for_wpli = 4.0`:
- DELTA (0.5 Hz × 1s = 0.5 cycles) → below threshold, `@warn`
- THETA (4.0 Hz × 1s = 4.0 cycles) → at threshold
- Option: `exclude_unreliable_bands = true` in `[connectivity]` to skip DELTA entirely

---

## 8. Strict Git rules

```bash
# Before ANY commit, verify no data files:
git ls-files | grep -E '(^data/|^results/|\.DS_Store$|^\.claude/|^\.vscode/)'
# → MUST return nothing
```

**Never commit:**
- `data/` — EEG signals, BIDS raw
- `results/` — per-subject derivatives
- `reports/`, `exports/`, `logs/`
- `.claude/`, `.vscode/`, `.cursor/`
- `.DS_Store`, `.env`, clinical logs, secrets

**Branch policy:** never commit directly to `main`. Use `feat/<name>` branches.

---

## 9. Common commands

```bash
# Install dependencies (once after cloning)
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Syntax check
julia --project=. -e 'include("src/NeuroMIND.jl"); println("OK")'

# Unit tests
julia --project=. test/runtests.jl
julia --project=. -e 'using Pkg; Pkg.test()'

# Single-subject pipeline
julia --project=. scripts/run_single_subject.jl

# Batch pipeline (all subjects)
julia --project=. scripts/run_batch_pipeline.jl --dry-run    # preview
julia --project=. scripts/run_batch_pipeline.jl --skip-done  # resume interrupted run

# Group analyses — estadística + TODAS las figuras (exploratorias + manuscrito EC +
# síntesis si el otro dominio ya tiene resultados) en una única pasada, siempre
julia --project=. scripts/run_transversal_analysis.jl
julia --project=. scripts/run_longitudinal_analysis.jl
# → results/{transversal,longitudinal}/{eyesclosed,eyesopen}/{tables,figures}/
# → results/summary/ (figura de síntesis, generada por el segundo análisis que se lance)
# Sync manual a Report_Pre/figures/plots/ cuando se quiera actualizar el informe:
#   cp results/{transversal,longitudinal}/eyesclosed/figures/*.png results/summary/*.png \
#      Report_Pre/figures/plots/

# Cohort interactive viewers (standalone; not Genie)
julia --project=. src/interactive/plot_transversal.jl      # → http://127.0.0.1:8781/
julia --project=. src/interactive/plot_longitudinal.jl    # → http://127.0.0.1:8780/

# Dashboard
julia --project=. scripts/launch_dashboard.jl
# → http://localhost:8080

# MNE-Python validation
cd mne_brain && python3 scripts/run_phase3_m05.py && python3 scripts/run_phase4_m05.py --apply-suggestions
```

---

## 10. Project status (2026-07-28)

### Dataset: MINDEM-IMIBIC
- **41 MS patients** (M4–M44) + **37 controls** (MC1–MC40) = 78 subjects
- **212 recordings** (.vhdr BrainVision), 206 valid for pipeline
- **27 complete longitudinal pairs** (T1+T2, EC+EO), plus M07 with a valid EC-only
  pair. Selection is now condition-specific: 28 EC / 27 EO candidates, of which
  16 (EC) / 18 (EO) pass per-session QC and are included in `results/longitudinal/`.
- Phases A+B complete; **Phase C (batch) executed 2026-07-25** — 201 OK, 5 SKIP, 0 ERR (~51 min, surrogates OFF)

> **Audit correction (2026-07-25).** The real roster in `inventory.csv` has 77 unique
> subjects (41 MS + 36 controls, not 37 — MC5/6/11/12 are missing from the declared
> range), and only 75 have data in `results/subjects/` (M4/M6 are ODDBALL recordings,
> correctly excluded as non-resting-state). An uncurated duplicate recording (M16 T1 EC,
> two `.vhdr` files, neither flagged `excluded`) explains the gap between 206 nominally
> "valid" recordings and 205 real rows in `qc_decision_table.csv`. Full detail:
> `docs/audits/batch_and_group_analysis_audit_2026-07-25.md` §3.3.

### Implemented and working
- [x] Full 8-step pipeline with custom ICA + surrogates
- [x] **wPLI multi-method** (Hilbert / FourierCSD / Multitaper) — configurable from TOML
- [x] Debiased dwPLI available in all 3 methods
- [x] Surrogates use the same estimator as the observed run
- [x] Dashboard panels 0–15 fully implemented
- [x] Panel 15: NeuroMIND-Julia vs mne_brain-MNE cross-validation, 3 stages, band power comparison
- [x] Native BrainVisionLoader (IEEE_FLOAT_32 / INT_16)
- [x] Batch runner with CLI filters, CSV log, skip-done
- [x] Transversal and longitudinal analyses with FDR-BH (full cohort, config vigente)
- [x] Cohort viz: `Transversal.jl`/`Longitudinal.jl` (stats+figures unified) + interactive viewers in `src/interactive/` (`:8780` / `:8781`)
- [x] `erf` removed from all scripts (replaced with A&S polynomial approximation 7.1.26)
- [x] `longitudinal_pairs.csv` format fix in `run_longitudinal_analysis.jl`
- [x] Unit tests (241+ passing)
- [x] PSD normalization fixed: `Pseg ./= (n_samp * fs)`

### Active branches
- Work locally; merge policy: never commit directly to `main` — use `feat/<name>`

### Pending
- [ ] Push/merge outstanding feature branches to main
- [ ] Integration tests for full pipeline
- [ ] GitHub Actions CI (syntax check + tests)
- [ ] Optional: channel-intersection policy for group analyses (currently hard intersect; viewers annotate N channels)

---

## 11. Julia dependencies (`Project.toml`)

```toml
Base64, CSV, CairoMakie, DSP, DataFrames, Dates, FFTW,
Genie, LinearAlgebra, Printf, Random, Serialization, Statistics,
StatsBase, TOML
```

**Not present:** `MultivariateStats`, `Pluto`, `JSON3`, `SpecialFunctions`.
`DSP.dpss` and `DSP.dpsseig` are used for multitaper estimation (DSP already in deps).

---

## 12. Rules for all AI agents

```
❌ Do not push data/, results/, reports/, logs/ to repo
❌ Do not push .claude/, .cursor/, .vscode/, .DS_Store
❌ Do not commit directly to main
❌ Do not change scientific logic without cross-checking with EEG_Julia/
❌ Do not use MultivariateStats or JSON3
❌ Do not use erf() directly — use the A&S polynomial approximation defined in each script
✅ Run syntax check (julia --project=. -e 'include("src/NeuroMIND.jl")') before committing
✅ Run tests before any merge to main
✅ Read only the relevant section of large files (dashboard.html is ~13 500 lines)
✅ Interactive viewers (cohort + single-subject aux) live in `src/interactive/`, never in `scripts/`
✅ Stats + figures for transversal/longitudinal live together in `Transversal.jl`/`Longitudinal.jl`
   (one module per domain, generated in the same pass — no separate figures-only script)
```

---

## 13. Changelog summary

### 2026-07-28 (cont.) — Visor transversal con estadística única de producción

- Contrato transversal v2 obligatorio en `statistics_contract.toml`, con
  procedencia, ámbitos y tamaños FDR, bootstrap, cuantiles, Mann–Whitney,
  `effect_rrb` y `effect_d_pooled`; CSV/JSON guardan precisión completa y el
  redondeo queda exclusivamente en la presentación.
- `plot_transversal.jl` valida contrato, columnas, familias y artefactos antes
  de cargar. Los resultados antiguos o incompletos se muestran como
  incompatibles y deben regenerarse; no existen fallbacks estadísticos ni de
  matrices en el navegador.
- El productor persiste medias/IC, medianas/cuartiles, diferencia/IC, r_rb,
  probabilidad de superioridad y d pooled/IC. El visor separa aristas
  evaluadas, nominales, Top-N y FDR, y consume sus flags/rangos persistidos.
- Tras auditar consumidores longitudinales y transversales, se retiraron de
  `viewer_common.js` los estimadores estadísticos compartidos que ya no tenían
  llamadas. Potencia conserva BH entre canales dentro de cada banda y muestra
  cobertura por canal.

### 2026-07-28 (cont.) — Visor longitudinal con estadística única de producción

- Contrato longitudinal v2 obligatorio en `statistics_contract.toml`: registra
  fuente estadística, ámbitos/tamaños FDR, bootstrap, cuantiles, `effect_rrb` y
  `effect_dz`. Los artefactos se escriben con precisión completa y nombres
  canónicos (`effect_dz`, `effect_dz_ci_*`, `median_diff`, `q1_diff`, `q3_diff`);
  el redondeo queda exclusivamente en la presentación.
- `plot_longitudinal.jl` valida contrato, columnas y artefactos antes de cargar.
  Resultados antiguos/incompletos se muestran como incompatibles y deben
  regenerarse: no se reconstruyen matrices ni estadísticos como fallback.
- `Longitudinal.jl::paired_change_summary` centraliza cuantiles, `r_rb` pareado,
  Cohen dz e IC95% bootstrap. `network_global_statistics.csv` materializa todos
  esos valores, su expresión mean wPLI normalizada y la familia BH de 7 bandas;
  el visor longitudinal ya no recalcula estadísticos en JavaScript.
- Mean wPLI se presenta como normalización exacta de mean strength
  (`mean_strength=(n_channels−1)×mean_wPLI`), no como estimando independiente.
  ALPHA-EC queda fijado por regresión: `r_rb≈−0.1618`, mediana
  `≈−0.0051785`, 351 aristas, 20 nominales, Top-20 con solapamiento 17 y 0 FDR.
- La UI separa recuentos total/nominal/Top-N/FDR, explicita las tres familias
  de corrección, usa −log10(p) con color FDR en el volcano y representa en el
  topograma de potencia la cobertura por canal (BH entre 31 canales/banda).

### 2026-07-28 (cont.) — Suite de tests consolidada

- Eliminada la duplicidad `test/` + `tests/`: la suite completa vive ahora en
  la ubicación canónica Julia `test/runtests.jl`; ya no existe un wrapper que
  delegue a una segunda carpeta.
- `Project.toml` declara `Test` mediante `[extras]` y `[targets]`, por lo que
  tanto `Pkg.test()` como la ejecución directa usan el mismo archivo sin
  modificar manualmente `LOAD_PATH`. `Printf`, usado por los módulos de cohorte,
  queda además declarado explícitamente en `[deps]` (antes lo ocultaba ese
  `LOAD_PATH` global).
- Encabezado de la suite rehecho con el formato documental del proyecto,
  comandos de ejecución e inventario de los 36 `@testset` con su finalidad.

### 2026-07-28 — Auditoría y endurecimiento longitudinal

- Selección T1–T2 independiente por condición: un sujeto ya no queda fuera de
  EC por carecer de un par EO. La cohorte real pasa a 16 pares EC y conserva
  18 EO; M07 queda correctamente recuperado en EC.
- Wilcoxon signed-rank exacto condicional mediante programación dinámica para
  N≤30, con fallback asintótico documentado y método registrado en CSV, JSON y
  snapshot. La arista marginal que habría superado FDR al aplicar el exacto
  sobre el antiguo N=15 deja de hacerlo al recuperar M07: 0 aristas y 0 celdas
  de potencia superan FDR en EC/EO.
- `subject_band_means.csv` y el gráfico pareado usan ahora el mismo montaje
  común que el estimando C (`mean_wpli_equiv`), con 27 canales EC / 28 EO.
  `n_bands_ok` conserva el recuento real también para excluidos, y las tablas
  añaden los aliases explícitos `n_edges` y `mean_abs_dz`.
- Nueva figura `power_effect_heatmap_channel_band.png`, canales en orden
  anatómico y marcas naranjas para N<70 % del máximo de la condición; los
  topomapas incorporan la misma advertencia. Fp2 queda marcado con N=8/16 EC
  y N=4/18 EO.
- Forest, cambios individuales, matrices ALPHA y gráfico pareado ajustados para
  evitar recortes y hacer explícitos N, montaje común, prueba exacta y leyenda
  T1/T2. Suite completa y ejecución end-to-end verificadas.

### 2026-07-27 (cont.) — Figuras de síntesis nuevas en Transversal.jl

- `fig_raincloud_band`: generaliza el raincloud de wPLI-por-sujeto (antes solo
  ALPHA-EC de manuscrito, vía `fig_raincloud_alpha` que se conserva intacta)
  a las 7 bandas, en ambas condiciones — `figures/raincloud_{BAND}.png`.
- `fig_power_effect_heatmap`: heatmap canal×banda de Cohen d para potencia
  (EM−Control), con overlay de celdas FDR-significativas — resume los 14
  topomapas Δ-potencia individuales en una sola figura por condición —
  `figures/power_effect_heatmap_channel_band.png`.
- `fig_interaction_ec_eo`: interacción grupo×condición (EC−EO). Lee
  `subject_band_means.csv` de `eyesclosed/` y `eyesopen/` (ya escritos por el
  bucle principal — mismo principio CSV→figura del resto del módulo), los
  empareja por sujeto+banda y compara Δ(EC−EO) entre EM y Control
  (Mann–Whitney + FDR-BH entre las 7 bandas + bootstrap CI). Como no
  pertenece a ninguna condición sola, sale a una tercera carpeta
  `results/transversal/combined/{tables,figures}/` (forest del diff-of-diff +
  slope plots por sujeto), generada una vez tras procesar ambas condiciones.
- Manifiesto de trazabilidad (`figures_manifest.tsv`) ahora se escribe para
  **ambas** condiciones (antes solo `eyesclosed`, restringido a las 4
  figuras de manuscrito) — la variable `manifest` se inicializa al principio
  de cada iteración del bucle `for cond in [...]` en vez de dentro del
  bloque `if cond == "eyesclosed"`.
- Verificado contra datos reales: paridad numérica exacta con el estado
  previo (mismos p/q/d por banda en EC/EO); interacción ALPHA
  (ΔEC−EO=0.106 EM vs 0.063 Control) coincide con el hallazgo descriptivo
  que motivó la petición; ninguna banda sobrevive FDR en la interacción
  (q mínimo 0.786) — resultado honesto, no forzado.

### 2026-07-27 — Simplificación radical de la capa de cohorte (auditoría + rediseño)

Reestructuración en dos pasadas el mismo día: una relocalización inicial
preservando comportamiento (`src/viz/` → `src/visualization/`, flags
`--figures-only`/`--manuscript`), seguida de una auditoría completa
(evidencia por grep + 2 subagentes de exploración) y un rediseño desde
principios tras 3 ideas del usuario. Solo el estado final se documenta aquí
en detalle — ver historial de git para los pasos intermedios.

- **Estadística + figuras unificadas en un único módulo por dominio**,
  generadas siempre juntas en cada lanzamiento (sin modo `--figures-only`
  separado): `src/transversal/Transversal.jl` y `src/longitudinal/Longitudinal.jl`
  sustituyen `run_{transversal,longitudinal}_analysis.jl` (lógica) +
  `{Transversal,Longitudinal}Figures.jl` + `{Transversal,Longitudinal}Manuscript.jl`
  + `src/visualization/{GroupVizCommon,PublicationCommon,PublicationTheme,SummaryFigures}.jl`
  (8 archivos, ~1865 líneas con ~40 % de duplicación real, sustituidos por 2
  módulos de dominio de ~1250 líneas cada uno). `scripts/run_{transversal,longitudinal}_analysis.jl`
  quedan como lanzadores finos (`include` + `Transversal.run(config)`).
- **Duplicación aceptada a propósito entre los dos módulos de dominio**
  (primitivas de dibujo + constantes `CH_POS`/`BAND_ORDER`, ~150 líneas): es
  el precio de tener exactamente 2 archivos en vez de 2 + un módulo
  compartido — decisión explícita del usuario, documentada en el propio
  código como "nota de diseño".
- **`results/{transversal,longitudinal}/EC|EO/` → `{eyesclosed,eyesopen}/{tables,figures}/`**,
  alineado con la convención que ya usa `results/subjects/` (vía `norm_cond()`).
  Todos los CSV bajo `tables/` (antes sueltos en la raíz + `tables/spectral/`
  + `tables/network/` — ahora plano); todas las PNG bajo `figures/`
  (exploratorias + antiguas "manuscript", ya sin distinción de carpeta).
  Tocó también `src/webapp/App.jl` (paneles 13/14, más superficie de la
  esperada) y los dos visores interactivos.
- **`results/publication/` eliminada por completo** (Idea 3): sin PDF (solo
  PNG); las figuras de manuscrito viven en `results/{dominio}/eyesclosed/figures/`
  junto a las exploratorias; la única figura cross-dominio
  (`summary_transversal_longitudinal_effects`) vive en `results/summary/`
  (carpeta nueva y mínima), generada por el segundo de los dos análisis que
  se lance, comprobando en disco si el otro ya tiene resultados.
- **`reconstruct_longitudinal_C!` eliminada de raíz** (no solo relocalizada):
  `mean_strength` por sujeto×banda×T1/T2 se captura una única vez, dentro del
  propio bucle que ya calcula `network_global_statistics.csv`
  (`tables/mean_strength_scores.csv`) — ya no hay una segunda pasada que
  recalcule desde wPLI crudo ni asserts para verificar que dos caminos
  independientes coinciden, porque ahora solo hay un camino.
- **Visores interactivos reubicados a `src/interactive/`** (Idea 1): no son
  lanzadores, son una tercera categoría de artefacto (HTTP+JS en vivo) —
  `scripts/aux_viewers/` → `src/interactive/aux/`, `plot_{transversal,longitudinal}.jl`
  + `group_viewer_common.js` (→ `viewer_common.js`) desde `src/{transversal,longitudinal}/`
  y `src/visualization/`. Nuevo `src/interactive/viewer_support.jl` (ex
  `GroupVizCommon.jl`, sección IO/JSON) — detectado durante la migración que
  los visores dependían de ese archivo antes de que se borrara; corregido
  antes de completar el borrado.
- **Código muerto podado** (single-subject, confirmado por grep exhaustivo,
  0 llamadas): `src/visualization/ClinicalPlots.jl` (archivo entero),
  `Heatmaps.jl::plot_group_comparison`, `Spectra.jl::plot_spectrum`,
  `GraphPlots.jl::plot_graph_metrics`, `GraphPlots.jl::plot_surrogate_distribution`
  (duplicada por `SingleSubjectPipeline.jl::_plot_surrogate_null_band`, que sí
  se usa — se conserva esa). `src/report/HTMLReport.jl` (335 líneas, huérfano,
  0 llamadas, sin relación con `Report_Pre/`) y `src/dashboard/Dashboard.jl`
  (stub de 6 líneas nunca incluido) eliminados enteros.
- **Duplicación estadística de `_assign_ranks` (4 copias: `GroupStats.jl`
  testeada-pero-no-usada, `run_transversal_analysis.jl`, `run_longitudinal_analysis.jl`,
  `PublicationCommon.jl`) reducida a 2** (una por módulo de dominio) — no
  eliminada del todo a propósito: consolidar con `GroupStats.jl`/`FDR.jl`
  queda fuera de esta ronda (toca lógica estadística, `AGENTS.md` §12 pide no
  tocarla sin verificar contra `EEG_Julia/`).
- Limpieza de 7 PNG sin uso en `Report_Pre/figures/plots/` (6 exploratorios +
  1 QC legacy, confirmados por grep de `\includegraphics` en cada capítulo).
  Los 8 PNG de manuscrito que sí se usan (verificado 8/8 en LaTeX) quedan sin
  tocar pero **desactualizados** — la sincronización con `Report_Pre/` ya no
  es automática (`copy_to_report_plots` eliminado, Idea 3): es un paso manual
  (`cp`/`rsync` desde `results/{transversal,longitudinal}/eyesclosed/figures/`
  y `results/summary/` a `Report_Pre/figures/plots/`).
- Verificado extremo a extremo contra los datos reales del dataset (no solo
  sintéticos): ambos lanzadores, paridad numérica exacta con el estado previo
  a la sesión (ALPHA transversal ctrl=0.196/MS=0.260, 12 sig; longitudinal EC
  15 pares, asserts de reconstrucción implícitos por construcción), ambos
  visores cargando desde las nuevas rutas.

### 2026-07-25 (cont.) — Manuscript publication figures (EC)

- `scripts/export_publication_figures.jl` + `src/viz/publication/`: eight
  static manuscript figures (PNG+PDF) from existing transversal/longitudinal
  CSVs — forest plots, raincloud ALPHA, connectivity matrices with FDR,
  network+power, longitudinal individual Δ (estimand C), and a
  transversal–longitudinal summary. Bootstrap CI (B=5000, seed fixed);
  estimand A transversal / C longitudinal; reconstructs per-subject
  `mean_strength` with assert `ΔwPLI = Δstrength/(n_ch−1)`. Integrated in
  `Report_Pre/chapters/{12,13,10}_*.tex`.

### 2026-07-25 (cont.) — Audit, Report_Pre group chapters, aux viewers versioned

- Full audit of the Phase C batch and both group analyses:
  `docs/audits/batch_and_group_analysis_audit_2026-07-25.md` — channel hard-intersect
  quantified (24/31 transversal EC, 27–28/31 longitudinal), a longitudinal
  candidate-filter bug (`audit_full_dataset.jl:411` requires EC *and* EO complete before
  per-condition QC), `GroupStats.jl`/`FDR.jl` (the only tested statistics module) not
  actually used by `run_transversal_analysis.jl`/`run_longitudinal_analysis.jl`, an
  uncurated duplicate recording (M16 T1 EC), and cohort-count mismatches vs this file.
- Applied the audit's low-risk fixes: this file's §7 clarified (dynamic QC exclusion is
  independent of `exclude_fp2`), A&S citation corrected (7.1.26, not 26.2.17), `load_wpli`
  wrapped in try/catch in both group scripts, EO cohort figures regenerated
  (`regenerate_group_figures.jl EO` — previously only EC had been run).
- `figures/aux/` (12 single-subject viewers, previously inside `results/`, gitignored and
  at risk of being lost on regeneration) versioned to
  `scripts/aux_viewers/sub-M05_ses-T2_eyesclosed/`, with each script's path constants
  rewritten to resolve `results/` from the new location.
- `Report_Pre/chapters/12_transversal.tex` and `13_longitudinal.tex` added — the report's
  first group-level chapters (previous chapters covered only sub-M05).
  `10_conclusiones.tex` updated: dropped now-obsolete limitations ("batch not launched",
  "no group contrast yet citable") and added ones the audit surfaced (channel-intersection
  variability, the longitudinal filter bug, missing clinical metadata).

### 2026-07-25 — Phase C + cohort visualization overhaul
- Batch: 201 OK / 5 SKIP / 0 ERR (~51 min, surrogates OFF); transversal + longitudinal regenerated.
- `src/viz/GroupVizCommon.jl` + `group_viewer_common.js`: shared PNG/JS for cohort plots.
- Viewers (`plot_longitudinal.jl` :8780, `plot_transversal.jl` :8781): full `edge_stats` tooltips, shared color scales, FDR banner, volcano, SEM on strip/spaghetti, `Reverse(:RdBu)` aligned with UI.
- Analysis scripts write `heatmap_triplet_*`, topo `sig_network_*`, and `explore_network_topN_*` when FDR is empty; honest `best_band` (empty if `n_total_sig=0`).
- `scripts/regenerate_group_figures.jl` regenerates PNGs from CSV without re-running stats.

### 2026-05-26 — wPLI multi-method + script fixes
- `src/connectivity/wPLI.jl`: multi-method architecture (Hilbert/FourierCSD/Multitaper) using abstract types + multiple dispatch. Surrogates use the same estimator.
- `config/batch_pipeline.toml` and `config/single_subject.toml`: new `[connectivity]` block with `wpli_method`, `[connectivity.fourier_csd]`, `[connectivity.multitaper]`.
- `scripts/run_transversal_analysis.jl` and `run_longitudinal_analysis.jl`: fix `erf` → A&S polynomial approximation 7.1.26.
- `scripts/run_longitudinal_analysis.jl`: same `erf` fix.
- `scripts/run_longitudinal_analysis.jl`: fix `longitudinal_pairs.csv` reader (new boolean flag format).
- Segmentation: `profile = "default"`, `segment_length_seconds = 2.0` in batch config.

### 2026-05-25 — dwPLI, valid surrogates, GraphMetrics, QC v2
- Debiased dwPLI (Vinck 2011) added to `wPLI.jl`.
- `_circular_shift_surrogate` — methodologically correct surrogate (independent shift per channel).
- Monte Carlo correction (+1): p-values never 0.
- `qc_decision_table.csv` — automatic QC decision table.
- Fp2 exclusion across all subjects.

### 2026-05-24 — Full dataset, BrainVision, Batch pipeline
- `src/io/BrainVisionLoader.jl`: native reader IEEE_FLOAT_32 / INT_16.
- `scripts/audit_full_dataset.jl`: 212 vhdr audited, 0 errors.
- `scripts/build_bids_full.jl`: 205 metadata JSON created.
- `scripts/run_batch_pipeline.jl`: batch pipeline with CLI filters and CSV log.

### 2026-05-23 — Full dashboard (panels 0–14), PSD fix
- Panels 11–14 implemented (Results, Export, Transversal, Longitudinal).
- PSD normalization fixed: `Pseg ./= (n_samp * fs)`.
- Panel 15: MNE-Python validation — 3 stages, band power comparison, 6-metric per-channel table.
