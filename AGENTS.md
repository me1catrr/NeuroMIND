# AGENTS.md — NeuroMIND

> **Single source of truth for AI agents.**
> Read automatically by Claude Code, Cursor, OpenAI Codex CLI, and any other AI assistant
> on project open. **Keep this file up to date — it is the authoritative project context.**
> Last updated: 2026-10-15.
>
> **README.md is the authoritative source for the pipeline, config/pipeline.toml,
> outputs and results structure.** This file keeps agent-specific rules (Git, dashboard)
> and cross-references the README to avoid drift. Full changelog: `docs/CHANGELOG.md`
> (§13 points there).

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
[1/8] Load          BrainVisionLoader (BIDS-EEG tree) → EEGRecording (channels × samples)
[2/8] QC            flag_bad_channels by z-score (threshold 3.0) + amplitude_warning
[3/8] Filter        HP 0.5 Hz + LP 150 Hz + Notch 50 Hz + Bandreject (Butterworth filtfilt)
[4/8] ICA           Symmetric FastICA on CONTINUOUS filtered signal → ICAResult
[5/8] Segment       segment_recording(rec_ICA) → EpochSet + baseline + AR rejection
[6/8] Spectral      Hamming-taper PSD + band power → SpectralResult
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
| `scripts/audit_full_dataset.jl` | **Phase A** — scans `data/sourcedata/brainvision/`, applies `config/raw_selection.csv`, merges clinical sources → `data/tables/inventory.csv` + `data/bids/participants.tsv` |
| `scripts/build_bids_full.jl` | **Phase B** — builds the real BIDS-EEG dataset in `data/bids/` (renamed triplets, crops for M04/M05/M06 T1, sidecars); templates in `config/bids/` |
| `config/raw_selection.csv` | Raw-selection decisions from Javier (`segment|use|discard|exclude`); single source for which raw file/window enters |
| `scripts/build_clinical_master.jl` | Documento maestro para Javier/Alejandro: une `participants.tsv` + `inventory.csv` + observaciones de Javier → `data/sourcedata/clinical/datos_clinicos_unificados.{xlsx,csv}` (valores EDSS/D-FIS/BDI-II, categorías, selección EEG y control de calidad). Regenerar tras `audit_full_dataset.jl` y `build_bids_full.jl`; con `--check` verifica además la consistencia de todo el conjunto y las cifras del cap. 4 |
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

# Dataset preparation (after adding recordings or updating the clinical sheets)
julia --project=. scripts/audit_full_dataset.jl      # → data/tables/inventory.csv + data/bids/participants.tsv
julia --project=. scripts/build_bids_full.jl         # → data/bids/sub-*/ses-*/eeg/ (BIDS-EEG)
julia --project=. scripts/build_clinical_master.jl --check   # → data/sourcedata/clinical/datos_clinicos_unificados.{xlsx,csv}
#   --check cruza BIDS en disco ↔ inventory ↔ raw_selection ↔ participants ↔ libro/CSV ↔ base original ↔ cifras
#   del cap. 4 del informe, imprime las filas de los Cuadros 4.1/4.2 y sale con error si algo no coincide
bids-validator-deno data/bids                        # must report 0 errors (uv tool install bids-validator-deno)

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

# Knowledge graph (graphify 0.9.79, versionado en graphify-out/) — consultar ANTES de leer código
graphify query "<pregunta>"          # subgrafo acotado (BFS)
graphify explain "GroupStats.jl"     # nodo + vecinos
graphify path "A" "B"                # camino más corto entre dos conceptos
graphify update .                    # re-extrae código (AST, sin LLM); el hook post-commit ya lo hace
open graphify-out/graph.html         # grafo interactivo
# Setup por equipo (una vez): pipx install graphifyy==0.9.79 (o uv tool install),
# aplicar graphify-out/patches/graphify_julia_returntype.patch (comando en su cabecera)
# y `graphify hook install`. Limitación conocida: no resuelve llamadas a funciones
# exportadas vía `using NeuroMIND` (p. ej. desde test/) — "sin aristas" ≠ código muerto.
```

---

## 10. Project status (2026-10-15)

### Dataset: MINDEM-IMIBIC
- **41 MS patients** (M4–M44) + **36 controls** = 77 subjects (`data/bids/participants.tsv`)
- **212 recordings** (.vhdr BrainVision) → **210 valid** recordings in `data/tables/inventory.csv` (105 EC + 105 EO; M4/M5/M6 T1 single raws are split into 2 each; M7 T1 EO is the file without eye suffix, confirmed by Alejandro 2026-09-28)
- **28 complete longitudinal pairs** (T1+T2, EC+EO; M07 recovered its T1 EO on
  2026-10-15). Selection is condition-specific: 28 EC / 28 EO candidates. The 16 (EC) /
  18 (EO) pairs that passed per-session QC in `results/longitudinal/` come from the
  July run and must be recomputed on the new BIDS.
- Phases A+B rebuilt 2026-10-15 (BIDS-EEG validated with bids-validator-deno 3.0.2: 0 errors, 4 documented warnings). **Phase C (batch) last executed 2026-07-25** — 201 OK, 5 SKIP, 0 ERR (~51 min, surrogates OFF) on the old light-BIDS layout; M04/M05/M06 T1 and the group analyses still need reprocessing

> **Audit correction (2026-07-25).** The real roster in `inventory.csv` has 77 unique
> subjects (41 MS + 36 controls, not 37 — MC5/6/11/12 are missing from the declared
> range), and only 75 have data in `results/subjects/` (M4/M6 were excluded as ODDBALL;
> **corrected 2026-10-14**: per Javier their single T1 raw holds eyes-open + eyes-closed
> and is split by seconds — see `config/raw_selection.csv`; not yet reprocessed). An
> uncurated duplicate recording (M16 T1 EC,
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
- Merge policy: never commit directly to `main` — use `feat/<name>`
- `feat/graphify-setup` → PR me1catrr/NeuroMIND#5 (21 commits: trabajo de julio + graphify)
- `feat/ci` → PR me1catrr/NeuroMIND#6, apilado sobre #5 (fusionar #5 primero)

### Pending
- [ ] Review/merge PR #5 and then PR #6 into main
- [ ] Reprocess M04/M05/M06 T1 (now cropped in `data/bids/`) and relaunch the batch on the new BIDS tree,
      then rerun transversal/longitudinal analyses; M05 T2 EC regression check vs the old TSV-based run
- [ ] Port the dashboard raw-signal/PSD routes (still read `*_eeg_data.tsv`, absent since the light BIDS) to the BrainVision reader
- [ ] Electrode coordinate system is undocumented in the source (declared `Other`); confirm origin with the cap supplier
- [ ] Age discrepancies between sources: M13 (34 vs 33), M29 (43 vs 42), M35 (41 vs 43 by dates) — confirm with the clinical team
- [ ] `data/_archive/` (pre-reform copies) can be deleted once the batch is reprocessed on the new BIDS
- [ ] Integration tests for full pipeline
- [ ] GitHub Actions CI: workflow added in PR #6 (`.github/workflows/ci.yml`); confirm first green run
- [ ] Decide GroupStats.jl/FDR.jl: production (Transversal/Longitudinal) does not call them —
      only test/runtests.jl does (confirmed 2026-10-07 with graphify + grep)
- [ ] Report_Pre/figures/plots/ out of sync with results/ since 2026-07-27 (manual sync)
- [ ] Report graphify Julia return-type bug upstream (local patch in graphify-out/patches/)
- [ ] Optional: channel-intersection policy for group analyses (currently hard intersect; viewers annotate N channels)

---

## 11. Julia dependencies (`Project.toml`)

Dependency list: see `Project.toml` (`[deps]`).

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
✅ data/ layout: `sourcedata/` is read-only; `bids/` and `tables/` are generated by
   audit_full_dataset.jl + build_bids_full.jl. After regenerating, run bids-validator (0 errors)
✅ Never put identifying clinical data outside `data/` (participants.tsv holds birth/inclusion dates)
✅ Read only the relevant section of large files (dashboard.html is ~13 500 lines)
✅ Interactive viewers (cohort + single-subject aux) live in `src/interactive/`, never in `scripts/`
✅ Stats + figures for transversal/longitudinal live together in `Transversal.jl`/`Longitudinal.jl`
   (one module per domain, generated in the same pass — no separate figures-only script)
```

---

## 13. Changelog summary

Historial completo, con fechas y justificación de cada cambio, en
[`docs/CHANGELOG.md`](docs/CHANGELOG.md) — léelo cuando necesites saber por qué algo
es como es o cuándo cambió. Al cerrar un cambio relevante, añade la entrada arriba
de ese archivo, no aquí. Hitos citados desde `Report_Pre/`: corrección de la
normalización PSD (`Pseg ./= (n_samp * fs)`, 2026-05-23).
