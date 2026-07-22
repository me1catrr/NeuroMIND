# AGENTS.md — NeuroMIND

> **Single source of truth for AI agents.**
> Read automatically by Claude Code, Cursor, OpenAI Codex CLI, and any other AI assistant
> on project open. **Keep this file up to date — it is the authoritative project context.**
> Last updated: 2026-07-21.
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
| `src/webapp/App.jl` | Genie server + all API routes (`/api/phase*`) |
| `web/views/dashboard.html` | Full dashboard SPA (~13 500 lines; inline HTML+CSS+JS) |
| `src/ica/ICACore.jl` | Pure-Julia FastICA (PCA whitening + tanh); profiles `eeg_julia`/`default` |
| `src/ica/ICAClassification.jl` | `compute_ica_features` (7 features) + `evaluate_ica_components` |
| `src/preprocessing/Filtering.jl` | `filter_recording` (HP/LP/Notch/Bandreject, Butterworth) |
| `src/segmentation/Epochs.jl` | `segment_recording`, `apply_baseline`, `reject_artifacts` |
| `src/spectral/PowerSpectrum.jl` | `compute_psd`, `plot_spectrum_grid` |
| `src/connectivity/wPLI.jl` | `compute_wpli` — 3 estimators: Hilbert, FourierCSD, Multitaper |
| `src/connectivity/CSD.jl` | `apply_csd` (optional, `use_csd = false` by default) |
| `src/io/BrainVisionLoader.jl` | Native BrainVision binary reader (.vhdr + .eeg, IEEE_FLOAT_32) |
| `config/pipeline.toml` | **Single config** for all 7 active scripts (see README §5) |
| `scripts/run_batch_pipeline.jl` | **Phase C** — batch pipeline with CLI filters |
| `scripts/run_transversal_analysis.jl` | Group analysis MS vs controls |
| `scripts/run_longitudinal_analysis.jl` | Longitudinal analysis T1 → T2 |
| `mne_brain/` | MNE-Python cross-validation pipeline |
| `tests/runtests.jl` | Unit test suite |

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
`results/transversal/` and `results/longitudinal/` (siblings of subjects/). See README §7–8.

---

## 7. Signal quality policy (QC v2 — 2026-05-24)

### Main montage: 31 channels (current config)
`exclude_fp2 = false` → Fp2 is KEPT: wPLI matrices of **31 × 31** with **465 edges per band**.
`n_channels_used = 31` so Fp2 also passes the ±70 µV artifact rejection.
For sensitivity analysis without Fp2: `exclude_fp2 = true` + `n_channels_used = 30`. See README §9.

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
julia --project=. tests/runtests.jl

# Single-subject pipeline
julia --project=. scripts/run_single_subject.jl

# Batch pipeline (all subjects)
julia --project=. scripts/run_batch_pipeline.jl --dry-run    # preview
julia --project=. scripts/run_batch_pipeline.jl --skip-done  # resume interrupted run

# Group analyses
julia --project=. scripts/run_transversal_analysis.jl
julia --project=. scripts/run_longitudinal_analysis.jl

# Dashboard
julia --project=. scripts/launch_dashboard.jl
# → http://localhost:8080

# MNE-Python validation
cd mne_brain && python3 scripts/run_phase3_m05.py && python3 scripts/run_phase4_m05.py --apply-suggestions
```

---

## 10. Project status (2026-05-26)

### Dataset: MINDEM-IMIBIC
- **41 MS patients** (M4–M44) + **37 controls** (MC1–MC40) = 78 subjects
- **212 recordings** (.vhdr BrainVision), 206 valid for pipeline
- **27 complete longitudinal pairs** (T1+T2, EC+EO)
- Phases A+B complete; Phase C (batch) ready to execute

### Implemented and working
- [x] Full 8-step pipeline with custom ICA + surrogates
- [x] **wPLI multi-method** (Hilbert / FourierCSD / Multitaper) — configurable from TOML
- [x] Debiased dwPLI available in all 3 methods
- [x] Surrogates use the same estimator as the observed run
- [x] Dashboard panels 0–15 fully implemented
- [x] Panel 15: NeuroMIND-Julia vs mne_brain-MNE cross-validation, 3 stages, band power comparison
- [x] Native BrainVisionLoader (IEEE_FLOAT_32 / INT_16)
- [x] Batch runner with CLI filters, CSV log, skip-done
- [x] Transversal and longitudinal analyses with FDR-BH
- [x] `erf` removed from all scripts (replaced with A&S polynomial approximation 26.2.17)
- [x] `longitudinal_pairs.csv` format fix in `run_longitudinal_analysis.jl`
- [x] Unit tests (241+ passing)
- [x] PSD normalization fixed: `Pseg ./= (n_samp * fs)`

### Active branches
- `feat/phase10-surrogates` ← **current branch**
- Phases 5–9 have branches pending PR to main

### Pending
- [ ] Run Phase C full dataset (`run_batch_pipeline.jl`) — ~18–50 h
- [ ] Run transversal and longitudinal analyses with full cohort
- [ ] Push/merge feat/phase5 → feat/phase10 branches to main
- [ ] Integration tests for full pipeline
- [ ] GitHub Actions CI (syntax check + tests)

---

## 11. Julia dependencies (`Project.toml`)

```toml
Base64, CSV, CairoMakie, DSP, DataFrames, Dates, FFTW,
Genie, LinearAlgebra, Random, Serialization, Statistics,
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
```

---

## 13. Changelog summary

### 2026-05-26 — wPLI multi-method + script fixes
- `src/connectivity/wPLI.jl`: multi-method architecture (Hilbert/FourierCSD/Multitaper) using abstract types + multiple dispatch. Surrogates use the same estimator.
- `config/batch_pipeline.toml` and `config/single_subject.toml`: new `[connectivity]` block with `wpli_method`, `[connectivity.fourier_csd]`, `[connectivity.multitaper]`.
- `scripts/run_transversal_analysis.jl` and `run_longitudinal_analysis.jl`: fix `erf` → A&S polynomial approximation 26.2.17.
- `src/longitudinal/LongitudinalAnalysis.jl`: same `erf` fix.
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
