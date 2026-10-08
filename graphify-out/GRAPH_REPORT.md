# Graph Report - NeuroMIND  (2026-10-07)

## Corpus Check
- 98 files · ~189,682 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 9 file(s) not represented in the graph (top: (none) 6, .toml 3)

## Summary
- 1479 nodes · 3022 edges · 71 communities (62 shown, 9 thin omitted)
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 174 edges (avg confidence: 0.87)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `c8d8e1c2`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- compute_wpli
- run_full_pipeline.py
- Transversal.jl: estadística y figuras
- Longitudinal
- epochs.py
- run_phase8_validation.py
- Suite de tests Julia
- audit_full_dataset.jl
- load_config
- loader.py
- types.py
- Tipos centrales Julia (types.jl)
- compare_neuromind_mne_rejections.py
- Visor aux: surrogates
- plot_transversal.jl
- Visor aux: épocas
- Visor aux: componentes ICA
- Visor aux: baseline
- plot_longitudinal.jl
- mne_brain: filtrado y QC
- viewer_common.js (UI compartida)
- numpy
- SingleSubjectPipeline.jl
- Visor aux: conectividad
- Visor aux: espectro
- Visor aux: histograma raw
- plot_raw.jl
- Visor aux: ICA antes/después
- GroupStats.jl
- Visor aux: PSD espectral
- Visor aux: filtrado vs raw
- run_batch_pipeline.jl
- Visor aux: butterfly raw
- Módulo NeuroMIND.jl
- Estimadores wPLI (wPLI.jl)
- Dashboard Genie (App.jl)
- IncompatibleTransversalResultsError
- Reglas de agentes (AGENTS.md)
- Soporte de visores (ViewerSupport)
- Visión general (README)
- Lanzadores de análisis de grupo
- IncompatibleResultsError
- Config.jl (rutas)
- Pipeline de 8 pasos (doc)
- ICA FastICA (ICACore.jl)
- argparse
- Single-subject viewer flow raw->filter->ICA->epochs->spectral->connectivity->surrogate
- Clean-slate per leaf unit + manifest
- Longitudinal.jl
- ica_core.py
- Batch & group analysis audit 2026-07-25
- Plan de continuación — Report_Pre + revisión de código/resultados
- PowerSpectrum.jl
- Lanzador del dashboard
- Clasificación ICA (Julia)
- src/interactive README (viewer launch guide)
- Filtrado (Filtering.jl)
- Surrogates.jl
- GraphPlots.jl
- Topomaps.jl
- CSD (CSD.jl)
- Métricas de grafo (GraphMetrics.jl)
- BIDSLoader.jl
- Control de calidad (QualityControl.jl)
- Heatmaps.jl
- Spectra.jl
- mne_brain: pyproject
- Inspección ICA (ICAInspection.jl)
- Corrección FDR (FDR.jl)

## God Nodes (most connected - your core abstractions)
1. `Transversal` - 77 edges
2. `Longitudinal` - 74 edges
3. `load_config()` - 41 edges
4. `run()` - 35 edges
5. `run_single_subject_pipeline()` - 32 edges
6. `PipelineConfig` - 30 edges
7. `run()` - 28 edges
8. `EEGRecording` - 25 edges
9. `run_phase4()` - 19 edges
10. `make_epochs()` - 17 edges

## Surprising Connections (you probably didn't know these)
- `Dynamic per-subject QC channel exclusion (z>3)` --references--> `run_single_subject_pipeline()`  [INFERRED]
  AGENTS.md → src/SingleSubjectPipeline.jl
- `Circular-shift surrogate test` --references--> `compute_wpli()`  [INFERRED]
  README.md → src/connectivity/wPLI.jl
- `load_ss_config()` --shares_data_with--> `config/pipeline.toml (single config)`  [EXTRACTED]
  src/SingleSubjectPipeline.jl → README.md
- `8-step single-subject pipeline` --implements--> `run_single_subject_pipeline()`  [EXTRACTED]
  README.md → src/SingleSubjectPipeline.jl
- `paired_change_summary()` --implements--> `Bootstrap 95% CI (B=5000, fixed seed)`  [EXTRACTED]
  src/longitudinal/Longitudinal.jl → README.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **M05 traceability chain: validation prompt, July-14 audit, July-25 audit** — docs_prompt_claudecode_lanzar_m05, docs_audits_run_consistency_audit, docs_audits_batch_and_group_analysis_audit_2026_07_25, docs_prompt_claudecode_lanzar_m05_sub_m05_reference_run [EXTRACTED 1.00]
- **Single-subject 8-step pipeline flow** — src_singlesubjectpipeline_load_ss_config, src_io_brainvisionloader, src_preprocessing_filtering_filter_recording, src_ica_icacore, src_segmentation_epochs_segment_recording, src_spectral_powerspectrum_compute_psd, src_connectivity_wpli_compute_wpli, src_singlesubjectpipeline__save_all_results [EXTRACTED 1.00]
- **Proposed src/runtime/ module set around RunContext** — docs_tdd_ejecucion_rutinas_runtime_layer, docs_tdd_ejecucion_rutinas_runcontext, docs_tdd_ejecucion_rutinas_resultsmanager, docs_tdd_ejecucion_rutinas_runmanifest, docs_tdd_ejecucion_rutinas_runlogging, docs_tdd_ejecucion_rutinas_console, docs_tdd_ejecucion_rutinas_progress, docs_tdd_ejecucion_rutinas_preflight, docs_tdd_ejecucion_rutinas_summary [EXTRACTED 1.00]
- **Configurable wPLI estimators** — readme_hilbert_estimator, readme_fourier_csd_estimator, readme_multitaper_estimator, readme_dwpli [EXTRACTED 1.00]
- **Cohort statistics production → contract → viewer flow** — src_transversal_transversal, src_longitudinal_longitudinal, readme_statistics_contract_v2, src_interactive_plot_transversal, src_interactive_plot_longitudinal, agents_single_statistics_source [INFERRED 0.85]
- **Threats to group-analysis validity identified in Phase C audit** — docs_audits_batch_and_group_analysis_audit_2026_07_25_untested_group_statistics, docs_audits_batch_and_group_analysis_audit_2026_07_25_longitudinal_pairing_filter_bug, docs_audits_batch_and_group_analysis_audit_2026_07_25_hard_intersect_channel_policy, docs_audits_batch_and_group_analysis_audit_2026_07_25_dynamic_qc_channel_exclusion, docs_audits_batch_and_group_analysis_audit_2026_07_25_group_results_summary [INFERRED 0.85]

## Communities (71 total, 9 thin omitted)

### Community 0 - "compute_wpli"
Cohesion: 0.11
Nodes (14): main(), compute_wpli(), connectivity_matrices(), _plot_wpli_heatmap(), save_connectivity_results(), _make_coupled_epochs(), _make_epochs(), test_coupled_channels_have_higher_alpha_wpli() (+6 more)

### Community 1 - "run_full_pipeline.py"
Cohesion: 0.13
Nodes (23): append_log_row(), build_parser(), discover_recordings(), _elapsed(), _exclude_recording(), _exclusion_path(), filter_recordings(), is_batch_mode() (+15 more)

### Community 2 - "Transversal.jl: estadística y figuras"
Cohesion: 0.07
Nodes (77): anatomical_order(), apply_theme!(), _assign_ranks(), _band_row(), bh_qvalues(), bootstrap_cohen_d_ci(), bootstrap_mean_diff_ci(), _ch_xy() (+69 more)

### Community 3 - "Longitudinal"
Cohesion: 0.06
Nodes (74): Bootstrap 95% CI (B=5000, fixed seed), anatomical_order(), apply_theme!(), _as_bool(), _assign_ranks(), _band_scores(), bh_qvalues(), bootstrap_cohen_dz_ci() (+66 more)

### Community 4 - "epochs.py"
Cohesion: 0.09
Nodes (19): main(), apply_baseline(), epoch_summary(), load_cleaned_raw(), make_epochs(), _npz_to_raw(), reject_artifacts(), save_epoch_exclusion() (+11 more)

### Community 5 - "run_phase8_validation.py"
Cohesion: 0.06
Nodes (35): CLAUDE.md (@AGENTS.md include), mne_brain AGENTS.md, Batch errors 'Epochs-object is empty' (AR ±70 µV), mne_brain Phase 8 surrogates + FDR (pending), mne_brain pipeline_config.yaml, EEG frequency bands (ALPHA 7.8-11.7 Hz, 7 bands), surrogates config (200, phase_shuffle, BH), mne_brain README (+27 more)

### Community 6 - "Suite de tests Julia"
Cohesion: 0.15
Nodes (15): Test suite consolidated in test/runtests.jl, Test, _cfg_seg(), _cfg_with_ica(), _cfg_with_root(), _cfg_with_surrogates(), CSV, DataFrames (+7 more)

### Community 7 - "audit_full_dataset.jl"
Cohesion: 0.11
Nodes (21): FDR-BH applied per band (7 families), Group results: transversal EC 21 sig (ALPHA), EO 3, longitudinal 0, Longitudinal pairing filter-order bug, Dataset preparation (Phase A audit + Phase B BIDS), _bids_subject_id(), is_excluded(), Dates, load_demographics() (+13 more)

### Community 8 - "load_config"
Cohesion: 0.09
Nodes (17): default_config_path(), load_config(), compute_band_power(), compute_psd(), mean_spectrum(), save_spectral_results(), _make_alpha_epochs(), _make_epochs() (+9 more)

### Community 9 - "loader.py"
Cohesion: 0.11
Nodes (20): bids_root_dir(), BrainVisionHeader, find_vhdr_by_name(), _first_existing(), load_eeg_bids(), load_eeg_brainvision(), load_eeg_tsv(), load_electrode_positions() (+12 more)

### Community 10 - "types.py"
Cohesion: 0.08
Nodes (14): _normalize_bands(), ClinicalData, ConnectivityMatrix, EpochSet, GraphMetrics, GroupAnalysis, LongitudinalAnalysis, PipelineConfig (+6 more)

### Community 11 - "Tipos centrales Julia (types.jl)"
Cohesion: 0.16
Nodes (16): ClinicalData, ConnectivityMatrix, duration(), EEGRecording, EpochSet, GraphMetrics, GroupAnalysis, ICAResult (+8 more)

### Community 12 - "compare_neuromind_mne_rejections.py"
Cohesion: 0.20
Nodes (18): _as_int(), _as_pct(), build_parser(), compare_recordings(), _decision(), _diff(), _fmt_top(), load_julia_recordings() (+10 more)

### Community 13 - "Visor aux: surrogates"
Cohesion: 0.08
Nodes (47): _auto_interpret(), _band_json(), BandPack, _cache_get(), _ch_index(), _contrast_json(), _edge_row(), _extract_null_dist() (+39 more)

### Community 14 - "plot_transversal.jl"
Cohesion: 0.13
Nodes (23): _cond_dir(), handle_request(), html_page(), CairoMakie, CSV, DataFrames, Dates, Sockets (+15 more)

### Community 15 - "Visor aux: épocas"
Cohesion: 0.10
Nodes (32): _default_rejected_epoch(), _epoch_window(), EpochStore, handle_request(), html_page(), CairoMakie, CSV, DataFrames (+24 more)

### Community 16 - "Visor aux: componentes ICA"
Cohesion: 0.11
Nodes (30): _components_json(), compute_ic_psd(), _decision_explanation(), _feature_json(), handle_request(), html_page(), ICACompStore, CairoMakie (+22 more)

### Community 17 - "Visor aux: baseline"
Cohesion: 0.10
Nodes (30): _apply_bl(), _baseline_end_s(), BaselineStore, compute_offset_table(), _epoch_max_offset(), _epoch_raw(), handle_request(), html_page() (+22 more)

### Community 18 - "plot_longitudinal.jl"
Cohesion: 0.13
Nodes (24): _cond_dir(), handle_request(), _honest_best_from_summary(), html_page(), CairoMakie, CSV, DataFrames, Dates (+16 more)

### Community 19 - "mne_brain: filtrado y QC"
Cohesion: 0.11
Nodes (18): EEGRecording, apply_bandreject(), apply_highpass(), apply_lowpass(), apply_notch(), _apply_sos(), describe_filter_chain(), filter_recording() (+10 more)

### Community 20 - "viewer_common.js (UI compartida)"
Cohesion: 0.22
Nodes (15): barChart(), bindHeatmap(), drawColorbar(), drawGraph(), drawHeatmap(), drawStrength(), drawTopo(), drawVolcano() (+7 more)

### Community 21 - "numpy"
Cohesion: 0.13
Nodes (15): run_phase4(), main(), _save_matrix(), _save_topomaps(), _bandpower(), compute_ica_features(), evaluate_ica_components(), _kurtosis() (+7 more)

### Community 22 - "SingleSubjectPipeline.jl"
Cohesion: 0.13
Nodes (35): RunLogging (LoggingExtras TeeLogger), _bh_qvalues(), detect_first_subject(), _format_eta_s(), _ica_config_hash(), _ica_effective_params(), CairoMakie, Serialization (+27 more)

### Community 23 - "Visor aux: conectividad"
Cohesion: 0.09
Nodes (43): _band_strength(), _betweenness(), _binary_adj(), _clustering(), ConnStore, _degree_from_adj(), _density_label(), _edges_for_band() (+35 more)

### Community 24 - "Visor aux: espectro"
Cohesion: 0.12
Nodes (29): _band_payload(), _band_row(), _curve_series(), handle_request(), html_page(), _index_row(), CairoMakie, CSV (+21 more)

### Community 25 - "Visor aux: histograma raw"
Cohesion: 0.15
Nodes (25): _channel_stats(), _grid_dims(), handle_request(), _hist_json_channel(), _histogram(), html_page(), CairoMakie, CSV (+17 more)

### Community 26 - "plot_raw.jl"
Cohesion: 0.16
Nodes (16): handle_request(), html_page(), CairoMakie, CSV, DataFrames, Dates, Sockets, load_store() (+8 more)

### Community 27 - "Visor aux: ICA antes/después"
Cohesion: 0.14
Nodes (22): handle_request(), html_page(), ICASignalStore, CairoMakie, CSV, DataFrames, Dates, Sockets (+14 more)

### Community 28 - "GroupStats.jl"
Cohesion: 0.29
Nodes (13): Untested production group statistics (GroupStats.jl unused), Welch/paired t-tests use Z (normal) approximation, statistics section (mann_whitney, wilcoxon, fdr_q 0.05), _assign_ranks(), _collect_wpli(), compare_groups_stats(), _empty_stat(), _erf_approx() (+5 more)

### Community 29 - "Visor aux: PSD espectral"
Cohesion: 0.14
Nodes (22): _band_row(), handle_request(), html_page(), _index_row(), CairoMakie, CSV, DataFrames, Dates (+14 more)

### Community 30 - "Visor aux: filtrado vs raw"
Cohesion: 0.16
Nodes (18): handle_request(), html_page(), CairoMakie, CSV, DataFrames, Dates, Sockets, _load_stage_csv() (+10 more)

### Community 31 - "run_batch_pipeline.jl"
Cohesion: 0.09
Nodes (26): Phase C batch run 2026-07-25 (201 OK/5 SKIP/0 ERR), config/pipeline.toml (single config), Report_Pre citability precondition (reprocess with current config), Uncurated duplicate recording M16 T1 EC, Phase C batch run (201 OK / 5 SKIP / 0 ERR), Progress (ProgressMeter), Logging, 31-channel montage (exclude_fp2=false, 465 edges) (+18 more)

### Community 32 - "Visor aux: butterfly raw"
Cohesion: 0.16
Nodes (18): handle_request(), html_page(), CairoMakie, CSV, DataFrames, Dates, Sockets, Statistics (+10 more)

### Community 33 - "Módulo NeuroMIND.jl"
Cohesion: 0.14
Nodes (12): CSV, DataFrames, Dates, DSP, FFTW, LinearAlgebra, Random, Serialization (+4 more)

### Community 34 - "Estimadores wPLI (wPLI.jl)"
Cohesion: 0.07
Nodes (43): Condition-specific T1–T2 pair selection, MINDEM-IMIBIC dataset (41 MS + 36/37 controls), Hypothesis: MS alters alpha/beta connectivity, Multiple Sclerosis (MS), NeuroMIND EEG connectivity framework, Amplitude warning (σ̄ > 20 µV), Single BIDS output tree results/subjects/, config_snapshot.toml (per-recording provenance) (+35 more)

### Community 35 - "Dashboard Genie (App.jl)"
Cohesion: 0.18
Nodes (10): Base64, Genie, Genie.Renderer.Html, Genie.Renderer.Json, Genie.Requests, Genie.Router, Genie.Server, Dates (+2 more)

### Community 36 - "IncompatibleTransversalResultsError"
Cohesion: 0.26
Nodes (9): CondStore, IncompatibleTransversalResultsError, load_condition(), _require_columns(), _required_csv(), _required_matrix(), TransViewer, validate_condition_contract() (+1 more)

### Community 37 - "Reglas de agentes (AGENTS.md)"
Cohesion: 0.33
Nodes (3): EEG_Julia original project (reference), AGENTS.md — NeuroMIND agent guide, Current Source Density (CSD)

### Community 38 - "Soporte de visores (ViewerSupport)"
Cohesion: 0.13
Nodes (7): df_rows_json(), CairoMakie, CSV, DataFrames, TOML, json_escape(), ViewerSupport

### Community 39 - "Visión general (README)"
Cohesion: 0.47
Nodes (5): 12 single-subject aux viewers (:8765–:8775), Genie dashboard (panels 0–15), mne_brain MNE-Python cross-validation, README.md — NeuroMIND user guide, sub-M05/ses-T2/eyesclosed reference case

### Community 40 - "Lanzadores de análisis de grupo"
Cohesion: 0.29
Nodes (5): A&S 7.1.26 polynomial erf approximation, Longitudinal, TOML, TOML, Transversal

### Community 41 - "IncompatibleResultsError"
Cohesion: 0.26
Nodes (9): CondStore, IncompatibleResultsError, load_condition(), LongViewer, _require_columns(), _required_csv(), _required_matrix(), validate_condition_contract() (+1 more)

### Community 42 - "Config.jl (rutas)"
Cohesion: 0.27
Nodes (6): ensure_dirs(), _f(), load_subjects(), results_dir(), _s(), subject_results_dir()

### Community 43 - "Pipeline de 8 pasos (doc)"
Cohesion: 0.40
Nodes (3): 8-step single-subject pipeline, load_eeg_brainvision(), read_vhdr_header()

### Community 44 - "ICA FastICA (ICACore.jl)"
Cohesion: 0.53
Nodes (5): Independent Component Analysis (ICA), _fastica_attempt(), run_ica(), _sym_decorr(), _whiten_pca()

### Community 45 - "argparse"
Cohesion: 0.14
Nodes (6): main(), _load_ica_suggestions(), _resolve_dirs(), main(), _resolve_dirs(), _resolve_dirs()

### Community 46 - "Single-subject viewer flow raw->filter->ICA->epochs->spectral->connectivity->surrogate"
Cohesion: 0.15
Nodes (13): Dynamic QC channel exclusion (bad_ch always merged), Hard-intersect channel policy in group analyses, Prompt: Validar M05 tras la unificación, Lazy Genie loading (dashboard not compiled by pipeline), Single BIDS output tree (legacy results/{ID}/{SES} removed), sub-M05/ses-T2/eyesclosed reference run, M05 validation checklist (checks A-G), TDD — Sistema de ejecución de rutinas (+5 more)

### Community 47 - "Clean-slate per leaf unit + manifest"
Cohesion: 0.20
Nodes (10): Versioning of 12 single-subject aux viewers out of results/, Atomic staging directory ({unit}.tmp/ + mv swap), Console (Crayons, TTY/ASCII degradation), Preflight checks (fail early), ResultsManager, run_manifest.json (per execution unit), RunContext, RunManifest module (+2 more)

### Community 48 - "Longitudinal.jl"
Cohesion: 0.28
Nodes (5): Cohort layer radical simplification 2026-07-27, mean_strength = (n_channels−1) × mean_wPLI normalization, Manual sync of figures to Report_Pre/figures/plots, Cohort interactive viewers (:8780/:8781), statistics_contract.toml schema v2

### Community 49 - "ica_core.py"
Cohesion: 0.27
Nodes (5): ICAResult, _explained_variance(), MNEICAResult, recording_to_raw(), run_ica()

### Community 50 - "Batch & group analysis audit 2026-07-25"
Cohesion: 0.27
Nodes (10): Hard channel intersection across cohort, Dynamic per-subject QC channel exclusion (z>3), Batch & group analysis audit 2026-07-25, Traceability audit — M05 run (2026-07-14), channel_statistics_compare.csv without generator, ICA cache temporal mix (FastICA May vs pipeline July), ICA auto-rejection log vs load_ica_labels code mismatch, Missing git_commit / config_hash / run_id (+2 more)

### Community 51 - "Plan de continuación — Report_Pre + revisión de código/resultados"
Cohesion: 0.25
Nodes (7): 0. Estado al cerrar, 1. Cómo trabajar en Cursor, 2. Bucle por capítulo (Report_Pre ↔ código ↔ resultados), 3. Mapa capítulo → código → salidas, 4. Cola de trabajo (orden sugerido), Plan de continuación — Report_Pre + revisión de código/resultados, Prompt de arranque para Cursor

### Community 52 - "PowerSpectrum.jl"
Cohesion: 0.60
Nodes (3): _band_power_from_psd(), compute_psd(), _hamming_taper()

### Community 53 - "Lanzador del dashboard"
Cohesion: 0.50
Nodes (3): NeuroMIND, Pkg, TOML

### Community 54 - "Clasificación ICA (Julia)"
Cohesion: 0.53
Nodes (5): _bandpower(), compute_ica_features(), evaluate_ica_components(), _kurtosis_simple(), _zscore()

### Community 56 - "Filtrado (Filtering.jl)"
Cohesion: 0.36
Nodes (8): apply_bandpass(), apply_bandreject(), apply_highpass(), apply_lowpass(), apply_notch(), _build(), _filt(), filter_recording()

### Community 59 - "Topomaps.jl"
Cohesion: 0.67
Nodes (3): CairoMakie, plot_topomap(), _topomap_resolve_style()

### Community 61 - "CSD (CSD.jl)"
Cohesion: 0.83
Nodes (3): apply_csd(), _apply_surface_laplacian(), _spline_matrices()

### Community 62 - "Métricas de grafo (GraphMetrics.jl)"
Cohesion: 0.67
Nodes (5): compute_graph_metrics(), _get_threshold(), _path_length_efficiency(), _threshold_matrix(), _weighted_clustering()

### Community 63 - "BIDSLoader.jl"
Cohesion: 0.83
Nodes (3): load_eeg_bids(), _load_electrode_positions(), _parse_bids_json()

### Community 65 - "Control de calidad (QualityControl.jl)"
Cohesion: 0.48
Nodes (5): compute_channel_spectral_qc(), compute_channel_stats(), flag_bad_channels(), qc_report(), welch_psd_raw()

### Community 70 - "Inspección ICA (ICAInspection.jl)"
Cohesion: 0.47
Nodes (3): has_manual_ica_labels(), _ica_manual_label_candidates(), load_ica_labels()

## Ambiguous Edges - Review These
- `Untested production group statistics (GroupStats.jl unused)` → `statistics section (mann_whitney, wilcoxon, fdr_q 0.05)`  [AMBIGUOUS]
  mne_brain/config/pipeline_config.yaml · relation: conceptually_related_to
- `mne_brain pipeline_config.yaml` → `Config mirror sync policy (YAML <-> Julia TOML)`  [AMBIGUOUS]
  mne_brain/config/pipeline_config.yaml · relation: conceptually_related_to

## Knowledge Gaps
- **208 isolated node(s):** `0. Estado al cerrar`, `1. Cómo trabajar en Cursor`, `2. Bucle por capítulo (Report_Pre ↔ código ↔ resultados)`, `3. Mapa capítulo → código → salidas`, `Prompt de arranque para Cursor` (+203 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 451 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **9 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `Untested production group statistics (GroupStats.jl unused)` and `statistics section (mann_whitney, wilcoxon, fdr_q 0.05)`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `Single-subject viewer flow raw->filter->ICA->epochs->spectral->connectivity->surrogate` connect `Single-subject viewer flow raw->filter->ICA->epochs->spectral->connectivity->surrogate` to `Visor aux: surrogates`, `Visor aux: componentes ICA`, `src/interactive README (viewer launch guide)`, `Visor aux: conectividad`, `plot_raw.jl`?**
  _High betweenness centrality (0.210) - this node is a cross-community bridge._
- **Are the 8 inferred relationships involving `load_config()` (e.g. with `run_batch()` and `run_single()`) actually correct?**
  _`load_config()` has 8 INFERRED edges - model-reasoned connections that need verification._
- **What connects `0. Estado al cerrar`, `1. Cómo trabajar en Cursor`, `2. Bucle por capítulo (Report_Pre ↔ código ↔ resultados)` to the rest of the system?**
  _208 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `compute_wpli` be split into smaller, more focused modules?**
  _Cohesion score 0.11174242424242424 - nodes in this community are weakly interconnected._
- **What is the exact relationship between `mne_brain pipeline_config.yaml` and `Config mirror sync policy (YAML <-> Julia TOML)`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `AGENTS.md — NeuroMIND agent guide` connect `Reglas de agentes (AGENTS.md)` to `Estimadores wPLI (wPLI.jl)`, `Dashboard Genie (App.jl)`, `run_phase8_validation.py`, `Suite de tests Julia`, `Visión general (README)`, `ICA FastICA (ICACore.jl)`, `Longitudinal.jl`, `Batch & group analysis audit 2026-07-25`, `run_batch_pipeline.jl`?**
  _High betweenness centrality (0.171) - this node is a cross-community bridge._