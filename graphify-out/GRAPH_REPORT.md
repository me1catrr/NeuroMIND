# Graph Report - NeuroMIND  (2026-10-07)

## Corpus Check
- cluster-only mode — file stats not available

## Summary
- 1471 nodes · 3015 edges · 64 communities (55 shown, 9 thin omitted)
- Extraction: 94% EXTRACTED · 6% INFERRED · 0% AMBIGUOUS · INFERRED: 174 edges (avg confidence: 0.87)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `f884fa21`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- mne_brain: wPLI y tests de conectividad
- mne_brain: scripts de ejecución por fases
- Transversal.jl: estadística y figuras
- Longitudinal.jl: estadística y figuras
- mne_brain: épocas y preprocesado
- mne_brain: ICA y validación fase 8
- Suite de tests Julia
- Auditoría y construcción BIDS
- mne_brain: configuración y espectro
- mne_brain: carga BIDS/BrainVision
- mne_brain: tipos y config
- Tipos centrales Julia (types.jl)
- Comparación de rechazos NeuroMIND vs MNE
- Visor aux: surrogates
- Visor transversal interactivo
- Visor aux: épocas
- Visor aux: componentes ICA
- Visor aux: baseline
- Visor longitudinal interactivo
- mne_brain: filtrado y QC
- viewer_common.js (UI compartida)
- mne_brain: tests de procesado
- Pipeline single-subject y batch
- Visor aux: conectividad
- Visor aux: espectro
- Visor aux: histograma raw
- Visor aux: señal raw y auditorías
- Visor aux: ICA antes/después
- Estadística de grupo (GroupStats.jl)
- Visor aux: PSD espectral
- Visor aux: filtrado vs raw
- mne_brain: PSD
- Visor aux: butterfly raw
- Módulo NeuroMIND.jl
- Estimadores wPLI (wPLI.jl)
- Dashboard Genie (App.jl)
- Visor transversal: contrato y errores
- Reglas de agentes (AGENTS.md)
- Soporte de visores (ViewerSupport)
- Visión general (README)
- Lanzadores de análisis de grupo
- Visor longitudinal: contrato y errores
- Config.jl (rutas)
- Pipeline de 8 pasos (doc)
- ICA FastICA (ICACore.jl)
- Simplificación de la capa de cohorte
- mne_brain: EpochSet
- Lanzador del dashboard
- Clasificación ICA (Julia)
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
- `Circular-shift surrogate test` --references--> `compute_wpli()`  [INFERRED]
  README.md → src/connectivity/wPLI.jl
- `Dynamic per-subject QC channel exclusion (z>3)` --references--> `run_single_subject_pipeline()`  [INFERRED]
  AGENTS.md → src/SingleSubjectPipeline.jl
- `8-step single-subject pipeline` --implements--> `run_single_subject_pipeline()`  [EXTRACTED]
  README.md → src/SingleSubjectPipeline.jl
- `paired_change_summary()` --implements--> `Bootstrap 95% CI (B=5000, fixed seed)`  [EXTRACTED]
  src/longitudinal/Longitudinal.jl → README.md
- `compute_wpli()` --implements--> `FourierCSD wPLI estimator`  [EXTRACTED]
  src/connectivity/wPLI.jl → README.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **M05 traceability chain: validation prompt, July-14 audit, July-25 audit** — docs_prompt_claudecode_lanzar_m05, docs_audits_run_consistency_audit, docs_audits_batch_and_group_analysis_audit_2026_07_25, docs_prompt_claudecode_lanzar_m05_sub_m05_reference_run [EXTRACTED 1.00]
- **Single-subject 8-step pipeline flow** — src_singlesubjectpipeline_load_ss_config, src_io_brainvisionloader, src_preprocessing_filtering_filter_recording, src_ica_icacore, src_segmentation_epochs_segment_recording, src_spectral_powerspectrum_compute_psd, src_connectivity_wpli_compute_wpli, src_singlesubjectpipeline__save_all_results [EXTRACTED 1.00]
- **Proposed src/runtime/ module set around RunContext** — docs_tdd_ejecucion_rutinas_runtime_layer, docs_tdd_ejecucion_rutinas_runcontext, docs_tdd_ejecucion_rutinas_resultsmanager, docs_tdd_ejecucion_rutinas_runmanifest, docs_tdd_ejecucion_rutinas_runlogging, docs_tdd_ejecucion_rutinas_console, docs_tdd_ejecucion_rutinas_progress, docs_tdd_ejecucion_rutinas_preflight, docs_tdd_ejecucion_rutinas_summary [EXTRACTED 1.00]
- **Configurable wPLI estimators** — readme_hilbert_estimator, readme_fourier_csd_estimator, readme_multitaper_estimator, readme_dwpli [EXTRACTED 1.00]
- **Cohort statistics production → contract → viewer flow** — src_transversal_transversal, src_longitudinal_longitudinal, readme_statistics_contract_v2, src_interactive_plot_transversal, src_interactive_plot_longitudinal, agents_single_statistics_source [INFERRED 0.85]
- **Threats to group-analysis validity identified in Phase C audit** — docs_audits_batch_and_group_analysis_audit_2026_07_25_untested_group_statistics, docs_audits_batch_and_group_analysis_audit_2026_07_25_longitudinal_pairing_filter_bug, docs_audits_batch_and_group_analysis_audit_2026_07_25_hard_intersect_channel_policy, docs_audits_batch_and_group_analysis_audit_2026_07_25_dynamic_qc_channel_exclusion, docs_audits_batch_and_group_analysis_audit_2026_07_25_group_results_summary [INFERRED 0.85]

## Communities (64 total, 9 thin omitted)

### Community 0 - "mne_brain: wPLI y tests de conectividad"
Cohesion: 0.10
Nodes (15): main(), _resolve_dirs(), compute_wpli(), connectivity_matrices(), _plot_wpli_heatmap(), save_connectivity_results(), _make_coupled_epochs(), _make_epochs() (+7 more)

### Community 1 - "mne_brain: scripts de ejecución por fases"
Cohesion: 0.05
Nodes (41): build_all(), build_metadata_dict(), copy_support_files(), load_inventory(), load_participants(), main(), read_vhdr_with_mne(), resolve_vhdr_path() (+33 more)

### Community 2 - "Transversal.jl: estadística y figuras"
Cohesion: 0.07
Nodes (77): anatomical_order(), apply_theme!(), _assign_ranks(), _band_row(), bh_qvalues(), bootstrap_cohen_d_ci(), bootstrap_mean_diff_ci(), _ch_xy() (+69 more)

### Community 3 - "Longitudinal.jl: estadística y figuras"
Cohesion: 0.06
Nodes (75): mean_strength = (n_channels−1) × mean_wPLI normalization, Bootstrap 95% CI (B=5000, fixed seed), anatomical_order(), apply_theme!(), _as_bool(), _assign_ranks(), _band_scores(), bh_qvalues() (+67 more)

### Community 4 - "mne_brain: épocas y preprocesado"
Cohesion: 0.14
Nodes (10): main(), apply_baseline(), epoch_summary(), load_cleaned_raw(), make_epochs(), _npz_to_raw(), reject_artifacts(), save_epoch_results() (+2 more)

### Community 5 - "mne_brain: ICA y validación fase 8"
Cohesion: 0.06
Nodes (33): main(), _save_matrix(), _save_topomaps(), _align_channels(), _fig_band_power(), _fig_heatmaps(), _fig_psd(), _fig_wpli_scatter_ba() (+25 more)

### Community 6 - "Suite de tests Julia"
Cohesion: 0.15
Nodes (15): Test suite consolidated in test/runtests.jl, Test, _cfg_seg(), _cfg_with_ica(), _cfg_with_root(), _cfg_with_surrogates(), CSV, DataFrames (+7 more)

### Community 7 - "Auditoría y construcción BIDS"
Cohesion: 0.10
Nodes (23): Dynamic QC channel exclusion (bad_ch always merged), FDR-BH applied per band (7 families), Group results: transversal EC 21 sig (ALPHA), EO 3, longitudinal 0, Hard-intersect channel policy in group analyses, Longitudinal pairing filter-order bug, Dataset preparation (Phase A audit + Phase B BIDS), _bids_subject_id(), is_excluded() (+15 more)

### Community 8 - "mne_brain: configuración y espectro"
Cohesion: 0.12
Nodes (14): default_config_path(), load_config(), compute_psd(), _make_alpha_epochs(), _make_epochs(), test_band_power_keys_match_config(), test_band_power_shape_and_positive(), test_mean_spectrum_equals_epoch_mean() (+6 more)

### Community 9 - "mne_brain: carga BIDS/BrainVision"
Cohesion: 0.11
Nodes (20): bids_root_dir(), BrainVisionHeader, find_vhdr_by_name(), _first_existing(), load_eeg_bids(), load_eeg_brainvision(), load_eeg_tsv(), load_electrode_positions() (+12 more)

### Community 10 - "mne_brain: tipos y config"
Cohesion: 0.10
Nodes (13): _normalize_bands(), ClinicalData, ConnectivityMatrix, GraphMetrics, GroupAnalysis, LongitudinalAnalysis, PipelineConfig, Session (+5 more)

### Community 11 - "Tipos centrales Julia (types.jl)"
Cohesion: 0.16
Nodes (16): ClinicalData, ConnectivityMatrix, duration(), EEGRecording, EpochSet, GraphMetrics, GroupAnalysis, ICAResult (+8 more)

### Community 12 - "Comparación de rechazos NeuroMIND vs MNE"
Cohesion: 0.20
Nodes (18): _as_int(), _as_pct(), build_parser(), compare_recordings(), _decision(), _diff(), _fmt_top(), load_julia_recordings() (+10 more)

### Community 13 - "Visor aux: surrogates"
Cohesion: 0.08
Nodes (47): _auto_interpret(), _band_json(), BandPack, _cache_get(), _ch_index(), _contrast_json(), _edge_row(), _extract_null_dist() (+39 more)

### Community 14 - "Visor transversal interactivo"
Cohesion: 0.13
Nodes (22): handle_request(), html_page(), CairoMakie, CSV, DataFrames, Dates, Sockets, Statistics (+14 more)

### Community 15 - "Visor aux: épocas"
Cohesion: 0.10
Nodes (32): _default_rejected_epoch(), _epoch_window(), EpochStore, handle_request(), html_page(), CairoMakie, CSV, DataFrames (+24 more)

### Community 16 - "Visor aux: componentes ICA"
Cohesion: 0.11
Nodes (30): _components_json(), compute_ic_psd(), _decision_explanation(), _feature_json(), handle_request(), html_page(), ICACompStore, CairoMakie (+22 more)

### Community 17 - "Visor aux: baseline"
Cohesion: 0.10
Nodes (30): _apply_bl(), _baseline_end_s(), BaselineStore, compute_offset_table(), _epoch_max_offset(), _epoch_raw(), handle_request(), html_page() (+22 more)

### Community 18 - "Visor longitudinal interactivo"
Cohesion: 0.13
Nodes (23): handle_request(), _honest_best_from_summary(), html_page(), CairoMakie, CSV, DataFrames, Dates, Sockets (+15 more)

### Community 19 - "mne_brain: filtrado y QC"
Cohesion: 0.11
Nodes (18): EEGRecording, apply_bandreject(), apply_highpass(), apply_lowpass(), apply_notch(), _apply_sos(), describe_filter_chain(), filter_recording() (+10 more)

### Community 20 - "viewer_common.js (UI compartida)"
Cohesion: 0.22
Nodes (15): barChart(), bindHeatmap(), drawColorbar(), drawGraph(), drawHeatmap(), drawStrength(), drawTopo(), drawVolcano() (+7 more)

### Community 21 - "mne_brain: tests de procesado"
Cohesion: 0.16
Nodes (8): _inject_artifact(), _make_raw(), test_baseline_idempotent_on_zero_mean(), test_baseline_zeroes_epoch_mean(), test_clean_signal_survives_ar(), test_make_epochs_count_and_duration(), test_make_epochs_no_overlap(), test_reject_drops_contaminated_epoch()

### Community 22 - "Pipeline single-subject y batch"
Cohesion: 0.05
Nodes (63): Hard channel intersection across cohort, Dynamic per-subject QC channel exclusion (z>3), Phase C batch run 2026-07-25 (201 OK/5 SKIP/0 ERR), config/pipeline.toml (single config), Progress (ProgressMeter), RunLogging (LoggingExtras TeeLogger), Logging, 31-channel montage (exclude_fp2=false, 465 edges) (+55 more)

### Community 23 - "Visor aux: conectividad"
Cohesion: 0.09
Nodes (43): _band_strength(), _betweenness(), _binary_adj(), _clustering(), ConnStore, _degree_from_adj(), _density_label(), _edges_for_band() (+35 more)

### Community 24 - "Visor aux: espectro"
Cohesion: 0.12
Nodes (29): _band_payload(), _band_row(), _curve_series(), handle_request(), html_page(), _index_row(), CairoMakie, CSV (+21 more)

### Community 25 - "Visor aux: histograma raw"
Cohesion: 0.15
Nodes (25): _channel_stats(), _grid_dims(), handle_request(), _hist_json_channel(), _histogram(), html_page(), CairoMakie, CSV (+17 more)

### Community 26 - "Visor aux: señal raw y auditorías"
Cohesion: 0.05
Nodes (48): Batch & group analysis audit 2026-07-25, Versioning of 12 single-subject aux viewers out of results/, Report_Pre citability precondition (reprocess with current config), Uncurated duplicate recording M16 T1 EC, Phase C batch run (201 OK / 5 SKIP / 0 ERR), Traceability audit — M05 run (2026-07-14), channel_statistics_compare.csv without generator, ICA cache temporal mix (FastICA May vs pipeline July) (+40 more)

### Community 27 - "Visor aux: ICA antes/después"
Cohesion: 0.14
Nodes (22): handle_request(), html_page(), ICASignalStore, CairoMakie, CSV, DataFrames, Dates, Sockets (+14 more)

### Community 28 - "Estadística de grupo (GroupStats.jl)"
Cohesion: 0.12
Nodes (23): CLAUDE.md (@AGENTS.md include), Untested production group statistics (GroupStats.jl unused), Welch/paired t-tests use Z (normal) approximation, mne_brain AGENTS.md, Batch errors 'Epochs-object is empty' (AR ±70 µV), mne_brain Phase 8 surrogates + FDR (pending), mne_brain pipeline_config.yaml, EEG frequency bands (ALPHA 7.8-11.7 Hz, 7 bands) (+15 more)

### Community 29 - "Visor aux: PSD espectral"
Cohesion: 0.14
Nodes (22): _band_row(), handle_request(), html_page(), _index_row(), CairoMakie, CSV, DataFrames, Dates (+14 more)

### Community 30 - "Visor aux: filtrado vs raw"
Cohesion: 0.16
Nodes (18): handle_request(), html_page(), CairoMakie, CSV, DataFrames, Dates, Sockets, _load_stage_csv() (+10 more)

### Community 31 - "mne_brain: PSD"
Cohesion: 0.22
Nodes (3): compute_band_power(), mean_spectrum(), save_spectral_results()

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

### Community 36 - "Visor transversal: contrato y errores"
Cohesion: 0.22
Nodes (11): _cond_dir(), CondStore, IncompatibleTransversalResultsError, load_condition(), load_viewer(), _require_columns(), _required_csv(), _required_matrix() (+3 more)

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

### Community 41 - "Visor longitudinal: contrato y errores"
Cohesion: 0.22
Nodes (11): _cond_dir(), CondStore, IncompatibleResultsError, load_condition(), load_viewer(), LongViewer, _require_columns(), _required_csv() (+3 more)

### Community 42 - "Config.jl (rutas)"
Cohesion: 0.27
Nodes (6): ensure_dirs(), _f(), load_subjects(), results_dir(), _s(), subject_results_dir()

### Community 43 - "Pipeline de 8 pasos (doc)"
Cohesion: 0.40
Nodes (3): 8-step single-subject pipeline, load_eeg_brainvision(), read_vhdr_header()

### Community 44 - "ICA FastICA (ICACore.jl)"
Cohesion: 0.53
Nodes (5): Independent Component Analysis (ICA), _fastica_attempt(), run_ica(), _sym_decorr(), _whiten_pca()

### Community 48 - "Simplificación de la capa de cohorte"
Cohesion: 0.29
Nodes (4): Cohort layer radical simplification 2026-07-27, Manual sync of figures to Report_Pre/figures/plots, Cohort interactive viewers (:8780/:8781), statistics_contract.toml schema v2

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
- **203 isolated node(s):** `ClinicalData`, `GraphMetrics`, `GroupAnalysis`, `LongitudinalAnalysis`, `Session` (+198 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 445 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **9 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `Untested production group statistics (GroupStats.jl unused)` and `statistics section (mann_whitney, wilcoxon, fdr_q 0.05)`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `Single-subject viewer flow raw->filter->ICA->epochs->spectral->connectivity->surrogate` connect `Visor aux: señal raw y auditorías` to `Visor aux: componentes ICA`, `Visor aux: surrogates`, `Visor transversal interactivo`, `Visor aux: conectividad`?**
  _High betweenness centrality (0.213) - this node is a cross-community bridge._
- **Are the 8 inferred relationships involving `load_config()` (e.g. with `run_batch()` and `run_single()`) actually correct?**
  _`load_config()` has 8 INFERRED edges - model-reasoned connections that need verification._
- **What connects `ClinicalData`, `GraphMetrics`, `GroupAnalysis` to the rest of the system?**
  _203 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `mne_brain: wPLI y tests de conectividad` be split into smaller, more focused modules?**
  _Cohesion score 0.10252100840336134 - nodes in this community are weakly interconnected._
- **What is the exact relationship between `mne_brain pipeline_config.yaml` and `Config mirror sync policy (YAML <-> Julia TOML)`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `AGENTS.md — NeuroMIND agent guide` connect `Reglas de agentes (AGENTS.md)` to `Estimadores wPLI (wPLI.jl)`, `Longitudinal.jl: estadística y figuras`, `Dashboard Genie (App.jl)`, `Suite de tests Julia`, `Visión general (README)`, `ICA FastICA (ICACore.jl)`, `Simplificación de la capa de cohorte`, `Pipeline single-subject y batch`, `Visor aux: señal raw y auditorías`, `Estadística de grupo (GroupStats.jl)`?**
  _High betweenness centrality (0.157) - this node is a cross-community bridge._