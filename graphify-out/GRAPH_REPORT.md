# Graph Report - NeuroMIND  (2026-10-07)

## Corpus Check
- 98 files · ~188,307 words
- Verdict: corpus is large enough that graph structure adds value.
- Unclassified: 8 file(s) not represented in the graph (top: (none) 5, .toml 3)

## Summary
- 1151 nodes · 2215 edges · 74 communities (53 shown, 21 thin omitted)
- Extraction: 92% EXTRACTED · 8% INFERRED · 0% AMBIGUOUS · INFERRED: 174 edges (avg confidence: 0.87)
- Token cost: 303,068 input · 0 output

## Community Hubs (Navigation)
- mne_brain: espectro, conectividad y tests
- mne_brain: scripts de ejecución por fases
- Transversal.jl: estadística y figuras
- Longitudinal.jl: estadística y figuras
- mne_brain: épocas y preprocesado
- mne_brain: ICA (fase 4)
- Reglas de agentes y suite de tests
- Construcción BIDS (Python y Julia)
- mne_brain: validación fase 8
- mne_brain: carga BIDS/BrainVision
- mne_brain: configuración y tipos
- Tipos centrales Julia (types.jl)
- Comparación de rechazos NeuroMIND vs MNE
- Visor aux: surrogates
- Visor transversal interactivo
- Visor aux: épocas
- Visor aux: componentes ICA
- Visor aux: baseline
- Visor longitudinal y contrato estadístico
- mne_brain: filtrado
- viewer_common.js (UI compartida)
- Batch fase C y su ejecución
- Pipeline single-subject de 8 pasos
- Visor aux: conectividad
- Visor aux: espectro
- Visor aux: histograma raw
- Auditorías de trazabilidad M05
- Visor aux: ICA antes/después
- Validez estadística y espejo de config MNE
- Visor aux: PSD espectral
- Visor aux: filtrado vs raw
- Visor aux: señal raw
- Visor aux: butterfly raw
- Módulo NeuroMIND.jl
- Estimadores wPLI (wPLI.jl)
- Dashboard Genie (App.jl)
- Diseño TDD de ejecución (runtime)
- mne_brain: control de calidad
- Soporte de visores interactivos
- Visión general y pipeline (README)
- Proyecto, dataset y estimadores wPLI
- Visor longitudinal: almacén de datos
- Config.jl (rutas)
- Métodos estadísticos de grupo
- mne_brain: tests de BIDS loader
- Auditoría y preparación del dataset
- Política de canales y QC
- mne_brain: tests de preprocesado
- Simplificación de la capa de cohorte
- Hallazgos de la auditoría 2026-07-25
- Salidas y tabla de decisión QC
- Hipótesis EM y bandas de frecuencia
- mne_brain: EpochSet
- Lanzador del dashboard
- Lanzador single-subject
- GraphPlots.jl
- Topomaps.jl
- Heatmaps.jl
- Spectra.jl
- mne_brain: pyproject
- apply_baseline (doc)

## God Nodes (most connected - your core abstractions)
1. `Transversal` - 54 edges
2. `Longitudinal` - 50 edges
3. `load_config()` - 41 edges
4. `PipelineConfig` - 30 edges
5. `EEGRecording` - 25 edges
6. `run()` - 24 edges
7. `run_single_subject_pipeline()` - 22 edges
8. `run()` - 21 edges
9. `run_phase4()` - 19 edges
10. `load_eeg_bids()` - 17 edges

## Surprising Connections (you probably didn't know these)
- `Dynamic per-subject QC channel exclusion (z>3)` --references--> `run_single_subject_pipeline()`  [INFERRED]
  AGENTS.md → src/SingleSubjectPipeline.jl
- `run_single_subject_pipeline()` --calls--> `load_ss_config`  [EXTRACTED]
  src/SingleSubjectPipeline.jl → README.md
- `8-step single-subject pipeline` --implements--> `run_single_subject_pipeline()`  [EXTRACTED]
  README.md → src/SingleSubjectPipeline.jl
- `paired_change_summary()` --implements--> `Bootstrap 95% CI (B=5000, fixed seed)`  [EXTRACTED]
  src/longitudinal/Longitudinal.jl → README.md
- `statistics section (mann_whitney, wilcoxon, fdr_q 0.05)` --conceptually_related_to--> `Untested production group statistics (GroupStats.jl unused)`  [AMBIGUOUS]
  mne_brain/config/pipeline_config.yaml → docs/audits/batch_and_group_analysis_audit_2026-07-25.md

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **Single-subject 8-step pipeline flow** — src_singlesubjectpipeline_load_ss_config, src_io_brainvisionloader, src_preprocessing_filtering_filter_recording, src_ica_icacore, src_segmentation_epochs_segment_recording, src_spectral_powerspectrum_compute_psd, src_connectivity_wpli_compute_wpli, src_singlesubjectpipeline__save_all_results [EXTRACTED 1.00]
- **Configurable wPLI estimators** — readme_hilbert_estimator, readme_fourier_csd_estimator, readme_multitaper_estimator, readme_dwpli [EXTRACTED 1.00]
- **Cohort statistics production → contract → viewer flow** — src_transversal_transversal, src_longitudinal_longitudinal, readme_statistics_contract_v2, src_interactive_plot_transversal, src_interactive_plot_longitudinal, agents_single_statistics_source [INFERRED 0.85]
- **Proposed src/runtime/ module set around RunContext** — docs_tdd_ejecucion_rutinas_runtime_layer, docs_tdd_ejecucion_rutinas_runcontext, docs_tdd_ejecucion_rutinas_resultsmanager, docs_tdd_ejecucion_rutinas_runmanifest, docs_tdd_ejecucion_rutinas_runlogging, docs_tdd_ejecucion_rutinas_console, docs_tdd_ejecucion_rutinas_progress, docs_tdd_ejecucion_rutinas_preflight, docs_tdd_ejecucion_rutinas_summary [EXTRACTED 1.00]
- **M05 traceability chain: validation prompt, July-14 audit, July-25 audit** — docs_prompt_claudecode_lanzar_m05, docs_audits_run_consistency_audit, docs_audits_batch_and_group_analysis_audit_2026_07_25, docs_prompt_claudecode_lanzar_m05_sub_m05_reference_run [EXTRACTED 1.00]
- **Threats to group-analysis validity identified in Phase C audit** — docs_audits_batch_and_group_analysis_audit_2026_07_25_untested_group_statistics, docs_audits_batch_and_group_analysis_audit_2026_07_25_longitudinal_pairing_filter_bug, docs_audits_batch_and_group_analysis_audit_2026_07_25_hard_intersect_channel_policy, docs_audits_batch_and_group_analysis_audit_2026_07_25_dynamic_qc_channel_exclusion, docs_audits_batch_and_group_analysis_audit_2026_07_25_group_results_summary [INFERRED 0.85]

## Communities (74 total, 21 thin omitted)

### Community 0 - "mne_brain: espectro, conectividad y tests"
Cohesion: 0.05
Nodes (29): default_config_path(), load_config(), compute_wpli(), connectivity_matrices(), _plot_wpli_heatmap(), save_connectivity_results(), compute_band_power(), compute_psd() (+21 more)

### Community 1 - "mne_brain: scripts de ejecución por fases"
Cohesion: 0.07
Nodes (32): append_log_row(), build_parser(), discover_recordings(), _elapsed(), _exclude_recording(), _exclusion_path(), filter_recordings(), is_batch_mode() (+24 more)

### Community 2 - "Transversal.jl: estadística y figuras"
Cohesion: 0.09
Nodes (54): apply_theme!(), _band_row(), bootstrap_cohen_d_ci(), bootstrap_mean_diff_ci(), _ch_xy(), fig_forest_by_band(), fig_interaction_common_channels(), fig_interaction_ec_eo() (+46 more)

### Community 3 - "Longitudinal.jl: estadística y figuras"
Cohesion: 0.08
Nodes (50): Bootstrap 95% CI (B=5000, fixed seed), apply_theme!(), _band_scores(), bootstrap_cohen_dz_ci(), _ch_xy(), evaluate_inclusion(), fig_alpha_delta_pareado(), fig_forest_by_band_dz() (+42 more)

### Community 4 - "mne_brain: épocas y preprocesado"
Cohesion: 0.09
Nodes (19): main(), apply_baseline(), epoch_summary(), load_cleaned_raw(), make_epochs(), _npz_to_raw(), reject_artifacts(), save_epoch_exclusion() (+11 more)

### Community 5 - "mne_brain: ICA (fase 4)"
Cohesion: 0.10
Nodes (19): main(), _save_matrix(), _save_topomaps(), ICAResult, _bandpower(), compute_ica_features(), evaluate_ica_components(), _kurtosis() (+11 more)

### Community 6 - "Reglas de agentes y suite de tests"
Cohesion: 0.06
Nodes (28): A&S 7.1.26 polynomial erf approximation, EEG_Julia original project (reference), AGENTS.md — NeuroMIND agent guide, Test suite consolidated in test/runtests.jl, CLAUDE.md (@AGENTS.md include), Longitudinal, Current Source Density (CSD), Independent Component Analysis (ICA) (+20 more)

### Community 7 - "Construcción BIDS (Python y Julia)"
Cohesion: 0.10
Nodes (18): build_all(), build_metadata_dict(), copy_support_files(), load_inventory(), load_participants(), main(), read_vhdr_with_mne(), resolve_vhdr_path() (+10 more)

### Community 8 - "mne_brain: validación fase 8"
Cohesion: 0.14
Nodes (14): _align_channels(), _fig_band_power(), _fig_heatmaps(), _fig_psd(), _fig_wpli_scatter_ba(), _get_paths(), _load_band_power(), _load_mb_psd() (+6 more)

### Community 9 - "mne_brain: carga BIDS/BrainVision"
Cohesion: 0.14
Nodes (14): bids_root_dir(), BrainVisionHeader, _first_existing(), load_eeg_bids(), load_eeg_brainvision(), load_eeg_tsv(), load_electrode_positions(), normalize_task() (+6 more)

### Community 10 - "mne_brain: configuración y tipos"
Cohesion: 0.10
Nodes (13): _normalize_bands(), ClinicalData, ConnectivityMatrix, GraphMetrics, GroupAnalysis, LongitudinalAnalysis, PipelineConfig, Session (+5 more)

### Community 11 - "Tipos centrales Julia (types.jl)"
Cohesion: 0.16
Nodes (16): ClinicalData, ConnectivityMatrix, duration(), EEGRecording, EpochSet, GraphMetrics, GroupAnalysis, ICAResult (+8 more)

### Community 12 - "Comparación de rechazos NeuroMIND vs MNE"
Cohesion: 0.20
Nodes (18): _as_int(), _as_pct(), build_parser(), compare_recordings(), _decision(), _diff(), _fmt_top(), load_julia_recordings() (+10 more)

### Community 13 - "Visor aux: surrogates"
Cohesion: 0.11
Nodes (15): BandPack, handle_request(), CairoMakie, CSV, DataFrames, Dates, Serialization, Sockets (+7 more)

### Community 14 - "Visor transversal interactivo"
Cohesion: 0.11
Nodes (21): CondStore, handle_request(), IncompatibleTransversalResultsError, CairoMakie, CSV, DataFrames, Dates, Sockets (+13 more)

### Community 15 - "Visor aux: épocas"
Cohesion: 0.11
Nodes (16): EpochStore, handle_request(), CairoMakie, CSV, DataFrames, Dates, NeuroMIND, Serialization (+8 more)

### Community 16 - "Visor aux: componentes ICA"
Cohesion: 0.11
Nodes (17): compute_ic_psd(), handle_request(), ICACompStore, CairoMakie, CSV, DataFrames, Dates, FFTW (+9 more)

### Community 17 - "Visor aux: baseline"
Cohesion: 0.12
Nodes (15): BaselineStore, handle_request(), CairoMakie, CSV, DataFrames, Dates, NeuroMIND, Serialization (+7 more)

### Community 18 - "Visor longitudinal y contrato estadístico"
Cohesion: 0.14
Nodes (17): mean_strength = (n_channels−1) × mean_wPLI normalization, statistics_contract.toml schema v2, handle_request(), CairoMakie, CSV, DataFrames, Dates, Sockets (+9 more)

### Community 19 - "mne_brain: filtrado"
Cohesion: 0.19
Nodes (9): EEGRecording, apply_bandreject(), apply_highpass(), apply_lowpass(), apply_notch(), _apply_sos(), describe_filter_chain(), filter_recording() (+1 more)

### Community 20 - "viewer_common.js (UI compartida)"
Cohesion: 0.22
Nodes (15): barChart(), bindHeatmap(), drawColorbar(), drawGraph(), drawHeatmap(), drawStrength(), drawTopo(), drawVolcano() (+7 more)

### Community 21 - "Batch fase C y su ejecución"
Cohesion: 0.13
Nodes (16): Phase C batch run 2026-07-25 (201 OK/5 SKIP/0 ERR), Report_Pre citability precondition (reprocess with current config), Uncurated duplicate recording M16 T1 EC, Phase C batch run (201 OK / 5 SKIP / 0 ERR), Progress (ProgressMeter), Logging, BatchJob, Dates (+8 more)

### Community 22 - "Pipeline single-subject de 8 pasos"
Cohesion: 0.19
Nodes (18): RunLogging (LoggingExtras TeeLogger), filter_recording, segment_recording, _ica_effective_params(), CairoMakie, Serialization, _log(), _print_cols_table() (+10 more)

### Community 23 - "Visor aux: conectividad"
Cohesion: 0.14
Nodes (14): ConnStore, handle_request(), CairoMakie, CSV, DataFrames, Dates, Sockets, Statistics (+6 more)

### Community 24 - "Visor aux: espectro"
Cohesion: 0.13
Nodes (13): handle_request(), CairoMakie, CSV, DataFrames, Dates, Sockets, Statistics, _json_escape() (+5 more)

### Community 25 - "Visor aux: histograma raw"
Cohesion: 0.16
Nodes (16): _channel_stats(), handle_request(), _histogram(), CairoMakie, CSV, DataFrames, Dates, Sockets (+8 more)

### Community 26 - "Auditorías de trazabilidad M05"
Cohesion: 0.16
Nodes (17): Batch & group analysis audit 2026-07-25, Traceability audit — M05 run (2026-07-14), channel_statistics_compare.csv without generator, ICA cache temporal mix (FastICA May vs pipeline July), ICA auto-rejection log vs load_ica_labels code mismatch, Missing git_commit / config_hash / run_id, Report_Pre / Slidev figure sync drift, Stale *_EC.png figures / dual layout (+9 more)

### Community 27 - "Visor aux: ICA antes/después"
Cohesion: 0.15
Nodes (13): handle_request(), ICASignalStore, CairoMakie, CSV, DataFrames, Dates, Sockets, _json_int() (+5 more)

### Community 28 - "Validez estadística y espejo de config MNE"
Cohesion: 0.14
Nodes (12): Untested production group statistics (GroupStats.jl unused), Welch/paired t-tests use Z (normal) approximation, mne_brain AGENTS.md, Batch errors 'Epochs-object is empty' (AR ±70 µV), mne_brain Phase 8 surrogates + FDR (pending), mne_brain pipeline_config.yaml, EEG frequency bands (ALPHA 7.8-11.7 Hz, 7 bands), statistics section (mann_whitney, wilcoxon, fdr_q 0.05) (+4 more)

### Community 29 - "Visor aux: PSD espectral"
Cohesion: 0.17
Nodes (13): handle_request(), CairoMakie, CSV, DataFrames, Dates, Sockets, Statistics, _json_escape() (+5 more)

### Community 30 - "Visor aux: filtrado vs raw"
Cohesion: 0.19
Nodes (11): handle_request(), CairoMakie, CSV, DataFrames, Dates, Sockets, main(), MultiStore (+3 more)

### Community 31 - "Visor aux: señal raw"
Cohesion: 0.19
Nodes (11): handle_request(), CairoMakie, CSV, DataFrames, Dates, Sockets, main(), open_browser() (+3 more)

### Community 32 - "Visor aux: butterfly raw"
Cohesion: 0.19
Nodes (12): handle_request(), CairoMakie, CSV, DataFrames, Dates, Sockets, Statistics, main() (+4 more)

### Community 33 - "Módulo NeuroMIND.jl"
Cohesion: 0.14
Nodes (12): CSV, DataFrames, Dates, DSP, FFTW, LinearAlgebra, Random, Serialization (+4 more)

### Community 34 - "Estimadores wPLI (wPLI.jl)"
Cohesion: 0.27
Nodes (4): AbstractWPLIEstimator, FourierCSDEstimator, HilbertEstimator, MultitaperEstimator

### Community 35 - "Dashboard Genie (App.jl)"
Cohesion: 0.18
Nodes (10): Base64, Genie, Genie.Renderer.Html, Genie.Renderer.Json, Genie.Requests, Genie.Router, Genie.Server, Dates (+2 more)

### Community 36 - "Diseño TDD de ejecución (runtime)"
Cohesion: 0.20
Nodes (10): Versioning of 12 single-subject aux viewers out of results/, Atomic staging directory ({unit}.tmp/ + mv swap), Console (Crayons, TTY/ASCII degradation), Preflight checks (fail early), ResultsManager, run_manifest.json (per execution unit), RunContext, RunManifest module (+2 more)

### Community 37 - "mne_brain: control de calidad"
Cohesion: 0.27
Nodes (4): compute_channel_stats(), flag_bad_channels(), qc_report(), _task_from_condition()

### Community 38 - "Soporte de visores interactivos"
Cohesion: 0.18
Nodes (7): src/interactive README (viewer launch guide), Port :8767 shared by plot_spectral and plot_spectral_PSD, CairoMakie, CSV, DataFrames, TOML, ViewerSupport

### Community 39 - "Visión general y pipeline (README)"
Cohesion: 0.24
Nodes (6): 12 single-subject aux viewers (:8765–:8775), Genie dashboard (panels 0–15), 8-step single-subject pipeline, mne_brain MNE-Python cross-validation, README.md — NeuroMIND user guide, sub-M05/ses-T2/eyesclosed reference case

### Community 40 - "Proyecto, dataset y estimadores wPLI"
Cohesion: 0.25
Nodes (9): MINDEM-IMIBIC dataset (41 MS + 36/37 controls), NeuroMIND EEG connectivity framework, Debiased wPLI² (dwPLI), FourierCSD wPLI estimator, Hilbert wPLI estimator (active), Multitaper (DPSS) wPLI estimator, Vinck et al. 2011, Weighted Phase Lag Index (wPLI) (+1 more)

### Community 41 - "Visor longitudinal: almacén de datos"
Cohesion: 0.22
Nodes (6): CondStore, IncompatibleResultsError, LongViewer, _require_columns(), _required_matrix(), _validate_provenance()

### Community 42 - "Config.jl (rutas)"
Cohesion: 0.31
Nodes (3): ensure_dirs(), results_dir(), subject_results_dir()

### Community 43 - "Métodos estadísticos de grupo"
Cohesion: 0.29
Nodes (8): Condition-specific T1–T2 pair selection, Group×condition EC×EO interaction (diff-of-diff), Benjamini–Hochberg FDR correction, Longitudinal analysis (MS T1→T2), Mann–Whitney U test (transversal), Circular-shift surrogate test, Transversal analysis (MS T1 vs Control), Exact conditional Wilcoxon signed-rank (N≤30)

### Community 44 - "mne_brain: tests de BIDS loader"
Cohesion: 0.32
Nodes (6): find_vhdr_by_name(), _cfg(), test_find_vhdr_by_name_skips_excluded(), test_load_eeg_bids_from_tsv(), test_read_vhdr_header(), test_resolve_vhdr_path_prefers_raw_data_root()

### Community 45 - "Auditoría y preparación del dataset"
Cohesion: 0.29
Nodes (5): Dataset preparation (Phase A audit + Phase B BIDS), Dates, main(), parse_vhdr_filename(), SubjectDemog

### Community 46 - "Política de canales y QC"
Cohesion: 0.29
Nodes (7): Hard channel intersection across cohort, Dynamic per-subject QC channel exclusion (z>3), config/pipeline.toml (single config), 31-channel montage (exclude_fp2=false, 465 edges), ±70 µV artifact rejection, reject_artifacts, load_ss_config

### Community 47 - "mne_brain: tests de preprocesado"
Cohesion: 0.53
Nodes (5): _cfg(), _meta(), test_channel_stats_and_bad_channel_flagging(), test_describe_filter_chain_eeg_julia_order(), test_filter_recording_reduces_50hz_component()

### Community 48 - "Simplificación de la capa de cohorte"
Cohesion: 0.40
Nodes (3): Cohort layer radical simplification 2026-07-27, Manual sync of figures to Report_Pre/figures/plots, Cohort interactive viewers (:8780/:8781)

### Community 49 - "Hallazgos de la auditoría 2026-07-25"
Cohesion: 0.40
Nodes (5): Dynamic QC channel exclusion (bad_ch always merged), FDR-BH applied per band (7 families), Group results: transversal EC 21 sig (ALPHA), EO 3, longitudinal 0, Hard-intersect channel policy in group analyses, Longitudinal pairing filter-order bug

### Community 50 - "Salidas y tabla de decisión QC"
Cohesion: 0.40
Nodes (5): Amplitude warning (σ̄ > 20 µV), Single BIDS output tree results/subjects/, config_snapshot.toml (per-recording provenance), QC decision table (qc_decision_table.csv), _save_all_results

### Community 51 - "Hipótesis EM y bandas de frecuencia"
Cohesion: 0.50
Nodes (4): Hypothesis: MS alters alpha/beta connectivity, Multiple Sclerosis (MS), Delta band below min_cycles_for_wpli, 7 frequency bands (δ θ α β_low β_mid β_high γ)

### Community 53 - "Lanzador del dashboard"
Cohesion: 0.50
Nodes (3): NeuroMIND, Pkg, TOML

## Ambiguous Edges - Review These
- `Untested production group statistics (GroupStats.jl unused)` → `statistics section (mann_whitney, wilcoxon, fdr_q 0.05)`  [AMBIGUOUS]
  mne_brain/config/pipeline_config.yaml · relation: conceptually_related_to
- `Config mirror sync policy (YAML <-> Julia TOML)` → `mne_brain pipeline_config.yaml`  [AMBIGUOUS]
  mne_brain/config/pipeline_config.yaml · relation: conceptually_related_to

## Knowledge Gaps
- **210 isolated node(s):** `mne-brain`, `ClinicalData`, `SpectralResult`, `SurrogateResult`, `GraphMetrics` (+205 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 478 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)
- **21 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **What is the exact relationship between `Untested production group statistics (GroupStats.jl unused)` and `statistics section (mann_whitney, wilcoxon, fdr_q 0.05)`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `AGENTS.md — NeuroMIND agent guide` connect `Reglas de agentes y suite de tests` to `Dashboard Genie (App.jl)`, `Visión general y pipeline (README)`, `Proyecto, dataset y estimadores wPLI`, `Política de canales y QC`, `Simplificación de la capa de cohorte`, `Visor longitudinal y contrato estadístico`, `Auditorías de trazabilidad M05`?**
  _High betweenness centrality (0.149) - this node is a cross-community bridge._
- **Are the 8 inferred relationships involving `load_config()` (e.g. with `run_batch()` and `run_single()`) actually correct?**
  _`load_config()` has 8 INFERRED edges - model-reasoned connections that need verification._
- **What connects `mne-brain`, `ClinicalData`, `SpectralResult` to the rest of the system?**
  _210 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `mne_brain: espectro, conectividad y tests` be split into smaller, more focused modules?**
  _Cohesion score 0.05403348554033485 - nodes in this community are weakly interconnected._
- **What is the exact relationship between `Config mirror sync policy (YAML <-> Julia TOML)` and `mne_brain pipeline_config.yaml`?**
  _Edge tagged AMBIGUOUS (relation: conceptually_related_to) - confidence is low._
- **Why does `Single-subject viewer flow raw->filter->ICA->epochs->spectral->connectivity->surrogate` connect `Auditorías de trazabilidad M05` to `Soporte de visores interactivos`, `Visor aux: surrogates`, `Visor aux: componentes ICA`, `Visor aux: conectividad`, `Visor aux: señal raw`?**
  _High betweenness centrality (0.125) - this node is a cross-community bridge._