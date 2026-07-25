"""
# NeuroMIND — Framework de conectividad funcional EEG

Estudio BRAIN: conectividad funcional en Esclerosis Múltiple
basada en wPLI (weighted Phase Lag Index).

Fichero:    src/NeuroMIND.jl
Autor:      Rafael Castro Triguero <me1catrr@uco.es>
Modificado: 22-07-2026

## Configuración

Única fuente de parámetros: `config/pipeline.toml`
(leída con `load_ss_config`).

## Uso básico

```julia
using NeuroMIND

# Pipeline de un sujeto (ruta al TOML)
run_single_subject_pipeline("config/pipeline.toml")

# O cargar la config y lanzar el dashboard
cfg = load_ss_config("config/pipeline.toml")
launch_webapp(cfg; port=8080)
```

Desde terminal (recomendado):

```bash
julia --project=. scripts/run_single_subject.jl
julia --project=. scripts/launch_dashboard.jl
```

Ver `README.md` y `AGENTS.md` para el resto de scripts
(batch, transversal, longitudinal, auditoría BIDS).
"""
module NeuroMIND

using CSV, DataFrames, Dates, Statistics, LinearAlgebra, Random
using Serialization, TOML, FFTW, DSP, StatsBase

# ═══════════════════════════════════════════════════════════════
#  Núcleo
# ═══════════════════════════════════════════════════════════════

# ─── Tipos ─────────────────────────────────────────────────────
#  Structs inmutables del framework (EEGRecording, ICAResult, …).
include("types.jl")
export PipelineConfig, RecordingMeta, EEGRecording, EpochSet
export ICAResult, SpectralResult, ConnectivityMatrix, SurrogateResult
export GraphMetrics, Session, Subject, GroupAnalysis
export ClinicalData, StatResult
export n_channels, n_samples, n_epochs, n_samples_epoch, duration

# ─── I/O y configuración ───────────────────────────────────────
#  TOML, BIDS, BrainVision nativo y serialización de resultados.
include("io/Config.jl")
include("io/BIDSLoader.jl")
include("io/BrainVisionLoader.jl")
include("io/Serializer.jl")
export load_config, load_subjects, load_eeg_bids, save_result, load_result, result_exists
export results_dir, subject_results_dir, ensure_dirs
export read_vhdr_header, load_eeg_brainvision, bv_electrode_positions


# ═══════════════════════════════════════════════════════════════
#  Pipeline científico (orden fijo: filtrado → ICA → segmentación)
# ═══════════════════════════════════════════════════════════════

# ─── Control de calidad ────────────────────────────────────────
#  Estadísticas por canal y marcado de bad_ch (z-score).
include("qc/QualityControl.jl")
export compute_channel_stats, flag_bad_channels, qc_report
export welch_psd_raw, compute_channel_spectral_qc, compute_correlation_summary

# ─── Preprocesado / filtrado ───────────────────────────────────
#  HP / LP / Notch / Bandreject (Butterworth, filtfilt).
include("preprocessing/Filtering.jl")
export apply_highpass, apply_lowpass, apply_notch, apply_bandpass
export apply_bandreject, filter_recording, describe_filter_chain

# ─── ICA ───────────────────────────────────────────────────────
#  FastICA puro Julia + clasificación / rechazo de componentes.
#  Siempre sobre señal continua filtrada, ANTES de segmentar.
include("ica/ICACore.jl")
include("ica/ICAClassification.jl")
include("ica/ICAInspection.jl")
export run_ica, apply_ica_rejection, load_ica_labels, write_ica_labels_auto, has_manual_ica_labels
export compute_ica_features, evaluate_ica_components

# ─── Segmentación ──────────────────────────────────────────────
#  Épocas, baseline y rechazo de artefactos por amplitud.
include("segmentation/Epochs.jl")
export segment_recording, apply_baseline, reject_artifacts
export compute_epoch_quality_report, compute_channel_coverage

# ─── Análisis espectral ────────────────────────────────────────
#  PSD (Hamming) y potencia por banda.
include("spectral/PowerSpectrum.jl")
export compute_psd, band_power

# ─── Conectividad ──────────────────────────────────────────────
#  wPLI (Hilbert / FourierCSD / Multitaper), CSD opcional, grafos.
include("connectivity/wPLI.jl")
include("connectivity/CSD.jl")
include("connectivity/GraphMetrics.jl")
export compute_wpli, apply_csd, compute_graph_metrics

# ─── Inferencia estadística ────────────────────────────────────
#  Surrogates, FDR y tests de grupo.
include("statistics/Surrogates.jl")
include("statistics/FDR.jl")
include("statistics/GroupStats.jl")
export surrogate_test, fdr_correction, threshold_connectivity
export validate_connectivity_matrix, validate_surrogate_result
export mann_whitney_test, wilcoxon_signed_rank, spearman_correlation
export compare_groups_stats


# ═══════════════════════════════════════════════════════════════
#  Análisis de cohorte y salida
# ═══════════════════════════════════════════════════════════════

# ─── Longitudinal ──────────────────────────────────────────────
#  Análisis T1→T2: scripts/run_longitudinal_analysis.jl
#  Visor interactivo: src/longitudinal/plot_longitudinal.jl (CLI :8780)
#  Helpers figuras: src/viz/GroupVizCommon.jl + group_viewer_common.js
#  Regenerar PNG: scripts/regenerate_group_figures.jl

# ─── Transversal ───────────────────────────────────────────────
#  Análisis MS vs Control: scripts/run_transversal_analysis.jl
#  Visor interactivo: src/transversal/plot_transversal.jl (CLI :8781)

# ─── Visualización ─────────────────────────────────────────────
#  Figuras PNG/SVG (topomapas, heatmaps, espectros, grafos).
include("visualization/Topomaps.jl")
include("visualization/Heatmaps.jl")
include("visualization/Spectra.jl")
include("visualization/GraphPlots.jl")
include("visualization/ClinicalPlots.jl")
export plot_topomap, plot_connectivity_heatmap, plot_spectrum
export plot_spectrum_grid, plot_group_comparison
export plot_graph_metrics, plot_clinical_correlation, plot_surrogate_distribution
export save_figure

# ─── Informes HTML ─────────────────────────────────────────────
#  Reportes por sujeto, grupo y longitudinal.
include("report/HTMLReport.jl")
export generate_report, generate_group_report, generate_longitudinal_report


# ═══════════════════════════════════════════════════════════════
#  Orquestación
# ═══════════════════════════════════════════════════════════════

# ─── Pipeline de un solo sujeto ────────────────────────────────
#  Cadena activa: load_ss_config + run_single_subject_pipeline.
include("SingleSubjectPipeline.jl")
export load_ss_config, detect_first_subject, load_single_subject, validate_channels
export run_single_subject_pipeline

# ─── Dashboard web (carga DIFERIDA) ────────────────────────────
#  webapp/App.jl hace `using Genie` (framework web pesado). Para que
#  `using NeuroMIND` y el pipeline NO compilen Genie, App.jl NO se incluye
#  aquí: se carga la primera vez que se llama a launch_webapp (normalmente
#  desde scripts/launch_dashboard.jl). Así el arranque del pipeline es rápido.
const _DASHBOARD_LOADED = Ref(false)

"""
    launch_webapp(cfg; port=8080, open_browser=true)

Lanza el dashboard web (Genie). Carga perezosa: la primera invocación incluye
`webapp/App.jl` (compila Genie, puede tardar); las siguientes son directas.
"""
function launch_webapp(cfg; kwargs...)
    if !_DASHBOARD_LOADED[]
        @info "Cargando el dashboard (Genie) por primera vez — puede tardar…"
        Base.include(@__MODULE__, joinpath(@__DIR__, "webapp", "App.jl"))
        _DASHBOARD_LOADED[] = true
    end
    # Tras el include existe el método específico launch_webapp(::PipelineConfig; …),
    # más específico que este stub; invokelatest resuelve el world-age.
    return Base.invokelatest(launch_webapp, cfg; kwargs...)
end
export launch_webapp

end # module NeuroMIND
