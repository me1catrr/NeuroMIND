"""
NeuroMIND — EEG Functional Connectivity Framework
Estudio BRAIN: Conectividad funcional en Esclerosis Múltiple basada en wPLI
Rafael Castro Triguero, 2026

Uso básico:
    using NeuroMIND
    cfg      = load_config("config/pipeline.toml")
    subjects = load_subjects(cfg)
    run_pipeline!(subjects, cfg)

    # Lanzar visor web interactivo
    launch_webapp(cfg; port=8080)

    # Generar informe HTML de un sujeto
    generate_report(subjects[1], cfg)
"""
module NeuroMIND

using CSV, DataFrames, Dates, Statistics, LinearAlgebra, Random
using Serialization, TOML, FFTW, DSP, StatsBase

# ─── Tipos centrales ──────────────────────────────────────────
include("types.jl")
export PipelineConfig, RecordingMeta, EEGRecording, EpochSet
export ICAResult, SpectralResult, ConnectivityMatrix, SurrogateResult
export GraphMetrics, Session, Subject, GroupAnalysis, LongitudinalAnalysis
export ClinicalData, StatResult
export n_channels, n_samples, n_epochs, n_samples_epoch, duration

# ─── I/O y configuración ──────────────────────────────────────
include("io/Config.jl")
include("io/BIDSLoader.jl")
include("io/Serializer.jl")
export load_config, load_subjects, load_eeg_bids, save_result, load_result, result_exists
export results_dir, subject_results_dir, ensure_dirs

# ─── Quality Control ──────────────────────────────────────────
include("qc/QualityControl.jl")
export compute_channel_stats, flag_bad_channels, qc_report

# ─── Preprocessing ────────────────────────────────────────────
include("preprocessing/Filtering.jl")
export apply_highpass, apply_lowpass, apply_notch, apply_bandpass
export apply_bandreject, filter_recording, describe_filter_chain

# ─── ICA ──────────────────────────────────────────────────────
include("ica/ICACore.jl")
include("ica/ICAInspection.jl")
export run_ica, apply_ica_rejection, load_ica_labels

# ─── Segmentación ─────────────────────────────────────────────
include("segmentation/Epochs.jl")
export segment_recording, apply_baseline, reject_artifacts

# ─── Spectral ─────────────────────────────────────────────────
include("spectral/PowerSpectrum.jl")
export compute_psd, band_power

# ─── Conectividad ─────────────────────────────────────────────
include("connectivity/wPLI.jl")
include("connectivity/CSD.jl")
include("connectivity/GraphMetrics.jl")
export compute_wpli, apply_csd, compute_graph_metrics

# ─── Estadística ──────────────────────────────────────────────
include("statistics/Surrogates.jl")
include("statistics/FDR.jl")
include("statistics/GroupStats.jl")
export surrogate_test, fdr_correction, threshold_connectivity
export mann_whitney_test, wilcoxon_signed_rank, spearman_correlation
export compare_groups_stats

# ─── Análisis longitudinal ────────────────────────────────────
include("longitudinal/LongitudinalAnalysis.jl")
export compute_longitudinal, group_mean_connectivity, compare_groups

# ─── Visualización (figuras PNG/SVG) ─────────────────────────
include("visualization/Topomaps.jl")
include("visualization/Heatmaps.jl")
include("visualization/Spectra.jl")
include("visualization/GraphPlots.jl")
include("visualization/ClinicalPlots.jl")
export plot_topomap, plot_connectivity_heatmap, plot_spectrum
export plot_spectrum_grid, plot_group_comparison, plot_longitudinal_evolution
export plot_graph_metrics, plot_clinical_correlation, plot_surrogate_distribution
export save_figure

# ─── Generación de informes HTML ─────────────────────────────
include("report/HTMLReport.jl")
export generate_report, generate_group_report, generate_longitudinal_report

# ─── Pipeline de alto nivel ───────────────────────────────────
include("Pipeline.jl")
export run_pipeline!, run_subject!, run_session!

# ─── Pipeline de un solo sujeto ──────────────────────────────
include("SingleSubjectPipeline.jl")
export load_ss_config, detect_first_subject, load_single_subject, validate_channels
export run_single_subject_pipeline, load_dashboard_data

# ─── Aplicación web (Genie + Stipple) ────────────────────────
include("webapp/App.jl")
export launch_webapp

end # module NeuroMIND
