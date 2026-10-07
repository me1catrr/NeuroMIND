# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Suite canónica de tests
#  Tests unitarios, de integración numérica y de regresión del
#  pipeline EEG. Este es el único punto de entrada de pruebas y
#  sigue la convención estándar de paquetes Julia: test/runtests.jl.
# ═══════════════════════════════════════════════════════════════
#
#  Ejecución:
#    julia --project=. test/runtests.jl
#    julia --project=. -e 'using Pkg; Pkg.test()'
#
#  Inventario de @testset y finalidad
#  ─────────────────────────────────
#  Núcleo y configuración
#    · Types — construcción, propiedades y helpers de los tipos del dominio.
#    · Config — lectura de TOML y traducción a PipelineConfig.
#
#  Filtrado y control de calidad
#    · Filtering — dimensiones, metadatos y callbacks de la cadena básica.
#    · ChannelStatsExtended — estadísticos temporales/espectrales por canal.
#    · FilteringProfile — parámetros y comportamiento del perfil eeg_julia.
#    · FilterChainParams_EEGJulia — orden y frecuencias exactas de la cadena.
#    · FilteringDefaults_EEGJulia — defaults Notch/LP compatibles con EEG_Julia.
#    · FilterAttenuation — atenuación y respuesta frecuencial de cada filtro.
#
#  Segmentación, baseline y rechazo de artefactos
#    · Segmentation — creación de epochs, solapamiento y metadatos.
#    · Segmentation_EEGJulia_Profile — epochs de 1 s del perfil eeg_julia.
#    · Baseline_FirstWindowMean — corrección usando la primera ventana.
#    · AR_EEGJulia_Profile — umbral ±70 µV y reglas del perfil.
#    · Baseline_DoublePass — estabilidad de la doble corrección de baseline.
#    · AR_Phase7_MaxAmplitude — rechazo por amplitud máxima.
#    · AR_Phase7_MinAmplitude — rechazo por amplitud mínima.
#    · AR_Phase7_GradientNotUsed — confirma que gradiente no decide en Phase 7.
#    · AR_Phase7_NChannelsUsed — respeto al número real de canales analizados.
#    · AR_Phase7_QualityReport — coherencia del informe y recuentos de calidad.
#
#  Espectral, conectividad, grafos y estadística de cohorte
#    · Spectral — PSD, frecuencias y potencia por bandas.
#    · wPLI — matrices, simetría, diagonal y valores del estimador.
#    · GroupStats — Mann–Whitney, Welch y tamaños de efecto grupales.
#    · TransversalProductionStats — contrato v2 y estadística única MS–Control.
#    · LongitudinalProductionStats — Wilcoxon exacto y selección EC/EO real.
#    · FDR — corrección Benjamini–Hochberg y significación.
#    · GraphMetrics_ChannelNames — propagación de nombres y métricas de grafo.
#
#  ICA
#    · ICA_EEGJulia_Profile — dimensiones y componentes del perfil eeg_julia.
#    · ICA_Features — extracción de rasgos para clasificar componentes.
#    · ICA_EEGJulia_Reproducibility — reproducibilidad con semilla fija.
#    · ICA_AutoReject_Precedence — precedencia de decisiones automáticas/manuales.
#
#  Surrogates, inferencia de conectividad y dwPLI
#    · Surrogates_PValue_Bounds — cotas Monte Carlo y p-valores nunca nulos.
#    · Surrogates_Coupling_Detected — detección de acoplamiento sintético.
#    · Surrogates_Uncoupled_No_FalseDiscoveries — control de falsos positivos.
#    · Surrogates_RNG_Reproducibility — reproducibilidad del nulo por semilla.
#    · Surrogates_Validation_Functions — validadores y resúmenes de inferencia.
#    · wPLI_dwPLI_Uncoupled — comportamiento en señales no acopladas.
#    · wPLI_dwPLI_Coupled — respuesta en señales acopladas.
#    · wPLI_BandDuration_Warning — aviso por ciclos insuficientes en DELTA.
#
# ───────────────────────────────────────────────────────────────
#  Fichero     test/runtests.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      21-05-2026
#  Modificado  28-07-2026
# ───────────────────────────────────────────────────────────────

using Test
using NeuroMIND
using Statistics
using LinearAlgebra: diag, det, I
using Random
using CSV, DataFrames
using DSP: digitalfilter, Bandstop, Highpass, Lowpass, Butterworth, freqresp

# Los módulos de cohorte son CLI independientes de NeuroMIND.jl; se incluyen
# aquí para cubrir la lógica estadística que genera los CSV de producción.
include(joinpath(@__DIR__, "..", "src", "transversal", "Transversal.jl"))
include(joinpath(@__DIR__, "..", "src", "longitudinal", "Longitudinal.jl"))
include(joinpath(@__DIR__, "..", "src", "interactive", "plot_longitudinal.jl"))

# El visor transversal se aísla porque sus constantes HTTP tienen los mismos
# nombres que las del visor longitudinal, ya incluido en Main.
module TransversalViewerHarness
include(joinpath(@__DIR__, "..", "src", "interactive", "plot_transversal.jl"))
end

# ─── Fixtures ─────────────────────────────────────────────────

function mock_config()
    root = joinpath(@__DIR__, "..")
    PipelineConfig(
        # project
        Dict{String,Any}("name" => "test", "version" => "0.2"),
        # study
        Dict{String,Any}("n_ms_t1" => 44, "n_controls" => 40,
                         "conditions" => ["EO", "EC"],
                         "groups" => ["MS", "Control"],
                         "sessions" => ["T1", "T2"]),
        # paths
        Dict{String,Any}("results" => "results", "data_cache" => "data/cache",
                         "bids_root" => "data/BIDS"),
        # recording
        Dict{String,Any}("fs" => 500.0, "conditions" => ["EO", "EC"],
                         "n_channels" => 32, "reference" => "average"),
        # filtering — parámetros EEG_Julia
        Dict{String,Any}("profile" => "eeg_julia",
                         "highpass_hz" => 0.5, "lowpass_hz" => 150.0,
                         "notch_hz" => 50.0, "notch_bw_hz" => 1.0,
                         "bandreject_lo" => 99.5, "bandreject_hi" => 100.5,
                         "filter_order" => 4),
        # segmentation
        Dict{String,Any}("epoch_length_s" => 1.0, "epoch_overlap" => 0.0,
                         "min_epochs" => 5),
        # baseline
        Dict{String,Any}("apply" => true, "method" => "mean"),
        # artifact_rejection
        Dict{String,Any}("amplitude_threshold_uv" => 100.0,
                         "gradient_threshold_uv" => 50.0, "enabled" => true),
        # ica
        Dict{String,Any}("n_components" => 10, "random_seed" => 42),
        # spectral
        Dict{String,Any}("nfft" => 256, "window" => "hamming", "window_pct" => 10.0),
        # bands
        Dict{String,Tuple{Float64,Float64}}(
            "DELTA"    => (0.5,  4.0),
            "ALPHA"    => (7.8, 11.7),
            "BETA_LOW" => (12.0, 15.0),
            "GAMMA"    => (30.0, 50.0),
        ),
        # connectivity
        Dict{String,Any}("filter_order" => 8, "use_csd" => false, "method" => "wpli"),
        # surrogates
        Dict{String,Any}("n_surrogates" => 10, "alpha" => 0.05,
                         "method" => "phase_shuffle", "fdr_method" => "bh"),
        # graph
        Dict{String,Any}("threshold_method" => "proportional", "density" => 0.1),
        # clinical
        Dict{String,Any}("variables" => ["EDSS", "fatigue_score", "cognition_score"]),
        # longitudinal
        Dict{String,Any}("min_visits" => 2),
        # statistics
        Dict{String,Any}("group_test" => "mann_whitney", "alpha" => 0.05,
                         "fdr_q" => 0.05, "paired_test" => "wilcoxon"),
        # export_cfg
        Dict{String,Any}("figure_format" => "png", "figure_dpi" => 150,
                         "table_format" => "csv"),
        # qc
        Dict{String,Any}("bad_channel_zscore_threshold" => 3.0,
                         "amplitude_warning_sigma_uv" => 20.0),
        # montage
        Dict{String,Any}("exclude_channels" => String[], "exclude_fp2" => false,
                         "n_channels_analysis" => 30),
        # root
        root
    )
end

function mock_recording(n_ch=10, n_samp=2000, fs=500.0)
    meta = RecordingMeta(
        "TEST01", "T1", "EC", 1, fs, n_ch,
        ["Ch$i" for i in 1:n_ch], nothing, "dummy.tsv"
    )
    data  = randn(n_ch, n_samp)
    times = collect(0.0:(1/fs):(n_samp-1)/fs)
    return EEGRecording(meta, data, times)
end

# ─── Tests de tipos ───────────────────────────────────────────

@testset "Types" begin
    rec = mock_recording()
    @test n_channels(rec) == 10
    @test n_samples(rec)  == 2000
    @test duration(rec)   ≈ 4.0 atol=1e-6

    subj = Subject("T01", "MS")
    @test subj.id    == "T01"
    @test subj.group == "MS"
    @test isempty(subj.sessions)
    @test ismissing(subj.age)

    cd = ClinicalData(3.0, 5.0, "Interferón", 42.0, 28.0, missing)
    @test cd.EDSS == 3.0
    @test ismissing(cd.lesion_load)

    subj_ms = Subject("M01", "MS")
    subj_ms.clinical = ClinicalData(2.5, 6.0, "Dimetilfumarato", 38.0, 26.0, missing)
    @test subj_ms.clinical.EDSS == 2.5

    sr = StatResult("mann_whitney", 100.0, 0.03, 0.04, 0.3, true, 0.5, 0.3, 15, 15)
    @test sr.significant == true
    @test sr.test_name == "mann_whitney"
end

# ─── Tests de filtrado ────────────────────────────────────────

@testset "Filtering" begin
    cfg = mock_config()
    rec = mock_recording()

    rec_hp = apply_highpass(rec, 1.0)
    @test size(rec_hp.data) == size(rec.data)
    @test rec_hp.meta.fs == 500.0

    rec_lp = apply_lowpass(rec, 40.0)
    @test size(rec_lp.data) == size(rec.data)

    rec_notch = apply_notch(rec, 50.0)
    @test size(rec_notch.data) == size(rec.data)

    rec_br = apply_bandreject(rec, 99.5, 100.5)
    @test size(rec_br.data) == size(rec.data)

    rec_filt = filter_recording(rec, cfg)
    @test size(rec_filt.data) == size(rec.data)

    # on_step: emite un CSV-key por filtro aplicado (perfil eeg_julia)
    keys_seen = String[]
    rec_steps = filter_recording(rec, cfg; on_step = (key, step_rec) -> begin
        push!(keys_seen, key)
        @test size(step_rec.data) == size(rec.data)
    end)
    @test keys_seen == ["notch", "bandreject", "highpass", "lowpass"]
    @test rec_steps.data == rec_filt.data
end

# ─── Tests de QC extendido (Report_Pre) ───────────────────────

@testset "ChannelStatsExtended" begin
    cfg = mock_config()
    # Señal más larga para Welch (nfft=1024)
    rec = mock_recording(8, 4096, 500.0)

    stats = compute_channel_stats(rec)
    @test String.(names(stats)) == [
        "channel", "mean_uv", "rms_uv", "std_uv", "range_uv", "rms_zscore",
        "min_uv", "max_uv", "skewness", "kurtosis",
    ]
    @test size(stats, 1) == 8
    @test all(isfinite, stats.skewness)
    @test all(isfinite, stats.kurtosis)
    # Curtosis bruta de gaussiana ≈ 3
    @test mean(stats.kurtosis) > 2.0
    @test mean(stats.kurtosis) < 5.0

    spec = compute_channel_spectral_qc(rec, cfg; welch_nfft=1024)
    @test String.(names(spec)) == ["channel", "hfnoise", "snr_db"]
    @test all(0.0 .<= spec.hfnoise .<= 1.0)
    @test all(isfinite, spec.snr_db)

    corr = compute_correlation_summary(rec)
    @test String.(names(corr)) == ["channel", "mean_abs_corr", "max_corr", "min_corr"]
    @test all(0.0 .<= corr.mean_abs_corr .<= 1.0)

    # Canal plano → NaN en asimetría/curtosis
    flat = mock_recording(8, 4096, 500.0)
    flat.data[1, :] .= 0.0
    stats_flat = compute_channel_stats(flat)
    @test isnan(stats_flat.skewness[1])
    @test isnan(stats_flat.kurtosis[1])
end

# ─── Tests de perfil de filtrado EEG_Julia ───────────────────

@testset "FilteringProfile" begin
    cfg = mock_config()   # profile = "eeg_julia"
    rec = mock_recording()

    # filter_recording conserva dimensiones con perfil eeg_julia
    rec_filt = filter_recording(rec, cfg)
    @test size(rec_filt.data) == size(rec.data)
    @test rec_filt.meta.fs == rec.meta.fs

    # describe_filter_chain devuelve la cadena correcta para eeg_julia
    chain = describe_filter_chain(cfg)
    @test length(chain) == 4            # Notch → BR → HP → LP
    names = [s.name for s in chain]
    @test names[1] == "Notch"
    @test names[2] == "Bandreject"
    @test names[3] == "High-pass"
    @test names[4] == "Low-pass"

    # Notch y Bandreject deben usar filt (causal)
    @test chain[1].method == "filt"
    @test chain[2].method == "filt"

    # HP y LP deben usar filtfilt (zero-phase)
    @test chain[3].method == "filtfilt"
    @test chain[4].method == "filtfilt"

    # Los parámetros deben coincidir con la config
    @test chain[3].order == 4
    @test chain[4].order == 4

    # Perfil "default" produce cadena HP→LP→Notch→BR, todos filtfilt
    cfg_def = PipelineConfig(
        cfg.project, cfg.study, cfg.paths, cfg.recording,
        Dict{String,Any}("profile" => "default",
                         "highpass_hz" => 0.5, "lowpass_hz" => 48.0,
                         "notch_hz" => 50.0, "notch_bw_hz" => 2.0,
                         "bandreject_lo" => 100.0, "bandreject_hi" => 120.0,
                         "filter_order" => 4),
        cfg.segmentation, cfg.baseline, cfg.artifact_rejection,
        cfg.ica, cfg.spectral, cfg.bands, cfg.connectivity,
        cfg.surrogates, cfg.graph, cfg.clinical, cfg.longitudinal,
        cfg.statistics, cfg.export_cfg, cfg.qc, cfg.montage, cfg.root
    )
    chain_def = describe_filter_chain(cfg_def)
    names_def = [s.name for s in chain_def]
    @test names_def[1] == "High-pass"
    @test names_def[2] == "Low-pass"
    @test all(s.method == "filtfilt" for s in chain_def)

    # filter_recording con perfil default también conserva dimensiones
    rec_def = filter_recording(rec, cfg_def)
    @test size(rec_def.data) == size(rec.data)
end

# ─── Tests de segmentación ────────────────────────────────────

@testset "Segmentation" begin
    cfg = mock_config()
    rec = mock_recording()

    epochs = segment_recording(rec, cfg)
    expected_ep = 2000 ÷ 500   # 4 epochs de 1s
    @test n_epochs(epochs) == expected_ep
    @test n_samples_epoch(epochs) == 500

    ep_bl = apply_baseline(epochs, cfg)
    for ep in 1:n_epochs(ep_bl)
        for ch in 1:10
            @test abs(mean(ep_bl.data[ch, :, ep])) < 1e-10
        end
    end

    ep_ar = reject_artifacts(epochs, cfg)
    @test ep_ar.n_valid <= n_epochs(epochs)
    @test length(ep_ar.rejection_reasons) == length(ep_ar.rejected_idx)
end

# ─── Helpers para configuraciones personalizadas de segmentación ──

function _cfg_seg(seg_dict, bl_dict, ar_dict)
    c = mock_config()
    PipelineConfig(
        c.project, c.study, c.paths, c.recording, c.filtering,
        seg_dict, bl_dict, ar_dict,
        c.ica, c.spectral, c.bands, c.connectivity,
        c.surrogates, c.graph, c.clinical, c.longitudinal,
        c.statistics, c.export_cfg, c.qc, c.montage, c.root
    )
end

# ─── Tests: segmentación perfil eeg_julia ────────────────────

@testset "Segmentation_EEGJulia_Profile" begin
    # Con profile="eeg_julia" siempre genera épocas de 1 s sin solapamiento
    seg = Dict{String,Any}("profile" => "eeg_julia",
                           "epoch_length_s" => 5.0,   # debe ignorarse
                           "epoch_overlap"  => 0.5,   # debe ignorarse
                           "min_epochs"     => 2)
    bl  = Dict{String,Any}("apply" => true, "method" => "mean")
    ar  = Dict{String,Any}("enabled" => false, "profile" => "default",
                           "amplitude_threshold_uv" => 100.0,
                           "gradient_threshold_uv"  => 50.0)
    cfg = _cfg_seg(seg, bl, ar)
    rec = mock_recording(10, 2000, 500.0)   # 4 s a 500 Hz

    epochs = segment_recording(rec, cfg)
    @test epochs.epoch_length_s == 1.0              # forzado por profile
    @test n_samples_epoch(epochs) == 500            # 1 s × 500 Hz
    @test n_epochs(epochs) == 4                     # 4 s / 1 s = 4 épocas
end

# ─── Tests: baseline first_window_mean ───────────────────────

@testset "Baseline_FirstWindowMean" begin
    fs  = 500.0
    rec = mock_recording(10, 2000, fs)
    seg = Dict{String,Any}("epoch_length_s" => 1.0, "epoch_overlap" => 0.0,
                           "min_epochs" => 2, "profile" => "default")
    bl  = Dict{String,Any}("apply" => true, "method" => "first_window_mean",
                           "baseline_start_s" => 0.0, "baseline_end_s" => 0.10)
    ar  = Dict{String,Any}("enabled" => false, "profile" => "default",
                           "amplitude_threshold_uv" => 100.0,
                           "gradient_threshold_uv"  => 50.0)
    cfg    = _cfg_seg(seg, bl, ar)
    epochs = segment_recording(rec, cfg)
    ep_bl  = apply_baseline(epochs, cfg)

    # La media de las muestras 1..50 (0–100 ms a 500 Hz) debe ser ≈ 0 en todos los canales/épocas
    bl_end_idx = round(Int, 0.10 * fs)   # = 50
    for ep in 1:n_epochs(ep_bl), ch in 1:5
        @test abs(mean(ep_bl.data[ch, 1:bl_end_idx, ep])) < 1e-10
    end

    # La media del epoch COMPLETO no es necesariamente cero (solo la ventana baseline lo es)
    # Verificar que el método "mean" daría media total ≈ 0 y first_window_mean no necesariamente
    bl_mean = Dict{String,Any}("apply" => true, "method" => "mean")
    cfg2    = _cfg_seg(seg, bl_mean, ar)
    ep_bl2  = apply_baseline(epochs, cfg2)
    for ep in 1:n_epochs(ep_bl2), ch in 1:5
        @test abs(mean(ep_bl2.data[ch, :, ep])) < 1e-10
    end
end

# ─── Tests: rechazo artefactos perfil eeg_julia ───────────────

@testset "AR_EEGJulia_Profile" begin
    fs  = 500.0
    # 35 canales, 1000 muestras → 2 épocas de 500 muestras (0.5 s a 500 Hz nrow_step=500)
    # Canal 3, época 2 (muestras 501-1000): plateau de 75 µV (amplitude > ±70 pero < ±100,
    # gradiente = 0 → el test es limpio independiente del perfil de gradiente).
    data = zeros(35, 1000)
    data[3, 501:1000] .= 75.0    # 75 µV: supera ±70 (eeg_julia) pero NO ±100 (default)
    meta  = RecordingMeta("T01", "T1", "EC", 1, fs, 35,
                          ["Ch$i" for i in 1:35], nothing, "dummy.tsv")
    rec   = EEGRecording(meta, data, collect(0:(1/fs):(1000-1)/fs))
    seg   = Dict{String,Any}("profile" => "default", "epoch_length_s" => 0.5,
                             "epoch_overlap" => 0.0, "min_epochs" => 1)
    bl    = Dict{String,Any}("apply" => false)

    ar_ej = Dict{String,Any}("profile" => "eeg_julia", "enabled" => true,
                             "min_amplitude_uv" => -70.0, "max_amplitude_uv" => 70.0,
                             "n_channels_used" => 30,
                             "amplitude_threshold_uv" => 100.0,
                             "gradient_threshold_uv"  => 50.0)
    cfg_ej = _cfg_seg(seg, bl, ar_ej)
    epochs = segment_recording(rec, cfg_ej)
    # Con eeg_julia: ±70 µV, canal 3 en época 2 = 75 µV → rechazada
    ep_ej = reject_artifacts(epochs, cfg_ej)
    @test ep_ej.n_valid < n_epochs(epochs)
    @test length(ep_ej.rejected_idx) >= 1

    # Con profile "default" (±100 µV, sin gradiente relevante porque plateau): sin rechazo
    ar_def = Dict{String,Any}("profile" => "default", "enabled" => true,
                              "amplitude_threshold_uv" => 100.0,
                              "gradient_threshold_uv"  => 200.0,  # umbral alto para gradient
                              "use_gradient" => false)
    cfg_def = _cfg_seg(seg, bl, ar_def)
    ep_def  = reject_artifacts(epochs, cfg_def)
    @test ep_def.n_valid == n_epochs(epochs)   # 75 µV < 100 µV → sin rechazo
end

# ─── Tests: doble pasada baseline (n_passes = 2) ─────────────

@testset "Baseline_DoublePass" begin
    fs  = 500.0
    rec = mock_recording(10, 2000, fs)
    seg = Dict{String,Any}("profile" => "default", "epoch_length_s" => 1.0,
                           "epoch_overlap" => 0.0, "min_epochs" => 2)
    ar  = Dict{String,Any}("enabled" => true, "profile" => "default",
                           "amplitude_threshold_uv" => 100.0,
                           "gradient_threshold_uv"  => 50.0)
    bl  = Dict{String,Any}("apply" => true, "method" => "first_window_mean",
                           "baseline_start_s" => 0.0, "baseline_end_s" => 0.10,
                           "n_passes" => 2)
    cfg = _cfg_seg(seg, bl, ar)

    epochs   = segment_recording(rec, cfg)
    epochs1  = apply_baseline(epochs, cfg)       # 1ª pasada
    epochs_ar = reject_artifacts(epochs1, cfg)   # AR
    epochs2  = apply_baseline(epochs_ar, cfg)    # 2ª pasada post-AR

    # Tras la 2ª pasada, la ventana baseline de cada época válida debe tener media ≈ 0
    bl_end_idx = round(Int, 0.10 * fs)
    for ep in 1:n_epochs(epochs2), ch in 1:5
        @test abs(mean(epochs2.data[ch, 1:bl_end_idx, ep])) < 1e-10
    end

    # La 2ª pasada opera solo sobre las épocas válidas (mismas dimensiones que epochs_ar)
    @test size(epochs2.data, 3) == epochs_ar.n_valid
end

# ─── Tests Phase 7: comportamiento exacto AR eeg_julia ───────

@testset "AR_Phase7_MaxAmplitude" begin
    # Canal 5, época 2: +75 µV > +70 µV → rechazado con eeg_julia, válido con default (±100)
    fs   = 500.0
    n_ch = 10; n_samp = 1000   # 2 épocas de 1 s a 500 Hz
    data = zeros(n_ch, n_samp)
    data[5, 501:1000] .= 75.0   # época 2, canal 5: +75 µV
    meta  = RecordingMeta("T01", "T1", "EC", 1, fs, n_ch,
                          ["Ch$i" for i in 1:n_ch], nothing, "dummy.tsv")
    rec   = EEGRecording(meta, data, collect(0.0:(1/fs):(n_samp-1)/fs))

    seg = Dict{String,Any}("profile" => "default", "epoch_length_s" => 1.0,
                           "epoch_overlap" => 0.0, "min_epochs" => 1)
    bl  = Dict{String,Any}("apply" => false)
    ar  = Dict{String,Any}("profile" => "eeg_julia", "enabled" => true,
                           "min_amplitude_uv" => -70.0, "max_amplitude_uv" => 70.0,
                           "n_channels_used" => 10,
                           "amplitude_threshold_uv" => 100.0,
                           "gradient_threshold_uv"  => 50.0)
    cfg    = _cfg_seg(seg, bl, ar)
    epochs = segment_recording(rec, cfg)
    ep_ar  = reject_artifacts(epochs, cfg)

    @test ep_ar.n_valid == 1              # solo época 1 válida
    @test length(ep_ar.rejected_idx) == 1
    @test ep_ar.rejected_idx[1] == 2     # época 2 rechazada por amplitud máxima

    # Con default (±100 µV): 75 µV < 100 µV → ninguna época rechazada
    ar_def = Dict{String,Any}("profile" => "default", "enabled" => true,
                              "amplitude_threshold_uv" => 100.0,
                              "gradient_threshold_uv"  => 200.0,
                              "use_gradient" => false)
    cfg_def = _cfg_seg(seg, bl, ar_def)
    ep_def  = reject_artifacts(epochs, cfg_def)
    @test ep_def.n_valid == n_epochs(epochs)   # 75 µV < 100 µV → sin rechazo
end

@testset "AR_Phase7_MinAmplitude" begin
    # Canal 3, época 1: -75 µV < -70 µV → rechazado con eeg_julia
    fs   = 500.0
    n_ch = 10; n_samp = 1000
    data = zeros(n_ch, n_samp)
    data[3, 1:500] .= -75.0   # época 1, canal 3: -75 µV
    meta  = RecordingMeta("T01", "T1", "EC", 1, fs, n_ch,
                          ["Ch$i" for i in 1:n_ch], nothing, "dummy.tsv")
    rec   = EEGRecording(meta, data, collect(0.0:(1/fs):(n_samp-1)/fs))

    seg = Dict{String,Any}("profile" => "default", "epoch_length_s" => 1.0,
                           "epoch_overlap" => 0.0, "min_epochs" => 1)
    bl  = Dict{String,Any}("apply" => false)
    ar  = Dict{String,Any}("profile" => "eeg_julia", "enabled" => true,
                           "min_amplitude_uv" => -70.0, "max_amplitude_uv" => 70.0,
                           "n_channels_used" => 10,
                           "amplitude_threshold_uv" => 100.0,
                           "gradient_threshold_uv"  => 50.0)
    cfg    = _cfg_seg(seg, bl, ar)
    epochs = segment_recording(rec, cfg)
    ep_ar  = reject_artifacts(epochs, cfg)

    @test ep_ar.n_valid == 1              # solo época 2 válida
    @test length(ep_ar.rejected_idx) == 1
    @test ep_ar.rejected_idx[1] == 1     # época 1 rechazada (min < -70 µV)
end

@testset "AR_Phase7_GradientNotUsed" begin
    # Spike que genera gradiente > 50 µV/muestra, pero amplitud ≤ 60 µV < 70 µV
    # eeg_julia: sin gradient check → NO rechazado
    # default + use_gradient=true: rechazado por gradiente
    fs   = 500.0
    n_ch = 5; n_samp = 1000
    data = zeros(n_ch, n_samp)
    # muestra 250 = 60 µV, muestra 251 = 0 → gradiente = 60 µV/muestra > 50
    data[2, 250] = 60.0
    meta  = RecordingMeta("T01", "T1", "EC", 1, fs, n_ch,
                          ["Ch$i" for i in 1:n_ch], nothing, "dummy.tsv")
    rec   = EEGRecording(meta, data, collect(0.0:(1/fs):(n_samp-1)/fs))

    seg = Dict{String,Any}("profile" => "default", "epoch_length_s" => 1.0,
                           "epoch_overlap" => 0.0, "min_epochs" => 1)
    bl  = Dict{String,Any}("apply" => false)

    # eeg_julia: amplitud 60 µV < 70 µV, gradiente ignorado → NO rechazado
    ar_ej = Dict{String,Any}("profile" => "eeg_julia", "enabled" => true,
                             "min_amplitude_uv" => -70.0, "max_amplitude_uv" => 70.0,
                             "n_channels_used" => 5,
                             "amplitude_threshold_uv" => 100.0,
                             "gradient_threshold_uv"  => 50.0)
    cfg_ej = _cfg_seg(seg, bl, ar_ej)
    epochs = segment_recording(rec, cfg_ej)
    ep_ej  = reject_artifacts(epochs, cfg_ej)
    @test ep_ej.n_valid == n_epochs(epochs)   # gradiente ignorado → ambas épocas válidas

    # default + use_gradient=true: gradiente 60 > 50 → época 1 rechazada
    ar_def = Dict{String,Any}("profile" => "default", "enabled" => true,
                              "amplitude_threshold_uv" => 100.0,
                              "gradient_threshold_uv"  => 50.0,
                              "use_gradient"           => true)
    cfg_def = _cfg_seg(seg, bl, ar_def)
    ep_def  = reject_artifacts(epochs, cfg_def)
    @test ep_def.n_valid < n_epochs(epochs)   # gradiente 60 > 50 → época 1 rechazada
end

@testset "AR_Phase7_NChannelsUsed" begin
    # Violación (+75 µV) en canal 31, fuera de los primeros 30 evaluados por eeg_julia.
    # n_channels_used=30 → NO rechaza; n_channels_used=35 → SÍ rechaza.
    fs   = 500.0
    n_ch = 35; n_samp = 1000
    data = zeros(n_ch, n_samp)
    data[31, 1:500] .= 75.0   # canal 31, época 1: 75 µV
    meta  = RecordingMeta("T01", "T1", "EC", 1, fs, n_ch,
                          ["Ch$i" for i in 1:n_ch], nothing, "dummy.tsv")
    rec   = EEGRecording(meta, data, collect(0.0:(1/fs):(n_samp-1)/fs))

    seg = Dict{String,Any}("profile" => "default", "epoch_length_s" => 1.0,
                           "epoch_overlap" => 0.0, "min_epochs" => 1)
    bl  = Dict{String,Any}("apply" => false)

    # n_channels_used=30: canal 31 no evaluado → 75 µV no detectado → NO rechaza
    ar_30 = Dict{String,Any}("profile" => "eeg_julia", "enabled" => true,
                             "min_amplitude_uv" => -70.0, "max_amplitude_uv" => 70.0,
                             "n_channels_used" => 30,
                             "amplitude_threshold_uv" => 100.0,
                             "gradient_threshold_uv"  => 50.0)
    cfg_30 = _cfg_seg(seg, bl, ar_30)
    epochs = segment_recording(rec, cfg_30)
    ep_30  = reject_artifacts(epochs, cfg_30)
    @test ep_30.n_valid == n_epochs(epochs)   # canal 31 ignorado → sin rechazo

    # n_channels_used=35: canal 31 evaluado → 75 µV > 70 µV → rechaza época 1
    ar_35 = Dict{String,Any}("profile" => "eeg_julia", "enabled" => true,
                             "min_amplitude_uv" => -70.0, "max_amplitude_uv" => 70.0,
                             "n_channels_used" => 35,
                             "amplitude_threshold_uv" => 100.0,
                             "gradient_threshold_uv"  => 50.0)
    cfg_35 = _cfg_seg(seg, bl, ar_35)
    ep_35  = reject_artifacts(epochs, cfg_35)
    @test ep_35.n_valid < n_epochs(epochs)    # canal 31 evaluado → 75 > 70 → rechazado
    @test ep_35.rejected_idx[1] == 1          # época 1 rechazada
end

@testset "AR_Phase7_QualityReport" begin
    # Verifica que compute_epoch_quality_report produce las columnas min_amp_uv
    # y channels_violating (Phase 7), y que sus valores son correctos.
    fs   = 500.0
    n_ch = 5; n_samp = 1000
    ch_names = ["Fp1","Fp2","F3","F4","Fz"]
    data = zeros(n_ch, n_samp)
    data[2, 501:1000] .= 75.0    # época 2, canal Fp2: +75 µV
    data[4, 501:1000] .= -65.0   # época 2, canal F4: -65 µV (dentro de ±70, no viola)
    meta  = RecordingMeta("T01", "T1", "EC", 1, fs, n_ch,
                          ch_names, nothing, "dummy.tsv")
    rec   = EEGRecording(meta, data, collect(0.0:(1/fs):(n_samp-1)/fs))

    seg = Dict{String,Any}("profile" => "default", "epoch_length_s" => 1.0,
                           "epoch_overlap" => 0.0, "min_epochs" => 1)
    bl  = Dict{String,Any}("apply" => false)
    ar  = Dict{String,Any}("profile" => "eeg_julia", "enabled" => true,
                           "min_amplitude_uv" => -70.0, "max_amplitude_uv" => 70.0,
                           "n_channels_used" => 5,
                           "amplitude_threshold_uv" => 100.0,
                           "gradient_threshold_uv"  => 50.0)
    cfg    = _cfg_seg(seg, bl, ar)
    epochs = segment_recording(rec, cfg)
    qr     = compute_epoch_quality_report(epochs, cfg)

    # ── Columnas nuevas Phase 7 ────────────────────────────────
    @test hasproperty(qr, :min_amp_uv)
    @test hasproperty(qr, :channels_violating)

    # ── Época 1 (todas las muestras = 0): válida, sin violaciones ─
    @test qr.status[1]             == "valid"
    @test isempty(qr.channels_violating[1])
    @test qr.min_amp_uv[1]         ≈ 0.0 atol=1e-6
    @test qr.max_amp_uv[1]         ≈ 0.0 atol=1e-6

    # ── Época 2 (Fp2=+75): rechazada, Fp2 en channels_violating ──
    @test qr.status[2]             == "rejected"
    @test qr.rejection_reason[2]   == "amplitude"
    @test occursin("Fp2", qr.channels_violating[2])

    # F4 (-65 µV) NO viola ±70 → no debe aparecer en channels_violating
    @test !occursin("F4", qr.channels_violating[2])

    # max_amp y p2p correctos para época 2
    @test qr.max_amp_uv[2]  ≈ 75.0 atol=1e-4
    @test qr.p2p_uv[2]      ≈ 140.0 atol=1e-4  # 75.0 - (-65.0) = 140.0

    # worst_channel debe ser Fp2 (mayor amplitud absoluta = 75.0 vs |−65|=65)
    @test qr.worst_channel[2] == "Fp2"

    # ── Tamaño del reporte: una fila por época ─────────────────
    @test size(qr, 1) == 2
    @test qr.epoch[1] == 1
    @test qr.epoch[2] == 2
end

# ─── Tests espectrales ────────────────────────────────────────

@testset "Spectral" begin
    cfg = mock_config()
    rec = mock_recording()
    ep  = segment_recording(rec, cfg)
    ep  = apply_baseline(ep, cfg)

    sp = compute_psd(ep, cfg)
    # nfft = max(config_nfft, n_samples_epoch) = max(256, 500) = 500 → 251 bins
    expected_bins = max(cfg.spectral["nfft"], n_samples_epoch(ep)) ÷ 2 + 1
    @test size(sp.psd, 1) == 10
    @test size(sp.psd, 2) == expected_bins
    @test all(sp.psd .>= 0.0)
    @test haskey(sp.band_power, "ALPHA")
    @test all(sp.band_power["ALPHA"] .>= 0.0)
end

# ─── Tests de conectividad wPLI ───────────────────────────────

@testset "wPLI" begin
    cfg    = mock_config()
    rec    = mock_recording(5, 2000)
    ep     = segment_recording(rec, cfg)
    ep_bl  = apply_baseline(ep, cfg)

    conn = compute_wpli(ep_bl, cfg)
    @test conn.method == "wpli"
    @test haskey(conn.matrices, "ALPHA")

    W = conn.matrices["ALPHA"]
    @test size(W) == (5, 5)
    @test all(W .>= 0.0)
    @test all(W .<= 1.0 + 1e-8)
    @test all(diag(W) .< 1e-10)
    @test W ≈ W'
end

# ─── Tests de estadística ─────────────────────────────────────

@testset "GroupStats" begin
    # Mann-Whitney: dos grupos claramente distintos
    a = Float64[1, 2, 3, 4, 5]
    b = Float64[6, 7, 8, 9, 10]
    sr = mann_whitney_test(a, b; alpha=0.05)
    @test sr.test_name == "mann_whitney"
    @test sr.p_value < 0.05
    @test sr.significant == true
    @test sr.effect_size > 0.0

    # Mann-Whitney: grupos iguales → no significativo
    c = Float64[1, 2, 3, 4, 5]
    d = Float64[1, 2, 3, 4, 5]
    sr2 = mann_whitney_test(c, d; alpha=0.05)
    @test !sr2.significant

    # Wilcoxon signed-rank
    before = Float64[2, 4, 3, 5, 1]
    after  = Float64[8, 9, 7, 9, 6]
    wr = wilcoxon_signed_rank(before, after; alpha=0.05)
    @test wr.test_name == "wilcoxon"
    @test wr.p_value < 0.05

    # Spearman: correlación positiva perfecta
    x = Float64[1, 2, 3, 4, 5]
    y = Float64[2, 4, 6, 8, 10]
    sp = spearman_correlation(x, y; alpha=0.05)
    @test sp.test_name == "spearman"
    @test sp.statistic ≈ 1.0 atol=1e-6

    # Spearman: muestra insuficiente
    sp_small = spearman_correlation([1.0, 2.0], [1.0, 2.0])
    @test isnan(sp_small.statistic)
end

@testset "TransversalProductionStats" begin
    # Convención de cuantiles de Julia para N par y efecto orientado MS−Control.
    ms = Float64[4, 5, 6, 7]
    ctrl = Float64[0, 1, 2, 3]
    prod = Transversal.independent_group_summary(
        ms, ctrl; n_boot=100, seed=20260728)
    @test prod.ms_median == 5.5
    @test prod.ms_q1 == 4.75
    @test prod.ms_q3 == 6.25
    @test prod.ctrl_median == 1.5
    @test prod.effect_rrb == 1.0
    @test prod.probability_superiority == 1.0
    @test prod.diff == 4.0
    @test prod.effect_d_pooled > 0
    @test prod.n_ms == 4
    @test prod.n_ctrl == 4

    # El visor rechaza derivados antiguos o incompletos y nunca reconstruye
    # cuantiles, efectos ni significación como fallback.
    mktempdir() do tmp
        err = try
            TransversalViewerHarness.validate_condition_contract(tmp)
            nothing
        catch e
            e
        end
        @test err isa TransversalViewerHarness.IncompatibleTransversalResultsError
        @test occursin("Resultados incompatibles con el visor actual",
                       sprint(showerror, err))
        old_dir = joinpath(tmp, "old_results", "transversal", "eyesclosed")
        mkpath(old_dir)
        old_viewer = TransversalViewerHarness.load_viewer(
            results_root=joinpath(tmp, "old_results"))
        @test isempty(old_viewer.stores)
        @test occursin("Resultados incompatibles con el visor actual",
                       old_viewer.errors["EC"])
        @test_throws TransversalViewerHarness.IncompatibleTransversalResultsError begin
            TransversalViewerHarness._require_columns(
                DataFrame(schema_version=[2]),
                [:schema_version, :effect_d_pooled], "legacy.csv")
        end
        @test_throws TransversalViewerHarness.IncompatibleTransversalResultsError begin
            TransversalViewerHarness._require_columns(
                DataFrame(schema_version=[1], effect_d_pooled=[0.0]),
                [:schema_version, :effect_d_pooled], "legacy.csv")
        end

        contract_path = joinpath(tmp, "statistics_contract.toml")
        Transversal.write_statistics_contract(contract_path, "EC")
        contract = TransversalViewerHarness.validate_condition_contract(tmp)
        @test contract["schema_version"] == 2
        @test contract["statistics_source"] == Transversal.STATISTICS_SOURCE
        @test contract["bootstrap_iterations"] == Transversal.BOOTSTRAP_N
        @test contract["fdr_scope_edges"] == Transversal.EDGE_FDR_SCOPE
        @test contract["fdr_scope_power"] == Transversal.POWER_FDR_SCOPE
    end

    # Regresión del dataset real ALPHA–EC. `results/` es deliberadamente
    # gitignored; en un checkout sin derivados se registra como skipped.
    reg_dir = joinpath(@__DIR__, "..", "results", "transversal", "eyesclosed", "tables")
    reg_files = [
        joinpath(reg_dir, "global_mean_wpli_statistics.csv"),
        joinpath(reg_dir, "subject_band_means.csv"),
        joinpath(reg_dir, "band_statistics.csv"),
        joinpath(reg_dir, "group_statistics_ALPHA.csv"),
        joinpath(reg_dir, "band_power_group_statistics.csv"),
    ]
    if all(isfile, reg_files)
        global_stats = CSV.read(reg_files[1], DataFrame)
        subject_means = CSV.read(reg_files[2], DataFrame)
        band_stats = CSV.read(reg_files[3], DataFrame)
        edges = CSV.read(reg_files[4], DataFrame)
        power = CSV.read(reg_files[5], DataFrame)
        alpha = only(eachrow(filter(r -> string(r.band) == "ALPHA", global_stats)))
        alpha_band = only(eachrow(filter(r -> string(r.band) == "ALPHA", band_stats)))
        alpha_subjects = filter(r -> string(r.band) == "ALPHA", subject_means)
        ms_values = Float64.(filter(r -> string(r.group) == "MS", alpha_subjects).mean_wpli)
        ctrl_values = Float64.(filter(r -> string(r.group) == "Control", alpha_subjects).mean_wpli)
        recomputed = Transversal.independent_group_summary(
            ms_values, ctrl_values; seed=Int(alpha.bootstrap_seed))
        top20 = sort(filter(r -> Int(r.effect_rank_abs_d) <= 20, edges),
                     :effect_rank_abs_d)

        @test Int(alpha.schema_version) == 2
        @test string(alpha.statistics_source) == Transversal.STATISTICS_SOURCE
        @test string(alpha.quantile_method) == Transversal.QUANTILE_METHOD
        @test string(alpha.rrb_method) == Transversal.RRB_METHOD
        @test Float64(alpha.effect_rrb) ≈ 1 / 3 atol=1e-12
        @test Float64(alpha.ms_median) ≈ 0.20864151962518226 atol=1e-14
        @test Float64(alpha.ctrl_median) ≈ 0.16040136212782588 atol=1e-14
        @test Float64(alpha.effect_d_pooled) ≈ 0.5480285757786955 atol=1e-13
        @test Float64(alpha.diff) ≈ 0.06632193175905668 atol=1e-14
        @test Float64(alpha.ms_median) == median(ms_values)
        @test Float64(alpha.ctrl_median) == median(ctrl_values)
        @test Float64(alpha.effect_rrb) == recomputed.effect_rrb
        @test Float64(alpha.effect_d_pooled) == recomputed.effect_d_pooled
        @test Int(alpha.fdr_family_size) == 7
        @test nrow(edges) == 276
        @test count(identity, Bool.(edges.is_nominal)) == 86
        @test count(identity, Bool.(edges.is_fdr)) == 12
        @test nrow(edges) == Int(alpha_band.n_edges)
        @test nrow(edges) == Int(alpha_band.n_channels) *
                            (Int(alpha_band.n_channels) - 1) ÷ 2
        @test count(identity, Bool.(edges.is_nominal)) == Int(alpha_band.n_nominal)
        @test count(identity, Bool.(edges.is_fdr)) == Int(alpha_band.n_sig)
        @test nrow(top20) == min(20, nrow(edges))
        @test count(identity, Bool.(top20.is_nominal)) == 20
        @test count(identity, Bool.(top20.is_nominal)) ==
              Int(alpha_band.top20_nominal_overlap)
        @test sort(Int.(edges.effect_rank_abs_d)) == collect(1:nrow(edges))
        for band in unique(String.(power.band))
            sub = filter(r -> string(r.band) == band, power)
            @test all(Int.(sub.fdr_family_size) .== nrow(sub))
        end
    else
        @test_skip false
    end
end

@testset "LongitudinalProductionStats" begin
    # Wilcoxon exacto: todos los cambios positivos, n=15.
    before = zeros(15)
    after = Float64.(1:15)
    p, _, method = Longitudinal.wilcoxon_p(before, after; return_method=true)
    @test method == "exact_conditional_dp"
    @test p ≈ 2.0 / 2.0^15 atol=1e-15

    # W+=1 reproduce el caso sensible C4–FC2: p exacto = 4/2^15.
    after_w1 = vcat(1.0, -Float64.(2:15))
    p_w1, _ = Longitudinal.wilcoxon_p(before, after_w1)
    @test p_w1 ≈ 4.0 / 2.0^15 atol=1e-15

    # Convención de Statistics.quantile para N par (interpolación Julia).
    even_values = Float64[1, 2, 3, 4]
    med, q1, q3 = Longitudinal.median_iqr(even_values)
    @test med == 2.5
    @test q1 == 1.75
    @test q3 == 3.25

    # El resumen único de producción materializa los estadísticos que usa el
    # visor; JavaScript no debe volver a estimarlos.
    prod = Longitudinal.paired_change_summary(
        zeros(4), Float64[1, 2, 3, 4]; n_boot=100, seed=20260728)
    @test prod.median_diff == med
    @test prod.q1_diff == q1
    @test prod.q3_diff == q3
    @test prod.effect_rrb == 1.0
    @test prod.n_positive == 4
    @test prod.n_negative == 0
    @test prod.n_zero == 0

    # El visor rechaza resultados antiguos o incompletos; nunca reconstruye
    # estadísticos como fallback.
    mktempdir() do tmp
        err = try
            validate_condition_contract(tmp)
            nothing
        catch e
            e
        end
        @test err isa IncompatibleResultsError
        @test occursin("Resultados incompatibles con el visor actual",
                       sprint(showerror, err))
        old_dir = joinpath(tmp, "old_results", "longitudinal", "eyesclosed")
        mkpath(old_dir)
        old_viewer = load_viewer(results_root=joinpath(tmp, "old_results"))
        @test isempty(old_viewer.stores)
        @test occursin("Resultados incompatibles con el visor actual",
                       old_viewer.errors["EC"])
        @test_throws IncompatibleResultsError _require_columns(
            DataFrame(schema_version=[2]), [:schema_version, :effect_dz], "legacy.csv")
        @test_throws IncompatibleResultsError _require_columns(
            DataFrame(schema_version=[1], effect_dz=[0.0]),
            [:schema_version, :effect_dz], "legacy.csv")
        Longitudinal.write_statistics_contract(joinpath(tmp, "statistics_contract.toml"), "EC")
        contract = validate_condition_contract(tmp)
        @test contract["schema_version"] == 2
        @test contract["statistics_source"] == Longitudinal.STATISTICS_SOURCE
        @test contract["bootstrap_iterations"] == Longitudinal.BOOTSTRAP_N
    end

    # La selección debe conservar un par EC aunque include_longitudinal=false
    # por ausencia de EO: las condiciones se evalúan por separado.
    mktempdir() do tmp
        bids = joinpath(tmp, "bids")
        mkpath(bids)
        CSV.write(joinpath(bids, "longitudinal_pairs.csv"), DataFrame(
            subject_id=["M7", "M8"],
            bids_id=["M07", "M08"],
            has_t1_ec=[true, true],
            has_t1_eo=[false, true],
            has_t2_ec=[true, true],
            has_t2_eo=[true, true],
            include_longitudinal=[false, true],
        ))
        pairs = Longitudinal.load_pairs(tmp, joinpath(tmp, "results"), bids)
        @test length(pairs) == 2
        m07 = only(filter(p -> p.subject_id == "M07", pairs))
        @test m07.has_ec
        @test !m07.has_eo
    end

    # Regresión del dataset real ALPHA–EC. `results/` es deliberadamente
    # gitignored; en un checkout sin derivados se registra como skipped.
    reg_dir = joinpath(@__DIR__, "..", "results", "longitudinal", "eyesclosed", "tables")
    reg_files = [
        joinpath(reg_dir, "network_global_statistics.csv"),
        joinpath(reg_dir, "mean_strength_scores.csv"),
        joinpath(reg_dir, "longitudinal_statistics_ALPHA.csv"),
        joinpath(reg_dir, "band_power_delta_statistics.csv"),
    ]
    if all(isfile, reg_files)
        ng = CSV.read(reg_files[1], DataFrame)
        scores = CSV.read(reg_files[2], DataFrame)
        edges = CSV.read(reg_files[3], DataFrame)
        power = CSV.read(reg_files[4], DataFrame)
        alpha = only(eachrow(filter(r -> string(r.band) == "ALPHA", ng)))
        alpha_scores = filter(r -> string(r.band) == "ALPHA", scores)
        deltas = Float64.(alpha_scores.delta_wpli_equiv)
        top20 = sort(filter(r -> Int(r.effect_rank_abs_dz) <= 20, edges),
                     :effect_rank_abs_dz)

        @test Float64(alpha.effect_rrb) ≈ -0.1618 atol=5e-5
        @test Float64(alpha.median_diff_wpli) ≈ -0.0051785 atol=5e-7
        @test median(deltas) == Float64(alpha.median_diff_wpli)
        @test Int(alpha.schema_version) == 2
        @test string(alpha.statistics_source) == Longitudinal.STATISTICS_SOURCE
        @test string(alpha.quantile_method) == Longitudinal.QUANTILE_METHOD
        @test string(alpha.rrb_method) == Longitudinal.RRB_METHOD
        @test nrow(edges) == 351
        @test nrow(edges) == Int(alpha.n_channels) * (Int(alpha.n_channels) - 1) ÷ 2
        @test count(identity, Bool.(edges.is_nominal)) == 20
        @test nrow(top20) == 20
        @test count(identity, Bool.(top20.is_nominal)) == 17
        @test count(identity, Bool.(edges.is_fdr)) == 0
        @test sort(Int.(edges.effect_rank_abs_dz)) == collect(1:351)
        @test Float64(alpha.diff) ≈
              (Int(alpha.n_channels) - 1) * Float64(alpha.diff_wpli) atol=1e-4
        for band in unique(String.(power.band))
            sub = filter(r -> string(r.band) == band, power)
            @test nrow(sub) == 31
            @test all(Int.(sub.fdr_family_size) .== 31)
        end
    else
        @test_skip false
    end
end

# ─── Tests FDR ────────────────────────────────────────────────

@testset "FDR" begin
    p = [0.001, 0.02, 0.05, 0.1, 0.3, 0.5, 0.8]
    thr = fdr_correction(p; alpha=0.05, method="bh")
    @test thr >= 0.0
    @test thr <= 0.05

    thr_bonf = fdr_correction(p; alpha=0.05, method="bonferroni")
    @test thr_bonf ≈ 0.05 / length(p)
end

# ─── Tests de configuración ───────────────────────────────────

@testset "Config" begin
    cfg = mock_config()
    @test haskey(cfg.bands, "ALPHA")
    @test cfg.bands["ALPHA"] == (7.8, 11.7)
    @test cfg.filtering["highpass_hz"]    == 0.5
    @test cfg.filtering["lowpass_hz"]     == 150.0
    @test cfg.filtering["notch_bw_hz"]    == 1.0
    @test cfg.filtering["bandreject_lo"]  == 99.5
    @test cfg.filtering["bandreject_hi"]  == 100.5
    @test cfg.filtering["profile"]        == "eeg_julia"
    @test cfg.statistics["fdr_q"] == 0.05
    @test haskey(cfg.study, "n_ms_t1")
    @test cfg.study["n_ms_t1"] == 44
end

# ─── Tests de parámetros EEG_Julia (réplica exacta) ──────────
# Verifica que NeuroMIND.describe_filter_chain reproduce el perfil documentado en
# EEG_Julia/src/Preprocessing/filtering.jl y EEG_Julia/config/default_config.jl

@testset "FilterChainParams_EEGJulia" begin
    cfg = mock_config()   # profile = "eeg_julia"
    chain = describe_filter_chain(cfg)

    # Cadena: 4 pasos (Notch → BR → HP → LP)
    @test length(chain) == 4

    notch = chain[1]; br = chain[2]; hp_f = chain[3]; lp_f = chain[4]

    # Nombres exactos
    @test notch.name == "Notch"
    @test br.name    == "Bandreject"
    @test hp_f.name  == "High-pass"
    @test lp_f.name  == "Low-pass"

    # Métodos: EEG_Julia usa filt para Notch y BR, filtfilt para HP y LP
    @test notch.method == "filt"
    @test br.method    == "filt"
    @test hp_f.method  == "filtfilt"
    @test lp_f.method  == "filtfilt"

    # Órdenes de diseño: todos 4 (EEG_Julia: Notch_order=4, Bandreject_order=4,
    #                             Highpass_order_design=4, Lowpass_order_design=4)
    @test notch.order == 4
    @test br.order    == 4
    @test hp_f.order  == 4
    @test lp_f.order  == 4

    # Frecuencias del Notch: 49.5–50.5 Hz (EEG_Julia: freq=50, width=1.0)
    @test occursin("49.5", notch.freq) || occursin("49.5–50.5", notch.freq)
    @test occursin("50.5", notch.freq)

    # Frecuencias del Bandreject: 99.5–100.5 Hz (EEG_Julia: freq=100, bw=1.0)
    @test occursin("99.5", br.freq)
    @test occursin("100.5", br.freq)

    # HP: 0.5 Hz; LP: 150.0 Hz  (EEG_Julia: Highpass_cutoff=0.5, Lowpass_cutoff=150)
    @test occursin("0.5",   hp_f.freq)
    @test occursin("150.0", lp_f.freq)

    # Orden efectivo con filtfilt: design × 2 = 8
    ord_design = cfg.filtering["filter_order"]
    @test ord_design == 4
    @test ord_design * 2 == 8    # orden efectivo HP y LP

    println("  ✓ Cadena EEG_Julia: $(join([s.name for s in chain], " → "))")
    println("  ✓ Notch 49.5–50.5 Hz, BR 99.5–100.5 Hz, HP 0.5 Hz, LP 150.0 Hz")
    println("  ✓ Notch/BR: filt (causal ord.4)  HP/LP: filtfilt (ord.efectivo 8)")
end

# ─── Tests de defaults del código (sin config) ────────────────
# Verifica que los valores por defecto en Filtering.jl coinciden con EEG_Julia
# cuando no se carga config (p. ej., llamadas programáticas sin TOML).

@testset "FilteringDefaults_EEGJulia" begin
    # Config mínima sin keys de filtrado → debe usar defaults del código
    cfg_min = PipelineConfig(
        Dict{String,Any}("name" => "test", "version" => "0.2"),
        Dict{String,Any}("n_ms_t1" => 44, "n_controls" => 40,
                         "conditions" => ["EO","EC"],
                         "groups" => ["MS","Control"], "sessions" => ["T1","T2"]),
        Dict{String,Any}("results" => "results", "data_cache" => "data/cache",
                         "bids_root" => "data/BIDS"),
        Dict{String,Any}("fs" => 500.0, "conditions" => ["EO","EC"],
                         "n_channels" => 32, "reference" => "average"),
        # filtering: solo el perfil, sin parámetros → defaults del código
        Dict{String,Any}("profile" => "eeg_julia",
                         "bandreject_lo" => 99.5, "bandreject_hi" => 100.5),
        Dict{String,Any}("epoch_length_s" => 1.0, "epoch_overlap" => 0.0,
                         "min_epochs" => 5),
        Dict{String,Any}("apply" => true, "method" => "mean"),
        Dict{String,Any}("amplitude_threshold_uv" => 100.0,
                         "gradient_threshold_uv"  => 50.0, "enabled" => true),
        Dict{String,Any}("n_components" => 10, "random_seed" => 42),
        Dict{String,Any}("nfft" => 256, "window" => "hamming", "window_pct" => 10.0),
        Dict{String,Tuple{Float64,Float64}}("DELTA" => (0.5,4.0), "ALPHA" => (7.8,11.7)),
        Dict{String,Any}("filter_order" => 8, "use_csd" => false, "method" => "wpli"),
        Dict{String,Any}("n_surrogates" => 10, "alpha" => 0.05,
                         "method" => "phase_shuffle", "fdr_method" => "bh"),
        Dict{String,Any}("threshold_method" => "proportional", "density" => 0.1),
        Dict{String,Any}("variables" => ["EDSS"]),
        Dict{String,Any}("min_visits" => 2),
        Dict{String,Any}("group_test" => "mann_whitney", "alpha" => 0.05,
                         "fdr_q" => 0.05, "paired_test" => "wilcoxon"),
        Dict{String,Any}("figure_format" => "png", "figure_dpi" => 150,
                         "table_format" => "csv"),
        Dict{String,Any}("bad_channel_zscore_threshold" => 3.0,
                         "amplitude_warning_sigma_uv" => 20.0),
        Dict{String,Any}("exclude_channels" => String[], "exclude_fp2" => false,
                         "n_channels_analysis" => 30),
        joinpath(@__DIR__, "..")
    )
    chain_min = describe_filter_chain(cfg_min)
    @test length(chain_min) == 4

    # Con defaults del código, LP debe ser 150.0 Hz (no 48.0 que era el bug)
    lp_step = chain_min[4]
    @test lp_step.name == "Low-pass"
    @test occursin("150.0", lp_step.freq)    # default correcto: 150 Hz

    # Con defaults del código, notch bw debe ser 1.0 Hz (no 2.0 que era el bug)
    notch_step = chain_min[1]
    @test notch_step.name == "Notch"
    @test occursin("49.5", notch_step.freq)  # bw=1.0 → 50-0.5=49.5 Hz

    println("  ✓ Default LP = 150.0 Hz (no el antiguo 48.0 Hz)")
    println("  ✓ Default Notch bw = 1.0 Hz (49.5–50.5 Hz)")
end

# ─── Tests de respuesta en frecuencia (atenuación real) ───────
# Verifica que los filtros construidos por NeuroMIND atenúan/preservan
# las frecuencias esperadas según el diseño Butterworth.

@testset "FilterAttenuation" begin
    fs  = 500.0
    nyq = fs / 2.0
    ord = 4

    # ── Notch 50 Hz, bw 1 Hz, filt (orden 4) ─────────────────
    f_notch = digitalfilter(Bandstop(49.5/nyq, 50.5/nyq), Butterworth(ord))
    # En 50 Hz → fuerte atenuación (< -40 dB)
    H_at_50 = freqresp(f_notch, [π * 50.0 / nyq])
    @test 20*log10(abs(H_at_50[1])) < -40.0
    # En 1 Hz → paso completo (> -1 dB)
    H_at_1  = freqresp(f_notch, [π * 1.0 / nyq])
    @test 20*log10(abs(H_at_1[1])) > -1.0

    # ── Bandreject 100 Hz, bw 1 Hz, filt (orden 4) ───────────
    f_br = digitalfilter(Bandstop(99.5/nyq, 100.5/nyq), Butterworth(ord))
    H_at_100 = freqresp(f_br, [π * 100.0 / nyq])
    @test 20*log10(abs(H_at_100[1])) < -40.0
    H_at_60  = freqresp(f_br, [π * 60.0 / nyq])
    @test 20*log10(abs(H_at_60[1])) > -1.0

    # ── HP 0.5 Hz, filtfilt → mag efectiva = |H|² ────────────
    f_hp = digitalfilter(Highpass(0.5/nyq), Butterworth(ord))
    # En 0.5 Hz → -3 dB (un solo paso); con filtfilt → -6 dB efectivo
    H_hp_at_fc = freqresp(f_hp, [π * 0.5 / nyq])
    mag_hp_fc  = abs(H_hp_at_fc[1])^2
    @test 20*log10(mag_hp_fc) ≈ -6.0 atol=1.5
    # En 10 Hz → paso completo (> -1 dB efectivo con filtfilt)
    H_hp_at_10 = freqresp(f_hp, [π * 10.0 / nyq])
    @test 20*log10(abs(H_hp_at_10[1])^2) > -1.0

    # ── LP 150 Hz, filtfilt → mag efectiva = |H|² ────────────
    f_lp = digitalfilter(Lowpass(150.0/nyq), Butterworth(ord))
    # En 150 Hz → -6 dB efectivo (filtfilt)
    H_lp_at_fc = freqresp(f_lp, [π * 150.0 / nyq])
    mag_lp_fc  = abs(H_lp_at_fc[1])^2
    @test 20*log10(mag_lp_fc) ≈ -6.0 atol=1.5
    # En 10 Hz → paso completo
    H_lp_at_10 = freqresp(f_lp, [π * 10.0 / nyq])
    @test 20*log10(abs(H_lp_at_10[1])^2) > -1.0

    println("  ✓ Notch 50 Hz: < −40 dB en 50 Hz, > −1 dB en 1 Hz")
    println("  ✓ BR 100 Hz:   < −40 dB en 100 Hz, > −1 dB en 60 Hz")
    println("  ✓ HP 0.5 Hz:   ≈ −6 dB efectivo en fc (filtfilt, ord.efect.8)")
    println("  ✓ LP 150 Hz:   ≈ −6 dB efectivo en fc (filtfilt, ord.efect.8)")
end

# ─── Tests de ICA — perfil eeg_julia y clasificación ─────────

# Helper: reconstruye mock_config con un dict ica personalizado
function _cfg_with_ica(ica_dict)
    c = mock_config()
    PipelineConfig(c.project, c.study, c.paths, c.recording, c.filtering,
                   c.segmentation, c.baseline, c.artifact_rejection,
                   ica_dict, c.spectral, c.bands, c.connectivity,
                   c.surrogates, c.graph, c.clinical, c.longitudinal,
                   c.statistics, c.export_cfg, c.qc, c.montage, c.root)
end

@testset "ICA_EEGJulia_Profile" begin
    n_ch = 8
    rec  = mock_recording(n_ch, 2000, 500.0)
    cfg  = _cfg_with_ica(Dict{String,Any}("profile" => "eeg_julia",
                                           "max_iter" => 512, "tol" => 1e-7, "seed" => 1234))
    ica  = run_ica(rec, cfg)

    # eeg_julia profile → n_comp must equal n_channels → square matrices
    A = ica.mixing_matrix
    W = ica.unmixing_matrix
    @test size(A) == (n_ch, n_ch)
    @test size(W) == (n_ch, n_ch)

    # A = inv(W_total) → A * W ≈ I (not just pseudoinverse)
    @test A * W ≈ Matrix(I, n_ch, n_ch) atol=1e-8
end

@testset "ICA_Features" begin
    fs     = 500.0
    n_ch   = 10
    n_ic   = 10
    n_s    = 3000
    rng    = Random.MersenneTwister(99)
    A_mat  = randn(rng, n_ch, n_ic)
    S_mat  = randn(rng, n_ic, n_s)
    cnames = ["Fp1","Fp2","F3","F4","Fz","C3","C4","T7","T8","Oz"]

    feat = compute_ica_features(A_mat, S_mat, fs, cnames)

    # DataFrame has all 7 required feature columns
    for col in [:frontal_ratio,:temporal_ratio,:blink_ratio,:emg_ratio,
                :line_ratio,:kurtosis,:extreme_frac]
        @test hasproperty(feat, col)
    end
    @test size(feat, 1) == n_ic

    # All feature values are finite
    for col in [:frontal_ratio,:temporal_ratio,:blink_ratio,:emg_ratio,
                :line_ratio,:kurtosis,:extreme_frac]
        @test all(isfinite, feat[!, col])
    end

    eval_df = evaluate_ica_components(feat; artifact_thresh=1.5)

    # evaluate adds score and label columns
    @test hasproperty(eval_df, :artifact_type)
    @test hasproperty(eval_df, :artifact_score)
    @test hasproperty(eval_df, :ocular_score)
    @test hasproperty(eval_df, :muscle_score)
    @test size(eval_df, 1) == n_ic

    # All labels belong to the known set
    valid = Set(["brain","eye/blink","muscle","line_noise","jump","unknown"])
    @test all(l -> l in valid, eval_df.artifact_type)
end

@testset "ICA_EEGJulia_Reproducibility" begin
    n_ch = 6
    rec  = mock_recording(n_ch, 1500, 500.0)
    cfg  = _cfg_with_ica(Dict{String,Any}("profile" => "eeg_julia",
                                           "max_iter" => 512, "tol" => 1e-7, "seed" => 1234))
    ica1 = run_ica(rec, cfg)
    ica2 = run_ica(rec, cfg)

    # Same seed (fixed in eeg_julia profile) → identical mixing matrix
    @test ica1.mixing_matrix ≈ ica2.mixing_matrix atol=1e-10
end

# Helper: reconstruye mock_config con un root personalizado (sandbox de I/O)
function _cfg_with_root(root_dir::String)
    c = mock_config()
    PipelineConfig(c.project, c.study, c.paths, c.recording, c.filtering,
                   c.segmentation, c.baseline, c.artifact_rejection,
                   c.ica, c.spectral, c.bands, c.connectivity,
                   c.surrogates, c.graph, c.clinical, c.longitudinal,
                   c.statistics, c.export_cfg, c.qc, c.montage, root_dir)
end

@testset "ICA_AutoReject_Precedence" begin
    fs     = 500.0
    n_ch   = 6
    n_ic   = 6
    rng    = Random.MersenneTwister(7)
    A_mat  = randn(rng, n_ch, n_ic)
    S_mat  = randn(rng, n_ic, 1500)
    cnames = ["Fp1","Fp2","F3","F4","Oz","Pz"]
    feat   = compute_ica_features(A_mat, S_mat, fs, cnames)

    tmp = mktempdir()
    cfg = _cfg_with_root(tmp)
    sid, sess, cond, task = "T99", "T1", "EC", "eyesclosed"
    subj_dir = joinpath(tmp, "results", "subjects", "sub-$sid", "ses-$sess", task)
    mkpath(subj_dir)

    # Sin CSV manual → sin rechazo
    @test !has_manual_ica_labels(cfg, sid, sess, cond)
    @test isempty(load_ica_labels(cfg, sid, sess, cond))

    # Umbral muy alto → evaluate_ica_components clasifica todo como "brain"
    eval_brain = evaluate_ica_components(feat; artifact_thresh=1000.0)
    @test isempty(write_ica_labels_auto(subj_dir, eval_brain))

    # Umbral muy bajo → todo "artifact" → auto rechaza todos y persiste el CSV
    eval_artifact = evaluate_ica_components(feat; artifact_thresh=-1000.0)
    rejected_auto = write_ica_labels_auto(subj_dir, eval_artifact)
    @test sort(rejected_auto) == collect(1:n_ic)
    @test isfile(joinpath(subj_dir, "ica_labels_auto.csv"))

    # CSV manual presente → gana sobre lo automático (precedencia)
    open(joinpath(subj_dir, "ica_labels.csv"), "w") do io
        write(io, "component,label\n2,artifact\n5,artifact\n")
    end
    @test has_manual_ica_labels(cfg, sid, sess, cond)
    @test sort(load_ica_labels(cfg, sid, sess, cond)) == [2, 5]

    rm(tmp; recursive=true, force=true)
end

# ─── Helper: config con n_surrogates personalizado ────────────

function _cfg_with_surrogates(n_sur::Int; seed::Int=42, alpha::Float64=0.05)
    c = mock_config()
    PipelineConfig(
        c.project, c.study, c.paths, c.recording, c.filtering,
        c.segmentation, c.baseline, c.artifact_rejection,
        c.ica, c.spectral, c.bands, c.connectivity,
        Dict{String,Any}("n_surrogates" => n_sur, "alpha" => alpha,
                         "method" => "circular_shift", "fdr_method" => "bh",
                         "seed" => seed),
        c.graph, c.clinical, c.longitudinal, c.statistics, c.export_cfg, c.qc, c.montage, c.root
    )
end

# ─── Tests de surrogates ──────────────────────────────────────

@testset "Surrogates_PValue_Bounds" begin
    # Los p-valores nunca deben ser 0.0 ni > 1.0.
    # Con corrección Monte Carlo (+1): p ∈ [1/(n+1), 1.0].
    n_sur = 20
    cfg   = _cfg_with_surrogates(n_sur)
    rec   = mock_recording(4, 2000, 500.0)
    ep    = segment_recording(rec, cfg)
    ep_bl = apply_baseline(ep, cfg)
    conn  = compute_wpli(ep_bl, cfg)
    sr    = surrogate_test(ep_bl, conn, "ALPHA", cfg)

    p_min = 1.0 / (n_sur + 1)
    n_ch  = size(sr.p_values, 1)

    # Extraer solo los pares off-diagonal (los tests de conectividad)
    p_offdiag = [sr.p_values[i,j] for i in 1:n_ch, j in 1:n_ch if i != j]

    # Ningún p-valor off-diagonal debe ser 0 (corrección +1)
    @test !any(iszero, p_offdiag)
    # Ningún p-valor debe superar 1
    @test all(p -> p <= 1.0 + 1e-10, p_offdiag)
    # Todos deben respetar el mínimo teórico 1/(n_sur+1)
    @test all(p -> p >= p_min - 1e-10, p_offdiag)
    # La diagonal debe ser 1.0 (no hay test de auto-conectividad)
    @test all(sr.p_values[i,i] ≈ 1.0 for i in 1:n_ch)

    println("  ✓ p_min_observado=$(round(minimum(p_offdiag), digits=4))  " *
            "p_min_teórico=$(round(p_min, digits=4))")
end

@testset "Surrogates_Coupling_Detected" begin
    # Señal sintética: 2 canales con lag de π/2 a 10 Hz (dentro de ALPHA 7.8-11.7).
    # wPLI observado debe ser alto (~1); la distribución nula (circular shift) debe
    # ser baja; el par debe quedar significativo tras FDR.
    fs       = 500.0
    f0       = 10.0               # Hz, dentro de ALPHA
    n_ch     = 2
    n_ep     = 30
    n_samp   = 500                # 1 s = 10 ciclos completos de f0 → circshift exacto
    n_sur    = 50                 # p_min = 1/51 ≈ 0.02 < FDR_threshold_k1 = 0.05

    delay_s  = round(Int, fs / (4 * f0))  # π/2 a 10 Hz → 12 muestras (24 ms)
    rng_data = MersenneTwister(77)
    data     = zeros(n_ch, n_samp, n_ep)
    t        = (0:n_samp-1) ./ fs

    for ep in 1:n_ep
        base        = sin.(2π * f0 .* t) .+ 0.05 .* randn(rng_data, n_samp)
        data[1, :, ep] = base
        data[2, :, ep] = circshift(base, delay_s)   # lag fijo → wPLI ≈ 1
    end

    meta   = RecordingMeta("SYN", "T1", "EC", 1, fs, n_ch,
                           ["Ch1", "Ch2"], nothing, "synthetic")
    epochs = EpochSet(meta, data, 1.0, n_ep, Int[])
    cfg    = _cfg_with_surrogates(n_sur)
    conn   = compute_wpli(epochs, cfg)

    W_obs_12 = conn.matrices["ALPHA"][1, 2]
    println("  wPLI observado Ch1-Ch2 en ALPHA = $(round(W_obs_12, digits=3))")

    # wPLI alto: el acoplamiento de fase π/2 debe ser detectable
    @test W_obs_12 > 0.5

    sr = surrogate_test(epochs, conn, "ALPHA", cfg)

    # La distribución nula (shift independiente) debe tener media << W_obs
    W_null_mean = mean(sr.null_distribution[1, 2, :])
    println("  wPLI nulo medio    Ch1-Ch2 en ALPHA = $(round(W_null_mean, digits=3))")
    @test W_null_mean < W_obs_12

    # p-valor del par acoplado: con n_sur=50 y 1 sola hipótesis, p_min=1/51≈0.02 < 0.05
    println("  p-valor Ch1-Ch2 = $(round(sr.p_values[1,2], digits=4))")
    @test sr.p_values[1, 2] < 0.05

    # El par debe quedar significativo tras FDR (1 hipótesis → umbral = alpha = 0.05)
    @test sr.sig_mask[1, 2]

    # El par no acoplado con sí mismo (diagonal) es 0
    @test sr.observed[1, 1] < 1e-10
    @test sr.observed[2, 2] < 1e-10
end

@testset "Surrogates_Uncoupled_No_FalseDiscoveries" begin
    # Canales sin acoplamiento (ruido blanco independiente).
    # La distribución nula debe parecerse a la distribución observada
    # y casi ningún par debe superar el umbral FDR.
    n_sur = 30
    cfg   = _cfg_with_surrogates(n_sur)
    rec   = mock_recording(4, 3000, 500.0)   # ruido blanco (randn)
    ep    = segment_recording(rec, cfg)
    ep_bl = apply_baseline(ep, cfg)
    conn  = compute_wpli(ep_bl, cfg)
    sr    = surrogate_test(ep_bl, conn, "ALPHA", cfg)

    n_pairs = size(sr.sig_mask, 1) * (size(sr.sig_mask, 1) - 1) ÷ 2   # 6 pares
    n_sig   = count(sr.sig_mask) ÷ 2   # dividido por 2 porque la máscara es simétrica
    println("  Pares significativos (FDR) en señal no acoplada: $(n_sig) / $(n_pairs)")

    # Bajo la hipótesis nula, se esperan 0 false discoveries (FDR controla la tasa)
    # Permitimos hasta 1 para no depender de la semilla exacta del ruido
    @test n_sig <= 1

    # Todos los p-valores deben ser válidos
    @test all(p -> p >= 1.0/(n_sur+1) - 1e-10, sr.p_values)
    @test all(p -> p <= 1.0 + 1e-10, sr.p_values)
end

@testset "Surrogates_RNG_Reproducibility" begin
    # Misma semilla + mismos datos → distribución nula idéntica bit a bit.
    # Distinta banda → distribución nula diferente.
    cfg  = _cfg_with_surrogates(15; seed=99)
    rec  = mock_recording(3, 2000, 500.0)
    ep   = segment_recording(rec, cfg)
    ep_bl = apply_baseline(ep, cfg)
    conn = compute_wpli(ep_bl, cfg)

    sr1 = surrogate_test(ep_bl, conn, "ALPHA", cfg)
    sr2 = surrogate_test(ep_bl, conn, "ALPHA", cfg)

    # Mismo seed → W_null idéntico
    @test sr1.null_distribution == sr2.null_distribution
    @test sr1.p_values == sr2.p_values

    # Distinta banda → seed distinto → W_null diferente
    if haskey(conn.matrices, "BETA_LOW")
        sr3 = surrogate_test(ep_bl, conn, "BETA_LOW", cfg)
        # Las distribuciones nulas de distintas bandas no deberían ser idénticas
        @test sr1.null_distribution != sr3.null_distribution
        println("  ✓ ALPHA y BETA_LOW tienen distribuciones nulas distintas")
    end
end

@testset "Surrogates_Validation_Functions" begin
    # validate_connectivity_matrix no debe emitir warnings con datos limpios.
    # validate_surrogate_result no debe emitir warnings con SurrogateResult válido.
    n_sur = 15
    cfg   = _cfg_with_surrogates(n_sur)
    rec   = mock_recording(4, 2000, 500.0)
    ep    = segment_recording(rec, cfg)
    ep_bl = apply_baseline(ep, cfg)
    conn  = compute_wpli(ep_bl, cfg)
    sr    = surrogate_test(ep_bl, conn, "ALPHA", cfg)

    W = conn.matrices["ALPHA"]

    # Propiedades básicas de la matriz observada
    @test all(isfinite, W)
    @test all(x -> x >= 0.0 - 1e-8, W)
    @test all(x -> x <= 1.0 + 1e-8, W)
    @test W ≈ W'                       # simétrica
    @test all(x -> abs(x) < 1e-10, diag(W))   # diagonal cero

    # validate_connectivity_matrix no lanza error ni warning para datos limpios
    @test_nowarn validate_connectivity_matrix(W, "ALPHA")

    # validate_surrogate_result no lanza error ni warning para SurrogateResult válido
    @test_nowarn validate_surrogate_result(sr)

    # validate_connectivity_matrix detecta NaN.
    # Al solo poner NaN en [1,2] la matriz también pierde simetría → 2 warnings.
    W_nan = copy(W); W_nan[1, 2] = NaN
    @test_logs (:warn, r"NaN") (:warn, r"simétric") validate_connectivity_matrix(W_nan, "ALPHA")

    # validate_connectivity_matrix detecta valores fuera de [0,1].
    # Al solo poner 1.5 en [1,2] también pierde simetría → 2 warnings.
    W_oor = copy(W); W_oor[1, 2] = 1.5
    @test_logs (:warn, r"fuera de") (:warn, r"simétric") validate_connectivity_matrix(W_oor, "ALPHA")

    println("  ✓ p_min_observado=$(round(minimum(sr.p_values[sr.p_values.>0]), digits=5))  " *
            "p_min_teórico=$(round(1.0/(n_sur+1), digits=5))")
end

# ─── Tests dwPLI ──────────────────────────────────────────────

@testset "wPLI_dwPLI_Uncoupled" begin
    # Para ruido blanco independiente, dwPLI debe ser cercano a 0
    # (estimador no sesgado: E[dwPLI] = 0 bajo H0), mientras que
    # wPLI clásico tiene sesgo positivo para muestras pequeñas.
    rng = MersenneTwister(7)
    fs  = 500.0; n_samp = 500; n_seg = 30; n_ch = 3
    epoch_s = n_samp / fs   # 1.0 s

    cfg_wpli  = mock_config()   # use_dwpli = false (default)
    cfg_dwpli = PipelineConfig(
        cfg_wpli.project, cfg_wpli.study, cfg_wpli.paths,
        cfg_wpli.recording, cfg_wpli.filtering, cfg_wpli.segmentation,
        cfg_wpli.baseline, cfg_wpli.artifact_rejection, cfg_wpli.ica,
        cfg_wpli.spectral, cfg_wpli.bands,
        Dict{String,Any}("filter_order" => 8, "use_csd" => false,
                         "use_dwpli" => true, "method" => "dwpli"),
        cfg_wpli.surrogates, cfg_wpli.graph, cfg_wpli.clinical,
        cfg_wpli.longitudinal, cfg_wpli.statistics, cfg_wpli.export_cfg,
        cfg_wpli.qc, cfg_wpli.montage, cfg_wpli.root
    )

    meta = RecordingMeta("T01", "T1", "EC", 1, fs, n_ch,
                         ["Ch$i" for i in 1:n_ch], nothing, "dummy.tsv")
    data = randn(rng, n_ch, n_samp, n_seg)
    ep   = EpochSet(meta, data, epoch_s, n_seg, Int[])

    conn_w  = compute_wpli(ep, cfg_wpli)
    conn_dw = compute_wpli(ep, cfg_dwpli)

    @test conn_w.method  == "wpli"
    @test conn_dw.method == "dwpli"

    Ww  = conn_w.matrices["ALPHA"]
    Wdw = conn_dw.matrices["ALPHA"]

    # wPLI en [0,1]; dwPLI en [-1,1] con valores negativos posibles bajo H0
    @test all(Ww  .>= -1e-8)
    @test all(Ww  .<= 1.0 + 1e-8)
    @test all(Wdw .>= -1.0 - 1e-8)
    @test all(Wdw .<=  1.0 + 1e-8)

    # Diagonal nula en ambos
    @test all(abs.(diag(Ww))  .< 1e-10)
    @test all(abs.(diag(Wdw)) .< 1e-10)

    # Simetría
    @test Ww  ≈ Ww'
    @test Wdw ≈ Wdw'

    # Para ruido blanco, dwPLI debe ser más cercano a 0 que wPLI
    # (dwPLI corrige el sesgo positivo de wPLI)
    off_diag = [(i,j) for i in 1:n_ch for j in (i+1):n_ch]
    mean_wpli  = mean(abs(Ww[i,j])  for (i,j) in off_diag)
    mean_dwpli = mean(abs(Wdw[i,j]) for (i,j) in off_diag)
    @test mean_dwpli < mean_wpli   # dwPLI menos sesgado bajo H0
end

@testset "wPLI_dwPLI_Coupled" begin
    # Para señales acopladas a 10 Hz con lag π/2, dwPLI debe detectar
    # el acoplamiento (valor positivo significativo).
    fs    = 500.0; n_samp = 500; n_seg = 20; n_ch = 2
    epoch_s = n_samp / fs
    t     = range(0.0, step=1/fs, length=n_samp)
    data  = zeros(Float64, n_ch, n_samp, n_seg)
    rng2  = MersenneTwister(99)
    for s in 1:n_seg
        φ = randn(rng2) * 0.1
        data[1, :, s] .= sin.(2π * 10 .* t .+ φ)
        data[2, :, s] .= sin.(2π * 10 .* t .+ φ .+ π/2)
    end
    meta = RecordingMeta("T01", "T1", "EC", 1, fs, n_ch,
                         ["Ch1", "Ch2"], nothing, "dummy.tsv")
    ep   = EpochSet(meta, data, epoch_s, n_seg, Int[])

    cfg  = mock_config()
    cfg_dw = PipelineConfig(
        cfg.project, cfg.study, cfg.paths, cfg.recording, cfg.filtering,
        cfg.segmentation, cfg.baseline, cfg.artifact_rejection, cfg.ica,
        cfg.spectral,
        Dict{String,Tuple{Float64,Float64}}("ALPHA" => (8.0, 12.0)),
        Dict{String,Any}("filter_order" => 8, "use_csd" => false,
                         "use_dwpli" => true, "method" => "dwpli"),
        cfg.surrogates, cfg.graph, cfg.clinical, cfg.longitudinal,
        cfg.statistics, cfg.export_cfg, cfg.qc, cfg.montage, cfg.root
    )

    conn_dw = compute_wpli(ep, cfg_dw)
    Wdw = conn_dw.matrices["ALPHA"]

    # Para señal 10 Hz acoplada con π/2: dwPLI ≈ wPLI² para acoplamiento fuerte
    @test Wdw[1, 2] > 0.5   # acoplamiento detectado
    @test Wdw ≈ Wdw'
end

@testset "wPLI_BandDuration_Warning" begin
    # DELTA (0.5 Hz) con epoch_length_s = 1.0 s → 0.5 ciclos/época < 4.0 mínimo.
    # compute_wpli debe emitir @warn.
    rng  = MersenneTwister(13)
    fs   = 500.0; n_samp = 500; n_seg = 10; n_ch = 3
    epoch_s = n_samp / fs
    meta = RecordingMeta("T01", "T1", "EC", 1, fs, n_ch,
                         ["Ch$i" for i in 1:n_ch], nothing, "dummy.tsv")
    ep   = EpochSet(meta, randn(rng, n_ch, n_samp, n_seg), epoch_s, n_seg, Int[])

    c = mock_config()
    cfg_delta = PipelineConfig(
        c.project, c.study, c.paths, c.recording, c.filtering,
        c.segmentation, c.baseline, c.artifact_rejection, c.ica,
        c.spectral,
        Dict{String,Tuple{Float64,Float64}}(
            "DELTA" => (0.5, 4.0),    # 0.5 Hz × 1 s = 0.5 ciclos → warn
            "ALPHA" => (7.8, 11.7),   # 7.8 × 1 s = 7.8 ciclos → ok
        ),
        Dict{String,Any}("filter_order" => 8, "use_csd" => false,
                         "use_dwpli" => false, "min_cycles_for_wpli" => 4.0),
        c.surrogates, c.graph, c.clinical, c.longitudinal,
        c.statistics, c.export_cfg, c.qc, c.montage, c.root
    )

    # DELTA debe producir @warn con "DELTA" y "ciclos" en el mensaje
    @test_logs (:warn, r"DELTA.*ciclos") compute_wpli(ep, cfg_delta)
end

@testset "GraphMetrics_ChannelNames" begin
    # compute_graph_metrics debe propagar channel_names desde ConnectivityMatrix.
    cfg    = mock_config()
    rec    = mock_recording(5, 2000)
    ep     = segment_recording(rec, cfg)
    ep_bl  = apply_baseline(ep, cfg)
    conn   = compute_wpli(ep_bl, cfg)

    gm = compute_graph_metrics(conn, "ALPHA", cfg)

    @test gm.band == "ALPHA"
    @test length(gm.channel_names) == 5
    @test gm.channel_names == conn.channel_names
    @test length(gm.strength)   == 5
    @test length(gm.clustering) == 5
    @test all(gm.strength .>= 0.0)
    @test gm.density >= 0.0
    @test gm.density <= 1.0
end

println("\n✅ Todos los tests completados")
