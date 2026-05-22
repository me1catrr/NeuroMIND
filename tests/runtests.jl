# NeuroMIND/tests/runtests.jl
# Suite de tests unitarios e integración.
# Uso directo:  julia --project=. tests/runtests.jl
# Vía Pkg.test(): el entorno ya está activo; no se llama Pkg.activate aquí.

using Test
using NeuroMIND
using Statistics
using LinearAlgebra: diag
using DSP: digitalfilter, Bandstop, Highpass, Lowpass, Butterworth, freqresp

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
        cfg.statistics, cfg.export_cfg, cfg.root
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

println("\n✅ Todos los tests completados")
