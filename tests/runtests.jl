# NeuroMIND/tests/runtests.jl
# Suite de tests unitarios e integración.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Test
using NeuroMIND
using Statistics

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
        # filtering
        Dict{String,Any}("highpass_hz" => 0.5, "lowpass_hz" => 48.0,
                         "notch_hz" => 50.0, "notch_bw_hz" => 2.0,
                         "bandreject_lo" => 100.0, "bandreject_hi" => 120.0,
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

    rec_br = apply_bandreject(rec, 100.0, 120.0)
    @test size(rec_br.data) == size(rec.data)

    rec_filt = filter_recording(rec, cfg)
    @test size(rec_filt.data) == size(rec.data)
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
    @test size(sp.psd, 1) == 10
    @test size(sp.psd, 2) == 256÷2+1
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
    @test cfg.filtering["highpass_hz"] == 0.5
    @test cfg.filtering["bandreject_lo"] == 100.0
    @test cfg.statistics["fdr_q"] == 0.05
    @test haskey(cfg.study, "n_ms_t1")
    @test cfg.study["n_ms_t1"] == 44
end

println("\n✅ Todos los tests completados")
