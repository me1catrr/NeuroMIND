# NeuroMIND/src/SingleSubjectPipeline.jl
# Pipeline mínimo para análisis de un solo sujeto EEG.
# Fase inicial antes de escalar a análisis transversal y longitudinal.

using CairoMakie

# ─── Configuración ────────────────────────────────────────────

"""
    load_ss_config(path) -> PipelineConfig

Carga single_subject.toml y construye un PipelineConfig compatible
con todas las funciones de algoritmos existentes.
"""
function load_ss_config(path::String)::PipelineConfig
    abs_path = isabspath(path) ? path : abspath(path)
    isfile(abs_path) || error("Config no encontrada: $abs_path")
    raw  = TOML.parsefile(abs_path)
    root = dirname(dirname(abs_path))   # sube config/ → NeuroMIND/

    get_s(sec, k, def) = get(get(raw, sec, Dict()), k, def)

    filt = get(raw, "filtering",          Dict{String,Any}())
    seg  = get(raw, "segmentation",       Dict{String,Any}())
    ar   = get(raw, "artifact_rejection", Dict{String,Any}())
    sp   = get(raw, "spectral",           Dict{String,Any}())
    conn = get(raw, "connectivity",       Dict{String,Any}())
    paths_r = get(raw, "paths",           Dict{String,Any}())
    out  = get(raw, "output",             Dict{String,Any}())
    rec  = get(raw, "recording",          Dict{String,Any}())

    bands_raw = get(raw, "bands", Dict{String,Any}())
    bands = Dict{String,Tuple{Float64,Float64}}(
        k => (Float64(v[1]), Float64(v[2])) for (k, v) in bands_raw
    )

    overlap_s  = Float64(get(seg, "overlap_seconds", 0.0))
    seg_len_s  = Float64(get(seg, "segment_length_seconds", 2.0))
    overlap_frac = seg_len_s > 0.0 ? overlap_s / seg_len_s : 0.0

    PipelineConfig(
        Dict{String,Any}("name" => "SingleSubject"),
        Dict{String,Any}("conditions" => ["EC"]),
        Dict{String,Any}(
            "bids_root"  => String(get(paths_r, "bids_root", "data/BIDS")),
            "results"    => String(get(paths_r, "results",   "results")),
            "data_cache" => "data/cache",
            "web_public" => "web/public",
        ),
        Dict{String,Any}("fs" => Float64(get(rec, "sampling_rate", 500.0))),
        Dict{String,Any}(
            "highpass_hz"   => Float64(get(filt, "highpass_hz",   0.5)),
            "lowpass_hz"    => Float64(get(filt, "lowpass_hz",   48.0)),
            "notch_hz"      => Float64(get(filt, "notch_hz",     50.0)),
            "notch_bw_hz"   => Float64(get(filt, "notch_bw_hz",   2.0)),
            "bandreject_lo" => Float64(get(filt, "bandreject_lo", 100.0)),
            "bandreject_hi" => Float64(get(filt, "bandreject_hi", 120.0)),
            "filter_order"  => Int(get(filt,    "filter_order",     4)),
        ),
        Dict{String,Any}(
            "epoch_length_s" => seg_len_s,
            "epoch_overlap"  => overlap_frac,
            "min_epochs"     => Int(get(seg, "min_segments", 10)),
        ),
        Dict{String,Any}("apply" => true, "method" => "mean"),
        Dict{String,Any}(
            "amplitude_threshold_uv" => Float64(get(ar, "amplitude_threshold_uv", 100.0)),
            "gradient_threshold_uv"  => Float64(get(ar, "gradient_threshold_uv",   50.0)),
            "enabled"                => Bool(get(ar,  "enabled", true)),
        ),
        Dict{String,Any}("n_components" => 30, "method" => "fastica"),
        Dict{String,Any}(
            "nfft"       => Int(get(sp, "nfft", 512)),
            "window_pct" => Float64(get(sp, "window_pct", 10.0)),
        ),
        bands,
        Dict{String,Any}(
            "filter_order" => Int(get(conn, "filter_order", 8)),
            "use_csd"      => Bool(get(conn, "use_csd", false)),
        ),
        Dict{String,Any}(),
        Dict{String,Any}(),
        Dict{String,Any}(),
        Dict{String,Any}(),
        Dict{String,Any}(),
        Dict{String,Any}(
            "figure_format" => String(get(out, "figure_format", "png")),
            "figure_dpi"    => Int(get(out, "figure_dpi", 150)),
        ),
        root,
    )
end

# ─── Detección automática de sujeto ───────────────────────────

"""
    detect_first_subject(bids_raw_dir) -> (subject_id, session_id, task, run)

Detecta el primer sujeto disponible en data/BIDS/raw/ por orden alfabético.
Patrón esperado: sub-{id}_ses-{sess}_task-{task}_run-{run}_eeg_data.tsv
"""
function detect_first_subject(bids_raw_dir::String)::Tuple{String,String,String,Int}
    isdir(bids_raw_dir) || error("Directorio BIDS/raw no encontrado: $bids_raw_dir")
    files = filter(f -> endswith(f, "_eeg_data.tsv"), readdir(bids_raw_dir))
    isempty(files) && error("No se encontraron archivos _eeg_data.tsv en: $bids_raw_dir")
    sort!(files)
    f = files[1]
    m = match(r"^sub-([^_]+)_ses-([^_]+)_task-([^_]+)_run-(\d+)_eeg_data\.tsv$", f)
    m === nothing && error("Nombre no reconoce patrón BIDS: $f")
    return String(m[1]), String(m[2]), String(m[3]), parse(Int, m[4])
end

# ─── Carga del sujeto ─────────────────────────────────────────

"""
    load_single_subject(cfg, subj_id, sess_id, task, run) -> EEGRecording

Carga señal EEG, metadatos y posiciones de electrodos.
El TSV tiene orientación filas=canales, columnas=muestras.
"""
function load_single_subject(
    cfg::PipelineConfig,
    subj_id::String,
    sess_id::String,
    task::String,
    run::Int
)::EEGRecording

    bids_dir     = bids_root_dir(cfg)
    raw_dir      = joinpath(bids_dir, "raw")
    elec_dir     = joinpath(bids_dir, "electrodes")
    base_prefix  = "sub-$(subj_id)_ses-$(sess_id)_task-$(task)_run-$(lpad(run,2,'0'))"
    data_path    = joinpath(raw_dir, "$(base_prefix)_eeg_data.tsv")
    # metadata puede existir con o sin sufijo _eeg
    meta_path    = let p1 = joinpath(raw_dir, "$(base_prefix)_eeg_metadata.json"),
                            p2 = joinpath(raw_dir, "$(base_prefix)_metadata.json")
        isfile(p1) ? p1 : p2
    end
    elec_path = joinpath(elec_dir, "sub-$(subj_id)_ses-$(sess_id)_electrodes.tsv")

    isfile(data_path) || error("EEG no encontrado: $data_path")
    isfile(meta_path) || error("Metadata no encontrado: $meta_path")

    # Metadata JSON (campos reales: "fs", "channel_names", etc.)
    meta_raw = _parse_bids_json(meta_path)
    fs = Float64(get(meta_raw, "fs", cfg.recording["fs"]))

    # TSV: filas=canales, primera col="Channel", resto=muestras
    df       = CSV.read(data_path, DataFrame; delim='\t')
    ch_names = string.(df[:, 1])
    data     = Float64.(Matrix(df[:, 2:end]))   # (channels × samples)

    n_ch, n_samp = size(data)
    @info "Cargado: $(subj_id)/$(sess_id)/$(task) — $(n_ch) canales × $(n_samp) muestras @ $(fs) Hz"

    ch_pos = isfile(elec_path) ? _load_electrode_positions(elec_path) : nothing

    condition = task == "eyesclosed" ? "EC" : (task == "eyesopen" ? "EO" : task)
    meta  = RecordingMeta(subj_id, sess_id, condition, run,
                          fs, n_ch, ch_names, ch_pos, data_path)
    times = collect(0.0:(1/fs):(n_samp - 1) / fs)

    return EEGRecording(meta, data, times)
end

# ─── Validación de electrodos ─────────────────────────────────

"""
    validate_channels(rec, elec_path) -> NamedTuple

Comprueba que los canales de la señal estén en el archivo de electrodos.
"""
function validate_channels(rec::EEGRecording, elec_path::String)
    if !isfile(elec_path)
        return (ok=true, matched=rec.meta.channel_names, missing_in_elec=String[],
                extra_in_elec=String[], msg="Sin archivo de electrodos")
    end
    df_e = CSV.read(elec_path, DataFrame; delim='\t')
    elec_names = uppercase.(string.(df_e[:, :name]))
    sig_names  = uppercase.(rec.meta.channel_names)
    missing_in_elec = setdiff(sig_names, elec_names)
    extra_in_elec   = setdiff(elec_names, sig_names)
    ok = isempty(missing_in_elec)
    ok || @warn "Canales en señal sin entrada en electrodos: $missing_in_elec"
    return (ok=ok, matched=intersect(sig_names, elec_names),
            missing_in_elec=missing_in_elec, extra_in_elec=extra_in_elec,
            msg=ok ? "Validación OK" : "Canales sin posición: $(join(missing_in_elec,','))")
end

# ─── Pipeline principal ───────────────────────────────────────

"""
    run_single_subject_pipeline(config_path)

Ejecuta el pipeline completo para el primer sujeto disponible.
Genera resultados en results/{subj}/{sess}/ y results/subjects/.
"""
function run_single_subject_pipeline(config_path::String)
    t0 = now()
    println("=" ^ 62)
    println(" NeuroMIND — Pipeline Sujeto Individual")
    println(" $(t0)")
    println("=" ^ 62)

    # ── 1. Configuración ──────────────────────────────────────
    cfg = load_ss_config(config_path)
    raw_cfg = TOML.parsefile(isabspath(config_path) ? config_path : abspath(config_path))
    sub_cfg = get(raw_cfg, "subject", Dict{String,Any}())
    out_cfg = get(raw_cfg, "output",  Dict{String,Any}())
    qc_cfg  = get(raw_cfg, "qc",      Dict{String,Any}())
    dpi     = Int(get(out_cfg, "figure_dpi", 150))

    # ── 2. Detección / resolución del sujeto ──────────────────
    raw_dir = joinpath(bids_root_dir(cfg), "raw")
    sid     = String(get(sub_cfg, "subject_id", "auto"))
    sess_id = String(get(sub_cfg, "session_id", "T2"))
    task    = String(get(sub_cfg, "task", "eyesclosed"))
    run_n   = Int(get(sub_cfg, "run", 1))

    if sid == "auto"
        sid, sess_id, task, run_n = detect_first_subject(raw_dir)
        println("▶ Sujeto detectado automáticamente: sub-$(sid) / ses-$(sess_id) / $(task)")
    else
        println("▶ Sujeto configurado: sub-$(sid) / ses-$(sess_id) / $(task)")
    end

    condition = task == "eyesclosed" ? "EC" : (task == "eyesopen" ? "EO" : task)

    # ── 3. Directorios de salida ───────────────────────────────
    subj_dir   = joinpath(results_dir(cfg), sid, sess_id)
    export_dir = joinpath(results_dir(cfg), "subjects",
                          "sub-$(sid)", "ses-$(sess_id)", task)
    for d in [joinpath(subj_dir, "figures"),
              joinpath(subj_dir, "tables"),
              joinpath(subj_dir, "cache"),
              joinpath(subj_dir, "logs"),
              export_dir]
        mkpath(d)
    end

    log_path = joinpath(export_dir, "pipeline_log.txt")
    log_io   = open(log_path, "w")
    _log(log_io, "NeuroMIND pipeline — $(t0)")
    _log(log_io, "Sujeto: sub-$(sid) | Sesión: ses-$(sess_id) | Tarea: $(task)")
    _log(log_io, "Config: $(config_path)")

    # ── 4. Carga de datos ─────────────────────────────────────
    println("\n[1/8] Cargando datos EEG...")
    _log(log_io, "\n[1/8] Carga de datos")
    rec = load_single_subject(cfg, sid, sess_id, task, run_n)
    _log(log_io, "  Canales: $(rec.meta.n_channels) | Muestras: $(n_samples(rec)) | fs: $(rec.meta.fs) Hz")
    _log(log_io, "  Duración: $(round(duration(rec), digits=2)) s")

    # Validación de electrodos
    elec_path = joinpath(bids_root_dir(cfg), "electrodes",
                         "sub-$(sid)_ses-$(sess_id)_electrodes.tsv")
    val = validate_channels(rec, elec_path)
    _log(log_io, "  Electrodos: $(val.msg)")

    # ── 5. QC básico ──────────────────────────────────────────
    println("[2/8] Control de calidad...")
    _log(log_io, "\n[2/8] QC de canales")
    qc_stats  = compute_channel_stats(rec)
    z_thresh  = Float64(get(qc_cfg, "bad_channel_zscore_threshold", 3.0))
    bad_ch    = flag_bad_channels(rec; z_threshold=z_thresh)
    qc_stats[!, :is_bad] = [ch in bad_ch for ch in qc_stats.channel]
    if isempty(bad_ch)
        _log(log_io, "  Sin canales sospechosos (umbral z=$(z_thresh))")
    else
        _log(log_io, "  Canales sospechosos (z>$(z_thresh)): " * join(bad_ch, ", "))
    end

    # ── 6. Filtrado ───────────────────────────────────────────
    println("[3/8] Filtrado...")
    _log(log_io, "\n[3/8] Filtrado")
    rec_filt = filter_recording(rec, cfg)
    _log(log_io, "  HP=$(cfg.filtering["highpass_hz"])Hz | LP=$(cfg.filtering["lowpass_hz"])Hz | Notch=$(cfg.filtering["notch_hz"])Hz")

    # ── ICA: antes de segmentar (señal continua filtrada) ─────
    println("[4/8] ICA (separación de fuentes)...")
    _log(log_io, "\n[4/8] ICA")
    t_ica = now()
    ica_result = try
        run_ica(rec_filt, cfg)
    catch e
        @warn "ICA falló: $e · continuando sin ICA"
        nothing
    end
    rec_ica = rec_filt   # señal que va a segmentación (filtrada o limpiada)
    if ica_result !== nothing
        n_comp_ica = size(ica_result.activations, 1)
        _log(log_io, "  Componentes: $(n_comp_ica) | Método: FastICA + PCA whitening")
        rej_labels = load_ica_labels(cfg, sid, sess_id, condition)
        ica_result = ICAResult(
            ica_result.meta, ica_result.mixing_matrix, ica_result.unmixing_matrix,
            ica_result.activations, rej_labels, ica_result.variance_explained
        )
        if !isempty(rej_labels)
            _log(log_io, "  Componentes rechazados: $(join(rej_labels, ", "))")
            rec_ica = apply_ica_rejection(rec_filt, ica_result, rej_labels)
            _log(log_io, "  Señal limpiada · $(length(rej_labels)) componente(s) eliminado(s)")
        else
            _log(log_io, "  Sin labels de rechazo · Se conservan todos los componentes")
        end
        ica_dur = round(Dates.value(now() - t_ica) / 1000, digits=1)
        _log(log_io, "  Duración ICA: $(ica_dur) s")
        _save_ica_results(ica_result, rec_filt, rec_ica, export_dir, ica_dur, log_io)
    else
        _log(log_io, "  ICA omitido por error · continuando con señal filtrada")
    end

    # ── Segmentación sobre señal limpiada (o filtrada si sin ICA) ──
    println("[5/8] Segmentación y rechazo de artefactos...")
    _log(log_io, "\n[5/8] Segmentación")
    t_seg      = now()
    epochs_raw = segment_recording(rec_ica, cfg)
    epochs_bl  = apply_baseline(epochs_raw, cfg)   # todos los epochs, baseline corregido
    epochs     = reject_artifacts(epochs_bl, cfg)  # solo epochs válidos
    n_total    = n_epochs(epochs_raw)
    n_valid    = epochs.n_valid
    n_rejected = length(epochs.rejected_idx)
    seg_dur    = round(Dates.value(now() - t_seg) / 1000, digits=1)
    _log(log_io, "  Segmentos totales: $(n_total) | Válidos: $(n_valid) | Rechazados: $(n_rejected) ($(round(100*n_rejected/n_total, digits=1))%)")
    _log(log_io, "  Duración segmentación: $(seg_dur) s")
    seg_signal_label = (ica_result !== nothing && !isempty(ica_result.rejected_components)) ? "ICA-limpiada" : "filtrada"
    _save_segmentation_results(epochs_bl, epochs, rec_ica, cfg, export_dir, t_seg, log_io, seg_signal_label)

    # ── 8. Análisis espectral ─────────────────────────────────
    println("[6/8] Análisis espectral (PSD)...")
    _log(log_io, "\n[6/8] Espectral")
    spectra = compute_psd(epochs, cfg)
    _log(log_io, "  nfft=$(spectra.params["nfft"]) | Bins: $(length(spectra.freqs)) | Resolución: $(round(spectra.freqs[2]-spectra.freqs[1], digits=3)) Hz")
    for (band, _) in sort(collect(cfg.bands))
        bp_mean = mean(spectra.band_power[band])
        _log(log_io, "  Potencia $(lpad(band,9)): $(round(bp_mean, digits=4)) μV² (media canales)")
    end

    # ── 9. Conectividad wPLI ──────────────────────────────────
    println("[7/8] Conectividad wPLI...")
    _log(log_io, "\n[7/8] Conectividad wPLI")
    use_csd = Bool(get(cfg.connectivity, "use_csd", false))
    epochs_conn = use_csd ? apply_csd(epochs, cfg) : epochs
    conn = compute_wpli(epochs_conn, cfg)
    _log(log_io, "  Espacio: $(conn.space) | Bandas: " * join(sort(collect(keys(conn.matrices))), ", "))

    # ── 10. Guardar resultados ────────────────────────────────
    println("[8/8] Guardando resultados...")
    _log(log_io, "\n[8/8] Guardando resultados")
    _save_all_results(rec, rec_filt, qc_stats, bad_ch, spectra, conn, val,
                      subj_dir, export_dir, condition, cfg, dpi, log_io)

    # Snapshot de config
    _save_config_snapshot(config_path, export_dir)

    # subjects_index.csv
    _update_subjects_index(results_dir(cfg), sid, sess_id, task, rec, n_valid, n_rejected)

    elapsed = round(Dates.value(now() - t0) / 1000, digits=1)
    _log(log_io, "\n✓ Pipeline completado en $(elapsed) s")
    close(log_io)

    println("\n" * "=" ^ 62)
    println("✅ Pipeline completado en $(elapsed) s")
    println()
    println("Resultados principales:")
    println("  Dashboard:  results/$(sid)/$(sess_id)/")
    println("  Exportación BIDS: results/subjects/sub-$(sid)/ses-$(sess_id)/$(task)/")
    println("=" ^ 62)

    return (subject_id=sid, session_id=sess_id, task=task,
            condition=condition, recording=rec, spectra=spectra, connectivity=conn)
end

# ─── Guardado de resultados ───────────────────────────────────

function _save_segmentation_results(
    epochs_bl::EpochSet,       # todos los epochs (pre-AR), baseline corregido
    epochs::EpochSet,          # solo válidos (post-AR)
    rec::EEGRecording,         # señal de entrada a la segmentación
    cfg::PipelineConfig,
    export_dir::String,
    t_start::DateTime,
    log_io::IO,
    signal_input::String = "filtrada"
)
    n_total   = size(epochs_bl.data, 3)
    n_valid   = epochs.n_valid
    n_rejected = length(epochs.rejected_idx)
    ret_pct   = round(100.0 * n_valid / max(n_total, 1), digits=1)
    fs        = rec.meta.fs
    epoch_s   = epochs_bl.epoch_length_s
    overlap_f = Float64(get(cfg.segmentation, "epoch_overlap", 0.0))
    overlap_s = round(epoch_s * overlap_f, digits=3)
    step_s    = round(epoch_s * (1.0 - overlap_f), digits=3)
    samp_ep   = round(Int, epoch_s * fs)
    sig_dur   = round(n_samples(rec) / fs, digits=2)
    amp_thr   = Float64(get(cfg.artifact_rejection, "amplitude_threshold_uv", 100.0))
    grad_thr  = Float64(get(cfg.artifact_rejection, "gradient_threshold_uv",  50.0))
    bl_method = String(get(cfg.baseline, "method", "mean"))
    ts        = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    dur_s     = round(Dates.value(now() - t_start) / 1000, digits=1)

    # ── Informe de calidad por época ──────────────────────────────────────────
    qr = try
        compute_epoch_quality_report(epochs_bl, cfg)
    catch e
        @warn "compute_epoch_quality_report falló: $e"
        DataFrame()
    end

    n_rej_amp  = isempty(qr) ? 0 : count(==("amplitude"), qr.rejection_reason)
    n_rej_grad = isempty(qr) ? 0 : count(==("gradient"),  qr.rejection_reason)
    q_vals     = isempty(qr) ? Float64[] : Float64.(qr.quality)
    q_mean     = isempty(q_vals) ? 0.0 : round(mean(q_vals), digits=3)
    q_median   = isempty(q_vals) ? 0.0 : round(median(q_vals), digits=3)

    # Histograma de calidad (10 bins × 0.1)
    hist_counts = zeros(Int, 10)
    for q in q_vals
        b = min(floor(Int, q * 10), 9)
        hist_counts[b+1] += 1
    end
    hist_json = join(["[$(round((i-1)*0.1, digits=1)),$(hist_counts[i])]" for i in 1:10], ",")

    # ── Cobertura por canal ───────────────────────────────────────────────────
    cov_df = try
        compute_channel_coverage(epochs_bl, cfg)
    catch e
        @warn "compute_channel_coverage falló: $e"
        DataFrame()
    end
    cov_mean = isempty(cov_df) ? 0.0 : round(mean(cov_df.coverage_pct), digits=1)
    cov_min  = isempty(cov_df) ? 0.0 : round(minimum(cov_df.coverage_pct), digits=1)
    cov_max  = isempty(cov_df) ? 0.0 : round(maximum(cov_df.coverage_pct), digits=1)

    # ── Guardado de ficheros ──────────────────────────────────────────────────

    # 1) segmentation_summary.json
    open(joinpath(export_dir, "segmentation_summary.json"), "w") do f
        write(f, """{
  "n_total": $(n_total),
  "n_valid": $(n_valid),
  "n_rejected": $(n_rejected),
  "retention_pct": $(ret_pct),
  "epoch_length_s": $(epoch_s),
  "overlap_s": $(overlap_s),
  "overlap_pct": $(round(overlap_f*100, digits=1)),
  "step_s": $(step_s),
  "n_channels": $(rec.meta.n_channels),
  "fs": $(fs),
  "samples_per_epoch": $(samp_ep),
  "signal_duration_s": $(sig_dur),
  "signal_input": "$(signal_input)",
  "amp_threshold_uv": $(amp_thr),
  "grad_threshold_uv": $(grad_thr),
  "baseline_method": "$(bl_method)",
  "n_rejected_amplitude": $(n_rej_amp),
  "n_rejected_gradient": $(n_rej_grad),
  "quality_mean": $(q_mean),
  "quality_median": $(q_median),
  "quality_threshold": 0.5,
  "quality_histogram": [$(hist_json)],
  "coverage_mean_pct": $(cov_mean),
  "coverage_min_pct": $(cov_min),
  "coverage_max_pct": $(cov_max),
  "timestamp": "$(ts)",
  "duration_s": $(dur_s)
}""")
    end
    _log(log_io, "  Guardado: segmentation_summary.json")

    # 2) segments_table.csv
    if !isempty(qr)
        CSV.write(joinpath(export_dir, "segments_table.csv"), qr)
        _log(log_io, "  Guardado: segments_table.csv ($(n_total) épocas)")
    end

    # 3) channel_coverage.csv
    if !isempty(cov_df)
        CSV.write(joinpath(export_dir, "channel_coverage.csv"), cov_df)
        _log(log_io, "  Guardado: channel_coverage.csv ($(rec.meta.n_channels) canales)")
    end

    # 4–6) Ficheros específicos de rechazo de artefactos
    if !isempty(qr)
        _save_ar_results(qr, cfg, export_dir, ts, dur_s,
                         n_total, n_valid, n_rejected, log_io)
    end
end

# ─── Guardado de resultados de rechazo de artefactos ──────────

function _save_ar_results(
    qr::DataFrame,              # informe completo pre-AR (todos los epochs)
    cfg::PipelineConfig,
    export_dir::String,
    ts::String,
    dur_s::Float64,
    n_total::Int,
    n_valid::Int,
    n_rejected::Int,
    log_io::IO
)
    ar_cfg     = cfg.artifact_rejection
    amp_thresh = Float64(get(ar_cfg, "amplitude_threshold_uv", 100.0))
    grad_thresh= Float64(get(ar_cfg, "gradient_threshold_uv",   50.0))
    ret_pct    = round(100.0 * n_valid / max(n_total, 1), digits=1)

    rej_mask   = qr.status .== "rejected"
    rej_df     = qr[rej_mask, :]
    n_rej_amp  = count(==("amplitude"), qr.rejection_reason)
    n_rej_grad = count(==("gradient"),  qr.rejection_reason)

    # ── Estadísticas P2P ─────────────────────────────────────────────────────
    p2p_vals = hasproperty(qr, :p2p_uv) ? Float64.(qr.p2p_uv) : Float64[]
    p2p_mean = isempty(p2p_vals) ? 0.0 : round(mean(p2p_vals),  digits=2)
    p2p_std  = isempty(p2p_vals) ? 0.0 : round(std(p2p_vals),   digits=2)
    p2p_max  = isempty(p2p_vals) ? 0.0 : round(maximum(p2p_vals), digits=2)
    p2p_thresh_2sd = round(p2p_mean + 2.0 * p2p_std, digits=2)

    # Histograma P2P (20 bins)
    p2p_hist_json = ""
    if !isempty(p2p_vals) && p2p_max > 0
        bin_w   = p2p_max / 20.0
        hist_c  = zeros(Int, 20)
        for v in p2p_vals
            b = min(floor(Int, v / bin_w), 19)
            hist_c[b+1] += 1
        end
        p2p_hist_json = join(
            ["[$(round((i-1)*bin_w, digits=1)),$(hist_c[i])]" for i in 1:20], ",")
    end

    # ── Tabla de canales más afectados ───────────────────────────────────────
    ch_bad = Dict{String, @NamedTuple{n::Int, amp::Int, grad::Int}}()
    if hasproperty(qr, :worst_channel)
        for row in eachrow(rej_df)
            ch   = string(row.worst_channel)
            prev = get(ch_bad, ch, (n=0, amp=0, grad=0))
            na   = string(row.rejection_reason) == "amplitude" ? prev.amp + 1 : prev.amp
            ng   = string(row.rejection_reason) == "gradient"  ? prev.grad + 1 : prev.grad
            ch_bad[ch] = (n=prev.n + 1, amp=na, grad=ng)
        end
    end
    ch_rows_sorted = sort(collect(ch_bad); by=x->x[2].n, rev=true)
    ca_df = isempty(ch_rows_sorted) ? DataFrame() :
        DataFrame(
            channel     = [r[1]  for r in ch_rows_sorted],
            n_bad       = [r[2].n   for r in ch_rows_sorted],
            pct_bad     = [round(100.0 * r[2].n / max(n_total,1), digits=1) for r in ch_rows_sorted],
            main_reason = [r[2].amp >= r[2].grad ? "amplitude" : "gradient" for r in ch_rows_sorted],
        )

    # ── 4) artifact_rejection_summary.json ───────────────────────────────────
    open(joinpath(export_dir, "artifact_rejection_summary.json"), "w") do f
        write(f, """{
  "n_total": $(n_total),
  "n_valid": $(n_valid),
  "n_rejected": $(n_rejected),
  "retention_pct": $(ret_pct),
  "n_rejected_amplitude": $(n_rej_amp),
  "n_rejected_gradient": $(n_rej_grad),
  "amp_threshold_uv": $(amp_thresh),
  "grad_threshold_uv": $(grad_thresh),
  "p2p_mean_uv": $(p2p_mean),
  "p2p_std_uv": $(p2p_std),
  "p2p_max_uv": $(p2p_max),
  "p2p_thresh_2sd": $(p2p_thresh_2sd),
  "p2p_histogram": [$(p2p_hist_json)],
  "timestamp": "$(ts)",
  "duration_s": $(dur_s)
}""")
    end
    _log(log_io, "  Guardado: artifact_rejection_summary.json")

    # ── 5) rejected_segments.csv ──────────────────────────────────────────────
    if !isempty(rej_df)
        CSV.write(joinpath(export_dir, "rejected_segments.csv"), rej_df)
        _log(log_io, "  Guardado: rejected_segments.csv ($(nrow(rej_df)) rechazados)")
    end

    # ── 6) channel_artifact_summary.csv ──────────────────────────────────────
    if !isempty(ca_df)
        CSV.write(joinpath(export_dir, "channel_artifact_summary.csv"), ca_df)
        _log(log_io, "  Guardado: channel_artifact_summary.csv")
    end
end

function _save_all_results(
    rec::EEGRecording, rec_filt::EEGRecording,
    qc_stats::DataFrame, bad_ch::Vector{String},
    spectra::SpectralResult, conn::ConnectivityMatrix,
    val, subj_dir::String, export_dir::String,
    condition::String, cfg::PipelineConfig, dpi::Int, log_io::IO
)
    tbl_dir = joinpath(subj_dir, "tables")
    fig_dir = joinpath(subj_dir, "figures")

    # ── Tablas QC ─────────────────────────────────────────────
    qc_path = joinpath(tbl_dir, "qc_channels_$(condition).csv")
    CSV.write(qc_path, qc_stats)
    cp(qc_path, joinpath(export_dir, "qc_summary.csv"); force=true)
    cp(qc_path, joinpath(export_dir, "channel_statistics.csv"); force=true)
    _log(log_io, "  Guardado: qc_channels_$(condition).csv")

    # ── Tabla overview ────────────────────────────────────────
    overview_df = DataFrame(
        subject_id   = [rec.meta.subject_id],
        session_id   = [rec.meta.session_id],
        task         = [condition == "EC" ? "eyesclosed" : "eyesopen"],
        sampling_hz  = [rec.meta.fs],
        n_channels   = [rec.meta.n_channels],
        n_samples    = [n_samples(rec)],
        duration_s   = [round(duration(rec), digits=2)],
        bad_channels = [join(bad_ch, ";")],
        n_bad_ch     = [length(bad_ch)],
        electrode_validation = [val.msg],
    )
    ov_path = joinpath(tbl_dir, "overview_$(condition).csv")
    CSV.write(ov_path, overview_df)
    cp(ov_path, joinpath(export_dir, "overview.csv"); force=true)

    # ── Tabla PSD ─────────────────────────────────────────────
    ch_names = rec.meta.channel_names
    n_freqs  = length(spectra.freqs)
    psd_rows = [(channel=ch_names[c], freq_hz=round(spectra.freqs[f], digits=3),
                 power_uv2=round(spectra.psd[c,f], digits=6))
                for c in 1:length(ch_names) for f in 1:n_freqs]
    psd_df = DataFrame(psd_rows)
    psd_path = joinpath(tbl_dir, "psd_by_channel_$(condition).csv")
    CSV.write(psd_path, psd_df)
    cp(psd_path, joinpath(export_dir, "psd_by_channel.csv"); force=true)
    _log(log_io, "  Guardado: psd_by_channel_$(condition).csv")

    # ── Tabla band power ──────────────────────────────────────
    band_names = sort(collect(keys(spectra.band_power)))
    bp_df = DataFrame(channel = ch_names)
    for b in band_names
        bp_df[!, b] = round.(spectra.band_power[b], digits=6)
    end
    bp_path = joinpath(tbl_dir, "band_power_$(condition).csv")
    CSV.write(bp_path, bp_df)
    cp(bp_path, joinpath(export_dir, "band_power_summary.csv"); force=true)
    _log(log_io, "  Guardado: band_power_$(condition).csv")

    # ── Tablas wPLI matrices + edges ──────────────────────────
    all_edges = DataFrame(ch_a=String[], ch_b=String[], band=String[],
                          wpli=Float64[], rank=Int[])
    for (band, W) in conn.matrices
        # Matriz como CSV
        mat_df = DataFrame(hcat(ch_names, W), vcat(["channel"], ch_names))
        mat_path = joinpath(tbl_dir, "wpli_matrix_$(band)_$(condition).csv")
        CSV.write(mat_path, mat_df)
        cp(mat_path, joinpath(export_dir, "wpli_$(band).csv"); force=true)

        # Edges (triángulo superior) — flat comprehension
        n_nodes = length(ch_names)
        edges = [(W[i,j], ch_names[i], ch_names[j])
                 for i in 1:n_nodes for j in (i+1):n_nodes]
        sort!(edges; rev=true)
        edge_df = DataFrame(
            ch_a  = [e[2] for e in edges],
            ch_b  = [e[3] for e in edges],
            band  = fill(band, length(edges)),
            wpli  = round.([e[1] for e in edges], digits=4),
            rank  = 1:length(edges)
        )
        edge_path = joinpath(tbl_dir, "wpli_edges_$(band)_$(condition).csv")
        CSV.write(edge_path, edge_df)
        append!(all_edges, edge_df)
    end
    edges_path = joinpath(export_dir, "connectivity_edges.csv")
    CSV.write(edges_path, all_edges)
    _log(log_io, "  Guardado: matrices y edges wPLI por banda")

    # ── Figura: señal preview ─────────────────────────────────
    try
        fig = _plot_signal_preview(rec, rec_filt)
        sig_path = joinpath(fig_dir, "signal_preview_$(condition).png")
        save_figure(fig, sig_path; dpi)
        cp(sig_path, joinpath(export_dir, "filtered_signal_preview.png"); force=true)
        _log(log_io, "  Guardado: signal_preview_$(condition).png")
    catch e
        @warn "No se pudo generar signal_preview: $e"
        _log(log_io, "  WARN: signal_preview fallido: $e")
    end

    # ── Figura: PSD grid de canales ───────────────────────────
    try
        fig = plot_spectrum_grid(spectra; xmax=50.0, cols=6)
        psd_fig_path = joinpath(fig_dir, "psd_all_channels_$(condition).png")
        save_figure(fig, psd_fig_path; dpi)
        cp(psd_fig_path, joinpath(export_dir, "psd_all_channels.png"); force=true)
        _log(log_io, "  Guardado: psd_all_channels_$(condition).png")
    catch e
        @warn "No se pudo generar psd_grid: $e"
        _log(log_io, "  WARN: psd_grid fallido: $e")
    end

    # ── Figuras: heatmap wPLI por banda ───────────────────────
    for (band, _) in conn.matrices
        try
            fig = plot_connectivity_heatmap(conn, band)
            wpli_fig = joinpath(fig_dir, "wpli_$(band)_$(condition).png")
            save_figure(fig, wpli_fig; dpi)
            cp(wpli_fig, joinpath(export_dir, "wpli_$(band).png"); force=true)
        catch e
            @warn "No se pudo generar heatmap wPLI $band: $e"
            _log(log_io, "  WARN: wpli_heatmap_$(band) fallido: $e")
        end
    end
    _log(log_io, "  Guardados: heatmaps wPLI por banda")

    # ── Figura: band power barplot ─────────────────────────────
    try
        fig = _plot_band_power_summary(spectra)
        bp_fig = joinpath(fig_dir, "band_power_summary_$(condition).png")
        save_figure(fig, bp_fig; dpi)
        cp(bp_fig, joinpath(export_dir, "band_power_summary.png"); force=true)
        _log(log_io, "  Guardado: band_power_summary_$(condition).png")
    catch e
        @warn "No se pudo generar band_power_summary: $e"
        _log(log_io, "  WARN: band_power_summary fallido: $e")
    end
end

# ─── Visualizaciones internas ─────────────────────────────────

function _plot_signal_preview(
    rec_raw::EEGRecording,
    rec_filt::EEGRecording;
    n_ch_show::Int = 5,
    t_window_s::Real = 10.0
)::CairoMakie.Figure

    fs          = rec_raw.meta.fs
    n_samp_win  = min(round(Int, t_window_s * fs), n_samples(rec_raw))
    n_ch_show   = min(n_ch_show, n_channels(rec_raw))
    t           = rec_raw.times[1:n_samp_win]

    fig = CairoMakie.Figure(size = (1200, 130 * n_ch_show + 70))
    CairoMakie.Label(fig[0, 1], "Vista previa de señal: raw (gris) vs filtrado (azul)";
                     fontsize=13, font=:bold)

    for i in 1:n_ch_show
        ax = CairoMakie.Axis(fig[i, 1];
            ylabel      = rec_raw.meta.channel_names[i],
            ylabelsize  = 10,
            yticklabelsize = 8,
            xticklabelsize = 8,
        )
        raw_sig  = rec_raw.data[i,  1:n_samp_win]
        filt_sig = rec_filt.data[i, 1:n_samp_win]
        CairoMakie.lines!(ax, t, raw_sig;  color=(:gray,   0.5), linewidth=0.8)
        CairoMakie.lines!(ax, t, filt_sig; color=:steelblue, linewidth=1.2)
        i < n_ch_show && CairoMakie.hidexdecorations!(ax; ticks=false)
    end
    CairoMakie.Label(fig[n_ch_show+1, 1], "Tiempo (s)"; fontsize=11)
    return fig
end

function _plot_band_power_summary(spectra::SpectralResult)::CairoMakie.Figure
    bands = sort(collect(keys(spectra.band_power)))
    n_bands = length(bands)
    ch_names = spectra.meta.channel_names
    n_ch = length(ch_names)

    cols = min(3, n_bands)
    rows = ceil(Int, n_bands / cols)

    fig = CairoMakie.Figure(size = (300 * cols, 220 * rows + 40))
    CairoMakie.Label(fig[0, 1:cols], "Potencia por banda (μV²)";
                     fontsize=13, font=:bold)

    palette = [:steelblue, :tomato, :green, :purple, :orange, :teal, :brown]

    for (k, band) in enumerate(bands)
        r = ((k-1) ÷ cols) + 1
        c = ((k-1) % cols) + 1
        ax = CairoMakie.Axis(fig[r, c];
            title = band, titlesize=10,
            xlabel="Canal", ylabel="μV²",
            xlabelsize=9, ylabelsize=9,
            xticklabelsize=7, yticklabelsize=7,
            xticks=(1:n_ch, ch_names),
            xticklabelrotation=π/3,
        )
        bp = spectra.band_power[band]
        col = palette[mod1(k, length(palette))]
        CairoMakie.barplot!(ax, 1:n_ch, bp; color=(col, 0.75))
    end
    return fig
end

# ─── Helpers ──────────────────────────────────────────────────

function _log(io::IO, msg::String)
    ts  = Dates.format(now(), "HH:MM:SS")
    println(io, "[$(ts)] $(msg)")
    flush(io)
end

function _save_ica_results(
    ica::ICAResult,
    rec_before::EEGRecording,
    rec_after::EEGRecording,
    export_dir::String,
    duration_s::Float64,
    log_io::IO
)
    n_comp  = size(ica.activations, 1)
    var_pct = clamp.(ica.variance_explained .* 100.0, 0.0, 100.0)

    # ── Tabla de componentes ──────────────────────────────────
    df = DataFrame(
        component     = 1:n_comp,
        variance_pct  = round.(var_pct, digits=2),
        rejected      = [i ∈ ica.rejected_components for i in 1:n_comp],
        label         = ["IC$(i)" for i in 1:n_comp],
        artifact_type = ["unknown" for _ in 1:n_comp],
    )
    CSV.write(joinpath(export_dir, "ica_components.csv"), df)
    _log(log_io, "  Guardado: ica_components.csv")

    # ── Resumen JSON ──────────────────────────────────────────
    n_rej   = length(ica.rejected_components)
    var_rej = sum(var_pct[i] for i in ica.rejected_components; init=0.0)
    var_ret = round(100.0 - var_rej, digits=1)
    ts      = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    open(joinpath(export_dir, "ica_summary.json"), "w") do f
        write(f, """{
  "n_components": $(n_comp),
  "n_rejected": $(n_rej),
  "n_accepted": $(n_comp - n_rej),
  "variance_retained": $(var_ret),
  "timestamp": "$(ts)",
  "duration_s": $(duration_s)
}""")
    end
    _log(log_io, "  Guardado: ica_summary.json")

    # ── Activaciones (primeros 10s) ───────────────────────────
    fs      = rec_before.meta.fs
    n_samp  = min(size(ica.activations, 2), Int(round(10.0 * fs)))
    t_vec   = collect(range(0.0, step=1.0/fs, length=n_samp))
    act_df  = DataFrame(:t_s => t_vec)
    for i in 1:n_comp
        act_df[!, Symbol("IC$(i)")] = round.(ica.activations[i, 1:n_samp], digits=4)
    end
    CSV.write(joinpath(export_dir, "ica_activations.csv"), act_df)
    _log(log_io, "  Guardado: ica_activations.csv ($(n_comp) componentes, $(round(n_samp/fs, digits=1))s)")

    # ── Señal antes y después (primeros 10s, todos los canales) ──
    ch_names = rec_before.meta.channel_names
    for (label, rec) in [("before", rec_before), ("after", rec_after)]
        ns  = min(size(rec.data, 2), Int(round(10.0 * fs)))
        tv  = collect(range(0.0, step=1.0/fs, length=ns))
        sdf = DataFrame(:t_s => tv)
        for (ci, ch) in enumerate(ch_names)
            sdf[!, Symbol(ch)] = round.(rec.data[ci, 1:ns], digits=4)
        end
        CSV.write(joinpath(export_dir, "ica_signal_$(label).csv"), sdf)
    end
    _log(log_io, "  Guardado: ica_signal_before.csv + ica_signal_after.csv")
end

function _save_config_snapshot(config_path::String, export_dir::String)
    dst = joinpath(export_dir, "config_snapshot.toml")
    try
        cp(isabspath(config_path) ? config_path : abspath(config_path), dst; force=true)
    catch e
        @warn "No se pudo copiar config snapshot: $e"
    end
end

function _update_subjects_index(results_dir::String, subj_id::String,
                                  sess_id::String, task::String,
                                  rec::EEGRecording, n_valid::Int, n_rejected::Int)
    idx_path = joinpath(results_dir, "subjects_index.csv")
    new_row = DataFrame(
        subject_id  = [subj_id],
        session_id  = [sess_id],
        task        = [task],
        n_channels  = [rec.meta.n_channels],
        n_samples   = [n_samples(rec)],
        duration_s  = [round(duration(rec), digits=2)],
        fs_hz       = [rec.meta.fs],
        n_segments  = [n_valid],
        n_rejected  = [n_rejected],
        processed_at = [string(now())],
    )
    if isfile(idx_path)
        old = CSV.read(idx_path, DataFrame)
        mask = .!(old.subject_id .== subj_id .&& old.session_id .== sess_id .&& old.task .== task)
        updated = vcat(old[mask, :], new_row)
        CSV.write(idx_path, updated)
    else
        CSV.write(idx_path, new_row)
    end
end

# ─── Carga de datos para dashboard ────────────────────────────

"""
    load_dashboard_data(results_dir, subj_id, sess_id, condition) -> Dict

Carga los resultados ya calculados para alimentar el dashboard.
"""
function load_dashboard_data(
    res_dir::String,
    subj_id::String,
    sess_id::String,
    condition::String
)::Dict{String,Any}

    tbl_dir = joinpath(res_dir, subj_id, sess_id, "tables")
    fig_dir = joinpath(res_dir, subj_id, sess_id, "figures")

    d = Dict{String,Any}()

    # Overview
    ov_path = joinpath(tbl_dir, "overview_$(condition).csv")
    isfile(ov_path) && (d["overview"] = CSV.read(ov_path, DataFrame))

    # QC
    qc_path = joinpath(tbl_dir, "qc_channels_$(condition).csv")
    isfile(qc_path) && (d["qc"] = CSV.read(qc_path, DataFrame))

    # Band power
    bp_path = joinpath(tbl_dir, "band_power_$(condition).csv")
    isfile(bp_path) && (d["band_power"] = CSV.read(bp_path, DataFrame))

    # Figuras disponibles
    d["figures"] = isdir(fig_dir) ?
        filter(f -> endswith(f, ".png"), readdir(fig_dir)) : String[]

    return d
end
