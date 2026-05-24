# NeuroMIND/src/SingleSubjectPipeline.jl
# Pipeline mínimo para análisis de un solo sujeto EEG.
# Fase inicial antes de escalar a análisis transversal y longitudinal.

using CairoMakie, Serialization

# ── Helper: hash de config ICA para invalidar caché ──────────

"""
    _ica_config_hash(cfg) -> String

Devuelve un string identificador del subconjunto de config que
afecta al resultado de ICA (perfil, n_components, seed, tol, max_iter).
El caché del ICAResult se invalida automáticamente si cambia alguno de estos.
"""
function _ica_config_hash(cfg::PipelineConfig)::String
    ica = get(cfg.ica, "", Dict{String,Any}())
    flt = get(cfg.filtering, "", Dict{String,Any}())
    # Incluir también params de filtrado porque ICA se calcula sobre la señal filtrada
    key = string(
        get(ica, "profile",      "default"),   "_",
        get(ica, "n_components", 0),           "_",
        get(ica, "seed",         42),          "_",
        get(ica, "tol",          1e-5),        "_",
        get(ica, "max_iter",     500),         "_",
        get(flt, "profile",      "default"),   "_",
        get(flt, "highpass_hz",  0.5),         "_",
        get(flt, "lowpass_hz",   150.0),       "_",
        get(flt, "filter_order", 4),
    )
    # Hash simple pero suficiente para detectar cambios de config
    h = zero(UInt32)
    for c in key
        h = xor(h * 31, UInt32(c))
    end
    return string(h, base=16)
end

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
    surr = get(raw, "surrogates",         Dict{String,Any}())
    paths_r = get(raw, "paths",           Dict{String,Any}())
    out  = get(raw, "output",             Dict{String,Any}())
    rec  = get(raw, "recording",          Dict{String,Any}())

    bands_raw = get(raw, "bands", Dict{String,Any}())
    bands = Dict{String,Tuple{Float64,Float64}}(
        k => (Float64(v[1]), Float64(v[2])) for (k, v) in bands_raw
    )

    overlap_s    = Float64(get(seg, "overlap_seconds", 0.0))
    seg_len_s    = Float64(get(seg, "segment_length_seconds", 2.0))
    overlap_frac = seg_len_s > 0.0 ? overlap_s / seg_len_s : 0.0
    bl_raw       = get(raw, "baseline", Dict{String,Any}())

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
            "profile"       => String(get(filt, "profile",        "default")),
            "highpass_hz"   => Float64(get(filt, "highpass_hz",     0.5)),
            "lowpass_hz"    => Float64(get(filt, "lowpass_hz",    150.0)),
            "notch_hz"      => Float64(get(filt, "notch_hz",       50.0)),
            "notch_bw_hz"   => Float64(get(filt, "notch_bw_hz",    1.0)),
            "bandreject_lo" => Float64(get(filt, "bandreject_lo",  99.5)),
            "bandreject_hi" => Float64(get(filt, "bandreject_hi", 100.5)),
            "filter_order"  => Int(get(filt,     "filter_order",     4)),
        ),
        Dict{String,Any}(
            "profile"        => String(get(seg, "profile",       "default")),
            "epoch_length_s" => seg_len_s,
            "epoch_overlap"  => overlap_frac,
            "min_epochs"     => Int(get(seg, "min_segments",       10)),
        ),
        Dict{String,Any}(
            "apply"            => true,
            "method"           => String(get(bl_raw, "method",           "mean")),
            "baseline_start_s" => Float64(get(bl_raw, "baseline_start_s", 0.0)),
            "baseline_end_s"   => Float64(get(bl_raw, "baseline_end_s",   0.10)),
            "n_passes"         => Int(get(bl_raw,     "n_passes",          1)),
        ),
        Dict{String,Any}(
            "profile"                => String(get(ar, "profile",                "default")),
            "amplitude_threshold_uv" => Float64(get(ar, "amplitude_threshold_uv", 100.0)),
            "gradient_threshold_uv"  => Float64(get(ar, "gradient_threshold_uv",   50.0)),
            "min_amplitude_uv"       => Float64(get(ar, "min_amplitude_uv",        -70.0)),
            "max_amplitude_uv"       => Float64(get(ar, "max_amplitude_uv",         70.0)),
            "n_channels_used"        => Int(get(ar,    "n_channels_used",           30)),
            "use_gradient"           => Bool(get(ar,   "use_gradient",             true)),
            "before_event_ms"        => Int(get(ar,    "before_event_ms",           200)),
            "after_event_ms"         => Int(get(ar,    "after_event_ms",            300)),
            "enabled"                => Bool(get(ar,   "enabled",                  true)),
        ),
        Dict{String,Any}(
            "profile"            => String(get(get(raw, "ica", Dict()), "profile", "default")),
            "n_components"       => get(get(raw, "ica", Dict()), "n_components", 30),
            "method"             => "fastica",
            "max_iter"           => Int(get(get(raw, "ica", Dict()), "max_iter", 500)),
            "tol"                => Float64(get(get(raw, "ica", Dict()), "tol", 1e-5)),
            "seed"               => Int(get(get(raw, "ica", Dict()), "seed", 42)),
            "artifact_threshold" => Float64(get(get(raw, "ica", Dict()), "artifact_threshold", 1.5)),
        ),
        Dict{String,Any}(
            "nfft"       => Int(get(sp, "nfft", 512)),
            "window_pct" => Float64(get(sp, "window_pct", 10.0)),
        ),
        bands,
        Dict{String,Any}(
            "filter_order" => Int(get(conn, "filter_order", 8)),
            "use_csd"      => Bool(get(conn, "use_csd", false)),
        ),
        Dict{String,Any}(
            "enabled"      => Bool(get(surr, "enabled",      false)),
            "n_surrogates" => Int( get(surr, "n_surrogates",  200)),
            "method"       => String(get(surr, "method",  "phase_shuffle")),
            "alpha"        => Float64(get(surr, "alpha",    0.05)),
            "fdr_method"   => String(get(surr, "fdr_method",  "bh")),
            "seed"         => Int(get(surr, "seed",            42)),
        ),
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

    isfile(meta_path) || error("Metadata no encontrado: $meta_path")

    # ── Metadata JSON ──────────────────────────────────────────
    meta_raw = _parse_bids_json(meta_path)

    # ── Formato BrainVision: carga desde binario nativo ───────
    # Cuando build_bids_full.jl genera el BIDS ligero (sin TSV),
    # pone data_format="brainvision" y vhdr_path en el JSON.
    data_fmt = _parse_json_string_field(meta_path, "data_format")
    if data_fmt == "brainvision" && !isfile(data_path)
        vhdr_path = _parse_json_string_field(meta_path, "vhdr_path")
        if !isempty(vhdr_path) && isfile(vhdr_path)
            ch_pos = isfile(elec_path) ? _load_electrode_positions(elec_path) : nothing
            return load_eeg_brainvision(vhdr_path, subj_id, sess_id, task;
                                        run=run, ch_pos=ch_pos)
        end
        error("data_format=brainvision pero vhdr_path no existe o no está en metadata: $meta_path")
    end

    isfile(data_path) || error("EEG no encontrado: $data_path")

    # ── Formato TSV (comportamiento original) ─────────────────
    fs = Float64(get(meta_raw, "fs", cfg.recording["fs"]))

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

# ── Helper: extrae un campo string de un JSON sin JSON3 ───────

function _parse_json_string_field(json_path::String, field::String)::String
    isfile(json_path) || return ""
    raw = read(json_path, String)
    pat = Regex("\"$(field)\"\\s*:\\s*\"([^\"]+)\"")
    m   = match(pat, raw)
    m === nothing ? "" : String(m[1])
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
    _filt_profile = get(cfg.filtering, "profile", "default")
    _log(log_io, "  Perfil: $(_filt_profile)")
    for step in describe_filter_chain(cfg)
        _log(log_io, "    [$(step.step)] $(step.name) $(step.freq)  ord=$(step.order)  método=$(step.method)")
    end

    # ── ICA: antes de segmentar (señal continua filtrada) ─────
    println("[4/8] ICA (separación de fuentes)...")
    _log(log_io, "\n[4/8] ICA")
    t_ica     = now()
    cache_dir = joinpath(subj_dir, "cache")
    ica_cache = joinpath(cache_dir, "ica_result.jls")
    ica_cfg_hash = _ica_config_hash(cfg)
    ica_hash_path = joinpath(cache_dir, "ica_config.hash")

    # Cache válido si: archivo existe Y el hash de config ICA coincide
    cache_valid = isfile(ica_cache) &&
                  isfile(ica_hash_path) &&
                  strip(read(ica_hash_path, String)) == ica_cfg_hash

    ica_result = if cache_valid
        println("  ↩ ICA cargado desde caché (config sin cambios)")
        _log(log_io, "  ICA cargado desde caché")
        try Serialization.deserialize(ica_cache) catch; nothing end
    else
        res = try
            run_ica(rec_filt, cfg)
        catch e
            @warn "ICA falló: $e · continuando sin ICA"
            nothing
        end
        if res !== nothing
            mkpath(cache_dir)
            Serialization.serialize(ica_cache, res)
            write(ica_hash_path, ica_cfg_hash)
        end
        res
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
    t_seg       = now()
    n_passes    = Int(get(cfg.baseline, "n_passes", 1))
    seg_profile = String(get(cfg.segmentation, "profile", "default"))
    ar_profile  = String(get(cfg.artifact_rejection, "profile", "default"))
    epochs_raw  = segment_recording(rec_ica, cfg)
    epochs_bl1  = apply_baseline(epochs_raw, cfg)   # 1ª pasada: todos los epochs
    epochs_ar   = reject_artifacts(epochs_bl1, cfg) # solo epochs válidos
    epochs      = n_passes >= 2 ? apply_baseline(epochs_ar, cfg) : epochs_ar   # 2ª pasada post-AR
    n_total    = n_epochs(epochs_raw)
    n_valid    = epochs.n_valid
    n_rejected = length(epochs_ar.rejected_idx)
    seg_dur    = round(Dates.value(now() - t_seg) / 1000, digits=1)
    _log(log_io, "  Perfil segmentación: $(seg_profile) | Perfil AR: $(ar_profile) | Baseline passes: $(n_passes)")
    _log(log_io, "  Segmentos totales: $(n_total) | Válidos: $(n_valid) | Rechazados: $(n_rejected) ($(round(100*n_rejected/n_total, digits=1))%)")
    _log(log_io, "  Duración segmentación: $(seg_dur) s")
    seg_signal_label = (ica_result !== nothing && !isempty(ica_result.rejected_components)) ? "ICA-limpiada" : "filtrada"
    _save_segmentation_results(epochs_bl1, epochs, rec_ica, cfg, export_dir, t_seg, log_io, seg_signal_label)

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

    # ── 9b. Surrogates (opcional — activar en config [surrogates] enabled=true) ─
    if Bool(get(cfg.surrogates, "enabled", false))
        println("[SUR] Inferencia por surrogates…")
        _log(log_io, "\n[SUR] Inferencia por surrogates")
        n_sur  = Int(get(cfg.surrogates, "n_surrogates", 200))
        _log(log_io, "  Método: $(get(cfg.surrogates,"method","phase_shuffle")) | N=$(n_sur) | FDR: $(get(cfg.surrogates,"fdr_method","bh"))")
        surr_results = SurrogateResult[]
        for band in sort(collect(keys(conn.matrices)))
            try
                sr = surrogate_test(epochs_conn, conn, band, cfg)
                push!(surr_results, sr)
                n_sig = count(sr.sig_mask) ÷ 2
                _log(log_io, "  Banda $(lpad(band,9)): $(n_sig) pares significativos  FDR-thr=$(round(sr.fdr_threshold,digits=4))")
            catch e
                @warn "Surrogate fallido para $band: $e"
                _log(log_io, "  WARN: surrogate $(band) fallido: $e")
            end
        end
        if !isempty(surr_results)
            _save_surrogate_results(surr_results, conn, export_dir, cfg, log_io)
        end
    end

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
    seg_profile  = String(get(cfg.segmentation,       "profile",  "default"))
    ar_profile   = String(get(cfg.artifact_rejection, "profile",  "default"))
    amp_thr      = Float64(get(cfg.artifact_rejection, "amplitude_threshold_uv", 100.0))
    grad_thr     = Float64(get(cfg.artifact_rejection, "gradient_threshold_uv",   50.0))
    min_amp_uv   = Float64(get(cfg.artifact_rejection, "min_amplitude_uv",        -70.0))
    max_amp_uv   = Float64(get(cfg.artifact_rejection, "max_amplitude_uv",          70.0))
    n_ch_ar      = Int(get(cfg.artifact_rejection,     "n_channels_used",           30))
    use_grad_ar  = Bool(get(cfg.artifact_rejection,    "use_gradient",             true))
    bl_method    = String(get(cfg.baseline, "method",           "mean"))
    bl_start_s   = Float64(get(cfg.baseline, "baseline_start_s", 0.0))
    bl_end_s     = Float64(get(cfg.baseline, "baseline_end_s",   0.10))
    n_bl_passes  = Int(get(cfg.baseline,     "n_passes",          1))
    ts           = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    dur_s        = round(Dates.value(now() - t_start) / 1000, digits=1)

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
  "profile": "$(seg_profile)",
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
  "baseline_method": "$(bl_method)",
  "baseline_start_s": $(bl_start_s),
  "baseline_end_s": $(bl_end_s),
  "n_baseline_passes": $(n_bl_passes),
  "artifact_profile": "$(ar_profile)",
  "min_amplitude_uv": $(min_amp_uv),
  "max_amplitude_uv": $(max_amp_uv),
  "n_channels_used_for_rejection": $(n_ch_ar),
  "use_gradient": $(use_grad_ar),
  "amp_threshold_uv": $(amp_thr),
  "grad_threshold_uv": $(grad_thr),
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
                         n_total, n_valid, n_rejected, log_io,
                         rec.meta.n_channels)
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
    log_io::IO,
    n_channels_total::Int = 0   # número total de canales (de la señal de entrada)
)
    ar_cfg      = cfg.artifact_rejection
    ar_profile  = String(get(ar_cfg, "profile",                "default"))
    amp_thresh  = Float64(get(ar_cfg, "amplitude_threshold_uv", 100.0))
    grad_thresh = Float64(get(ar_cfg, "gradient_threshold_uv",   50.0))
    min_amp_uv  = Float64(get(ar_cfg, "min_amplitude_uv",        -70.0))
    max_amp_uv  = Float64(get(ar_cfg, "max_amplitude_uv",         70.0))
    n_ch_used   = Int(get(ar_cfg,    "n_channels_used",           30))
    use_grad    = Bool(get(ar_cfg,   "use_gradient",             true))
    bef_ms      = Int(get(ar_cfg,    "before_event_ms",           200))
    aft_ms      = Int(get(ar_cfg,    "after_event_ms",            300))
    ret_pct     = round(100.0 * n_valid / max(n_total, 1), digits=1)

    # For default profile, min/max are symmetric around amp_thresh
    if ar_profile != "eeg_julia"
        min_amp_uv = -amp_thresh
        max_amp_uv =  amp_thresh
    end
    # n_channels_used capped to actual total if known
    n_ch_reported = (n_channels_total > 0 && ar_profile == "eeg_julia") ?
                    min(n_ch_used, n_channels_total) : n_ch_used

    rej_mask   = qr.status .== "rejected"
    rej_df     = qr[rej_mask, :]
    n_rej_amp  = count(==("amplitude"), qr.rejection_reason)
    n_rej_grad = count(==("gradient"),  qr.rejection_reason)

    # ── Estadísticas P2P ─────────────────────────────────────────────────────
    p2p_vals = hasproperty(qr, :p2p_uv) ? Float64.(qr.p2p_uv) : Float64[]
    p2p_mean = isempty(p2p_vals) ? 0.0 : round(mean(p2p_vals),    digits=2)
    p2p_std  = isempty(p2p_vals) ? 0.0 : round(std(p2p_vals),     digits=2)
    p2p_max  = isempty(p2p_vals) ? 0.0 : round(maximum(p2p_vals), digits=2)
    p2p_thresh_2sd = round(p2p_mean + 2.0 * p2p_std, digits=2)

    # Histograma P2P (20 bins)
    p2p_hist_json = ""
    if !isempty(p2p_vals) && p2p_max > 0
        bin_w  = p2p_max / 20.0
        hist_c = zeros(Int, 20)
        for v in p2p_vals
            b = min(floor(Int, v / bin_w), 19)
            hist_c[b+1] += 1
        end
        p2p_hist_json = join(
            ["[$(round((i-1)*bin_w, digits=1)),$(hist_c[i])]" for i in 1:20], ",")
    end

    # ── Tabla de canales más afectados ───────────────────────────────────────
    # Preferir channels_violating (lista real de canales), si no worst_channel
    ch_bad = Dict{String, @NamedTuple{n::Int, amp::Int, grad::Int}}()
    if hasproperty(rej_df, :channels_violating)
        for row in eachrow(rej_df)
            viols = filter(!isempty, split(string(row.channels_violating), ";"))
            for ch in viols
                prev = get(ch_bad, ch, (n=0, amp=0, grad=0))
                na   = string(row.rejection_reason) == "amplitude" ? prev.amp + 1 : prev.amp
                ng   = string(row.rejection_reason) == "gradient"  ? prev.grad + 1 : prev.grad
                ch_bad[ch] = (n=prev.n + 1, amp=na, grad=ng)
            end
        end
    elseif hasproperty(rej_df, :worst_channel)
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
            channel     = [r[1] for r in ch_rows_sorted],
            n_bad       = [r[2].n   for r in ch_rows_sorted],
            pct_bad     = [round(100.0 * r[2].n / max(n_total,1), digits=1) for r in ch_rows_sorted],
            main_reason = [r[2].amp >= r[2].grad ? "amplitude" : "gradient" for r in ch_rows_sorted],
        )

    # ── 4) artifact_rejection_summary.json ───────────────────────────────────
    open(joinpath(export_dir, "artifact_rejection_summary.json"), "w") do f
        write(f, """{
  "profile": "$(ar_profile)",
  "n_total": $(n_total),
  "n_valid": $(n_valid),
  "n_rejected": $(n_rejected),
  "retention_pct": $(ret_pct),
  "n_rejected_amplitude": $(n_rej_amp),
  "n_rejected_gradient": $(n_rej_grad),
  "min_amplitude_uv": $(min_amp_uv),
  "max_amplitude_uv": $(max_amp_uv),
  "amp_threshold_uv": $(amp_thresh),
  "grad_threshold_uv": $(grad_thresh),
  "use_gradient": $(use_grad),
  "n_channels_used": $(n_ch_reported),
  "n_channels_total": $(n_channels_total),
  "before_event_ms": $(bef_ms),
  "after_event_ms": $(aft_ms),
  "before_after_applied": false,
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
    # Seleccionar columnas canónicas (orden limpio para CSV)
    if !isempty(rej_df)
        wanted = [:epoch, :start_s, :end_s, :duration_s, :quality,
                  :status, :rejection_reason,
                  :max_amp_uv, :min_amp_uv, :p2p_uv,
                  :worst_channel, :channels_violating, :max_grad_uv]
        present = [c for c in wanted if hasproperty(rej_df, c)]
        CSV.write(joinpath(export_dir, "rejected_segments.csv"), rej_df[:, present])
        _log(log_io, "  Guardado: rejected_segments.csv ($(nrow(rej_df)) rechazados)")
    end

    # ── 6) channel_artifact_summary.csv ──────────────────────────────────────
    if !isempty(ca_df)
        CSV.write(joinpath(export_dir, "channel_artifact_summary.csv"), ca_df)
        _log(log_io, "  Guardado: channel_artifact_summary.csv")
    end
end

# ─── Extras espectrales: summary.json, regional_psd.csv, spectral_indices.csv ───

function _save_spectral_extras(
    spectra::SpectralResult,
    cfg::PipelineConfig,
    export_dir::String,
    log_io::IO
)
    ch_names  = spectra.meta.channel_names
    n_ch      = length(ch_names)
    fs_val    = Float64(spectra.meta.fs)
    nfft_used = Int(get(spectra.params, "nfft", length(spectra.freqs) * 2 - 2))
    win_pct   = Float64(get(spectra.params, "window_pct", 10.0))
    n_freqs   = length(spectra.freqs)
    delta_f   = n_freqs > 1 ? round(spectra.freqs[2] - spectra.freqs[1], digits=4) : 1.0
    epoch_s   = Float64(get(cfg.segmentation, "epoch_length_s",
                    get(cfg.segmentation, "segment_length_seconds", 1.0)))
    ts        = string(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
    ref_val   = String(get(cfg.recording, "reference", "average"))
    prof_val  = String(get(cfg.filtering, "profile", "default"))

    # ── 1) spectral_summary.json ──────────────────────────────────────────────
    # Total power = integral of global mean PSD
    psd_global_mean = vec(mean(spectra.psd, dims=1))   # (n_bins,)
    total_pw        = round(sum(psd_global_mean) * delta_f, digits=6)

    # Band relative power (fraction of band means)
    bnames_sorted    = sort(collect(keys(spectra.band_power)))
    band_mean_dict   = Dict{String,Float64}()
    for bname in bnames_sorted
        valid = filter(!isnan, spectra.band_power[bname])
        isempty(valid) || (band_mean_dict[bname] = mean(valid))
    end
    total_band_sum = max(sum(values(band_mean_dict)), 1e-12)

    band_json_parts = String[]
    for bname in bnames_sorted
        bvals = spectra.band_power[bname]
        valid = filter(!isnan, bvals)
        isempty(valid) && continue
        bpct  = round(get(band_mean_dict, bname, 0.0) / total_band_sum * 100, digits=2)
        push!(band_json_parts,
            "  \"$(bname)_mean\":$(round(mean(valid),digits=6))," *
            "\"$(bname)_std\":$(round(length(valid)>1 ? std(valid) : 0.0, digits=6))," *
            "\"$(bname)_median\":$(round(median(valid),digits=6))," *
            "\"$(bname)_min\":$(round(minimum(valid),digits=6))," *
            "\"$(bname)_max\":$(round(maximum(valid),digits=6))," *
            "\"$(bname)_pct\":$(bpct)"
        )
    end

    open(joinpath(export_dir, "spectral_summary.json"), "w") do f
        write(f, """{
  "timestamp": "$(ts)",
  "method": "fft_hamming_taper",
  "window": "hamming_taper",
  "window_pct": $(win_pct),
  "nfft": $(nfft_used),
  "epoch_length_s": $(epoch_s),
  "n_epochs": $(spectra.n_epochs_used),
  "fs": $(fs_val),
  "delta_f": $(delta_f),
  "freq_range_lo": $(round(spectra.freqs[1], digits=3)),
  "freq_range_hi": $(round(spectra.freqs[end], digits=3)),
  "n_freq_bins": $(n_freqs),
  "n_channels": $(n_ch),
  "reference": "$(ref_val)",
  "profile": "$(prof_val)",
  "total_power_uv2": $(total_pw),
$(join(band_json_parts, ",\n"))
}""")
    end
    _log(log_io, "  Guardado: spectral_summary.json")

    # ── 2) regional_psd.csv ───────────────────────────────────────────────────
    regions_map = Dict{String,Vector{String}}(
        "frontal"   => ["Fp1","Fp2","F3","F4","Fz","F7","F8","AF3","AF4","AF7","AF8"],
        "central"   => ["C3","C4","Cz","FC1","FC2","FC5","FC6"],
        "parietal"  => ["P3","P4","Pz","P7","P8","CP1","CP2","CP5","CP6"],
        "occipital" => ["O1","O2","Oz","PO3","PO4","PO7","PO8"],
    )
    reg_region = String[]; reg_band = String[]
    reg_mean   = Float64[]; reg_std = Float64[]; reg_n = Int[]

    for reg_ord in ["frontal","central","parietal","occipital"]
        reg_chs = regions_map[reg_ord]
        idxs    = [i for (i,c) in enumerate(ch_names) if c in reg_chs]
        isempty(idxs) && continue
        for bname in bnames_sorted
            bvals = spectra.band_power[bname]
            valid = [bvals[i] for i in idxs if !isnan(bvals[i])]
            isempty(valid) && continue
            push!(reg_region, reg_ord);  push!(reg_band, bname)
            push!(reg_mean, round(mean(valid), digits=6))
            push!(reg_std,  round(length(valid) > 1 ? std(valid) : 0.0, digits=6))
            push!(reg_n,    length(valid))
        end
    end

    if !isempty(reg_region)
        CSV.write(joinpath(export_dir, "regional_psd.csv"),
                  DataFrame(region=reg_region, band=reg_band,
                            mean_power=reg_mean, std_power=reg_std, n_channels=reg_n))
        _log(log_io, "  Guardado: regional_psd.csv")
    end

    # ── 3) spectral_indices.csv ───────────────────────────────────────────────
    get_bp(name) = get(spectra.band_power, name, fill(NaN, n_ch))

    alpha   = get_bp("ALPHA")
    theta   = get_bp("THETA")
    gamma   = get_bp("GAMMA")
    bl      = get_bp("BETA_LOW")
    bm      = get_bp("BETA_MID")
    bh_     = get_bp("BETA_HIGH")

    beta = map(1:n_ch) do i
        v = filter(!isnan, [bl[i], bm[i], bh_[i]])
        isempty(v) ? NaN : mean(v)
    end

    a_range = get(cfg.bands, "ALPHA", (7.8, 11.7))
    a_idx   = findall(f -> f >= a_range[1] && f <= a_range[2], spectra.freqs)
    safe(x, y) = (!isnan(x) && !isnan(y) && abs(y) > 1e-12) ? round(x/y, digits=3) : 0.0

    idx_ch = String[]; idx_at = Float64[]; idx_ba = Float64[]
    idx_tb = Float64[]; idx_ga = Float64[]
    idx_ph = Float64[]; idx_pu = Float64[]

    for i in 1:n_ch
        pk_hz = 0.0; pk_uv = 0.0
        if !isempty(a_idx)
            psd_a = spectra.psd[i, a_idx]
            if !any(isnan, psd_a)
                pi2   = argmax(psd_a)
                pk_hz = round(spectra.freqs[a_idx[pi2]], digits=2)
                pk_uv = round(psd_a[pi2], digits=6)
            end
        end
        push!(idx_ch, ch_names[i])
        push!(idx_at, safe(alpha[i], theta[i]))
        push!(idx_ba, safe(beta[i],  alpha[i]))
        push!(idx_tb, safe(theta[i], beta[i]))
        push!(idx_ga, safe(gamma[i], alpha[i]))
        push!(idx_ph, pk_hz); push!(idx_pu, pk_uv)
    end

    CSV.write(joinpath(export_dir, "spectral_indices.csv"),
              DataFrame(channel=idx_ch,
                        alpha_theta=idx_at, beta_alpha=idx_ba,
                        theta_beta=idx_tb, gamma_alpha=idx_ga,
                        peak_alpha_hz=idx_ph, peak_alpha_uv2=idx_pu))
    _log(log_io, "  Guardado: spectral_indices.csv")
end

# ─── Extras de conectividad: connectivity_summary.json + network_metrics.csv ───

function _save_connectivity_extras(
    conn::ConnectivityMatrix,
    export_dir::String,
    log_io::IO
)
    bands   = sort(collect(keys(conn.matrices)))
    ch_names = conn.channel_names
    n        = length(ch_names)
    thr      = 0.1   # densidad a umbral 0.1

    # ── Flat JSON summary ─────────────────────────────────────
    pairs = String[
        "\"method\": \"$(conn.method)\"",
        "\"space\": \"$(conn.space)\"",
        "\"n_channels\": $n",
        "\"n_epochs_used\": $(conn.n_epochs_used)",
        "\"timestamp\": \"$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))\"",
        "\"threshold\": $thr",
    ]
    for band in bands
        W = conn.matrices[band]
        vals = Float64[]
        for i in 1:n, j in (i+1):n
            push!(vals, W[i,j])
        end
        isempty(vals) && continue
        n_edges = length(vals)
        n_above = count(v -> v > thr, vals)
        push!(pairs, "\"$(band)_mean\": $(round(mean(vals),    digits=4))")
        push!(pairs, "\"$(band)_std\": $(round(std(vals),     digits=4))")
        push!(pairs, "\"$(band)_median\": $(round(median(vals), digits=4))")
        push!(pairs, "\"$(band)_max\": $(round(maximum(vals),  digits=4))")
        push!(pairs, "\"$(band)_n_edges\": $n_edges")
        push!(pairs, "\"$(band)_n_above\": $n_above")
        push!(pairs, "\"$(band)_density\": $(round(n_above / max(1, n_edges), digits=4))")
    end
    open(joinpath(export_dir, "connectivity_summary.json"), "w") do io
        write(io, "{\n" * join(["  " * p for p in pairs], ",\n") * "\n}")
    end

    # ── network_metrics.csv (strength + degree agregados por banda) ──
    strength_agg = zeros(Float64, n)
    degree_agg   = zeros(Int,     n)
    n_counted    = 0
    for (_, W) in conn.matrices
        for i in 1:n
            for j in 1:n
                i == j && continue
                strength_agg[i] += W[i,j]
                W[i,j] > thr && (degree_agg[i] += 1)
            end
        end
        n_counted += 1
    end
    if n_counted > 0
        strength_agg ./= n_counted
        degree_agg    = round.(Int, degree_agg ./ n_counted)
    end
    max_str = maximum(abs.(strength_agg))
    nm_df = DataFrame(
        channel       = ch_names,
        strength      = round.(strength_agg, digits=4),
        degree        = degree_agg,
        norm_strength = round.(max_str > 0 ? strength_agg ./ max_str : zeros(n), digits=4),
    )
    CSV.write(joinpath(export_dir, "network_metrics.csv"), nm_df)
    _log(log_io, "  Guardado: connectivity_summary.json + network_metrics.csv")
end

# ─── BH q-valores individuales (para surrogates) ─────────────

function _bh_qvalues(p_vec::Vector{Float64})::Vector{Float64}
    m = length(p_vec)
    m == 0 && return Float64[]
    order  = sortperm(p_vec)
    rank   = invperm(order)           # rank[i] = posición de p_vec[i] al ordenar
    q_vec  = p_vec .* m ./ rank      # BH q-value
    # Monotonicidad de derecha a izquierda sobre p ordenado
    q_sorted = q_vec[order]
    for i in (m-1):-1:1
        q_sorted[i] = min(q_sorted[i], q_sorted[i+1])
    end
    q_out         = zeros(Float64, m)
    q_out[order]  = q_sorted
    return min.(q_out, 1.0)
end

# ─── Guardar resultados de surrogates ──────────────────────────

function _save_surrogate_results(
    surr_results::Vector{SurrogateResult},
    conn::ConnectivityMatrix,
    export_dir::String,
    cfg::PipelineConfig,
    log_io::IO
)
    ch_names  = conn.channel_names
    n         = length(ch_names)
    alpha     = Float64(get(cfg.surrogates, "alpha", 0.05))
    method    = String(get(cfg.surrogates, "method", "phase_shuffle"))
    n_sur_cfg = Int(get(cfg.surrogates, "n_surrogates", 200))
    fdr_meth  = String(get(cfg.surrogates, "fdr_method", "bh"))
    seed_v    = Int(get(cfg.surrogates, "seed", 42))

    upper_idx = [(i,j) for i in 1:n for j in (i+1):n]

    # Acumular conexiones significativas globales
    all_sig_rows = NamedTuple[]

    # Métricas QC + summary por banda
    qc_bands = Dict{String,Any}[]

    for sr in surr_results
        band  = sr.band
        W_obs = sr.observed
        p_mat = sr.p_values
        null_d = sr.null_distribution   # (n_ch, n_ch, n_sur)

        p_vec = [p_mat[i,j] for (i,j) in upper_idx]
        q_vec = _bh_qvalues(p_vec)

        # ── Observado ─────────────────────────────────────────
        obs_df = DataFrame(hcat(ch_names, W_obs), vcat(["channel"], ch_names))
        CSV.write(joinpath(export_dir, "wpli_observed_$(band).csv"), obs_df)

        # ── p-values ──────────────────────────────────────────
        p_df = DataFrame(hcat(ch_names, p_mat), vcat(["channel"], ch_names))
        CSV.write(joinpath(export_dir, "wpli_pvalues_$(band).csv"), p_df)

        # ── q-values (matriz simétrica) ───────────────────────
        q_mat = zeros(Float64, n, n)
        for (k,(i,j)) in enumerate(upper_idx)
            q_mat[i,j] = q_vec[k]; q_mat[j,i] = q_vec[k]
        end
        q_df = DataFrame(hcat(ch_names, q_mat), vcat(["channel"], ch_names))
        CSV.write(joinpath(export_dir, "wpli_qvalues_$(band).csv"), q_df)

        # ── Máscara significativa ─────────────────────────────
        sig_int = Int.(sr.sig_mask)
        mask_df = DataFrame(hcat(ch_names, sig_int), vcat(["channel"], ch_names))
        CSV.write(joinpath(export_dir, "wpli_significant_$(band).csv"), mask_df)

        # ── Estadísticas de la distribución nula por par ──────
        null_mean_mat = dropdims(mean(null_d, dims=3), dims=3)
        null_std_mat  = dropdims(std( null_d, dims=3), dims=3)

        null_rows = [(
            ch_a      = ch_names[i],
            ch_b      = ch_names[j],
            wpli_obs  = round(W_obs[i,j], digits=4),
            null_mean = round(null_mean_mat[i,j], digits=4),
            null_std  = round(null_std_mat[i,j],  digits=4),
            p_value   = round(p_vec[k], digits=4),
            q_value   = round(q_vec[k], digits=4),
            z_score   = round(null_std_mat[i,j] > 0 ?
                            (W_obs[i,j] - null_mean_mat[i,j]) / null_std_mat[i,j] : 0.0,
                            digits=3)
        ) for (k,(i,j)) in enumerate(upper_idx)]
        CSV.write(joinpath(export_dir, "surrogate_null_stats_$(band).csv"),
                  DataFrame(null_rows))

        # ── Conexiones significativas (q < alpha) ─────────────
        for (k,(i,j)) in enumerate(upper_idx)
            q_vec[k] < alpha || continue
            zs = null_std_mat[i,j] > 0 ?
                 (W_obs[i,j] - null_mean_mat[i,j]) / null_std_mat[i,j] : 0.0
            push!(all_sig_rows, (
                ch_a     = ch_names[i],
                ch_b     = ch_names[j],
                band     = band,
                wpli_obs = round(W_obs[i,j], digits=4),
                p_value  = round(p_vec[k], digits=4),
                q_value  = round(q_vec[k], digits=4),
                z_score  = round(zs, digits=3)
            ))
        end

        # ── QC de esta banda ──────────────────────────────────
        n_sig   = count(q_vec .< alpha)
        n_total = length(p_vec)
        nzm_vals = [null_mean_mat[i,j] for (i,j) in upper_idx if !isnan(null_mean_mat[i,j])]
        nzs_vals = [null_std_mat[i,j]  for (i,j) in upper_idx if !isnan(null_std_mat[i,j])]
        push!(qc_bands, Dict{String,Any}(
            "band"           => band,
            "n_surrogates"   => sr.n_surrogates,
            "n_sig"          => n_sig,
            "n_total"        => n_total,
            "pct_sig"        => round(100.0*n_sig/max(1,n_total), digits=2),
            "mean_p"         => round(mean(p_vec), digits=4),
            "fdr_threshold"  => round(sr.fdr_threshold, digits=4),
            "mean_null_mean" => isempty(nzm_vals) ? 0.0 : round(mean(nzm_vals), digits=4),
            "mean_null_std"  => isempty(nzs_vals) ? 0.0 : round(mean(nzs_vals), digits=4),
            "obs_mean"       => round(mean([W_obs[i,j] for (i,j) in upper_idx]), digits=4),
            "obs_max"        => round(maximum(W_obs), digits=4),
        ))
    end

    # ── Guardar conexiones significativas ─────────────────────
    if !isempty(all_sig_rows)
        sig_df = DataFrame(all_sig_rows)
        sort!(sig_df, :q_value)
        CSV.write(joinpath(export_dir, "significant_connections.csv"), sig_df)
    else
        # CSV vacío con columnas correctas
        CSV.write(joinpath(export_dir, "significant_connections.csv"),
            DataFrame(ch_a=String[], ch_b=String[], band=String[],
                      wpli_obs=Float64[], p_value=Float64[],
                      q_value=Float64[], z_score=Float64[]))
    end

    # ── Guardar QC por banda ──────────────────────────────────
    qc_df = DataFrame(qc_bands)
    CSV.write(joinpath(export_dir, "surrogate_quality.csv"), qc_df)

    # ── Generar surrogate_summary.json ───────────────────────
    n_sig_total = length(all_sig_rows)
    n_pairs     = n * (n-1) ÷ 2
    json_pairs  = String[
        "\"method\": \"$(method)\"",
        "\"n_surrogates\": $(n_sur_cfg)",
        "\"alpha\": $(alpha)",
        "\"fdr_method\": \"$(fdr_meth)\"",
        "\"seed\": $(seed_v)",
        "\"n_channels\": $n",
        "\"n_total_pairs\": $(n_pairs)",
        "\"n_sig_total\": $(n_sig_total)",
        "\"timestamp\": \"$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))\"",
    ]
    for qc in qc_bands
        band = qc["band"]
        push!(json_pairs, "\"$(band)_n_sig\": $(qc["n_sig"])")
        push!(json_pairs, "\"$(band)_pct_sig\": $(qc["pct_sig"])")
        push!(json_pairs, "\"$(band)_mean_p\": $(qc["mean_p"])")
        push!(json_pairs, "\"$(band)_fdr_thr\": $(qc["fdr_threshold"])")
    end
    open(joinpath(export_dir, "surrogate_summary.json"), "w") do io
        write(io, "{\n" * join(["  " * p for p in json_pairs], ",\n") * "\n}")
    end

    _log(log_io, "  Surrogates guardados: $(n_sig_total) conexiones significativas en $(length(surr_results)) bandas")
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

    # ── Extras espectrales (summary + regional + indices) ─────
    _save_spectral_extras(spectra, cfg, export_dir, log_io)

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

    # ── Extras de conectividad (summary JSON + network metrics) ──
    _save_connectivity_extras(conn, export_dir, log_io)

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
    fs      = rec_before.meta.fs

    # ── Clasificación automática de componentes ───────────────
    artifact_thresh = 1.5
    feat_df = DataFrame()
    eval_df = DataFrame()
    try
        feat_df = compute_ica_features(
            Float64.(ica.mixing_matrix),
            Float64.(ica.activations),
            fs,
            rec_before.meta.channel_names
        )
        eval_df = evaluate_ica_components(feat_df; artifact_thresh=artifact_thresh)
        CSV.write(joinpath(export_dir, "ica_component_features.csv"), eval_df)
        _log(log_io, "  Guardado: ica_component_features.csv")
    catch e
        @warn "Clasificación ICA falló: $e"
    end

    # Resuelve artifact_type por componente (desde clasificación o fallback)
    artifact_types = if !isempty(eval_df) && hasproperty(eval_df, :artifact_type)
        String.(eval_df.artifact_type)
    else
        ["unknown" for _ in 1:n_comp]
    end
    # Componentes ya marcados en rejected_components → forzar etiqueta coherente
    for i in ica.rejected_components
        if 1 ≤ i ≤ length(artifact_types) && artifact_types[i] == "brain"
            artifact_types[i] = "unknown"
        end
    end

    # ── Tabla de componentes ──────────────────────────────────
    df = DataFrame(
        component     = 1:n_comp,
        variance_pct  = round.(var_pct, digits=2),
        rejected      = [i ∈ ica.rejected_components for i in 1:n_comp],
        label         = ["IC$(i)" for i in 1:n_comp],
        artifact_type = artifact_types,
    )
    CSV.write(joinpath(export_dir, "ica_components.csv"), df)
    _log(log_io, "  Guardado: ica_components.csv")

    # ── Matrices de mezcla / desmezcla ───────────────────────
    ch_names = rec_before.meta.channel_names
    mix_df = DataFrame(ica.mixing_matrix,
                       [Symbol("IC$(i)") for i in 1:n_comp])
    mix_df[!, :channel] = ch_names
    select!(mix_df, :channel, Not(:channel))
    CSV.write(joinpath(export_dir, "ica_mixing_matrix.csv"), mix_df)

    unmix_df = DataFrame(ica.unmixing_matrix,
                         [Symbol(c) for c in ch_names])
    unmix_df[!, :IC] = ["IC$(i)" for i in 1:n_comp]
    select!(unmix_df, :IC, Not(:IC))
    CSV.write(joinpath(export_dir, "ica_unmixing_matrix.csv"), unmix_df)
    _log(log_io, "  Guardado: ica_mixing_matrix.csv + ica_unmixing_matrix.csv")

    # ── Topomaps por componente ───────────────────────────────
    n_topo = 0
    figs_dir = joinpath(export_dir, "figures")
    mkpath(figs_dir)
    if rec_before.meta.channel_positions !== nothing
        n_topo = _save_ica_topomaps(ica, rec_before, figs_dir, log_io)
    else
        _log(log_io, "  Topomaps: sin posiciones de electrodos, omitido")
    end

    # ── Resumen JSON ──────────────────────────────────────────
    n_rej   = length(ica.rejected_components)
    var_rej = sum(var_pct[i] for i in ica.rejected_components; init=0.0)
    var_ret = round(100.0 - var_rej, digits=1)
    ts      = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")

    # Conteo de tipos de artefacto
    type_counts = Dict{String,Int}()
    for t in artifact_types; type_counts[t] = get(type_counts, t, 0) + 1; end
    type_json = join(["\"$(k)\": $(v)" for (k,v) in type_counts], ", ")

    open(joinpath(export_dir, "ica_summary.json"), "w") do f
        write(f, """{
  "n_components": $(n_comp),
  "n_rejected": $(n_rej),
  "n_accepted": $(n_comp - n_rej),
  "variance_retained": $(var_ret),
  "timestamp": "$(ts)",
  "duration_s": $(duration_s),
  "artifact_types": {$(type_json)},
  "n_topomaps": $(n_topo),
  "has_features": $((!isempty(eval_df)))
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

function _save_ica_topomaps(
    ica::ICAResult,
    rec::EEGRecording,
    figs_dir::String,
    log_io::IO
)::Int
    ch_pos   = rec.meta.channel_positions
    ch_names = rec.meta.channel_names
    n_comp   = size(ica.mixing_matrix, 2)
    n_saved  = 0

    for ic in 1:n_comp
        weights = Float64.(ica.mixing_matrix[:, ic])
        fname   = "ica_topomap_$(lpad(ic, 3, '0')).png"
        fpath   = joinpath(figs_dir, fname)
        try
            clim_val = maximum(abs.(weights))
            clim_val == 0.0 && (clim_val = 1.0)
            fig = plot_topomap(weights, ch_names, ch_pos;
                               title    = "IC $(lpad(ic, 3, '0'))",
                               colormap = CairoMakie.Reverse(:RdBu),
                               clims    = (-clim_val, clim_val))
            CairoMakie.save(fpath, fig)
            n_saved += 1
        catch e
            @warn "Topomap IC$ic falló: $e"
        end
    end
    _log(log_io, "  Topomaps: $(n_saved)/$(n_comp) guardados en figures/")
    return n_saved
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
