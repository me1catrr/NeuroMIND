# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Pipeline de un solo sujeto EEG
# ═══════════════════════════════════════════════════════════════
#
#  Orquestador canónico de 8 pasos (incluido desde NeuroMIND.jl):
#    [1/8] Carga        BIDS / BrainVision → EEGRecording
#    [2/8] QC           bad_ch (z-score) + amplitude_warning
#    [3/8] Filtrado     HP / LP / Notch / Bandreject
#    [4/8] ICA          FastICA en señal continua (ANTES de segmentar)
#    [5/8] Segmentación épocas + baseline + rechazo de artefactos
#    [6/8] Espectral    PSD + potencia por banda
#    [7/8] Conectividad wPLI (+ surrogates opcionales)
#    [8/8] Guardado     CSV / JSON / PNG → results/subjects/…
#
#  Config:     config/pipeline.toml
#  Invocado:   scripts/run_single_subject.jl
#              scripts/run_batch_pipeline.jl  (TOML temporal por sujeto)
#
#  API pública (exportada por NeuroMIND)
#  ─────────────────────────────────────
#    load_ss_config(path)               → PipelineConfig
#    detect_first_subject(bids_raw_dir) → (id, sess, task, run)
#    load_single_subject(…)             → EEGRecording
#    validate_channels(rec, elec_path)
#    run_single_subject_pipeline(path)  → ejecuta los 8 pasos
#
#  Helpers internos
#  ────────────────
#    Config / ICA
#      _ica_config_hash
#      _parse_json_string_field
#    Persistencia por etapa
#      _save_segmentation_results
#      _save_ar_results
#      _save_spectral_phase_results
#      _save_spectral_extras
#      _save_connectivity_phase_results
#      _save_connectivity_extras
#      _save_surrogate_results
#      _save_ica_results
#      _save_ica_topomaps
#      _save_raw_signal
#      _save_filtered_signal
#      _save_all_results
#      _save_config_snapshot
#    Figuras
#      _plot_signal_preview
#      _plot_band_power_summary
#    Índices / QC / utilidades
#      _update_subjects_index
#      _update_qc_decision_table
#      _bh_qvalues
#      _log
#
# ───────────────────────────────────────────────────────────────
#  Fichero    src/SingleSubjectPipeline.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 22-07-2026
# ───────────────────────────────────────────────────────────────

using CairoMakie, Serialization

# ── Helper: hash de config ICA para invalidar caché ──────────

"""
    _ica_config_hash(cfg) -> String

Devuelve un string identificador del subconjunto de config que
afecta al resultado de ICA (perfil, n_components, seed, tol, max_iter).
El caché del ICAResult se invalida automáticamente si cambia alguno de estos.
"""
function _ica_config_hash(cfg::PipelineConfig)::String
    ica = cfg.ica
    flt = cfg.filtering
    # Incluir también params de filtrado porque ICA se calcula sobre la señal filtrada
    key = string(
        get(ica, "profile",      "default"),   "_",
        get(ica, "n_components", 0),           "_",
        get(ica, "seed",         42),          "_",
        get(ica, "tol",          1e-5),        "_",
        get(ica, "max_iter",     500),         "_",
        get(ica, "max_attempts", 1),           "_",
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

Carga `config/pipeline.toml` y construye un `PipelineConfig`.

Lee todas las secciones del TOML unificado relevantes al pipeline:
  subject (vía raw en el runner), recording, qc, paths, output,
  filtering, ica, segmentation, baseline, artifact_rejection,
  spectral, bands, connectivity (+ fourier_csd / multitaper),
  montage, graph, surrogates.

Notas de mapeo de claves (compatibilidad):
  · `segment_length_seconds` → `epoch_length_s`
  · `overlap_seconds` (s)    → `epoch_overlap` (fracción)
  · `min_segments`           → `min_epochs`
  También acepta ya las claves canónicas `epoch_length_s` /
  `epoch_overlap` / `min_epochs` si están en el TOML.
"""
function load_ss_config(path::String)::PipelineConfig
    abs_path = isabspath(path) ? path : abspath(path)
    isfile(abs_path) || error("Config no encontrada: $abs_path")
    raw  = TOML.parsefile(abs_path)
    root = dirname(dirname(abs_path))   # sube config/ → NeuroMIND/

    filt    = get(raw, "filtering",          Dict{String,Any}())
    seg     = get(raw, "segmentation",       Dict{String,Any}())
    ar      = get(raw, "artifact_rejection", Dict{String,Any}())
    sp      = get(raw, "spectral",           Dict{String,Any}())
    conn    = get(raw, "connectivity",       Dict{String,Any}())
    surr    = get(raw, "surrogates",         Dict{String,Any}())
    paths_r = get(raw, "paths",              Dict{String,Any}())
    out     = get(raw, "output",             Dict{String,Any}())
    rec     = get(raw, "recording",          Dict{String,Any}())
    bl_raw  = get(raw, "baseline",           Dict{String,Any}())
    ica_raw = get(raw, "ica",                Dict{String,Any}())
    qc_raw  = get(raw, "qc",                 Dict{String,Any}())
    mont    = get(raw, "montage",            Dict{String,Any}())
    graph_r = get(raw, "graph",              Dict{String,Any}())

    bands_raw = get(raw, "bands", Dict{String,Any}())
    bands = Dict{String,Tuple{Float64,Float64}}(
        k => (Float64(v[1]), Float64(v[2])) for (k, v) in bands_raw
    )

    # Segmentación: aceptar claves canónicas o legacy
    seg_len_s = Float64(get(seg, "epoch_length_s",
                        get(seg, "segment_length_seconds", 2.0)))
    if haskey(seg, "epoch_overlap")
        overlap_frac = Float64(seg["epoch_overlap"])
    else
        overlap_s    = Float64(get(seg, "overlap_seconds", 0.0))
        overlap_frac = seg_len_s > 0.0 ? overlap_s / seg_len_s : 0.0
    end
    min_ep = Int(get(seg, "min_epochs", get(seg, "min_segments", 10)))

    # Conectividad: pasar wpli_method + subtablas; no inventar "method"="wpli"
    fourier_csd = Dict{String,Any}(get(conn, "fourier_csd", Dict{String,Any}()))
    multitaper  = Dict{String,Any}(get(conn, "multitaper",  Dict{String,Any}()))
    conn_dict = Dict{String,Any}(
        "wpli_method"              => String(get(conn, "wpli_method", "hilbert")),
        "use_dwpli"                => Bool(get(conn, "use_dwpli", false)),
        "min_cycles_for_wpli"      => Float64(get(conn, "min_cycles_for_wpli", 4.0)),
        "exclude_unreliable_bands" => Bool(get(conn, "exclude_unreliable_bands", false)),
        "use_csd"                  => Bool(get(conn, "use_csd", false)),
        "filter_order"             => Int(get(conn, "filter_order", 8)),
        "fourier_csd"              => fourier_csd,
        "multitaper"               => multitaper,
    )

    PipelineConfig(
        Dict{String,Any}("name" => "SingleSubject"),
        Dict{String,Any}("conditions" => ["EC"]),
        Dict{String,Any}(
            "bids_root"  => String(get(paths_r, "bids_root", "data/bids")),
            "results"    => String(get(paths_r, "results",   "results")),
            "data_cache" => "data/cache",
            "web_public" => "web/public",
        ),
        Dict{String,Any}(
            "fs"        => Float64(get(rec, "sampling_rate", 500.0)),
            "reference" => String(get(rec, "reference", "average")),
        ),
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
            "profile"        => String(get(seg, "profile", "default")),
            "epoch_length_s" => seg_len_s,
            "epoch_overlap"  => overlap_frac,
            "min_epochs"     => min_ep,
        ),
        Dict{String,Any}(
            "apply"            => Bool(get(bl_raw, "apply", true)),
            "method"           => String(get(bl_raw, "method", "mean")),
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
            "profile"            => String(get(ica_raw, "profile", "default")),
            "n_components"       => get(ica_raw, "n_components", 30),
            "method"             => "fastica",
            "max_iter"           => Int(get(ica_raw, "max_iter", 500)),
            "tol"                => Float64(get(ica_raw, "tol", 1e-5)),
            "seed"               => Int(get(ica_raw, "seed", 42)),
            "max_attempts"       => max(1, Int(get(ica_raw, "max_attempts", 1))),
            "verbose"            => Bool(get(ica_raw, "verbose", false)),
            "artifact_threshold" => Float64(get(ica_raw, "artifact_threshold", 1.5)),
            "auto_reject"        => Bool(get(ica_raw, "auto_reject", true)),
        ),
        Dict{String,Any}(
            "nfft"       => Int(get(sp, "nfft", 512)),
            "window_pct" => Float64(get(sp, "window_pct", 10.0)),
        ),
        bands,
        conn_dict,
        Dict{String,Any}(
            "enabled"      => Bool(get(surr, "enabled", false)),
            "n_surrogates" => Int(get(surr, "n_surrogates", 200)),
            "method"       => String(get(surr, "method", "circular_shift")),
            "alpha"        => Float64(get(surr, "alpha", 0.05)),
            "fdr_method"   => String(get(surr, "fdr_method", "bh")),
            "seed"         => Int(get(surr, "seed", 42)),
        ),
        Dict{String,Any}(
            "density"          => Float64(get(graph_r, "density", 0.1)),
            "threshold_method" => String(get(graph_r, "threshold_method", "proportional")),
        ),
        Dict{String,Any}(),   # clinical
        Dict{String,Any}(),   # longitudinal
        Dict{String,Any}(),   # statistics
        Dict{String,Any}(
            "figure_format" => String(get(out, "figure_format", "png")),
            "figure_dpi"    => Int(get(out, "figure_dpi", 150)),
        ),
        Dict{String,Any}(
            "bad_channel_zscore_threshold" => Float64(get(qc_raw, "bad_channel_zscore_threshold", 3.0)),
            "amplitude_warning_sigma_uv"   => Float64(get(qc_raw, "amplitude_warning_sigma_uv", 20.0)),
            # QC extendido (Report_Pre §Asimetría/Curtosis/HFNoise/SNR/Corr.)
            "welch_nfft"       => Int(get(qc_raw, "welch_nfft", 1024)),
            "hfnoise_band"     => Float64.(get(qc_raw, "hfnoise_band",     [25.0, 45.0])),
            "hfnoise_ref_band" => Float64.(get(qc_raw, "hfnoise_ref_band", [1.0,  45.0])),
            "snr_signal_band"  => String(get(qc_raw, "snr_signal_band", "ALPHA")),
            "snr_noise_band"   => String(get(qc_raw, "snr_noise_band",  "GAMMA")),
        ),
        Dict{String,Any}(
            "exclude_channels"    => String.(get(mont, "exclude_channels", String[])),
            "exclude_fp2"         => Bool(get(mont, "exclude_fp2", true)),
            "n_channels_analysis" => Int(get(mont, "n_channels_analysis", 30)),
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
    println("  " * "═"^60)
    println("  NeuroMIND · Pipeline individual · $(Dates.format(t0, "yyyy-mm-dd HH:MM:SS"))")
    println("  " * "═"^60)

    # ── 1. Configuración ──────────────────────────────────────
    # load_ss_config lee config/pipeline.toml (o la ruta pasada) y
    # construye un PipelineConfig: objeto único que consumen filtrado,
    # ICA, segmentación, PSD, wPLI, surrogates y el resto de fases.
    # [subject] no entra en el struct: se lee del TOML crudo más abajo.
    cfg = load_ss_config(config_path)
    raw_cfg     = TOML.parsefile(isabspath(config_path) ? config_path : abspath(config_path))
    sub_cfg     = get(raw_cfg, "subject",  Dict{String,Any}())
    dpi         = Int(get(cfg.export_cfg, "figure_dpi", 150))
    qc_cfg      = cfg.qc
    montage_cfg = cfg.montage

    # Resumen compacto en terminal (detalle fino → pipeline_log)
    _cfg_abs  = isabspath(config_path) ? config_path : abspath(config_path)
    _cfg_name = basename(_cfg_abs)
    _filt_p   = String(get(cfg.filtering, "profile", "default"))
    _ica_p    = String(get(cfg.ica, "profile", "default"))
    _seg_p    = String(get(cfg.segmentation, "profile", "default"))
    _ep_s     = Float64(get(cfg.segmentation, "epoch_length_s", 2.0))
    _ov       = Float64(get(cfg.segmentation, "epoch_overlap", 0.0))
    _wpli     = String(get(cfg.connectivity, "wpli_method", "hilbert"))
    _dwpli    = Bool(get(cfg.connectivity, "use_dwpli", false))
    _surr_on  = Bool(get(cfg.surrogates, "enabled", false))
    _nsurr    = Int(get(cfg.surrogates, "n_surrogates", 0))
    _n_an     = Int(get(montage_cfg, "n_channels_analysis", 30))
    _excl_cfg = String.(get(montage_cfg, "exclude_channels", String[]))
    _excl_fp2 = Bool(get(montage_cfg, "exclude_fp2", false))
    # Misma mecánica que el paso de segmentación: exclude_fp2=false
    # saca Fp2 de la lista efectiva → Fp2 se conserva.
    _eff_excl = _excl_fp2 ? _excl_cfg : filter(c -> c != "Fp2", _excl_cfg)
    _fp2_drop = _excl_fp2 && ("Fp2" ∈ _excl_cfg)
    _mont_excl = if _fp2_drop
        _other = filter(c -> c != "Fp2", _eff_excl)
        isempty(_other) ? "−Fp2" : "−Fp2," * join(_other, ",")
    elseif isempty(_eff_excl)
        "Fp2 incluido"
    else
        "−" * join(_eff_excl, ",")
    end
    _ov_str   = _ov > 0 ? " overlap=$(Int(round(_ov * 100)))%" : ""
    _surr_str = _surr_on ? "on (n=$(_nsurr))" : "off"
    _dwpli_str = _dwpli ? "on" : "off"
    _print_kv_table([
        ("Config",       _cfg_name),
        ("Filtering",    _filt_p),
        ("ICA",          _ica_p),
        ("Segmentation", "$(_seg_p) · $(_ep_s)s$(_ov_str)"),
        ("wPLI",         _wpli),
        ("dwPLI",        _dwpli_str),
        ("Surrogates",   _surr_str),
        ("Montage",      "$(_n_an) ch ($(_mont_excl))"),
    ])

    # ── 2. Detección / resolución del sujeto ──────────────────
    raw_dir = joinpath(bids_root_dir(cfg), "raw")
    sid     = String(get(sub_cfg, "subject_id", "auto"))
    sess_id = String(get(sub_cfg, "session_id", "T2"))
    task    = String(get(sub_cfg, "task", "eyesclosed"))
    run_n   = Int(get(sub_cfg, "run", 1))

    _sid_mode = sid == "auto" ? "auto" : "manual"
    if sid == "auto"
        sid, sess_id, task, run_n = detect_first_subject(raw_dir)
    end

    condition = task == "eyesclosed" ? "EC" : (task == "eyesopen" ? "EO" : task)

    # ── 3. Directorio de salida ────────────────────────────────
    # results/subjects/sub-{id}/ses-{sess}/{task}/
    export_dir = joinpath(results_dir(cfg), "subjects",
                          "sub-$(sid)", "ses-$(sess_id)", task)
    mkpath(joinpath(export_dir, "cache"))
    _out_rel = joinpath("results", "subjects", "sub-$(sid)", "ses-$(sess_id)", task)

    _print_kv_table([
        ("Subject",   "sub-$(sid)"),
        ("Session",   "ses-$(sess_id)"),
        ("Task",      task),
        ("Run",       string(run_n)),
        ("Selection", _sid_mode),
        ("Output",    "$(_out_rel)/"),
    ])

    log_path = joinpath(export_dir, "pipeline_log.txt")
    log_io   = open(log_path, "w")
    _log(log_io, "NeuroMIND pipeline — $(t0)")
    _log(log_io, "Sujeto: sub-$(sid) | Sesión: ses-$(sess_id) | Tarea: $(task)")
    _log(log_io, "Config: $(_cfg_abs)")
    _log(log_io, "  filt=$(_filt_p) · ica=$(_ica_p) · seg=$(_seg_p) $(_ep_s)s$(_ov_str)")
    _log(log_io, "  wPLI=$(_wpli) (dwPLI=$(_dwpli_str)) · surrogates=$(_surr_str) · montage=$(_n_an)ch ($(_mont_excl))")
    _log(log_io, "Salida: $(_out_rel)/")

    # ── 4. Carga de datos ─────────────────────────────────────
    # load_single_subject: ruta BrainVision preferente (.vhdr + .eeg
    # nativo); aplica la resolución µV/bit del .vhdr y devuelve un
    # EEGRecording (canales × muestras) en µV, sin filtrar ni limpiar.
    # Fallback: TSV BIDS si existe. Electrodos/QC van en pasos aparte.
    _log(log_io, "\n[1/8] Carga de datos")
    t_step = now()
    print("  [1/8] Carga EEG... ")
    rec = load_single_subject(cfg, sid, sess_id, task, run_n)
    _fs_disp = isinteger(rec.meta.fs) ? string(Int(rec.meta.fs)) : string(rec.meta.fs)
    _load_dur = round(Dates.value(now() - t_step) / 1000, digits=1)
    println("✓  ($(_load_dur) s)")
    _log(log_io, "  Canales: $(rec.meta.n_channels) | Muestras: $(n_samples(rec)) | fs: $(rec.meta.fs) Hz")
    _log(log_io, "  Duración: $(round(duration(rec), digits=2)) s")

    # Validación de electrodos
    elec_path = joinpath(bids_root_dir(cfg), "electrodes",
                         "sub-$(sid)_ses-$(sess_id)_electrodes.tsv")
    val = validate_channels(rec, elec_path)
    _elec_str = if !isfile(elec_path)
        "sin archivo"
    elseif val.ok
        "OK"
    else
        val.msg
    end
    _log(log_io, "  Electrodos: $(val.msg)")

    _print_kv_table([
        ("Channels",      string(rec.meta.n_channels)),
        ("Samples",       string(n_samples(rec))),
        ("Sampling rate", "$(_fs_disp) Hz"),
        ("Duration",      "$(round(duration(rec), digits=1)) s"),
        ("Electrodes",    _elec_str),
    ]; indent="        ")

    # Vista rápida en tabla: todos los canales × primeras muestras
    _n_prev = min(5, n_samples(rec))
    _name_w = max(5, maximum(length.(rec.meta.channel_names); init=4))
    _val_w  = 8
    _cell(s, w; left=false) = left ? rpad(s, w) : lpad(s, w)
    _rule(l, m, r) = l * "─"^(_name_w + 2) * join(fill(m * "─"^(_val_w + 2), _n_prev)) * r
    _sep(cells) = "│ " * join(cells, " │ ") * " │"
    println("        señal (µV · primeras $(_n_prev) muestras):")
    _top = "        " * _rule("┌", "┬", "┐")
    _mid = "        " * _rule("├", "┼", "┤")
    _bot = "        " * _rule("└", "┴", "┘")
    _hcells = [_cell("canal", _name_w; left=true);
               [_cell("s$j", _val_w) for j in 0:(_n_prev - 1)]...]
    println(_top)
    println("        " * _sep(_hcells))
    println(_mid)
    _log(log_io, "  señal (µV · primeras $(_n_prev) muestras):")
    for i in 1:rec.meta.n_channels
        _cells = [_cell(rec.meta.channel_names[i], _name_w; left=true);
                  [_cell(string(round(rec.data[i, j]; digits=2)), _val_w)
                   for j in 1:_n_prev]...]
        _line = "        " * _sep(_cells)
        println(_line)
        _log(log_io, _line)
    end
    println(_bot)


    # ── 5. QC básico + métricas extendidas (Report_Pre) ───────
    # Estadísticas por canal (μ, RMS, σ, rango, asimetría, curtosis)
    # + HFNoise / SNR(dB) / corr. media; flag de canales anómalos
    # por z-score de RMS; aviso si σ̄ cruda sugiere filtro online OFF.
    # No elimina canales aquí: solo informa (exclusión efectiva en
    # el montaje previo a segmentación).
    _log(log_io, "\n[2/8] QC de canales")
    t_step    = now()
    print("  [2/8] QC canales... ")
    # mean / RMS / std / rango / min / max / skewness / kurtosis + z-RMS
    qc_stats  = compute_channel_stats(rec)

    # QC P2 — HFNoise, SNR y correlación intercanal (señal raw)
    welch_nfft       = Int(get(qc_cfg, "welch_nfft", 1024))
    _hf_b            = get(qc_cfg, "hfnoise_band",     [25.0, 45.0])
    _hf_ref          = get(qc_cfg, "hfnoise_ref_band", [1.0,  45.0])
    hfnoise_band     = (Float64(_hf_b[1]),   Float64(_hf_b[2]))
    hfnoise_ref_band = (Float64(_hf_ref[1]), Float64(_hf_ref[2]))
    snr_signal_band  = String(get(qc_cfg, "snr_signal_band", "ALPHA"))
    snr_noise_band   = String(get(qc_cfg, "snr_noise_band",  "GAMMA"))
    qc_spectral = compute_channel_spectral_qc(
        rec, cfg;
        welch_nfft=welch_nfft,
        hfnoise_band=hfnoise_band,
        hfnoise_ref_band=hfnoise_ref_band,
        snr_signal_band=snr_signal_band,
        snr_noise_band=snr_noise_band,
    )
    qc_corr = compute_correlation_summary(rec)
    qc_stats[!, :hfnoise]       = qc_spectral.hfnoise
    qc_stats[!, :snr_db]        = qc_spectral.snr_db
    qc_stats[!, :mean_abs_corr] = qc_corr.mean_abs_corr
    qc_stats[!, :max_corr]      = qc_corr.max_corr
    qc_stats[!, :min_corr]      = qc_corr.min_corr

    z_thresh  = Float64(get(qc_cfg, "bad_channel_zscore_threshold", 3.0))
    # |rms_z| > umbral → canal sospechoso
    bad_ch    = flag_bad_channels(rec; z_threshold=z_thresh)
    qc_stats[!, :is_bad] = [ch in bad_ch for ch in qc_stats.channel]
    # σ̄ de todos los canales vs umbral de amplitude_warning
    amplitude_warn_thr = Float64(get(qc_cfg, "amplitude_warning_sigma_uv", 20.0))
    raw_sigma_mean     = mean(std(rec.data[ch, :]) for ch in 1:rec.meta.n_channels)
    amplitude_warning  = raw_sigma_mean > amplitude_warn_thr
    _bad_str  = isempty(bad_ch) ? "ninguno" : join(bad_ch, ", ")
    _amp_str  = amplitude_warning ? "⚠ amplitude_warning" : "✓ filtro OK"
    _qc_dur   = round(Dates.value(now() - t_step) / 1000, digits=1)
    println("✓  ($(_qc_dur) s)")
    if amplitude_warning
        _log(log_io, "  ⚠ amplitude_warning: σ̄_raw=$(round(raw_sigma_mean, digits=1)) µV > $(amplitude_warn_thr) µV (posible filtro online OFF)")
    else
        _log(log_io, "  amplitude_warning: false (σ̄_raw=$(round(raw_sigma_mean, digits=1)) µV)")
    end
    if isempty(bad_ch)
        _log(log_io, "  Sin canales sospechosos (umbral z=$(z_thresh))")
    else
        _log(log_io, "  Canales sospechosos (z>$(z_thresh)): " * join(bad_ch, ", "))
    end

    _print_kv_table([
        ("σ̄ raw",              "$(round(raw_sigma_mean, digits=1)) µV"),
        ("Amplitude warning",  _amp_str),
        ("Bad channels (z>$z_thresh)", _bad_str),
        ("HFNoise / SNR",      "welch_nfft=$welch_nfft · $snr_signal_band/$snr_noise_band"),
    ]; indent="        ")

    # Salidas de este paso (señal cruda + stats QC). raw_signal.csv
    # se escribe aquí — no tras ICA — para que exista aunque ICA falle.
    tables_dir = joinpath(export_dir, "tables"); mkpath(tables_dir)
    CSV.write(joinpath(tables_dir, "channel_statistics.csv"), qc_stats)
    CSV.write(joinpath(tables_dir, "qc_summary.csv"), qc_stats)
    _save_raw_signal(rec, export_dir, log_io)
    _print_kv_table([
        ("channel_statistics.csv", "tables/channel_statistics.csv"),
        ("qc_summary.csv",         "tables/qc_summary.csv"),
        ("raw_signal.csv",         "tables/raw_signal.csv"),
    ]; indent="        ", headers=("Archivo", "Ruta"))
    _log(log_io, "  Guardado: channel_statistics.csv, qc_summary.csv, raw_signal.csv")
    _log(log_io, "    columnas extendidas: skewness, kurtosis, hfnoise, snr_db, mean_abs_corr")

    # ── 6. Filtrado ───────────────────────────────────────────
    # filter_recording: cadena Butterworth según cfg.filtering["profile"].
    #   eeg_julia → Notch(filt) → Bandreject(filt) → HP(filtfilt) → LP(filtfilt)
    #   default   → HP(filtfilt) → LP(filtfilt) → Notch(filtfilt) → BR(filtfilt)
    # Señal continua filtrada (sin segmentar) → paso ICA.
    # Salidas:
    #   tables/filtered_signal_{notch|bandreject|highpass|lowpass}.csv
    #     (señal completa tras cada paso aplicado; mismo formato que raw_signal.csv)
    #   figures/filtered_signal_preview.png (raw vs filtrada final)
    _log(log_io, "\n[3/8] Filtrado")
    t_step = now()
    print("  [3/8] Filtrado... ")
    _filt_saved = Tuple{String,String}[]   # (archivo, ruta relativa)
    rec_filt = filter_recording(rec, cfg; on_step = (key, step_rec) -> begin
        fname = _save_filtered_signal(step_rec, export_dir, key, log_io)
        push!(_filt_saved, (fname, "tables/$fname"))
    end)
    _filt_profile = String(get(cfg.filtering, "profile", "default"))
    _filt_ord     = Int(get(cfg.filtering, "filter_order", 4))
    _filt_dur     = round(Dates.value(now() - t_step) / 1000, digits=1)
    println("✓  ($(_filt_dur) s)")

    _chain = describe_filter_chain(cfg)
    _print_kv_table([
        ("Profile", _filt_profile),
        ("Order",   string(_filt_ord)),
        ("Steps",   string(length(_chain))),
    ]; indent="        ")
    _print_cols_table(
        ["#", "Filter", "Frequency", "Order", "Method"],
        [String[string(s.step), s.name, s.freq, string(s.order), s.method] for s in _chain];
        indent="        ",
    )
    _log(log_io, "  Perfil: $(_filt_profile) · order=$(_filt_ord)")
    for step in _chain
        _log(log_io, "    [$(step.step)] $(step.name) $(step.freq)  ord=$(step.order)  método=$(step.method)")
    end

    # Salida de esta fase: preview raw vs filtrada (CairoMakie).
    # Nombre canónico actual: figures/filtered_signal_preview.png
    # (antes Slidev usaba signal_preview_EC.png — deprecado).
    fig_dir = joinpath(export_dir, "figures"); mkpath(fig_dir)
    _filt_fig = "filtered_signal_preview.png"
    _filt_fig_ok = false
    try
        fig = _plot_signal_preview(rec, rec_filt)
        save_figure(fig, joinpath(fig_dir, _filt_fig); dpi)
        _filt_fig_ok = true
        _log(log_io, "  Guardado: $(_filt_fig)")
    catch e
        @warn "No se pudo generar $(_filt_fig): $e"
        _log(log_io, "  WARN: $(_filt_fig) fallido: $e")
    end
    _out_rows = Tuple{String,String}[(_filt_fig, _filt_fig_ok ? "figures/$(_filt_fig)" : "no generado")]
    append!(_out_rows, _filt_saved)
    _print_kv_table(_out_rows; indent="        ", headers=("Archivo", "Ruta"))

    # ── ICA: antes de segmentar (señal continua filtrada) ─────
    #
    # Regla científica crítica (AGENTS.md / .cursorrules):
    #   ICA SIEMPRE sobre la señal CONTINUA filtrada, ANTES de
    #   segmentar. Orden del pipeline: Filtrado [3/8] → ICA [4/8]
    #   → Segmentación [5/8]. No invertir.
    #
    # Flujo de este bloque:
    #   1) Caché  · Si existe ica_result.jls y el hash de config ICA
    #      (+ params de filtrado) coincide → deserializar (evita
    #      recompute ~minutos). El hash incluye filtering porque
    #      ICA se calcula sobre rec_filt; un cambio de Notch/HP
    #      invalidaría componentes antiguos.
    #   2) Cómputo · run_ica (ICACore.jl): PCA whitening + FastICA
    #      simétrico (tanh). Perfil eeg_julia → n_comp = n_channels;
    #      mixing A = pinv(W) (o inv si cuadrada). Sin MultivariateStats.
    #   3) Rechazo · precedencia manual > auto > ninguno:
    #      a) ica_labels.csv MANUAL (label=="artifact") → gana siempre.
    #      b) Sin CSV manual y cfg.ica["auto_reject"]=true → se
    #         rechazan los componentes con artifact_type != "brain"
    #         (evaluate_ica_components, mismo artifact_threshold que
    #         el resumen de terminal) y se persisten en
    #         ica_labels_auto.csv (write_ica_labels_auto).
    #      c) Si no → se conservan todos los componentes.
    #      apply_ica_rejection reconstruye la señal restando los IC
    #      resueltos por (a)/(b).
    #   4) Salida · rec_ica alimenta segmentación; si ICA falla se
    #      continúa con rec_filt (aviso, no aborta el pipeline).
    #   5) Archivos · _save_ica_results escribe tables/ica/*,
    #      json/ica_summary.json, figures/ica/topomaps, etc.
    #   6) Diagnósticos · ICAResult.diagnostics reporta convergencia
    #      real (n_iter, error final, intentos, seed, whitening…).
    #      max_attempts>1 (solo profile=default) reintenta con seed+k.
    #
    _log(log_io, "\n[4/8] ICA")
    t_ica     = now()
    cache_dir = joinpath(export_dir, "cache")
    ica_cache = joinpath(cache_dir, "ica_result.jls")
    ica_cfg_hash = _ica_config_hash(cfg)          # ICA + filtering
    ica_hash_path = joinpath(cache_dir, "ica_config.hash")

    # Cache válido solo si el artefacto y el hash coinciden a la vez
    cache_valid = isfile(ica_cache) &&
                  isfile(ica_hash_path) &&
                  strip(read(ica_hash_path, String)) == ica_cfg_hash

    ica_from_cache = false
    print("  [4/8] ICA... ")
    ica_result = if cache_valid
        # Reutilizar descomposición previa (mismos params)
        ica_from_cache = true
        _log(log_io, "  ICA cargado desde caché")
        try Serialization.deserialize(ica_cache) catch; nothing end
    else
        # FastICA sobre rec_filt (continua). Ante error → nothing
        res = try
            run_ica(rec_filt, cfg)
        catch e
            @warn "ICA falló: $e · continuando sin ICA"
            nothing
        end
        # Persistir para re-runs (solo si el cómputo tuvo éxito)
        if res !== nothing
            mkpath(cache_dir)
            Serialization.serialize(ica_cache, res)
            write(ica_hash_path, ica_cfg_hash)
        end
        res
    end

    # Por defecto la señal hacia segmentación es la filtrada;
    # si hay rechazo manual se sustituye por la reconstrucción limpia.
    rec_ica = rec_filt
    if ica_result !== nothing
        n_comp_ica = size(ica_result.activations, 1)
        _log(log_io, "  Componentes: $(n_comp_ica) | Método: FastICA + PCA whitening")

        # Clasificación automática (score/threshold + tipo) — se calcula una
        # sola vez y alimenta tanto la decisión de rechazo (si aplica) como
        # el resumen de terminal y ica_component_features.csv.
        _artifact_thresh = Float64(get(cfg.ica, "artifact_threshold", 1.5))
        _feat = DataFrame()
        _eval = DataFrame()
        try
            _feat = compute_ica_features(
                Float64.(ica_result.mixing_matrix),
                Float64.(ica_result.activations),
                rec_filt.meta.fs,
                rec_filt.meta.channel_names,
            )
            _eval = evaluate_ica_components(_feat; artifact_thresh=_artifact_thresh)
        catch e
            _log(log_io, "  WARN: clasificación automática ICA falló: $e")
        end

        # Precedencia de rechazo: manual > auto (auto_reject) > ninguno
        _has_manual = has_manual_ica_labels(cfg, sid, sess_id, condition)
        rejection_source, rej_labels = if _has_manual
            ("manual", load_ica_labels(cfg, sid, sess_id, condition))
        elseif Bool(get(cfg.ica, "auto_reject", true)) && !isempty(_eval)
            ("auto", write_ica_labels_auto(export_dir, _eval))
        else
            ("none", Int[])
        end

        ica_result = ICAResult(
            ica_result.meta, ica_result.mixing_matrix, ica_result.unmixing_matrix,
            ica_result.activations, rej_labels, ica_result.variance_explained,
            ica_result.diagnostics,
        )

        if !isempty(rej_labels)
            # Reconstrucción: X_clean = A[:,keep] * S[keep,:]
            _log(log_io, "  Componentes rechazados ($(rejection_source)): $(join(rej_labels, ", "))")
            rec_ica = apply_ica_rejection(rec_filt, ica_result, rej_labels)
            _log(log_io, "  Señal limpiada · $(length(rej_labels)) componente(s) eliminado(s)")
        else
            _log(log_io, "  Sin componentes rechazados (fuente: $(rejection_source)) · Se conservan todos")
        end

        ica_dur = round(Dates.value(now() - t_ica) / 1000, digits=1)
        _log(log_io, "  Duración ICA: $(ica_dur) s")
        _rej_str = isempty(rej_labels) ?
            "ninguno ($(rejection_source))" :
            "$(join(rej_labels, ", ")) ($(rejection_source))"
        _cache_str = ica_from_cache ? "caché" : "nuevo"

        # Diagnósticos de convergencia (vacíos si caché antiguo)
        _d = ica_result.diagnostics
        _conv     = Bool(get(_d, "converged", false))
        _n_iter   = Int(get(_d, "n_iter", 0))
        _max_iter = Int(get(_d, "max_iter", 0))
        _n_att    = Int(get(_d, "n_attempts", 1))
        _max_att  = Int(get(_d, "max_attempts", 1))
        _ferr     = Float64(get(_d, "final_error", NaN))
        _tol      = Float64(get(_d, "tol", NaN))
        _seed_e   = Int(get(_d, "seed", 0))
        _opt_s    = Float64(get(_d, "optimize_s", ica_dur))
        _restart  = Bool(get(_d, "restart_needed", false))
        _has_diag = !isempty(_d) && haskey(_d, "n_iter")

        # Línea de estado informativa (no solo ✓)
        if !_has_diag
            println("✓  desde caché sin diagnósticos  ($(ica_dur) s)")
        elseif _conv && _n_att == 1
            println("✓  convergida · $(_n_iter)/$(_max_iter) iteraciones · 1 intento · error=$(round(_ferr; sigdigits=3)) · $(_opt_s) s")
        elseif _conv && _n_att > 1
            println("⚠  convergida en intento $(_n_att) · $(_n_iter)/$(_max_iter) iteraciones · error=$(round(_ferr; sigdigits=3)) · $(_opt_s) s")
        else
            println("✗  no convergió tras $(_n_att)/$(_max_att) intentos · $(_n_iter)/$(_max_iter) iter · error=$(round(_ferr; sigdigits=3)) · $(_opt_s) s")
        end
        _log(log_io, "  Convergencia: conv=$(_conv) iter=$(_n_iter)/$(_max_iter) attempts=$(_n_att) err=$(_ferr) seed=$(_seed_e)")

        # Parámetros EFECTIVOS (perfil eeg_julia fija max_iter/tol/seed
        # en código; el TOML solo aplica si profile="default").
        _ica_p = _ica_effective_params(cfg, rec_filt.meta.n_channels)
        _print_kv_table([
            ("Method",       "FastICA + PCA whitening"),
            ("Profile",      _ica_p.profile),
            ("n_components", string(_ica_p.n_comp)),
            ("max_iter",     string(_ica_p.max_iter)),
            ("tol",          string(_ica_p.tol)),
            ("seed (cfg)",   string(_ica_p.seed)),
            ("max_attempts", string(_ica_p.max_attempts)),
            ("nonlinearity", "tanh (a=$(_ica_p.a))"),
            ("mixing",       _ica_p.mixing),
            ("Source",       _cache_str),
            ("Rejected ICs", _rej_str),
        ]; indent="        ")

        # Tabla de convergencia real (optimización FastICA)
        if _has_diag
            _estado = !_conv ? "no convergió" :
                      (_restart ? "convergido (reinicio)" : "convergido")
            _print_kv_table([
                ("Estado",                 _estado),
                ("Iteraciones realizadas", "$(_n_iter) / $(_max_iter)"),
                ("Intentos",               "$(_n_att) / $(_max_att)"),
                ("Tolerancia objetivo",    string(_tol)),
                ("Error final",            string(round(_ferr; sigdigits=4))),
                ("Tiempo de optimización", "$(_opt_s) s"),
                ("Seed efectiva",          string(_seed_e)),
                ("Reinicio necesario",     _restart ? "sí" : "no"),
            ]; indent="        ")

            # Detalle por intento (si hubo más de uno, o verbose)
            _atts = get(_d, "attempts", Dict{String,Any}[])
            _ica_verbose = Bool(get(cfg.ica, "verbose", false))
            if length(_atts) > 1 || (_ica_verbose && !isempty(_atts))
                _att_rows = Vector{Vector{String}}()
                for a in _atts
                    push!(_att_rows, [
                        string(get(a, "attempt", "?")),
                        string(get(a, "seed", "?")),
                        string(get(a, "n_iter", "?")),
                        string(round(Float64(get(a, "final_error", NaN)); sigdigits=3)),
                        Bool(get(a, "converged", false)) ? "converge" : "no converge",
                    ])
                end
                _print_cols_table(
                    ["Intento", "Seed", "Iteraciones", "Error final", "Estado"],
                    _att_rows;
                    indent="        ",
                )
            end

            # Calidad numérica PCA / whitening (siempre breve; detalle en verbose)
            _wcond = Float64(get(_d, "whitening_cond", NaN))
            _ortho = Float64(get(_d, "orthogonality", NaN))
            _prank = Int(get(_d, "pca_rank", n_comp_ica))
            _pfull = Int(get(_d, "pca_rank_full", n_comp_ica))
            _pvpct = Float64(get(_d, "pca_var_pct", NaN))
            if _ica_verbose
                _print_kv_table([
                    ("Condición whitening", string(round(_wcond; sigdigits=4))),
                    ("Ortogonalidad final", string(round(_ortho; sigdigits=3))),
                    ("Rango PCA retenido",  "$(_prank)/$(_pfull)"),
                    ("Varianza explicada PCA", string(round(_pvpct, digits=2)) * " %"),
                ]; indent="        ")
            end
        end

        # Resumen de calidad de componentes (clasificación calculada arriba,
        # antes de decidir el rechazo — ver _feat/_eval)
        if !isempty(_eval)
            _atypes = String.(_eval.artifact_type)
            _n_art  = count(!=("brain"), _atypes)
            _n_eye  = count(==("eye/blink"), _atypes)
            _n_mus  = count(==("muscle"), _atypes)
            _n_kurt = hasproperty(_feat, :kurtosis) ?
                      count(k -> k > 5.0, Float64.(_feat.kurtosis)) : 0
            _print_kv_table([
                ("ICs totales",                string(n_comp_ica)),
                ("ICs marcadas artefacto",     string(_n_art)),
                ("ICs rechazadas",             "$(length(rej_labels)) ($(rejection_source))"),
                ("ICs alta curtosis (>5)",     string(_n_kurt)),
                ("ICs tipo eye/blink",         string(_n_eye)),
                ("ICs tipo muscle",            string(_n_mus)),
            ]; indent="        ")
        else
            _log(log_io, "  Resumen calidad IC omitido: clasificación automática no disponible")
        end

        # Tablas/figuras ICA (features, mixing, topomaps, summary…)
        _ica_files = _save_ica_results(ica_result, rec_filt, rec_ica, export_dir, ica_dur, log_io;
            eval_df=_eval, rejection_source=rejection_source, artifact_thresh=_artifact_thresh)
        # Caché serializada (si se escribió o ya existía)
        push!(_ica_files, ("ica_result.jls", "cache/ica_result.jls"))
        push!(_ica_files, ("ica_config.hash", "cache/ica_config.hash"))
        if rejection_source == "auto"
            push!(_ica_files, ("ica_labels_auto.csv", "ica_labels_auto.csv"))
        end
        _print_kv_table(_ica_files; indent="        ", headers=("Archivo", "Ruta"))
    else
        # Fallo de cómputo/caché: no abortar; AR/wPLI siguen sobre filtrada
        ica_dur = round(Dates.value(now() - t_ica) / 1000, digits=1)
        println("✗  omitido por error — continuando con señal filtrada  ($(ica_dur) s)")
        _log(log_io, "  ICA omitido por error · continuando con señal filtrada")
        _ica_p = _ica_effective_params(cfg, rec_filt.meta.n_channels)
        _print_kv_table([
            ("Method",       "FastICA + PCA whitening"),
            ("Profile",      _ica_p.profile),
            ("n_components", string(_ica_p.n_comp)),
            ("max_iter",     string(_ica_p.max_iter)),
            ("tol",          string(_ica_p.tol)),
            ("seed",         string(_ica_p.seed)),
            ("Status",       "omitido (error)"),
        ]; indent="        ")
    end

    # ── Construir montaje de análisis ─────────────────────────
    # Quién elimina canales es `exclude_channels` (pipeline.toml
    # [montage]). `exclude_fp2` NO elimina por sí solo: solo filtra
    # Fp2 de esa lista cuando es false.
    #
    #   exclude_channels=["Fp2"] · exclude_fp2=true  → Fp2 SE ELIMINA
    #   exclude_channels=["Fp2"] · exclude_fp2=false → Fp2 SE CONSERVA
    #   exclude_channels=[]      · exclude_fp2=true  → Fp2 se conserva
    #     (pero fp2_excluded puede reportarse true si Fp2 está en señal)
    #
    # A la lista efectiva se unen siempre bad_ch del QC (z-score).
    # rec_ica queda intacto (export ICA / figuras); rec_for_seg es
    # la señal recortada que entra en segmentación + wPLI.
    exclude_fp2    = Bool(get(montage_cfg, "exclude_fp2", true))  # default código si falta clave
    montage_excl   = String.(get(montage_cfg, "exclude_channels", String[]))
    if !exclude_fp2
        # Conservar Fp2 aunque figure en exclude_channels
        montage_excl = filter(c -> c != "Fp2", montage_excl)
    end
    all_excl       = unique(vcat(montage_excl, bad_ch))
    fp2_in_signal  = "Fp2" ∈ rec_ica.meta.channel_names
    fp2_excluded   = exclude_fp2 && fp2_in_signal

    rec_for_seg = if !isempty(all_excl)
        excl_idx   = findall(ch -> ch ∈ Set(all_excl), rec_ica.meta.channel_names)
        good_idx   = setdiff(1:rec_ica.meta.n_channels, excl_idx)
        good_names = rec_ica.meta.channel_names[good_idx]
        new_meta   = RecordingMeta(
            rec_ica.meta.subject_id, rec_ica.meta.session_id,
            rec_ica.meta.condition, rec_ica.meta.run,
            rec_ica.meta.fs, length(good_idx), good_names,
            rec_ica.meta.channel_positions, rec_ica.meta.bids_path
        )
        montage_note = fp2_excluded ? "Fp2 (exclude_fp2=true)" : ""
        qc_note      = isempty(bad_ch) ? "" : join(bad_ch, ", ") * " (QC)"
        excl_note    = filter(!isempty, [montage_note, qc_note])
        _log(log_io, "  Canales excluidos del análisis: $(join(excl_note, "; ")) → $(length(good_idx)) canales activos")
        EEGRecording(new_meta, rec_ica.data[good_idx, :], rec_ica.times)
    else
        fp2_excluded = false
        rec_ica
    end
    n_ch_analysis  = rec_for_seg.meta.n_channels
    bad_ch_non_fp2 = filter(c -> c != "Fp2", bad_ch)   # bad QC distintos de Fp2

    # Resumen terminal: montaje efectivo hacia segmentación / wPLI
    _fp2_status = if fp2_excluded
        "excluido (exclude_fp2=true)"
    elseif fp2_in_signal
        "incluido"
    else
        "no está en la señal"
    end
    _excl_eff = isempty(all_excl) ? "ninguno" : join(all_excl, ", ")
    _print_kv_table([
        ("Channels (analysis)", string(n_ch_analysis)),
        ("Fp2",                 _fp2_status),
        ("exclude_fp2",         string(exclude_fp2)),
        ("Excluded (montage+QC)", _excl_eff),
        ("Bad QC (non-Fp2)",    isempty(bad_ch_non_fp2) ? "ninguno" : join(bad_ch_non_fp2, ", ")),
    ]; indent="        ")
    _log(log_io, "  Montaje análisis: $(n_ch_analysis) ch · Fp2=$(_fp2_status) · excluidos=$(_excl_eff)")

    # ── Segmentación + baseline + AR ──────────────────────────
    #
    # Orden (Epochs.jl / pipeline.toml [segmentation][baseline][artifact_rejection]):
    #   1) segment_recording  · corta rec_for_seg en epochs.
    #      profile="eeg_julia" → fuerza 1.0 s, overlap 0 (ignora
    #      epoch_length_s / epoch_overlap del config).
    #   2) apply_baseline     · 1ª pasada (siempre si apply=true).
    #      method="first_window_mean" → resta media de [0, baseline_end_s]
    #      por canal/epoch (EEG_Julia 0–100 ms). baseline_start_s es inerte.
    #   3) reject_artifacts   · AR ±70 µV (eeg_julia) o umbral+gradiente
    #      (default). Elimina epochs malos → epochs.n_valid.
    #   4) apply_baseline     · 2ª pasada solo si n_passes ≥ 2
    #      (recalcula baseline sobre epochs ya limpios).
    #
    # min_epochs (min_segments en TOML): umbral blando aquí (@warn);
    # en qc_decision_table es exclusión dura si n_valid < min_epochs.
    #
    _log(log_io, "\n[5/8] Segmentación")
    t_seg       = now()
    print("  [5/8] Segmentación + AR... ")
    n_passes    = Int(get(cfg.baseline, "n_passes", 1))
    seg_profile = String(get(cfg.segmentation, "profile", "default"))
    ar_profile  = String(get(cfg.artifact_rejection, "profile", "default"))

    # Epoch length / overlap EFECTIVOS (perfil eeg_julia los fija)
    _ep_cfg_s   = Float64(get(cfg.segmentation, "epoch_length_s", 1.0))
    _ov_cfg     = Float64(get(cfg.segmentation, "epoch_overlap", 0.0))
    _ep_eff_s   = seg_profile == "eeg_julia" ? 1.0 : _ep_cfg_s
    _ov_eff     = seg_profile == "eeg_julia" ? 0.0 : _ov_cfg
    _min_ep     = Int(get(cfg.segmentation, "min_epochs", 10))
    _bl_apply   = Bool(get(cfg.baseline, "apply", true))
    _bl_method  = String(get(cfg.baseline, "method", "mean"))
    _bl_end_s   = Float64(get(cfg.baseline, "baseline_end_s", 0.10))
    _bl_start_s = Float64(get(cfg.baseline, "baseline_start_s", 0.0))

    epochs_raw  = segment_recording(rec_for_seg, cfg)
    epochs_bl1  = apply_baseline(epochs_raw, cfg)          # pasada 1 (pre-AR)
    epochs_ar   = reject_artifacts(epochs_bl1, cfg)
    epochs      = n_passes >= 2 ? apply_baseline(epochs_ar, cfg) : epochs_ar  # pasada 2
    n_total    = n_epochs(epochs_raw)
    n_valid    = epochs.n_valid
    n_rejected = length(epochs_ar.rejected_idx)
    seg_dur    = round(Dates.value(now() - t_seg) / 1000, digits=1)
    _valid_pct = n_total > 0 ? round(100 * n_valid / n_total, digits=1) : 0.0
    _ar_thr    = Float64(get(cfg.artifact_rejection, "max_amplitude_uv", 70.0))
    _epoch_icon = n_valid == 0 ? "✗" : (n_valid < _min_ep ? "⚠" : "✓")
    println("$(_epoch_icon)  ($(seg_dur) s)")
    if n_valid == 0
        println("        ✗  Sin epochs válidos — grabación no utilizable con umbral ±$(_ar_thr) µV")
    end

    # Config tomada (efectiva) — segmentación + baseline
    _print_kv_table([
        ("Seg. profile",     seg_profile),
        ("Epoch length",     "$(_ep_eff_s) s" * (seg_profile == "eeg_julia" && _ep_cfg_s != 1.0 ?
                              " (forzado; TOML=$(_ep_cfg_s)s)" : "")),
        ("Overlap",          "$(Int(round(_ov_eff * 100)))%" *
                              (seg_profile == "eeg_julia" ? " (forzado)" : "")),
        ("min_epochs",       string(_min_ep)),
        ("Baseline apply",   string(_bl_apply)),
        ("Baseline method",  _bl_method),
        ("Baseline window",  "[$(_bl_start_s), $(_bl_end_s)] s" *
                              (_bl_method == "first_window_mean" ?
                               " (start inerte)" : "")),
        ("Baseline passes",  string(n_passes) *
                              (n_passes >= 2 ? " (pre-AR + post-AR)" : " (solo pre-AR)")),
        ("AR profile",       ar_profile),
        ("AR threshold",     "±$(_ar_thr) µV"),
    ]; indent="        ")

    # Resultado de la fase
    _sig_in = (ica_result !== nothing && !isempty(ica_result.rejected_components)) ?
              "ICA-limpiada" : "filtrada"
    _print_kv_table([
        ("Epochs total",    string(n_total)),
        ("Epochs valid",    "$(n_valid) ($(_valid_pct)%)"),
        ("Epochs rejected", string(n_rejected)),
        ("Channels",        string(n_ch_analysis)),
        ("Fp2",             fp2_excluded ? "excluido" : "incluido"),
        ("Signal input",    _sig_in),
    ]; indent="        ")

    _log(log_io, "  Perfil segmentación: $(seg_profile) | Perfil AR: $(ar_profile) | Baseline passes: $(n_passes)")
    _log(log_io, "  Epoch efectivo: $(_ep_eff_s) s · overlap=$(_ov_eff) · min_epochs=$(_min_ep)")
    _log(log_io, "  Baseline: apply=$(_bl_apply) method=$(_bl_method) window=[0,$(_bl_end_s)] s")
    _log(log_io, "  Segmentos totales: $(n_total) | Válidos: $(n_valid) | Rechazados: $(n_rejected) ($(_valid_pct)%)")
    _log(log_io, "  Duración segmentación: $(seg_dur) s")

    seg_signal_label = _sig_in
    _seg_files = _save_segmentation_results(epochs_bl1, epochs, rec_ica, cfg, export_dir, t_seg, log_io, seg_signal_label)
    _print_kv_table(_seg_files; indent="        ", headers=("Archivo", "Ruta"))

    # ── 8. Análisis espectral ─────────────────────────────────
    # compute_psd() — FFT Hamming-taper estilo BrainVision Analyzer
    # (PowerSpectrum.jl). PSD media sobre epochs válidos post-AR.
    #
    # nfft efectivo = max(cfg.spectral.nfft, n_samp). Con epoch 1.0 s
    # y fs=500 Hz (n_samp=500), nfft=1024 se mantiene vía zero-padding
    # → resolución ≈ fs/nfft ≈ 0.488 Hz/bin · n_bins = nfft/2+1.
    #
    # window_pct: % de la ventana con taper Hamming (10 → 5% cada lado).
    #
    # Bandas: intervalo semiabierto [flo, fhi) por banda (independientes;
    # no son un reparto mutuamente excluyente de bins). Con los límites
    # actuales hay SOLAPE THETA/ALPHA (7.8–8.0 Hz) y HUECO ALPHA/BETA_LOW
    # (11.7–12.0 Hz) — la suma de bandas ≠ potencia total; no es un bug.
    #
    _log(log_io, "\n[6/8] Espectral")
    t_step = now()
    print("  [6/8] PSD... ")
    spectra = compute_psd(epochs, cfg)
    _nfft_eff  = Int(spectra.params["nfft"])
    _n_bins    = length(spectra.freqs)
    _df_hz     = _n_bins > 1 ? round(spectra.freqs[2] - spectra.freqs[1], digits=4) : NaN
    _win_pct   = Float64(get(spectra.params, "window_pct",
                    get(cfg.spectral, "window_pct", 10.0)))
    _nfft_cfg  = Int(get(cfg.spectral, "nfft", 1024))
    _psd_dur   = round(Dates.value(now() - t_step) / 1000, digits=1)
    _log(log_io, "  nfft=$(_nfft_eff) (cfg=$(_nfft_cfg)) | Bins: $(_n_bins) | Δf: $(_df_hz) Hz | window_pct=$(_win_pct)")
    for (band, (flo, fhi)) in sort(collect(cfg.bands))
        bp_mean = mean(spectra.band_power[band])
        _log(log_io, "  Potencia $(lpad(band,9)) [$flo, $fhi): $(round(bp_mean, digits=4)) μV² (media canales)")
    end
    println("✓  ($(_psd_dur) s)")

    # Parámetros adoptados (efectivos tras max(nfft, n_samp))
    _print_kv_table([
        ("Método",       "FFT Hamming-taper (BrainVision)"),
        ("nfft (cfg)",   string(_nfft_cfg)),
        ("nfft (efect.)", string(_nfft_eff)),
        ("n_bins",       string(_n_bins)),
        ("Δf (Hz)",      string(_df_hz)),
        ("window_pct",   string(_win_pct)),
        ("fs (Hz)",      string(spectra.meta.fs)),
        ("n_epochs",     string(spectra.n_epochs_used)),
        ("n_channels",   string(length(spectra.meta.channel_names))),
    ]; indent="        ")

    # Potencia media por banda + intervalo semiabierto
    _band_rows = Vector{Vector{String}}()
    _bp_means  = Float64[]
    for (band, (flo, fhi)) in sort(collect(cfg.bands))
        _bp_valid = filter(!isnan, spectra.band_power[band])
        bp_m = isempty(_bp_valid) ? NaN : mean(_bp_valid)
        push!(_bp_means, bp_m)
        push!(_band_rows, [
            band,
            "[$flo, $fhi)",
            isnan(bp_m) ? "NaN" : string(round(bp_m, digits=4)),
        ])
    end
    _bp_sum = max(sum(filter(!isnan, _bp_means); init=0.0), 1e-12)
    for (i, row) in enumerate(_band_rows)
        pct = isnan(_bp_means[i]) ? "—" : string(round(100 * _bp_means[i] / _bp_sum, digits=1)) * "%"
        push!(row, pct)
    end
    _print_cols_table(
        ["Banda", "Intervalo", "Media µV²", "% rel."],
        _band_rows;
        indent="        ",
    )

    # CSV / JSON / PNG espectrales (antes del bloque final [8/8])
    _spec_files = _save_spectral_phase_results(spectra, cfg, export_dir, dpi, log_io)
    _print_kv_table(_spec_files; indent="        ", headers=("Archivo", "Ruta"))

    # ── 9. Conectividad wPLI ──────────────────────────────────
    # compute_wpli() — método vía cfg.connectivity["wpli_method"]
    # (_build_estimator, wPLI.jl). Surrogates (paso siguiente) usan
    # el mismo estimador que la matriz observada.
    #
    # Métodos: hilbert (Butterworth+Hilbert; default / M05) ·
    # fourier_csd · multitaper. use_dwpli=false → wPLI clásico [0,1];
    # true → dwPLI no sesgado [−1,1] (mejor para grupos, diverge de M05).
    #
    # use_csd: CSD opcional antes del wPLI (false en todos los
    # resultados existentes). min_cycles_for_wpli: con epochs 1.0 s
    # DELTA (0.5 Hz) tiene 0.5 ciclos → @warn; solo se omite si
    # exclude_unreliable_bands=true.
    #
    _log(log_io, "\n[7/8] Conectividad wPLI")
    t_step = now()
    print("  [7/8] wPLI... ")
    use_csd = Bool(get(cfg.connectivity, "use_csd", false))
    epochs_conn = use_csd ? apply_csd(epochs, cfg) : epochs
    conn = compute_wpli(epochs_conn, cfg)
    _n_edges   = n_ch_analysis * (n_ch_analysis - 1) ÷ 2
    _wpli_meth = String(get(conn.params, "wpli_method",
                    get(cfg.connectivity, "wpli_method", "hilbert")))
    _use_dwpli = Bool(get(conn.params, "use_dwpli",
                    get(cfg.connectivity, "use_dwpli", false)))
    _min_cyc   = Float64(get(cfg.connectivity, "min_cycles_for_wpli", 4.0))
    _excl_unr  = Bool(get(cfg.connectivity, "exclude_unreliable_bands", false))
    _filt_ord  = Int(get(cfg.connectivity, "filter_order", 8))
    _epoch_s   = size(epochs_conn.data, 2) / Float64(epochs_conn.meta.fs)
    _skipped   = String.(get(conn.params, "skipped_bands", String[]))
    _bands_ok  = sort(collect(keys(conn.matrices)))
    _wpli_dur  = round(Dates.value(now() - t_step) / 1000, digits=1)
    _log(log_io, "  Espacio: $(conn.space) | Método: $(_wpli_meth) | Estimador: $(conn.method)")
    _log(log_io, "  Bandas: " * join(_bands_ok, ", ") *
         (isempty(_skipped) ? "" : " | omitidas: " * join(_skipped, ", ")))
    println("✓  ($(_wpli_dur) s)")

    # Parámetros adoptados
    _conn_rows = Tuple{String,String}[
        ("Método",          _wpli_meth),
        ("Estimador",       _use_dwpli ? "dwPLI" : "wPLI"),
        ("use_dwpli",       string(_use_dwpli)),
        ("Espacio",         conn.space),
        ("use_csd",         string(use_csd)),
        ("n_channels",      string(n_ch_analysis)),
        ("n_edges/banda",   string(_n_edges)),
        ("n_epochs",        string(conn.n_epochs_used)),
        ("epoch_s",         string(round(_epoch_s, digits=3))),
        ("min_cycles",      string(_min_cyc)),
        ("excl. unreliable", string(_excl_unr)),
        ("Bandas OK",       join(_bands_ok, ", ")),
        ("Bandas omitidas", isempty(_skipped) ? "ninguna" : join(_skipped, ", ")),
    ]
    if _wpli_meth == "hilbert"
        push!(_conn_rows, ("filter_order", string(get(conn.params, "filter_order", _filt_ord))))
    elseif _wpli_meth == "fourier_csd"
        _fc = get(cfg.connectivity, "fourier_csd", Dict{String,Any}())
        push!(_conn_rows, ("window", string(get(conn.params, "window", get(_fc, "window", "hann")))))
        push!(_conn_rows, ("nfft",   string(get(conn.params, "nfft", get(_fc, "nfft", 0)))))
    elseif _wpli_meth == "multitaper"
        _mt = get(cfg.connectivity, "multitaper", Dict{String,Any}())
        push!(_conn_rows, ("nw",       string(get(conn.params, "nw", get(_mt, "nw", 4.0)))))
        push!(_conn_rows, ("n_tapers", string(get(conn.params, "n_tapers", get(_mt, "n_tapers", 0)))))
        push!(_conn_rows, ("low_bias", string(get(conn.params, "low_bias", get(_mt, "low_bias", true)))))
    end
    _print_kv_table(_conn_rows; indent="        ")

    # Media wPLI (triángulo superior) + ciclos/época por banda
    _wpli_band_rows = Vector{Vector{String}}()
    for band in sort(collect(keys(cfg.bands)))
        flo, fhi = cfg.bands[band]
        n_cyc = round(flo * _epoch_s, digits=2)
        if haskey(conn.matrices, band)
            W = conn.matrices[band]
            n = size(W, 1)
            vals = Float64[W[i,j] for i in 1:n for j in (i+1):n]
            μ = isempty(vals) ? NaN : mean(vals)
            push!(_wpli_band_rows, [
                band,
                "[$flo, $fhi)",
                string(n_cyc),
                isnan(μ) ? "NaN" : string(round(μ, digits=4)),
                "OK",
            ])
        else
            push!(_wpli_band_rows, [
                band,
                "[$flo, $fhi)",
                string(n_cyc),
                "—",
                "omitida",
            ])
        end
    end
    _print_cols_table(
        ["Banda", "Intervalo", "ciclos/ép", "μ wPLI", "Estado"],
        _wpli_band_rows;
        indent="        ",
    )

    # CSV / JSON / PNG de conectividad (antes del bloque final [8/8])
    _conn_files = _save_connectivity_phase_results(conn, cfg, export_dir, dpi, log_io)
    _print_kv_table(_conn_files; indent="        ", headers=("Archivo", "Ruta"))

    # ── Surrogates (opcional) ─────────────────────────────────
    # surrogate_test() (Surrogates.jl) — desplazamiento circular
    # independiente por canal/época, con el MISMO estimador wPLI
    # que la matriz observada (hilbert / fourier_csd / multitaper).
    #
    # method en TOML es metadato: el código siempre usa
    # circular_shift (phase_shuffle no está implementado).
    # p-valor Monte Carlo (+1): p ≥ 1/(N+1). Con N=200 → suelo ≈0.005;
    # con N≈20 el mínimo ≈0.048 queda al borde de α=0.05.
    # FDR: bh (Benjamini–Hochberg) o bonferroni; seed base + idx banda.
    #
    # enabled=false es el default (197/205 grabaciones del lote).
    #
    surr_results = SurrogateResult[]
    _surr_on   = Bool(get(cfg.surrogates, "enabled", false))
    _n_sur     = Int(get(cfg.surrogates, "n_surrogates", 200))
    _surr_meth = String(get(cfg.surrogates, "method", "circular_shift"))
    _surr_alpha = Float64(get(cfg.surrogates, "alpha", 0.05))
    _surr_fdr  = String(get(cfg.surrogates, "fdr_method", "bh"))
    _surr_seed = Int(get(cfg.surrogates, "seed", 42))
    _p_floor   = round(1.0 / (_n_sur + 1), digits=5)
    _wpli_est  = String(get(conn.params, "wpli_method",
                    get(cfg.connectivity, "wpli_method", "hilbert")))

    _log(log_io, "\n[SUR] Inferencia por surrogates")
    if !_surr_on
        print("  [SUR] Surrogates... ")
        println("⊘  omitido (enabled=false)")
        _print_kv_table([
            ("enabled",       "false"),
            ("n_surrogates",  string(_n_sur)),
            ("method (meta)", _surr_meth),
            ("método real",   "circular_shift"),
            ("estimador wPLI", _wpli_est),
            ("alpha",         string(_surr_alpha)),
            ("fdr_method",    _surr_fdr),
            ("seed",          string(_surr_seed)),
            ("p_min teórico", string(_p_floor)),
        ]; indent="        ")
        _log(log_io, "  omitido (enabled=false) · N=$(_n_sur) · FDR=$(_surr_fdr) · α=$(_surr_alpha)")
    else
        print("  [SUR] Surrogates... ")
        t_sur = now()
        _log(log_io, "  Método real: circular_shift (meta=$(_surr_meth)) | N=$(_n_sur) | FDR=$(_surr_fdr) | α=$(_surr_alpha) | seed=$(_surr_seed)")
        _log(log_io, "  Estimador wPLI: $(_wpli_est) | p_min ≈ $(_p_floor)")

        _band_sig_rows = Vector{Vector{String}}()
        println()  # progreso por banda debajo de la cabecera
        for band in sort(collect(keys(conn.matrices)))
            print("        $(lpad(band,9))... ")
            try
                sr = surrogate_test(epochs_conn, conn, band, cfg)
                push!(surr_results, sr)
                n_sig = count(sr.sig_mask) ÷ 2
                println("$(n_sig)/$(_n_edges) pares sig  (FDR thr=$(round(sr.fdr_threshold, digits=4)))")
                _log(log_io, "  Banda $(lpad(band,9)): $(n_sig)/$(_n_edges) pares sig  FDR-thr=$(round(sr.fdr_threshold, digits=4))")
                push!(_band_sig_rows, [
                    band,
                    string(sr.n_surrogates),
                    string(n_sig),
                    string(_n_edges),
                    string(round(100.0 * n_sig / max(1, _n_edges), digits=1)) * "%",
                    string(round(sr.fdr_threshold, digits=4)),
                ])
            catch e
                println("⚠ fallido: $e")
                @warn "Surrogate fallido para $band: $e"
                _log(log_io, "  WARN: surrogate $(band) fallido: $e")
                push!(_band_sig_rows, [band, "—", "—", string(_n_edges), "—", "error"])
            end
        end
        sur_dur = round(Dates.value(now() - t_sur) / 1000, digits=1)
        println("        ✓  total ($(sur_dur) s)")

        _print_kv_table([
            ("enabled",        "true"),
            ("n_surrogates",   string(_n_sur)),
            ("method (meta)",  _surr_meth),
            ("método real",    "circular_shift"),
            ("estimador wPLI", _wpli_est),
            ("alpha",          string(_surr_alpha)),
            ("fdr_method",     _surr_fdr),
            ("seed",           string(_surr_seed)),
            ("p_min teórico",  string(_p_floor)),
            ("n_edges/banda",  string(_n_edges)),
            ("n_bandas",       string(length(conn.matrices))),
        ]; indent="        ")

        if !isempty(_band_sig_rows)
            _print_cols_table(
                ["Banda", "N surr", "n_sig", "n_pares", "% sig", "FDR thr"],
                _band_sig_rows;
                indent="        ",
            )
        end

        if !isempty(surr_results)
            _surr_files = _save_surrogate_results(surr_results, conn, export_dir, cfg, log_io)
            _print_kv_table(_surr_files; indent="        ", headers=("Archivo", "Ruta"))
        end
    end

    # ── 10. Guardar resultados ────────────────────────────────
    print("  [8/8] Guardando resultados... ")
    _log(log_io, "\n[8/8] Guardando resultados")
    t_step = now()
    _save_all_results(rec, rec_filt, qc_stats, bad_ch, spectra, conn, val,
                      export_dir, condition, cfg, dpi, log_io)
    _save_config_snapshot(config_path, export_dir)
    _update_subjects_index(results_dir(cfg), sid, sess_id, task, rec, n_valid, n_rejected)

    # Tabla QC global (results/qc/qc_decision_table.csv)
    _update_qc_decision_table(
        results_dir(cfg), sid, sess_id, condition,
        amplitude_warning, raw_sigma_mean,
        bad_ch, bad_ch_non_fp2, fp2_excluded, n_ch_analysis,
        Float64(get(cfg.artifact_rejection, "max_amplitude_uv", 70.0)),
        n_total, n_rejected, n_valid,
        Int(get(cfg.segmentation, "min_epochs", 10)),
        log_io
    )
    println("✓  tablas · figuras · BIDS export · QC table  ($(round(Dates.value(now()-t_step)/1000,digits=1)) s)")

    elapsed = round(Dates.value(now() - t0) / 1000, digits=1)
    _log(log_io, "\n✓ Pipeline completado en $(elapsed) s")
    close(log_io)

    # ── Resumen final ──────────────────────────────────────────
    _valid_pct_f  = n_total > 0 ? round(100*n_valid/n_total, digits=1) : 0.0
    _qc_decision  = n_valid == 0 ? "exclude" :
                    (amplitude_warning || !isempty(bad_ch_non_fp2)) ?
                    (n_valid < 10 ? "manual_review" : "include_with_warning") : "include"
    _decision_icon = _qc_decision == "include" ? "✅" :
                     _qc_decision == "include_with_warning" ? "🟡" :
                     _qc_decision == "manual_review" ? "🟠" : "❌"
    _bad_summary  = isempty(bad_ch_non_fp2) ? "ninguno" : join(bad_ch_non_fp2, ", ")
    _warn_flag    = amplitude_warning ? " · ⚠ amplitude_warning" : ""
    println()
    println("  " * "─"^60)
    println("  $(_decision_icon)  sub-$(sid) / ses-$(sess_id) / $(condition)  →  $(_qc_decision)   ($(elapsed) s total)")
    println("     Epochs  : $(n_valid)/$(n_total) válidos ($(_valid_pct_f)%) · rechazados: $(n_rejected)")
    println("     Montaje : $(n_ch_analysis) canales$(fp2_excluded ? " (Fp2 excluido)" : "") · bad no-Fp2: $(_bad_summary)$(_warn_flag)")
    println("     PSD     : " * join(["$(b)=$(round(mean(spectra.band_power[b]),digits=2))" for (b,_) in sort(collect(cfg.bands))], " · ") * " µV²")
    if !isempty(surr_results)
        _sig_counts = ["$(sr.band): $(count(sr.sig_mask)÷2)" for sr in surr_results]
        println("     Sig wPLI: " * join(_sig_counts, " · ") * " pares (FDR 5%)")
    end
    println("  " * "─"^60)
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
)::Vector{Tuple{String,String}}
    written   = Tuple{String,String}[]
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
    tables_dir = joinpath(export_dir, "tables"); mkpath(tables_dir)
    json_dir   = joinpath(export_dir, "json");   mkpath(json_dir)

    # 1) segmentation_summary.json
    open(joinpath(json_dir, "segmentation_summary.json"), "w") do f
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
    push!(written, ("segmentation_summary.json", "json/segmentation_summary.json"))
    _log(log_io, "  Guardado: segmentation_summary.json")

    # 2) segments_table.csv
    if !isempty(qr)
        CSV.write(joinpath(tables_dir, "segments_table.csv"), qr)
        push!(written, ("segments_table.csv", "tables/segments_table.csv"))
        _log(log_io, "  Guardado: segments_table.csv ($(n_total) épocas)")
    end

    # 3) channel_coverage.csv
    if !isempty(cov_df)
        CSV.write(joinpath(tables_dir, "channel_coverage.csv"), cov_df)
        push!(written, ("channel_coverage.csv", "tables/channel_coverage.csv"))
        _log(log_io, "  Guardado: channel_coverage.csv ($(rec.meta.n_channels) canales)")
    end

    # 4–6) Ficheros específicos de rechazo de artefactos
    if !isempty(qr)
        append!(written, _save_ar_results(qr, cfg, export_dir, ts, dur_s,
                         n_total, n_valid, n_rejected, log_io,
                         rec.meta.n_channels))
    end
    return written
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
)::Vector{Tuple{String,String}}
    written     = Tuple{String,String}[]
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

    tables_dir = joinpath(export_dir, "tables"); mkpath(tables_dir)
    json_dir   = joinpath(export_dir, "json");   mkpath(json_dir)

    # ── 4) artifact_rejection_summary.json ───────────────────────────────────
    open(joinpath(json_dir, "artifact_rejection_summary.json"), "w") do f
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
    push!(written, ("artifact_rejection_summary.json", "json/artifact_rejection_summary.json"))
    _log(log_io, "  Guardado: artifact_rejection_summary.json")

    # ── 5) rejected_segments.csv ──────────────────────────────────────────────
    # Seleccionar columnas canónicas (orden limpio para CSV)
    if !isempty(rej_df)
        wanted = [:epoch, :start_s, :end_s, :duration_s, :quality,
                  :status, :rejection_reason,
                  :max_amp_uv, :min_amp_uv, :p2p_uv,
                  :worst_channel, :channels_violating, :max_grad_uv]
        present = [c for c in wanted if hasproperty(rej_df, c)]
        CSV.write(joinpath(tables_dir, "rejected_segments.csv"), rej_df[:, present])
        push!(written, ("rejected_segments.csv", "tables/rejected_segments.csv"))
        _log(log_io, "  Guardado: rejected_segments.csv ($(nrow(rej_df)) rechazados)")
    end

    # ── 6) channel_artifact_summary.csv ──────────────────────────────────────
    if !isempty(ca_df)
        CSV.write(joinpath(tables_dir, "channel_artifact_summary.csv"), ca_df)
        push!(written, ("channel_artifact_summary.csv", "tables/channel_artifact_summary.csv"))
        _log(log_io, "  Guardado: channel_artifact_summary.csv")
    end
    return written
end

# ─── Guardado fase espectral (CSV / JSON / PNG) ───────────────

"""
    _save_spectral_phase_results(spectra, cfg, export_dir, dpi, log_io)
        -> Vector{Tuple{String,String}}

Escribe en el paso [6/8] todos los artefactos espectrales:
  tables/psd_by_channel.csv, band_power_summary.csv,
  regional_psd.csv, spectral_indices.csv,
  json/spectral_summary.json,
  figures/psd_all_channels.png, band_power_summary.png
"""
function _save_spectral_phase_results(
    spectra::SpectralResult,
    cfg::PipelineConfig,
    export_dir::String,
    dpi::Int,
    log_io::IO,
)::Vector{Tuple{String,String}}
    written    = Tuple{String,String}[]
    tables_dir = joinpath(export_dir, "tables"); mkpath(tables_dir)
    fig_dir    = joinpath(export_dir, "figures"); mkpath(fig_dir)

    ch_names = spectra.meta.channel_names
    n_ch_psd, n_freqs_psd = size(spectra.psd)
    n_ch     = min(length(ch_names), n_ch_psd)
    n_freqs  = min(length(spectra.freqs), n_freqs_psd)
    if n_ch != length(ch_names) || n_freqs != length(spectra.freqs)
        _log(log_io, "  WARN: PSD dims ajustadas al guardar: psd=$(size(spectra.psd)), canales=$(length(ch_names)), freqs=$(length(spectra.freqs))")
    end

    psd_rows = [(channel=ch_names[c], freq_hz=round(spectra.freqs[f], digits=3),
                 power_uv2=round(spectra.psd[c,f], digits=6))
                for c in 1:n_ch for f in 1:n_freqs]
    CSV.write(joinpath(tables_dir, "psd_by_channel.csv"), DataFrame(psd_rows))
    push!(written, ("psd_by_channel.csv", "tables/psd_by_channel.csv"))
    _log(log_io, "  Guardado: psd_by_channel.csv")

    band_names = sort(collect(keys(spectra.band_power)))
    bp_df = DataFrame(channel = ch_names[1:n_ch])
    for b in band_names
        bp_df[!, b] = round.(spectra.band_power[b][1:n_ch], digits=6)
    end
    CSV.write(joinpath(tables_dir, "band_power_summary.csv"), bp_df)
    push!(written, ("band_power_summary.csv", "tables/band_power_summary.csv"))
    _log(log_io, "  Guardado: band_power_summary.csv")

    append!(written, _save_spectral_extras(spectra, cfg, export_dir, log_io))

    try
        fig = plot_spectrum_grid(spectra; xmax=50.0, cols=6)
        save_figure(fig, joinpath(fig_dir, "psd_all_channels.png"); dpi)
        push!(written, ("psd_all_channels.png", "figures/psd_all_channels.png"))
        _log(log_io, "  Guardado: psd_all_channels.png")
    catch e
        @warn "No se pudo generar psd_grid: $e"
        _log(log_io, "  WARN: psd_grid fallido: $e")
    end

    try
        fig = _plot_band_power_summary(spectra)
        save_figure(fig, joinpath(fig_dir, "band_power_summary.png"); dpi)
        push!(written, ("band_power_summary.png", "figures/band_power_summary.png"))
        _log(log_io, "  Guardado: band_power_summary.png")
    catch e
        @warn "No se pudo generar band_power_summary: $e"
        _log(log_io, "  WARN: band_power_summary fallido: $e")
    end

    return written
end

# ─── Extras espectrales: summary.json, regional_psd.csv, spectral_indices.csv ───

function _save_spectral_extras(
    spectra::SpectralResult,
    cfg::PipelineConfig,
    export_dir::String,
    log_io::IO,
)::Vector{Tuple{String,String}}
    written   = Tuple{String,String}[]
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
    tables_dir = joinpath(export_dir, "tables"); mkpath(tables_dir)
    json_dir   = joinpath(export_dir, "json");   mkpath(json_dir)

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

    open(joinpath(json_dir, "spectral_summary.json"), "w") do f
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
    push!(written, ("spectral_summary.json", "json/spectral_summary.json"))
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
        CSV.write(joinpath(tables_dir, "regional_psd.csv"),
                  DataFrame(region=reg_region, band=reg_band,
                            mean_power=reg_mean, std_power=reg_std, n_channels=reg_n))
        push!(written, ("regional_psd.csv", "tables/regional_psd.csv"))
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

    CSV.write(joinpath(tables_dir, "spectral_indices.csv"),
              DataFrame(channel=idx_ch,
                        alpha_theta=idx_at, beta_alpha=idx_ba,
                        theta_beta=idx_tb, gamma_alpha=idx_ga,
                        peak_alpha_hz=idx_ph, peak_alpha_uv2=idx_pu))
    push!(written, ("spectral_indices.csv", "tables/spectral_indices.csv"))
    _log(log_io, "  Guardado: spectral_indices.csv")
    return written
end

# ─── Guardado fase conectividad (CSV / JSON / PNG) ────────────

"""
    _save_connectivity_phase_results(conn, cfg, export_dir, dpi, log_io)
        -> Vector{Tuple{String,String}}

Escribe en el paso [7/8] los artefactos de conectividad observada:
  tables/connectivity/wpli_{band}.csv, connectivity_edges.csv,
  network_metrics.csv,
  json/connectivity_summary.json,
  figures/connectivity/wpli_{band}.png
"""
function _save_connectivity_phase_results(
    conn::ConnectivityMatrix,
    cfg::PipelineConfig,
    export_dir::String,
    dpi::Int,
    log_io::IO,
)::Vector{Tuple{String,String}}
    written         = Tuple{String,String}[]
    fig_conn_dir    = joinpath(export_dir, "figures", "connectivity"); mkpath(fig_conn_dir)
    tables_conn_dir = joinpath(export_dir, "tables", "connectivity");  mkpath(tables_conn_dir)
    ch_names        = conn.channel_names

    all_edges = DataFrame(ch_a=String[], ch_b=String[], band=String[],
                          wpli=Float64[], rank=Int[])
    for band in sort(collect(keys(conn.matrices)))
        W = conn.matrices[band]
        mat_df = DataFrame(hcat(ch_names, W), vcat(["channel"], ch_names))
        CSV.write(joinpath(tables_conn_dir, "wpli_$(band).csv"), mat_df)
        push!(written, ("wpli_$(band).csv", "tables/connectivity/wpli_$(band).csv"))

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
        append!(all_edges, edge_df)
    end
    CSV.write(joinpath(tables_conn_dir, "connectivity_edges.csv"), all_edges)
    push!(written, ("connectivity_edges.csv", "tables/connectivity/connectivity_edges.csv"))
    _log(log_io, "  Guardado: matrices y edges wPLI por banda")

    append!(written, _save_connectivity_extras(conn, export_dir, log_io, cfg))

    for band in sort(collect(keys(conn.matrices)))
        try
            fig = plot_connectivity_heatmap(conn, band)
            save_figure(fig, joinpath(fig_conn_dir, "wpli_$(band).png"); dpi)
            push!(written, ("wpli_$(band).png", "figures/connectivity/wpli_$(band).png"))
        catch e
            @warn "No se pudo generar heatmap wPLI $band: $e"
            _log(log_io, "  WARN: wpli_heatmap_$(band) fallido: $e")
        end
    end
    _log(log_io, "  Guardados: heatmaps wPLI por banda")
    return written
end

# ─── Extras de conectividad: connectivity_summary.json + network_metrics.csv ───

function _save_connectivity_extras(
    conn::ConnectivityMatrix,
    export_dir::String,
    log_io::IO,
    cfg::PipelineConfig,
)::Vector{Tuple{String,String}}
    written  = Tuple{String,String}[]
    bands    = sort(collect(keys(conn.matrices)))
    ch_names = conn.channel_names
    n        = length(ch_names)
    thr      = Float64(get(cfg.graph, "density", 0.1))
    json_dir        = joinpath(export_dir, "json");                mkpath(json_dir)
    tables_conn_dir = joinpath(export_dir, "tables", "connectivity"); mkpath(tables_conn_dir)

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
    open(joinpath(json_dir, "connectivity_summary.json"), "w") do io
        write(io, "{\n" * join(["  " * p for p in pairs], ",\n") * "\n}")
    end
    push!(written, ("connectivity_summary.json", "json/connectivity_summary.json"))

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
    CSV.write(joinpath(tables_conn_dir, "network_metrics.csv"), nm_df)
    push!(written, ("network_metrics.csv", "tables/connectivity/network_metrics.csv"))
    _log(log_io, "  Guardado: connectivity_summary.json + network_metrics.csv")
    return written
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
    log_io::IO,
)::Vector{Tuple{String,String}}
    written   = Tuple{String,String}[]
    ch_names  = conn.channel_names
    n         = length(ch_names)
    alpha     = Float64(get(cfg.surrogates, "alpha", 0.05))
    method    = String(get(cfg.surrogates, "method", "circular_shift"))
    n_sur_cfg = Int(get(cfg.surrogates, "n_surrogates", 200))
    fdr_meth  = String(get(cfg.surrogates, "fdr_method", "bh"))
    seed_v    = Int(get(cfg.surrogates, "seed", 42))

    upper_idx = [(i,j) for i in 1:n for j in (i+1):n]
    json_dir        = joinpath(export_dir, "json");                mkpath(json_dir)
    tables_conn_dir = joinpath(export_dir, "tables", "connectivity"); mkpath(tables_conn_dir)

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
        CSV.write(joinpath(tables_conn_dir, "wpli_observed_$(band).csv"), obs_df)
        push!(written, ("wpli_observed_$(band).csv", "tables/connectivity/wpli_observed_$(band).csv"))

        # ── p-values ──────────────────────────────────────────
        p_df = DataFrame(hcat(ch_names, p_mat), vcat(["channel"], ch_names))
        CSV.write(joinpath(tables_conn_dir, "wpli_pvalues_$(band).csv"), p_df)
        push!(written, ("wpli_pvalues_$(band).csv", "tables/connectivity/wpli_pvalues_$(band).csv"))

        # ── q-values (matriz simétrica) ───────────────────────
        q_mat = zeros(Float64, n, n)
        for (k,(i,j)) in enumerate(upper_idx)
            q_mat[i,j] = q_vec[k]; q_mat[j,i] = q_vec[k]
        end
        q_df = DataFrame(hcat(ch_names, q_mat), vcat(["channel"], ch_names))
        CSV.write(joinpath(tables_conn_dir, "wpli_qvalues_$(band).csv"), q_df)
        push!(written, ("wpli_qvalues_$(band).csv", "tables/connectivity/wpli_qvalues_$(band).csv"))

        # ── Máscara significativa ─────────────────────────────
        sig_int = Int.(sr.sig_mask)
        mask_df = DataFrame(hcat(ch_names, sig_int), vcat(["channel"], ch_names))
        CSV.write(joinpath(tables_conn_dir, "wpli_significant_$(band).csv"), mask_df)
        push!(written, ("wpli_significant_$(band).csv", "tables/connectivity/wpli_significant_$(band).csv"))

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
        CSV.write(joinpath(tables_conn_dir, "surrogate_null_stats_$(band).csv"),
                  DataFrame(null_rows))
        push!(written, ("surrogate_null_stats_$(band).csv",
                        "tables/connectivity/surrogate_null_stats_$(band).csv"))

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
        CSV.write(joinpath(tables_conn_dir, "significant_connections.csv"), sig_df)
    else
        CSV.write(joinpath(tables_conn_dir, "significant_connections.csv"),
            DataFrame(ch_a=String[], ch_b=String[], band=String[],
                      wpli_obs=Float64[], p_value=Float64[],
                      q_value=Float64[], z_score=Float64[]))
    end
    push!(written, ("significant_connections.csv",
                    "tables/connectivity/significant_connections.csv"))

    # ── Guardar QC por banda ──────────────────────────────────
    qc_df = DataFrame(qc_bands)
    CSV.write(joinpath(tables_conn_dir, "surrogate_quality.csv"), qc_df)
    push!(written, ("surrogate_quality.csv", "tables/connectivity/surrogate_quality.csv"))

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
    open(joinpath(json_dir, "surrogate_summary.json"), "w") do io
        write(io, "{\n" * join(["  " * p for p in json_pairs], ",\n") * "\n}")
    end
    push!(written, ("surrogate_summary.json", "json/surrogate_summary.json"))

    _log(log_io, "  Surrogates guardados: $(n_sig_total) conexiones significativas en $(length(surr_results)) bandas")
    return written
end

function _save_all_results(
    rec::EEGRecording, rec_filt::EEGRecording,
    qc_stats::DataFrame, bad_ch::Vector{String},
    spectra::SpectralResult, conn::ConnectivityMatrix,
    val, export_dir::String,
    condition::String, cfg::PipelineConfig, dpi::Int, log_io::IO
)
    # Salida ÚNICA: árbol BIDS (export_dir). El antiguo árbol heredado
    # results/{ID}/{SES}/{tables,figures} se eliminó el 2026-07-21.
    # Desde 2026-07-22: subdivisión interna tables/figures/json + connectivity/
    # (TDD_ejecucion_rutinas.md §3.1) — solo connectivity/ e ica/ tienen
    # subcarpeta propia, por volumen de ficheros.
    # PSD/wPLI/figuras de fase: ya escritos en [6/8] y [7/8].
    tables_dir = joinpath(export_dir, "tables"); mkpath(tables_dir)

    # ── Tablas QC ─────────────────────────────────────────────
    CSV.write(joinpath(tables_dir, "qc_summary.csv"), qc_stats)
    CSV.write(joinpath(tables_dir, "channel_statistics.csv"), qc_stats)
    _log(log_io, "  Guardado: qc_summary.csv, channel_statistics.csv")

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
    CSV.write(joinpath(tables_dir, "overview.csv"), overview_df)
    _log(log_io, "  Guardado: overview.csv")

    # PSD / band power / spectral extras / figuras espectrales:
    # ya escritos en [6/8] vía _save_spectral_phase_results.
    # Matrices wPLI / edges / heatmaps / connectivity extras:
    # ya escritos en [7/8] vía _save_connectivity_phase_results.
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

"""
    _print_kv_table(rows; indent="  ", headers=("Variable","Valor"))

Imprime una tabla Unicode de dos columnas (Variable | Valor).
`rows` es un vector de pares `(label, value)` (ambos se convierten a String).
"""
function _print_kv_table(
    rows::AbstractVector{<:Tuple};
    indent::String = "  ",
    headers::Tuple{String,String} = ("Variable", "Valor"),
)
    isempty(rows) && return
    labs = String[string(r[1]) for r in rows]
    vals = String[string(r[2]) for r in rows]
    w1 = max(length(headers[1]), maximum(length, labs))
    w2 = max(length(headers[2]), maximum(length, vals))
    rule(l, m, r) = l * "─"^(w1 + 2) * m * "─"^(w2 + 2) * r
    row(a, b) = "│ " * rpad(a, w1) * " │ " * rpad(b, w2) * " │"
    println(indent * rule("┌", "┬", "┐"))
    println(indent * row(headers[1], headers[2]))
    println(indent * rule("├", "┼", "┤"))
    for (a, b) in zip(labs, vals)
        println(indent * row(a, b))
    end
    println(indent * rule("└", "┴", "┘"))
end

"""
    _print_cols_table(headers, rows; indent="  ")

Tabla Unicode de N columnas. `headers` y cada fila de `rows` son
vectores de String de la misma longitud.
"""
function _print_cols_table(
    headers::Vector{String},
    rows::Vector{Vector{String}};
    indent::String = "  ",
)
    isempty(headers) && return
    ncols = length(headers)
    widths = [length(h) for h in headers]
    for row in rows
        @assert length(row) == ncols
        for i in 1:ncols
            widths[i] = max(widths[i], length(row[i]))
        end
    end
    rule(l, m, r) = l * join(["─"^(w + 2) for w in widths], m) * r
    fmt(cells) = "│ " * join([rpad(cells[i], widths[i]) for i in 1:ncols], " │ ") * " │"
    println(indent * rule("┌", "┬", "┐"))
    println(indent * fmt(headers))
    println(indent * rule("├", "┼", "┤"))
    for row in rows
        println(indent * fmt(row))
    end
    println(indent * rule("└", "┴", "┘"))
end

"""
    _ica_effective_params(cfg, n_ch) -> NamedTuple

Parámetros que `run_ica` aplicará de verdad. Con `profile="eeg_julia"`
max_iter/tol/seed/n_comp están fijados en ICACore.jl (el TOML no los
sobrescribe). Con `"default"` se leen de `cfg.ica`.
"""
function _ica_effective_params(cfg::PipelineConfig, n_ch::Int)
    ica_cfg = cfg.ica
    profile = String(get(ica_cfg, "profile", "default"))
    if profile == "eeg_julia"
        return (
            profile      = profile,
            n_comp       = n_ch,
            max_iter     = 512,
            tol          = 1e-7,
            seed         = 1234,
            max_attempts = 1,
            a            = 1.0,
            mixing       = "inv(W)  (n_comp = n_ch)",
        )
    else
        n_comp_raw = get(ica_cfg, "n_components", 30)
        n_comp = if n_comp_raw == 0 || n_comp_raw == "auto"
            n_ch
        else
            min(Int(n_comp_raw), n_ch)
        end
        return (
            profile      = profile,
            n_comp       = n_comp,
            max_iter     = Int(get(ica_cfg, "max_iter", 500)),
            tol          = Float64(get(ica_cfg, "tol", 1e-5)),
            seed         = Int(get(ica_cfg, "seed", 42)),
            max_attempts = max(1, Int(get(ica_cfg, "max_attempts", 1))),
            a            = 1.0,
            mixing       = n_comp == n_ch ? "inv(W)" : "pinv(W)",
        )
    end
end

function _save_ica_results(
    ica::ICAResult,
    rec_before::EEGRecording,
    rec_after::EEGRecording,
    export_dir::String,
    duration_s::Float64,
    log_io::IO;
    eval_df::DataFrame = DataFrame(),
    rejection_source::String = "none",
    artifact_thresh::Real = 1.5,
)::Vector{Tuple{String,String}}
    written = Tuple{String,String}[]
    n_comp  = size(ica.activations, 1)
    var_pct = clamp.(ica.variance_explained .* 100.0, 0.0, 100.0)
    fs      = rec_before.meta.fs
    tables_ica_dir = joinpath(export_dir, "tables", "ica"); mkpath(tables_ica_dir)
    json_dir       = joinpath(export_dir, "json");          mkpath(json_dir)

    # ── Clasificación automática de componentes (precomputada por el caller,
    #    mismo cálculo usado para decidir el rechazo — sin duplicar umbral) ──
    if !isempty(eval_df)
        try
            CSV.write(joinpath(tables_ica_dir, "ica_component_features.csv"), eval_df)
            push!(written, ("ica_component_features.csv", "tables/ica/ica_component_features.csv"))
            _log(log_io, "  Guardado: ica_component_features.csv")
        catch e
            @warn "No se pudo guardar ica_component_features.csv: $e"
        end
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
    CSV.write(joinpath(tables_ica_dir, "ica_components.csv"), df)
    push!(written, ("ica_components.csv", "tables/ica/ica_components.csv"))
    _log(log_io, "  Guardado: ica_components.csv")

    # ── Matrices de mezcla / desmezcla ───────────────────────
    ch_names = rec_before.meta.channel_names
    mix_df = DataFrame(ica.mixing_matrix,
                       [Symbol("IC$(i)") for i in 1:n_comp])
    mix_df[!, :channel] = ch_names
    select!(mix_df, :channel, Not(:channel))
    CSV.write(joinpath(tables_ica_dir, "ica_mixing_matrix.csv"), mix_df)

    unmix_df = DataFrame(ica.unmixing_matrix,
                         [Symbol(c) for c in ch_names])
    unmix_df[!, :IC] = ["IC$(i)" for i in 1:n_comp]
    select!(unmix_df, :IC, Not(:IC))
    CSV.write(joinpath(tables_ica_dir, "ica_unmixing_matrix.csv"), unmix_df)
    push!(written, ("ica_mixing_matrix.csv", "tables/ica/ica_mixing_matrix.csv"))
    push!(written, ("ica_unmixing_matrix.csv", "tables/ica/ica_unmixing_matrix.csv"))
    _log(log_io, "  Guardado: ica_mixing_matrix.csv + ica_unmixing_matrix.csv")

    # ── Topomaps por componente ───────────────────────────────
    n_topo = 0
    figs_dir = joinpath(export_dir, "figures", "ica")
    mkpath(figs_dir)
    if rec_before.meta.channel_positions !== nothing
        n_topo = _save_ica_topomaps(ica, rec_before, figs_dir, log_io; artifact_types=artifact_types)
        if n_topo > 0
            push!(written, ("ica_topomap_*.png ($n_topo)", "figures/ica/"))
        end
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

    open(joinpath(json_dir, "ica_summary.json"), "w") do f
        d = ica.diagnostics
        conv_json = if isempty(d)
            """
  "converged": null,
  "n_iter": null,
  "final_error": null,
  "n_attempts": null,
  "seed": null"""
        else
            """
  "converged": $(Bool(get(d, "converged", false))),
  "n_iter": $(Int(get(d, "n_iter", 0))),
  "max_iter": $(Int(get(d, "max_iter", 0))),
  "final_error": $(Float64(get(d, "final_error", NaN))),
  "tol": $(Float64(get(d, "tol", NaN))),
  "n_attempts": $(Int(get(d, "n_attempts", 1))),
  "seed": $(Int(get(d, "seed", 0))),
  "restart_needed": $(Bool(get(d, "restart_needed", false))),
  "optimize_s": $(Float64(get(d, "optimize_s", duration_s))),
  "whitening_cond": $(Float64(get(d, "whitening_cond", NaN))),
  "orthogonality": $(Float64(get(d, "orthogonality", NaN))),
  "pca_var_pct": $(Float64(get(d, "pca_var_pct", NaN)))"""
        end
        write(f, """{
  "n_components": $(n_comp),
  "n_rejected": $(n_rej),
  "n_accepted": $(n_comp - n_rej),
  "variance_retained": $(var_ret),
  "timestamp": "$(ts)",
  "duration_s": $(duration_s),
  "artifact_types": {$(type_json)},
  "n_topomaps": $(n_topo),
  "has_features": $((!isempty(eval_df))),
  "rejection_source": "$(rejection_source)",
  "artifact_threshold": $(Float64(artifact_thresh)),
$(conv_json)
}""")
    end
    push!(written, ("ica_summary.json", "json/ica_summary.json"))
    _log(log_io, "  Guardado: ica_summary.json")

    # ── Activaciones (primeros 10s) ───────────────────────────
    fs      = rec_before.meta.fs
    n_samp  = min(size(ica.activations, 2), Int(round(10.0 * fs)))
    t_vec   = collect(range(0.0, step=1.0/fs, length=n_samp))
    act_df  = DataFrame(:t_s => t_vec)
    for i in 1:n_comp
        act_df[!, Symbol("IC$(i)")] = round.(ica.activations[i, 1:n_samp], digits=4)
    end
    CSV.write(joinpath(tables_ica_dir, "ica_activations.csv"), act_df)
    push!(written, ("ica_activations.csv", "tables/ica/ica_activations.csv"))
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
        CSV.write(joinpath(tables_ica_dir, "ica_signal_$(label).csv"), sdf)
        push!(written, ("ica_signal_$(label).csv", "tables/ica/ica_signal_$(label).csv"))
    end
    _log(log_io, "  Guardado: ica_signal_before.csv + ica_signal_after.csv")
    return written
end

function _save_raw_signal(
    rec::EEGRecording,
    export_dir::String,
    log_io::IO
)
    fs      = rec.meta.fs
    n_total = size(rec.data, 2)
    tables_dir = joinpath(export_dir, "tables"); mkpath(tables_dir)
    # Save the full raw (pre-filter) signal for cross-pipeline validation
    # (no length cap — serves any window the dashboard requests)
    t_vec   = collect(range(0.0, step=1.0/fs, length=n_total))
    sdf     = DataFrame(:t_s => t_vec)
    for (ci, ch) in enumerate(rec.meta.channel_names)
        sdf[!, Symbol(ch)] = round.(rec.data[ci, :], digits=4)
    end
    CSV.write(joinpath(tables_dir, "raw_signal.csv"), sdf)
    _log(log_io, "  Guardado: raw_signal.csv ($(n_total) muestras = $(round(n_total/fs, digits=1))s)")
end

"""
    _save_filtered_signal(rec, export_dir, filter_key, log_io) -> String

Guarda la señal completa tras un paso de filtrado en
`tables/filtered_signal_<filter_key>.csv` (mismo esquema que `raw_signal.csv`:
`t_s` + un canal por columna).

`filter_key` ∈ `notch` | `bandreject` | `highpass` | `lowpass`.
Devuelve el nombre de archivo escrito.
"""
function _save_filtered_signal(
    rec::EEGRecording,
    export_dir::String,
    filter_key::AbstractString,
    log_io::IO,
)::String
    key = lowercase(replace(String(filter_key), r"[^a-z0-9_]+" => "_"))
    isempty(key) && (key = "step")
    fname = "filtered_signal_$(key).csv"
    fs      = rec.meta.fs
    n_total = size(rec.data, 2)
    tables_dir = joinpath(export_dir, "tables"); mkpath(tables_dir)
    t_vec = collect(range(0.0, step=1.0/fs, length=n_total))
    sdf   = DataFrame(:t_s => t_vec)
    for (ci, ch) in enumerate(rec.meta.channel_names)
        sdf[!, Symbol(ch)] = round.(rec.data[ci, :], digits=4)
    end
    CSV.write(joinpath(tables_dir, fname), sdf)
    _log(log_io, "  Guardado: $fname ($(n_total) muestras = $(round(n_total/fs, digits=1))s)")
    return fname
end

function _save_ica_topomaps(
    ica::ICAResult,
    rec::EEGRecording,
    figs_dir::String,
    log_io::IO;
    artifact_types::Vector{String} = String[],
)::Int
    ch_pos   = rec.meta.channel_positions
    ch_names = rec.meta.channel_names
    n_comp   = size(ica.mixing_matrix, 2)
    n_saved  = 0

    var_pct = clamp.(ica.variance_explained .* 100.0, 0.0, 100.0)
    for ic in 1:n_comp
        weights = Float64.(ica.mixing_matrix[:, ic])
        fname   = "ica_topomap_$(lpad(ic, 3, '0')).png"
        fpath   = joinpath(figs_dir, fname)
        try
            vp = ic ≤ length(var_pct) ? round(var_pct[ic]; digits=1) : NaN
            typ = (ic ≤ length(artifact_types)) ? artifact_types[ic] : ""
            title = isempty(typ) ?
                "IC$(lpad(ic, 3, '0')) · $(vp)% var" :
                "IC$(lpad(ic, 3, '0')) · $(typ) · $(vp)%"
            fig = plot_topomap(weights, ch_names, ch_pos;
                               title = title,
                               style = :ica,
                               contours = true,
                               show_channel_labels = true,
                               label_fontsize = 8.5)
            CairoMakie.save(fpath, fig)
            n_saved += 1
        catch e
            @warn "Topomap IC$ic falló: $e"
        end
    end
    _log(log_io, "  Topomaps: $(n_saved)/$(n_comp) guardados en figures/ica/")
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

# ─── Tabla QC global ──────────────────────────────────────────

"""
    _update_qc_decision_table(results_dir, ...)

Escribe o actualiza results/qc/qc_decision_table.csv con la decisión
de inclusión/exclusión de cada grabación procesada.

### Criterios de `final_decision`
- `include`              — válido sin alertas
- `include_with_warning` — válido pero con amplitude_warning o canales malos adicionales
- `manual_review`        — amplitude_warning Y epochs cerca del mínimo, o canales malos graves
- `exclude`              — 0 epochs válidos o duración insuficiente
"""
function _update_qc_decision_table(
    res_dir        ::String,
    subject_id     ::String,
    session_id     ::String,
    condition      ::String,
    amplitude_warning ::Bool,
    raw_sigma_uv   ::Float64,
    bad_ch_original::Vector{String},
    bad_ch_non_fp2 ::Vector{String},
    fp2_excluded   ::Bool,
    n_ch_analysis  ::Int,
    ar_threshold_uv::Float64,
    n_epochs_total ::Int,
    n_epochs_rejected::Int,
    n_epochs_valid ::Int,
    min_epochs     ::Int,
    log_io         ::IO
)
    qc_dir  = joinpath(res_dir, "qc")
    mkpath(qc_dir)
    out_path = joinpath(qc_dir, "qc_decision_table.csv")

    valid_pct = n_epochs_total > 0 ?
        round(100.0 * n_epochs_valid / n_epochs_total, digits=1) : 0.0

    # ── Decisión final ────────────────────────────────────────
    # Exclusión dura: sin epochs utilizables
    exclusion_reason = ""
    notes_list       = String[]

    final_decision = if n_epochs_valid == 0
        exclusion_reason = "0 epochs válidos tras AR ±$(ar_threshold_uv) µV"
        "exclude"
    elseif n_epochs_valid < min_epochs
        exclusion_reason = "$(n_epochs_valid) epochs < mínimo requerido ($(min_epochs))"
        "exclude"
    elseif amplitude_warning && length(bad_ch_non_fp2) >= 2
        push!(notes_list, "amplitude_warning + $(length(bad_ch_non_fp2)) canales malos adicionales")
        "manual_review"
    elseif amplitude_warning && valid_pct < 50.0
        push!(notes_list, "amplitude_warning + solo $(valid_pct)% epochs válidos")
        "manual_review"
    elseif amplitude_warning || !isempty(bad_ch_non_fp2)
        amplitude_warning && push!(notes_list, "amplitude_warning (σ̄=$(round(raw_sigma_uv,digits=1)) µV)")
        !isempty(bad_ch_non_fp2) && push!(notes_list, "canales malos adicionales: $(join(bad_ch_non_fp2, ", "))")
        "include_with_warning"
    else
        "include"
    end

    notes = join(notes_list, "; ")

    new_row = DataFrame(
        subject_id               = [subject_id],
        session_id               = [session_id],
        condition                = [condition],
        amplitude_warning        = [amplitude_warning],
        raw_sigma_mean_uv        = [round(raw_sigma_uv, digits=2)],
        bad_channels_original    = [isempty(bad_ch_original) ? "" : join(bad_ch_original, "; ")],
        fp2_removed              = [fp2_excluded],
        n_channels_analysis      = [n_ch_analysis],
        n_bad_channels_non_fp2   = [length(bad_ch_non_fp2)],
        bad_channels_non_fp2     = [isempty(bad_ch_non_fp2) ? "" : join(bad_ch_non_fp2, "; ")],
        ar_threshold_uv          = [ar_threshold_uv],
        n_epochs_initial         = [n_epochs_total],
        n_epochs_rejected        = [n_epochs_rejected],
        n_epochs_valid           = [n_epochs_valid],
        valid_epochs_pct         = [valid_pct],
        final_decision           = [final_decision],
        exclusion_reason         = [exclusion_reason],
        notes                    = [notes],
        processed_at             = [string(now())],
    )

    if isfile(out_path)
        old = CSV.read(out_path, DataFrame)
        mask = .!(old.subject_id .== subject_id .&&
                  old.session_id .== session_id .&&
                  old.condition  .== condition)
        CSV.write(out_path, vcat(old[mask, :], new_row))
    else
        CSV.write(out_path, new_row)
    end

    _log(log_io, "  QC decision: $(final_decision)" *
        (isempty(exclusion_reason) ? "" : " ($(exclusion_reason))") *
        (isempty(notes) ? "" : " | $(notes)"))
end
