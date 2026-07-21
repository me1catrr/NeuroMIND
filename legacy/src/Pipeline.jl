# NeuroMIND/src/Pipeline.jl
# Orquestador de alto nivel: ejecuta todo el pipeline para uno o varios sujetos.
# Cada etapa puede omitirse si el resultado ya está en cache.

"""
    run_pipeline!(subjects::Vector{Subject}, cfg::PipelineConfig;
                  force=false, verbose=true)

Ejecuta el pipeline completo para todos los sujetos y condiciones configuradas.
`force=true` re-ejecuta todas las etapas aunque exista cache.
"""
function run_pipeline!(
    subjects::Vector{Subject},
    cfg::PipelineConfig;
    force::Bool = false,
    verbose::Bool = true
)
    conditions = get(cfg.recording, "conditions", ["EO", "EC"])

    for subj in subjects
        verbose && println("▶ Sujeto: $(subj.id) [$(subj.group)]")
        for sess_id in sort(collect(keys(subj.sessions)))
            verbose && println("  ↳ Sesión: $sess_id")
            run_session!(subj, sess_id, conditions, cfg; force, verbose)
        end
    end
    verbose && println("✅ Pipeline completado")
end

"""
    run_session!(subj::Subject, sess_id::String,
                 conditions::Vector{String}, cfg::PipelineConfig;
                 force=false, verbose=true)

Ejecuta todas las etapas del pipeline para una sesión.
"""
function run_session!(
    subj::Subject,
    sess_id::String,
    conditions::Vector{String},
    cfg::PipelineConfig;
    force::Bool = false,
    verbose::Bool = true
)
    sess = subj.sessions[sess_id]
    ensure_dirs(cfg, subj.id, sess_id)

    for cond in conditions
        verbose && println("    ↳ Condición: $cond")
        run_subject!(subj, sess, cond, cfg; force, verbose)
    end
end

"""
    run_subject!(subj, sess, condition, cfg; force, verbose)

Ejecuta etapas individuales. Cachea cada resultado en disco.
"""
function run_subject!(
    subj::Subject,
    sess::Session,
    condition::String,
    cfg::PipelineConfig;
    force::Bool = false,
    verbose::Bool = true
)
    id, sess_id = subj.id, sess.id

    # ── 1. Carga BIDS ──────────────────────────────────────────
    rec = _cached(EEGRecording, cfg, id, sess_id, "raw", condition; force) do
        verbose && println("      [1/7] Cargando BIDS...")
        load_eeg_bids(cfg, id, sess_id, condition)
    end
    sess.recordings[condition] = rec

    # ── 2. QC inicial ──────────────────────────────────────────
    _cached(DataFrame, cfg, id, sess_id, "qc", condition; force) do
        verbose && println("      [2/7] QC de canales...")
        qc_report(rec, cfg)
    end

    # ── 3. Filtrado ────────────────────────────────────────────
    rec_filt = _cached(EEGRecording, cfg, id, sess_id, "filtered", condition; force) do
        verbose && println("      [3/7] Filtrado...")
        filter_recording(rec, cfg)
    end

    # ── 4. ICA ─────────────────────────────────────────────────
    ica = _cached(ICAResult, cfg, id, sess_id, "ica", condition; force) do
        verbose && println("      [4/7] ICA...")
        run_ica(rec_filt, cfg)
    end
    sess.ica[condition] = ica

    rejected = load_ica_labels(cfg, id, sess_id, condition)
    rec_clean = apply_ica_rejection(rec_filt, ica, rejected)

    # ── 5. Segmentación + Baseline + AR ────────────────────────
    epochs = _cached(EpochSet, cfg, id, sess_id, "epochs", condition; force) do
        verbose && println("      [5/7] Segmentación...")
        ep = segment_recording(rec_clean, cfg)
        ep = apply_baseline(ep, cfg)
        reject_artifacts(ep, cfg)
    end
    sess.epochs[condition] = epochs

    # ── 6. Spectral ────────────────────────────────────────────
    spectra = _cached(SpectralResult, cfg, id, sess_id, "spectra", condition; force) do
        verbose && println("      [6/7] PSD...")
        compute_psd(epochs, cfg)
    end
    sess.spectra[condition] = spectra

    # ── 7. Conectividad wPLI (+ CSD opcional) ──────────────────
    conn = _cached(ConnectivityMatrix, cfg, id, sess_id, "wpli", condition; force) do
        verbose && println("      [7/7] wPLI...")
        epochs_for_conn = apply_csd(epochs, cfg)  # no-op si use_csd=false
        compute_wpli(epochs_for_conn, cfg)
    end
    sess.connectivity[condition] = conn

    verbose && println("      ✓ Sesión completada")
end

# ─── Helper de cache transparente ────────────────────────────

function _cached(
    f::Function,
    ::Type{T},
    cfg::PipelineConfig,
    subj_id::String,
    sess_id::String,
    stage::String,
    condition::String;
    force::Bool = false
)::T where T

    if !force && result_exists(cfg, subj_id, sess_id, stage; condition)
        return load_result(T, cfg, subj_id, sess_id, stage; condition)
    end
    result = f()
    save_result(result, cfg, subj_id, sess_id, stage; condition)
    return result
end
