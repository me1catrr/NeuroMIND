# NeuroMIND/src/io/Config.jl

function load_config(path::String)::PipelineConfig
    abs_path = isabspath(path) ? path : abspath(path)
    isfile(abs_path) || error("Config no encontrada: $abs_path")
    raw  = TOML.parsefile(abs_path)
    root = dirname(dirname(abs_path))  # sube config/ → NeuroMIND/

    bands = Dict{String,Tuple{Float64,Float64}}()
    for (k, v) in get(raw, "bands", Dict())
        bands[k] = (Float64(v[1]), Float64(v[2]))
    end

    PipelineConfig(
        get(raw, "project",            Dict()),
        get(raw, "study",              Dict()),
        get(raw, "paths",              Dict()),
        get(raw, "recording",          Dict()),
        get(raw, "filtering",          Dict()),
        get(raw, "segmentation",       Dict()),
        get(raw, "baseline",           Dict()),
        get(raw, "artifact_rejection", Dict()),
        get(raw, "ica",                Dict()),
        get(raw, "spectral",           Dict()),
        bands,
        get(raw, "connectivity",       Dict()),
        get(raw, "surrogates",         Dict()),
        get(raw, "graph",              Dict()),
        get(raw, "clinical",           Dict()),
        get(raw, "longitudinal",       Dict()),
        get(raw, "statistics",         Dict()),
        get(raw, "export",             Dict()),
        root,
    )
end

function load_subjects(cfg::PipelineConfig)::Vector{Subject}
    path = joinpath(cfg.root, "config", "subjects.toml")
    isfile(path) || error("subjects.toml no encontrado: $path")
    raw = TOML.parsefile(path)

    subjects = Subject[]
    for s in get(raw, "subjects", [])
        clin_raw = get(s, "clinical", Dict())
        clin = ClinicalData(
            _f(clin_raw, "EDSS"),
            _f(clin_raw, "disease_duration_y"),
            _s(clin_raw, "medication"),
            _f(clin_raw, "fatigue_score"),
            _f(clin_raw, "cognition_score"),
            _f(clin_raw, "lesion_load"),
        )

        subj = Subject(
            s["id"], s["group"],
            get(s, "age", missing),
            get(s, "sex", missing),
            clin,
            Dict{String,Session}()
        )

        for (i, sess_id) in enumerate(get(s, "sessions", String[]))
            subj.sessions[sess_id] = Session(sess_id, i)
        end
        push!(subjects, subj)
    end
    return subjects
end

_f(d, k) = haskey(d, k) ? Float64(d[k]) : missing
_s(d, k) = haskey(d, k) ? String(d[k]) : missing

# ─── Accesores de rutas ───────────────────────────────────────

results_dir(cfg::PipelineConfig)    = joinpath(cfg.root, cfg.paths["results"])
data_cache_dir(cfg::PipelineConfig) = joinpath(cfg.root, cfg.paths["data_cache"])
bids_root_dir(cfg::PipelineConfig)  = joinpath(cfg.root, cfg.paths["bids_root"])
web_public_dir(cfg::PipelineConfig) = joinpath(cfg.root, get(cfg.paths, "web_public", "web/public"))

function subject_results_dir(cfg::PipelineConfig, subject_id::String, session_id::String)
    # Árbol ÚNICO: BIDS. Antes devolvía results/{ID}/{SES} (árbol heredado,
    # eliminado el 2026-07-21). Ahora apunta al nivel sub-/ses- de BIDS para
    # que ningún consumidor pueda recrear el árbol paralelo.
    joinpath(results_dir(cfg), "subjects", "sub-$(subject_id)", "ses-$(session_id)")
end

function ensure_dirs(cfg::PipelineConfig, subject_id::String, session_id::String)
    base = subject_results_dir(cfg, subject_id, session_id)
    for sub in ["figures", "tables", "logs", "cache", "reports"]
        mkpath(joinpath(base, sub))
    end
    # Directorios compartidos de grupo
    # Nivel grupo: transversal/ y longitudinal/ son hermanos de subjects/
    for grp in ("transversal", "longitudinal")
        mkpath(joinpath(results_dir(cfg), grp, "figures"))
        mkpath(joinpath(results_dir(cfg), grp, "tables"))
    end
    base
end
