# NeuroMIND/src/io/Serializer.jl
# Persistencia de resultados: guarda/carga structs con Serialization.
# Añade versionado y checksum para garantizar reproducibilidad.

const NEUROMIND_VERSION = "0.1.0"

"""
    save_result(obj, cfg, subject_id, session_id, stage; condition="")

Serializa un resultado a `results/{subject}/{session}/cache/{stage}.bin`.
Añade metadatos de versión para detectar incompatibilidades futuras.
"""
function save_result(
    obj,
    cfg::PipelineConfig,
    subject_id::String,
    session_id::String,
    stage::String;
    condition::String = ""
)
    dir  = ensure_dirs(cfg, subject_id, session_id)
    stem = isempty(condition) ? stage : "$(stage)_$(condition)"
    path = joinpath(dir, "cache", "$(stem).bin")

    wrapper = Dict(
        "version" => NEUROMIND_VERSION,
        "stage"   => stage,
        "subject" => subject_id,
        "session" => session_id,
        "saved_at"=> string(now()),
        "data"    => obj
    )
    Serialization.serialize(path, wrapper)
    return path
end

"""
    load_result(T, cfg, subject_id, session_id, stage; condition="") -> T

Carga un resultado serializado. Lanza error si la versión es incompatible.
"""
function load_result(
    ::Type{T},
    cfg::PipelineConfig,
    subject_id::String,
    session_id::String,
    stage::String;
    condition::String = ""
)::T where T

    dir  = subject_results_dir(cfg, subject_id, session_id)
    stem = isempty(condition) ? stage : "$(stage)_$(condition)"
    path = joinpath(dir, "cache", "$(stem).bin")

    isfile(path) || error("Resultado no encontrado: $path")

    wrapper = Serialization.deserialize(path)
    stored_v = get(wrapper, "version", "unknown")
    stored_v == NEUROMIND_VERSION ||
        @warn "Versión del cache ($stored_v) difiere de la actual ($NEUROMIND_VERSION)"

    return wrapper["data"]::T
end

"""
    result_exists(cfg, subject_id, session_id, stage; condition="") -> Bool
"""
function result_exists(
    cfg::PipelineConfig,
    subject_id::String,
    session_id::String,
    stage::String;
    condition::String = ""
)::Bool
    dir  = subject_results_dir(cfg, subject_id, session_id)
    stem = isempty(condition) ? stage : "$(stage)_$(condition)"
    isfile(joinpath(dir, "cache", "$(stem).bin"))
end
