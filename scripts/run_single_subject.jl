#!/usr/bin/env julia
# NeuroMIND — Pipeline para un solo sujeto EEG
#
# Uso:
#   julia --project=. scripts/run_single_subject.jl
#   julia --project=. scripts/run_single_subject.jl --config config/pipeline.toml
#   julia --project=. scripts/run_single_subject.jl --force
#
# Config: config/pipeline.toml — fichero unificado (2026-07-21) que sustituye
# a single_subject.toml y batch_pipeline.toml, archivados en deprecated/code/config/.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()   # instala dependencias si no están presentes

using NeuroMIND

# ─── Argumentos de línea de comandos ──────────────────────────
args = Dict{String,String}()
let i = 1
    while i <= length(ARGS)
        if startswith(ARGS[i], "--")
            k = ARGS[i][3:end]
            if i < length(ARGS) && !startswith(ARGS[i+1], "--")
                i += 1
                args[k] = ARGS[i]
            else
                args[k] = "true"
            end
        end
        i += 1
    end
end

config_path = get(args, "config",
    joinpath(@__DIR__, "..", "config", "pipeline.toml"))

# ─── Ejecución ────────────────────────────────────────────────
run_single_subject_pipeline(config_path)
