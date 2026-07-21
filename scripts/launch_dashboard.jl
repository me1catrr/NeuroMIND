#!/usr/bin/env julia
# NeuroMIND — Dashboard interactivo
#
# Uso:
#   julia --project=. scripts/launch_dashboard.jl
#   julia --project=. scripts/launch_dashboard.jl --port 9090
#   julia --project=. scripts/launch_dashboard.jl --no-browser

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()   # instala dependencias si no están presentes

using NeuroMIND
using TOML

# ─── Argumentos ───────────────────────────────────────────────
port      = 8080
open_brsr = true
let args = ARGS
    for i in eachindex(args)
        if args[i] == "--port" && i < length(args)
            port = parse(Int, args[i+1])
        elseif args[i] == "--no-browser"
            open_brsr = false
        end
    end
end

# ─── Configuración ────────────────────────────────────────────
# Fuente única: config/pipeline.toml (unificado 2026-07-21; sustituye a
# single_subject.toml y batch_pipeline.toml, archivados en deprecated/code/config/).
root        = joinpath(@__DIR__, "..")
config_path = joinpath(root, "config", "pipeline.toml")

isfile(config_path) ||
    error("Config no encontrada: $(config_path)\n" *
          "→ config/pipeline.toml es la configuración única del proyecto.")

cfg = load_ss_config(config_path)

# [dashboard] del TOML — --port en CLI tiene prioridad
if !any(a -> a == "--port", ARGS)
    raw       = TOML.parsefile(config_path)
    dash      = get(raw, "dashboard", Dict())
    port      = Int(get(dash, "port", port))
    open_brsr = Bool(get(dash, "open_browser", open_brsr))
end

println("Iniciando NeuroMIND Dashboard → http://localhost:$(port)")
launch_webapp(cfg; port, open_browser=open_brsr)

try
    while true; sleep(3600); end
catch e
    e isa InterruptException || rethrow(e)
end
