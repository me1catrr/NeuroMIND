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
root      = joinpath(@__DIR__, "..")
ss_config = joinpath(root, "config", "single_subject.toml")
# Fallback legacy: solo se usa si single_subject.toml no existe.
# pipeline.toml se archivó en legacy/config/ el 2026-07-21.
pl_config = joinpath(root, "legacy", "config", "pipeline.toml")

cfg = if isfile(ss_config)
    load_ss_config(ss_config)
else
    load_config(pl_config)
end

# Puerto desde single_subject.toml si no se pasó --port
if !any(a -> a == "--port", ARGS) && isfile(ss_config)
    raw = TOML.parsefile(ss_config)
    port      = Int(get(get(raw, "dashboard", Dict()), "port", port))
    open_brsr = Bool(get(get(raw, "dashboard", Dict()), "open_browser", open_brsr))
end

println("Iniciando NeuroMIND Dashboard → http://localhost:$(port)")
launch_webapp(cfg; port, open_browser=open_brsr)

try
    while true; sleep(3600); end
catch e
    e isa InterruptException || rethrow(e)
end
