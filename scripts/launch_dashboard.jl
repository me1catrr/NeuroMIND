# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Dashboard interactivo (Genie.jl)
# ═══════════════════════════════════════════════════════════════
#
#  Arranca el servidor local del dashboard.
#  Lee [dashboard] de config/pipeline.toml; --port en CLI tiene
#  prioridad sobre el TOML.
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/launch_dashboard.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 22-07-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/<este-script>.jl …
#  Sin shebang: #!/usr/bin/env julia no activaría --project=.
#
#  Config
#  ──────
#    config/pipeline.toml  →  sección [dashboard]
#
#  Uso
#  ───
#    julia --project=. scripts/launch_dashboard.jl
#    julia --project=. scripts/launch_dashboard.jl --port 9090
#    julia --project=. scripts/launch_dashboard.jl --no-browser
#
#  Opciones
#  ────────
#    --port N        puerto HTTP (defecto: 8080 o [dashboard].port)
#    --no-browser    no abrir el navegador al arrancar
#
#  Salida
#  ──────
#    http://localhost:<port>

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
