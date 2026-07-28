# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Lanzador: análisis transversal (EM vs Control, T1)
# ═══════════════════════════════════════════════════════════════
#
#  Lanzador fino — toda la lógica (estadística + figuras, en una
#  única pasada) vive en src/transversal/Transversal.jl.
#
#  Antes de lanzar, borra results/transversal/ por completo para
#  no dejar residuos (p. ej. EC/EO antiguos u otros huérfanos).
#  Transversal.run vuelve a crear eyesclosed/eyesopen.
#
#  Diseño experimental (Fig. 3.1):
#    · Caso-control en T1: EM (N diseño=44) vs Control (N diseño=40)
#    · eyesclosed y eyesopen en paralelo, sin pooling
#
#  Salida:
#    results/transversal/{eyesclosed|eyesopen}/
#      config_snapshot.toml · transversal_summary.json
#      tables/ (CSV + figures_manifest.tsv) · figures/ (PNG)
#    results/transversal/combined/ — interacción grupo×condición (EC−EO),
#      tables/ + figures/, generada tras procesar ambas condiciones
#    results/summary/ — figura de síntesis con el longitudinal, si
#      ese análisis ya tiene resultados en disco
#
#  Uso
#  ───
#    julia --project=. scripts/run_transversal_analysis.jl
#    julia --project=. scripts/run_transversal_analysis.jl config/pipeline.toml
#
# ───────────────────────────────────────────────────────────────
#  Fichero     scripts/run_transversal_analysis.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      22-05-2026
#  Modificado  27-07-2026
# ───────────────────────────────────────────────────────────────

using TOML

include(joinpath(@__DIR__, "..", "src", "transversal", "Transversal.jl"))
using .Transversal

const PROJ     = dirname(@__DIR__)
const CONFIG_P = length(ARGS) > 0 ? ARGS[1] : joinpath(PROJ, "config", "pipeline.toml")

# Vaciar salida previa (carpeta completa → no quedan EC/EO ni huérfanos)
let
    paths_raw = get(TOML.parsefile(CONFIG_P), "paths", Dict{String,Any}())
    r = get(paths_raw, "results", "results")
    res_root = isabspath(r) ? r : joinpath(PROJ, r)
    trans_dir = joinpath(res_root, "transversal")
    if isdir(trans_dir)
        println("  Limpiando $trans_dir …")
        rm(trans_dir; recursive=true, force=true)
    end
end

Transversal.run(CONFIG_P)
