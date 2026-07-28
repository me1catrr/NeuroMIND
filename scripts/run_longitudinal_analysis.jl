# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Lanzador: análisis longitudinal (T1 → T2, solo EM)
# ═══════════════════════════════════════════════════════════════
#
#  Lanzador fino — toda la lógica (estadística + figuras, en una
#  única pasada) vive en src/longitudinal/Longitudinal.jl.
#
#  Antes de lanzar, borra results/longitudinal/ por completo para
#  no dejar residuos (p. ej. EC/EO antiguos u otros huérfanos).
#  Longitudinal.run vuelve a crear eyesclosed/eyesopen.
#
#  Diseño experimental (Fig. 3.1):
#    · Solo pacientes EM con par completo T1+T2 (N diseño = 30)
#    · Controles NO entran; eyesclosed y eyesopen en paralelo, sin pooling
#
#  Salida:
#    results/longitudinal/{eyesclosed|eyesopen}/
#      config_snapshot.toml · longitudinal_summary.json
#      tables/ (CSV + figures_manifest.tsv) · figures/ (PNG)
#    results/summary/ — figura de síntesis con el transversal, si
#      ese análisis ya tiene resultados en disco
#
#  Uso
#  ───
#    julia --project=. scripts/run_longitudinal_analysis.jl
#    julia --project=. scripts/run_longitudinal_analysis.jl config/pipeline.toml
#
# ───────────────────────────────────────────────────────────────
#  Fichero     scripts/run_longitudinal_analysis.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      22-05-2026
#  Modificado  27-07-2026
# ───────────────────────────────────────────────────────────────

using TOML

include(joinpath(@__DIR__, "..", "src", "longitudinal", "Longitudinal.jl"))
using .Longitudinal

const PROJ     = dirname(@__DIR__)
const CONFIG_P = length(ARGS) > 0 ? ARGS[1] : joinpath(PROJ, "config", "pipeline.toml")

# Vaciar salida previa (carpeta completa → no quedan EC/EO ni huérfanos)
let
    paths_raw = get(TOML.parsefile(CONFIG_P), "paths", Dict{String,Any}())
    r = get(paths_raw, "results", "results")
    res_root = isabspath(r) ? r : joinpath(PROJ, r)
    long_dir = joinpath(res_root, "longitudinal")
    if isdir(long_dir)
        println("  Limpiando $long_dir …")
        rm(long_dir; recursive=true, force=true)
    end
end

Longitudinal.run(CONFIG_P)
