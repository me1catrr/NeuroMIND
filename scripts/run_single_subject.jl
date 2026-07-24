# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Pipeline para un solo sujeto EEG
# ═══════════════════════════════════════════════════════════════
#
#  Ejecuta el pipeline de 8 pasos sobre una grabación
#  (run_single_subject_pipeline → SingleSubjectPipeline.jl).
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/run_single_subject.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 22-07-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/<este-script>.jl …
#  Sin shebang: #!/usr/bin/env julia no activaría --project=.
#
#  Config
#  ──────
#  Por defecto: config/pipeline.toml
#
#  El bloque [subject] define sujeto / sesión / tarea.
#  subject_id = "auto" → primer sujeto en {paths.bids_root}/raw/
#
#  Uso
#  ───
#    julia --project=. scripts/run_single_subject.jl
#    julia --project=. scripts/run_single_subject.jl --config ruta.toml
#
#  Opciones
#  ────────
#    --config  PATH   TOML de parámetros (defecto: config/pipeline.toml)
#
#  Salida
#  ──────
#    results/subjects/sub-{id}/ses-{sess}/{task}/

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
