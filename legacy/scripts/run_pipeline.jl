#!/usr/bin/env julia
# NeuroMIND/scripts/run_pipeline.jl
#
# Script de entrada para ejecutar el pipeline completo desde la línea de comandos.
#
# Uso:
#   julia scripts/run_pipeline.jl                      # todos los sujetos
#   julia scripts/run_pipeline.jl --subject M05        # sujeto específico
#   julia scripts/run_pipeline.jl --force              # fuerza re-cómputo
#   julia scripts/run_pipeline.jl --subject M05 --session T2 --condition EC

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using NeuroMIND

# ─── Parseo de argumentos ─────────────────────────────────────
args   = Dict{String,String}()
i = 1
while i <= length(ARGS)
    if startswith(ARGS[i], "--")
        k = ARGS[i][3:end]
        args[k] = i < length(ARGS) && !startswith(ARGS[i+1], "--") ? (i += 1; ARGS[i]) : "true"
    end
    i += 1
end

force     = get(args, "force", "false") == "true"
subj_filt = get(args, "subject", "")
sess_filt = get(args, "session", "")
cond_filt = get(args, "condition", "")

# ─── Configuración ────────────────────────────────────────────
cfg_path = joinpath(@__DIR__, "..", "config", "pipeline.toml")
cfg      = load_config(cfg_path)
subjects = load_subjects(cfg)

# Filtrar sujetos/sesiones/condiciones si se especificaron
if !isempty(subj_filt)
    subjects = filter(s -> s.id == subj_filt, subjects)
    isempty(subjects) && error("Sujeto '$subj_filt' no encontrado en subjects.toml")
end

conditions = isempty(cond_filt) ?
    get(cfg.recording, "conditions", ["EO", "EC"]) : [cond_filt]

# ─── Ejecución ────────────────────────────────────────────────
println("=" ^ 60)
println("NeuroMIND — Pipeline EEG Conectividad")
println("=" ^ 60)
println("Sujetos: $(length(subjects))")
println("Condiciones: $(join(conditions, ", "))")
println("Force: $force")
println("=" ^ 60)
println()

for subj in subjects
    println("▶ $(subj.id) [$(subj.group)]")
    sess_ids = isempty(sess_filt) ? sort(collect(keys(subj.sessions))) : [sess_filt]
    for sess_id in sess_ids
        haskey(subj.sessions, sess_id) || continue
        run_session!(subj, sess_id, conditions, cfg; force, verbose=true)
    end
end

println()
println("=" ^ 60)
println("✅ Pipeline completado")
println("=" ^ 60)
