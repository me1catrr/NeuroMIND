#!/usr/bin/env julia
# NeuroMIND — scripts/run_batch_pipeline.jl
#
# FASE C: Ejecución del pipeline individual para todos los sujetos del
# dataset completo MINDEM-IMIBIC.
#
# Itera sobre el inventory.csv generado por audit_full_dataset.jl,
# crea un config temporal por sujeto/sesión/condición y llama a
# run_single_subject_pipeline().  Los resultados se guardan en:
#   results/subjects/sub-{id}/ses-{sess}/{task}/
#
# Uso:
#   julia --project=. scripts/run_batch_pipeline.jl
#   julia --project=. scripts/run_batch_pipeline.jl --condition EC
#   julia --project=. scripts/run_batch_pipeline.jl --condition EO
#   julia --project=. scripts/run_batch_pipeline.jl --group MS --session T1
#   julia --project=. scripts/run_batch_pipeline.jl --subjects M11,M12,M13
#   julia --project=. scripts/run_batch_pipeline.jl --dry-run
#
# Opciones de filtrado (combinables):
#   --condition  EC | EO | ALL  (defecto: ALL)
#   --group      MS | HC | ALL  (defecto: ALL)
#   --session    T1 | T2 | ALL  (defecto: ALL)
#   --subjects   lista separada por comas de subject_id (ej. M11,M12,MC01)
#   --skip-done  no reprocesar si ya existe overview.csv en results/
#   --dry-run    solo listar lo que se procesaría, sin ejecutar
#   --max-subjects  N  procesar solo los primeros N sujetos
#
# Salida adicional:
#   logs/batch_run_YYYY-MM-DD_HH-MM.csv

using Dates, TOML

const PROJ_ROOT  = dirname(@__DIR__)
const BIDS_DIR   = joinpath(PROJ_ROOT, "data", "bids")
const RESULTS    = joinpath(PROJ_ROOT, "results")
const INVENTORY  = joinpath(PROJ_ROOT, "data", "full_data", "inventory.csv")
const BASE_CFG   = joinpath(PROJ_ROOT, "config", "batch_pipeline.toml")
const LOGS_DIR   = joinpath(PROJ_ROOT, "logs")

# ─── Carga de NeuroMIND ───────────────────────────────────────

# Asegurar que el módulo está en el load path
push!(LOAD_PATH, PROJ_ROOT)
include(joinpath(PROJ_ROOT, "src", "NeuroMIND.jl"))
using .NeuroMIND

# ─── Parseo de argumentos ─────────────────────────────────────

function parse_cli_args()
    args = Dict{String,Any}(
        "condition"    => "ALL",
        "group"        => "ALL",
        "session"      => "ALL",
        "subjects"     => String[],
        "skip_done"    => false,
        "dry_run"      => false,
        "max_subjects" => typemax(Int),
    )
    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if arg == "--condition" && i+1 <= length(ARGS)
            args["condition"] = uppercase(ARGS[i+1]); i += 2
        elseif arg == "--group" && i+1 <= length(ARGS)
            args["group"] = uppercase(ARGS[i+1]); i += 2
        elseif arg == "--session" && i+1 <= length(ARGS)
            args["session"] = uppercase(ARGS[i+1]); i += 2
        elseif arg == "--subjects" && i+1 <= length(ARGS)
            args["subjects"] = String.(split(ARGS[i+1], ",")); i += 2
        elseif arg == "--skip-done"
            args["skip_done"] = true; i += 1
        elseif arg == "--dry-run"
            args["dry_run"] = true; i += 1
        elseif arg == "--max-subjects" && i+1 <= length(ARGS)
            args["max_subjects"] = parse(Int, ARGS[i+1]); i += 2
        else
            i += 1
        end
    end
    return args
end

# ─── Lectura del inventario ───────────────────────────────────

struct BatchJob
    subject_id  :: String
    bids_id     :: String
    group       :: String
    session     :: String
    bids_session:: String
    condition   :: String
    vhdr_path   :: String
end

function load_jobs(inventory_path::String, cli::Dict)::Vector{BatchJob}
    jobs = BatchJob[]
    isfile(inventory_path) ||
        error("Inventario no encontrado: $inventory_path\n→ Ejecuta primero audit_full_dataset.jl")

    filter_cond    = cli["condition"]
    filter_group   = cli["group"]
    filter_session = cli["session"]
    filter_subjs   = Set(cli["subjects"])

    for line in readlines(inventory_path)[2:end]
        isempty(strip(line)) && continue
        parts = split(line, ",")
        length(parts) < 16 && continue

        excluded = strip(String(parts[12])) == "true"
        excluded && continue

        condition = strip(String(parts[9]))
        condition ∈ ("EC","EO") || continue

        subject_id  = strip(String(parts[2]))
        bids_id     = strip(String(parts[3]))
        group       = strip(String(parts[4]))
        session     = strip(String(parts[5]))
        bids_session= strip(String(parts[7]))
        filepath    = strip(String(parts[16]))

        # Filtros
        filter_cond    != "ALL" && condition != filter_cond    && continue
        filter_group   != "ALL" && group     != filter_group   && continue
        filter_session != "ALL" && session   != filter_session && continue
        !isempty(filter_subjs) && subject_id ∉ filter_subjs   && continue

        push!(jobs, BatchJob(subject_id, bids_id, group,
                             session, bids_session, condition, filepath))
    end

    return jobs
end

# ─── Check si ya está procesado ──────────────────────────────

function already_done(bids_id::String, sess::String, cond::String)::Bool
    task = cond == "EC" ? "eyesclosed" : "eyesopen"
    overview = joinpath(RESULTS, "subjects",
                        "sub-$(bids_id)", "ses-$(sess)", task, "overview.csv")
    isfile(overview)
end

# ─── Generación de config temporal ───────────────────────────

function write_temp_config(job::BatchJob, base_cfg_raw::Dict)::String
    cfg = deepcopy(base_cfg_raw)

    # Sobreescribir sujeto/sesión/tarea
    task = job.condition == "EC" ? "eyesclosed" : "eyesopen"
    cfg["subject"] = Dict{String,Any}(
        "subject_id" => job.bids_id,
        "session_id" => job.bids_session,
        "task"       => task,
        "run"        => 1,
    )

    # Asegura que bids_root y results apuntan al BIDS del proyecto.
    # IMPORTANTE: se usan rutas ABSOLUTAS porque load_ss_config calcula
    # cfg.root = dirname(dirname(config_path)) y los configs temporales
    # viven en config/.batch_tmp/ (un nivel más profundo que config/).
    # En Julia, joinpath(cualquier_root, "/ruta/absoluta") devuelve la ruta
    # absoluta, por lo que los paths absolutos ignorarán el root relativo.
    if !haskey(cfg, "paths")
        cfg["paths"] = Dict{String,Any}()
    end
    cfg["paths"]["bids_root"] = BIDS_DIR    # ruta absoluta
    cfg["paths"]["results"]   = RESULTS     # ruta absoluta

    # Guardar en directorio temporal del proyecto
    tmp_dir = joinpath(PROJ_ROOT, "config", ".batch_tmp")
    mkpath(tmp_dir)
    tmp_path = joinpath(tmp_dir,
        "$(job.bids_id)_$(job.bids_session)_$(task).toml")
    open(tmp_path, "w") do io
        TOML.print(io, cfg)
    end
    return tmp_path
end

# ─── Runner ──────────────────────────────────────────────────

function run_batch(cli::Dict)
    t_start = now()
    println("=" ^ 62)
    println(" NeuroMIND — Fase C: Pipeline en Lote")
    println(" $(t_start)")
    println("=" ^ 62)

    jobs = load_jobs(INVENTORY, cli)
    println("✓ Trabajos cargados del inventario: $(length(jobs))")

    if cli["dry_run"]
        println("\n── Modo DRY-RUN ─────────────────────────────────────────")
        for (i, j) in enumerate(jobs)
            done = already_done(j.bids_id, j.bids_session, j.condition)
            status = done ? "DONE" : "PENDIENTE"
            println("  [$i] sub-$(j.bids_id) ses-$(j.bids_session) $(j.condition) [$status]")
        end
        println("Total: $(length(jobs)) trabajos")
        return
    end

    # Cargar config base
    isfile(BASE_CFG) || error("Config base no encontrada: $BASE_CFG\n→ Crea config/batch_pipeline.toml")
    base_cfg_raw = TOML.parsefile(BASE_CFG)

    # Log
    mkpath(LOGS_DIR)
    log_fname = "batch_run_$(Dates.format(t_start, "yyyy-mm-dd_HH-MM")).csv"
    log_path  = joinpath(LOGS_DIR, log_fname)
    log_io    = open(log_path, "w")
    println(log_io, "subject_id,bids_id,group,session,condition,status,duration_s,error,timestamp")

    n_done = 0; n_skip = 0; n_err = 0
    max_n  = min(length(jobs), cli["max_subjects"])

    for (idx, job) in enumerate(jobs[1:max_n])
        tag = "sub-$(job.bids_id) ses-$(job.bids_session) $(job.condition)"
        println("\n─── [$idx/$max_n] $tag ─────────────────────────────")

        # Skip si ya está procesado
        if cli["skip_done"] && already_done(job.bids_id, job.bids_session, job.condition)
            println("  ⏭ Ya procesado — saltando")
            n_skip += 1
            println(log_io, "$(job.subject_id),$(job.bids_id),$(job.group),$(job.bids_session),$(job.condition),SKIPPED,0,,$(now())")
            continue
        end

        t_job = now()
        try
            tmp_cfg = write_temp_config(job, base_cfg_raw)
            run_single_subject_pipeline(tmp_cfg)
            dur = round(Dates.value(now() - t_job) / 1000, digits=1)
            println("  ✓ Completado en $(dur) s")
            n_done += 1
            println(log_io, "$(job.subject_id),$(job.bids_id),$(job.group),$(job.bids_session),$(job.condition),OK,$(dur),,$(now())")
        catch e
            dur = round(Dates.value(now() - t_job) / 1000, digits=1)
            errmsg = replace(string(e), "," => ";", "\n" => " ")
            @warn "Error en $tag: $e"
            n_err += 1
            println(log_io, "$(job.subject_id),$(job.bids_id),$(job.group),$(job.bids_session),$(job.condition),ERROR,$(dur),\"$(errmsg)\",$(now())")
        end
        flush(log_io)
    end

    close(log_io)

    t_total = round(Dates.value(now() - t_start) / 1000 / 60, digits=1)
    println("\n═══ Resumen batch ═══════════════════════════════════════")
    println("  Completados : $n_done")
    println("  Saltados    : $n_skip")
    println("  Errores     : $n_err")
    println("  Tiempo total: $(t_total) min")
    println("  Log guardado: $log_path")
    println("═════════════════════════════════════════════════════════")

    if n_done + n_skip == max_n && n_err == 0
        println("\n✅ Fase C completada")
        println("\nPróximos pasos:")
        println("  julia --project=. scripts/run_transversal_analysis.jl")
        println("  julia --project=. scripts/run_longitudinal_analysis.jl")
    end
end

# ─── Punto de entrada ─────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    cli = parse_cli_args()
    run_batch(cli)
end
