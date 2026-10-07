# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Pipeline en lote (dataset completo)
# ═══════════════════════════════════════════════════════════════
#
#  Fase C — Ejecuta run_single_subject_pipeline() para cada
#  grabación del inventory.csv (MINDEM-IMIBIC).
#
#  Por cada trabajo crea un TOML temporal: copia config/pipeline.toml
#  y sobreescribe [subject], [paths] y [surrogates].enabled
#  (write_temp_config). Surrogates quedan OFF salvo --with-surrogates.
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/run_batch_pipeline.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 25-07-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/<este-script>.jl …
#  Sin shebang: #!/usr/bin/env julia no activaría --project=.
#
#  Prerrequisito
#  ─────────────
#    data/full_data/inventory.csv  (audit_full_dataset.jl)
#    data/bids/                    (build_bids_full.jl)
#
#  Config base
#  ───────────
#    config/pipeline.toml
#
#  Uso
#  ───
#    julia --project=. scripts/run_batch_pipeline.jl
#    julia --project=. scripts/run_batch_pipeline.jl --condition EC
#    julia --project=. scripts/run_batch_pipeline.jl --group MS --session T1
#    julia --project=. scripts/run_batch_pipeline.jl --subjects M11,M12,MC1
#    julia --project=. scripts/run_batch_pipeline.jl --dry-run --skip-done
#    julia --project=. scripts/run_batch_pipeline.jl --skip-done
#
#  Opciones (combinables)
#  ──────────────────────
#    --condition       EC | EO | ALL          (defecto: ALL)
#    --group           MS | HC | ALL          (defecto: ALL)
#    --session         T1 | T2 | ALL          (defecto: ALL)
#    --subjects        id1,id2,…  (subject_id o bids_id: M5 o M05)
#    --skip-done       omitir si existe tables/overview.csv
#    --dry-run         listar sin ejecutar (+ tabla plan)
#    --max-subjects N  limitar a los primeros N jobs (grabaciones)
#    --with-surrogates activar surrogates (por defecto OFF en lote)
#    --verbose         mostrar tablas del pipeline por sujeto
#                      (por defecto: silencioso; solo plan + dashboard)
#
#  Salida
#  ──────
#    results/subjects/sub-{id}/ses-{sess}/{task}/   (pipeline_log.txt por sujeto)
#    results/logs/batch_run_YYYY-MM-DD_HH-MM.csv
#
#  Nota: en modo lote (sin --verbose) se suprime el stdout del pipeline
#  individual. El detalle sigue en cada pipeline_log.txt.

using Dates, TOML, Logging

const PROJ_ROOT  = dirname(@__DIR__)
const BIDS_DIR   = joinpath(PROJ_ROOT, "data", "bids")
const RESULTS    = joinpath(PROJ_ROOT, "results")
const INVENTORY  = joinpath(PROJ_ROOT, "data", "full_data", "inventory.csv")
const BASE_CFG   = joinpath(PROJ_ROOT, "config", "pipeline.toml")
const LOGS_DIR   = joinpath(RESULTS,   "logs")

# ─── Carga de NeuroMIND ───────────────────────────────────────

# Asegurar que el módulo está en el load path
push!(LOAD_PATH, PROJ_ROOT)
include(joinpath(PROJ_ROOT, "src", "NeuroMIND.jl"))
using .NeuroMIND

# ─── Parseo de argumentos ─────────────────────────────────────

function parse_cli_args()
    args = Dict{String,Any}(
        "condition"       => "ALL",
        "group"           => "ALL",
        "session"         => "ALL",
        "subjects"        => String[],
        "skip_done"       => false,
        "dry_run"         => false,
        "max_subjects"    => typemax(Int),
        "with_surrogates" => false,
        "verbose"         => false,
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
        elseif arg == "--with-surrogates"
            args["with_surrogates"] = true; i += 1
        elseif arg == "--verbose"
            args["verbose"] = true; i += 1
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
        # Aceptar subject_id (M5) o bids_id (M05)
        if !isempty(filter_subjs) && subject_id ∉ filter_subjs && bids_id ∉ filter_subjs
            continue
        end

        push!(jobs, BatchJob(subject_id, bids_id, group,
                             session, bids_session, condition, filepath))
    end

    return jobs
end

# ─── Check si ya está procesado ──────────────────────────────

function already_done(bids_id::String, sess::String, cond::String)::Bool
    task = cond == "EC" ? "eyesclosed" : "eyesopen"
    overview = joinpath(RESULTS, "subjects",
                        "sub-$(bids_id)", "ses-$(sess)", task, "tables", "overview.csv")
    isfile(overview)
end

# ─── Helpers de progreso / tablas ─────────────────────────────

"""
    _format_eta_s(seconds) -> String

Formatea segundos restantes como `Xs`, `Xm Ys` o `Xh Ym`.
"""
function _format_eta_s(seconds::Real)::String
    s = max(0.0, Float64(seconds))
    if s < 60
        return string(round(Int, s), "s")
    elseif s < 3600
        m = floor(Int, s / 60)
        r = round(Int, s - 60 * m)
        return "$(m)m $(r)s"
    else
        h = floor(Int, s / 3600)
        m = floor(Int, (s - 3600 * h) / 60)
        return "$(h)h $(m)m"
    end
end

"""
    _progress_bar(done, total; width=28) -> String

Barra Unicode `[████░░░░]`.
"""
function _progress_bar(done::Int, total::Int; width::Int=28)::String
    pct = total > 0 ? done / total : 1.0
    filled = clamp(round(Int, pct * width), 0, width)
    return "[" * "█"^filled * "░"^(width - filled) * "]"
end

"""
    job_plan_stats(jobs) -> NamedTuple

Contadores DONE/PEND por grupo, condición y sesión.
"""
function job_plan_stats(jobs::Vector{BatchJob})
    groups = ["MS", "HC"]
    by_group = Dict(g => (total=0, done=0, pend=0) for g in groups)
    by_cond  = Dict("EC" => 0, "EO" => 0)
    by_sess  = Dict("T1" => 0, "T2" => 0)
    done_tags = String[]

    for j in jobs
        is_done = already_done(j.bids_id, j.bids_session, j.condition)
        if !haskey(by_group, j.group)
            by_group[j.group] = (total=0, done=0, pend=0)
        end
        st = by_group[j.group]
        by_group[j.group] = (
            total = st.total + 1,
            done  = st.done + (is_done ? 1 : 0),
            pend  = st.pend + (is_done ? 0 : 1),
        )
        by_cond[j.condition] = get(by_cond, j.condition, 0) + 1
        by_sess[j.bids_session] = get(by_sess, j.bids_session, 0) + 1
        if is_done
            push!(done_tags, "sub-$(j.bids_id) $(j.bids_session) $(j.condition)")
        end
    end

    n_total = length(jobs)
    n_done  = count(j -> already_done(j.bids_id, j.bids_session, j.condition), jobs)
    n_pend  = n_total - n_done
    return (; by_group, by_cond, by_sess, n_total, n_done, n_pend, done_tags)
end

"""
    print_job_plan(jobs; title="Plan de trabajos")

Tabla pre-run / dry-run: MS/HC × DONE/PEND + condición/sesión.
"""
function print_job_plan(jobs::Vector{BatchJob}; title::String="Plan de trabajos")
    st = job_plan_stats(jobs)
    println("\n── $(title) ──────────────────────────────────")
    println("  Grupo   Total  DONE  PEND")
    for g in ("MS", "HC")
        if haskey(st.by_group, g)
            s = st.by_group[g]
            println("  $(rpad(g, 6))  $(lpad(string(s.total), 5))  $(lpad(string(s.done), 4))  $(lpad(string(s.pend), 4))")
        end
    end
    # Grupos inesperados
    for g in sort(collect(keys(st.by_group)))
        g ∈ ("MS", "HC") && continue
        s = st.by_group[g]
        println("  $(rpad(g, 6))  $(lpad(string(s.total), 5))  $(lpad(string(s.done), 4))  $(lpad(string(s.pend), 4))")
    end
    println("  ─────────────────────────")
    note = if st.n_done > 0 && length(st.done_tags) <= 3
        "  (" * join(st.done_tags, ", ") * " ya DONE)"
    elseif st.n_done > 0
        "  ($(st.n_done) ya DONE)"
    else
        ""
    end
    println("  Total   $(lpad(string(st.n_total), 5))  $(lpad(string(st.n_done), 4))  $(lpad(string(st.n_pend), 4))$(note)")

    cond_parts = ["$(k)=$(st.by_cond[k])" for k in sort(collect(keys(st.by_cond)))]
    sess_parts = ["$(k)=$(st.by_sess[k])" for k in sort(collect(keys(st.by_sess)))]
    println("  Por condición: $(join(cond_parts, "  "))")
    println("  Por sesión:    $(join(sess_parts, "  "))")
    println("──────────────────────────────────────────────────")
    return st
end

"""
    _run_pipeline_quiet(cfg_path)

Ejecuta el pipeline silenciando stdout/stderr y @warn/@info.
El detalle permanece en `pipeline_log.txt` del sujeto.
"""
function _run_pipeline_quiet(cfg_path::String)
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            with_logger(NullLogger()) do
                run_single_subject_pipeline(cfg_path)
            end
        end
    end
    return nothing
end

"""
    print_batch_dashboard(...; dash_lines, current, last_status, in_place)

Dashboard de evolución del lote. En TTY reescribe el bloque in-place
para no inundar la terminal; si no es TTY, imprime líneas nuevas.
"""
function print_batch_dashboard(
    idx::Int, max_n::Int,
    n_ok::Int, n_skip::Int, n_err::Int,
    durations::Vector{Float64},
    group_ok::Dict{String,Int}, group_pend::Dict{String,Int};
    dash_lines::Ref{Int} = Ref(0),
    current::String = "",
    last_status::String = "",
    in_place::Bool = true,
)
    processed = n_ok + n_skip + n_err
    bar = _progress_bar(processed, max_n)
    pct = max_n > 0 ? round(Int, 100 * processed / max_n) : 100
    remaining_slots = max(0, max_n - idx)
    mean_dur = isempty(durations) ? 0.0 : sum(durations) / length(durations)
    eta_s = mean_dur * remaining_slots
    eta_str = isempty(durations) ? "ETA —" : "ETA $(_format_eta_s(eta_s))"
    mean_str = isempty(durations) ? "—" : _format_eta_s(mean_dur)

    ms_done = get(group_ok, "MS", 0); ms_pend = get(group_pend, "MS", 0)
    hc_done = get(group_ok, "HC", 0); hc_pend = get(group_pend, "HC", 0)

    lines = String[
        "── Lote en curso ──────────────────────────────────",
        "  $bar  $processed/$max_n  $(lpad(string(pct), 3))%  · $eta_str  · media $mean_str",
        "  OK=$n_ok  SKIP=$n_skip  ERR=$n_err",
        "  Grupo   done  pend",
        "  MS    $(lpad(string(ms_done), 5))  $(lpad(string(ms_pend), 4))",
        "  HC    $(lpad(string(hc_done), 5))  $(lpad(string(hc_pend), 4))",
        "  Actual: $(isempty(current) ? "—" : current)",
        "  Último: $(isempty(last_status) ? "—" : last_status)",
        "──────────────────────────────────────────────────",
    ]

    use_tty = in_place && isa(stdout, Base.TTY)
    if use_tty && dash_lines[] > 0
        # Subir N líneas y borrar hasta el final de pantalla
        print("\e[$(dash_lines[])A\e[J")
    elseif !use_tty && dash_lines[] > 0
        println()  # separador si no podemos reescribir
    end

    for ln in lines
        println(ln)
    end
    flush(stdout)
    dash_lines[] = length(lines)
    return nothing
end

# ─── Generación de config temporal ───────────────────────────

function write_temp_config(job::BatchJob, base_cfg_raw::Dict;
                           with_surrogates::Bool=false)::String
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

    # Surrogates: OFF por defecto en lote (costoso; verificación solo M05).
    # Opt-in con --with-surrogates.
    if !haskey(cfg, "surrogates")
        cfg["surrogates"] = Dict{String,Any}()
    end
    cfg["surrogates"]["enabled"] = with_surrogates

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
    surr_on = Bool(cli["with_surrogates"])
    verbose = Bool(cli["verbose"])
    println("=" ^ 62)
    println(" NeuroMIND — Fase C: Pipeline en Lote")
    println(" $(t_start)")
    println(" surrogates=$(surr_on ? "ON (--with-surrogates)" : "OFF")")
    println(" salida=$(verbose ? "verbose (tablas por sujeto)" : "silenciosa (solo dashboard; detalle en pipeline_log.txt)")")
    println("=" ^ 62)

    jobs = load_jobs(INVENTORY, cli)
    println("✓ Trabajos cargados del inventario: $(length(jobs))")
    if cli["max_subjects"] < typemax(Int)
        println("  (límite --max-subjects = $(cli["max_subjects"]) jobs)")
    end

    max_n = min(length(jobs), cli["max_subjects"])
    jobs_run = jobs[1:max_n]
    plan = print_job_plan(jobs_run)

    if cli["dry_run"]
        println("\n── Modo DRY-RUN ─────────────────────────────────────────")
        for (i, j) in enumerate(jobs_run)
            done = already_done(j.bids_id, j.bids_session, j.condition)
            status = done ? "DONE" : "PENDIENTE"
            println("  [$i] sub-$(j.bids_id) ses-$(j.bids_session) $(j.condition) $(j.group) [$status]")
        end
        println("Total: $(length(jobs_run)) jobs  ·  surrogates=$(surr_on ? "ON" : "OFF")")
        return
    end

    # Cargar config base
    isfile(BASE_CFG) || error("Config base no encontrada: $BASE_CFG\n→ config/pipeline.toml es la configuración única del proyecto.")
    base_cfg_raw = TOML.parsefile(BASE_CFG)

    # Log
    mkpath(LOGS_DIR)
    log_fname = "batch_run_$(Dates.format(t_start, "yyyy-mm-dd_HH-MM")).csv"
    log_path  = joinpath(LOGS_DIR, log_fname)
    log_io    = open(log_path, "w")
    println(log_io, "subject_id,bids_id,group,session,condition,status,duration_s,error,timestamp")

    n_ok = 0; n_skip = 0; n_err = 0
    durations = Float64[]
    # Desglose vivo: seen = resueltos en este run; pend = aún no tocados
    n_by_group = Dict{String,Int}()
    for j in jobs_run
        n_by_group[j.group] = get(n_by_group, j.group, 0) + 1
    end
    group_ok   = Dict(g => 0 for g in keys(n_by_group))
    group_seen = Dict(g => 0 for g in keys(n_by_group))

    n_pending_at_start = plan.n_pend
    dash_lines = Ref(0)
    last_status = ""
    # Dashboard in-place solo en modo silencioso (en verbose el pipeline
    # empuja el cursor y rompería el redraw).
    in_place = !verbose && isa(stdout, Base.TTY)

    for (idx, job) in enumerate(jobs_run)
        tag = "sub-$(job.bids_id) ses-$(job.bids_session) $(job.condition)"
        tag_g = "$tag ($(job.group))"

        # Skip si ya está procesado
        if cli["skip_done"] && already_done(job.bids_id, job.bids_session, job.condition)
            n_skip += 1
            group_seen[job.group] = get(group_seen, job.group, 0) + 1
            group_ok[job.group]   = get(group_ok, job.group, 0) + 1
            last_status = "$tag · SKIP"
            println(log_io, "$(job.subject_id),$(job.bids_id),$(job.group),$(job.bids_session),$(job.condition),SKIPPED,0,,$(now())")
            flush(log_io)
            group_pend_live = Dict(g => get(n_by_group, g, 0) - get(group_seen, g, 0)
                                   for g in keys(n_by_group))
            # En skips masivos: actualizar dashboard cada 10 o al final del bloque
            # (cada skip) — in-place es barato; siempre actualizar.
            print_batch_dashboard(idx, max_n, n_ok, n_skip, n_err, durations,
                                  group_ok, group_pend_live;
                                  dash_lines=dash_lines, current="—",
                                  last_status=last_status, in_place=in_place)
            continue
        end

        group_pend_live = Dict(g => get(n_by_group, g, 0) - get(group_seen, g, 0)
                               for g in keys(n_by_group))
        print_batch_dashboard(idx - 1, max_n, n_ok, n_skip, n_err, durations,
                              group_ok, group_pend_live;
                              dash_lines=dash_lines,
                              current="$tag_g · ejecutando…",
                              last_status=last_status, in_place=in_place)

        if verbose
            # Salir del bloque dashboard antes del dump del pipeline
            dash_lines[] = 0
            println("\n─── [$idx/$max_n] $tag_g ─────────────────────────────")
        end

        t_job = now()
        try
            tmp_cfg = write_temp_config(job, base_cfg_raw; with_surrogates=surr_on)
            if verbose
                run_single_subject_pipeline(tmp_cfg)
            else
                _run_pipeline_quiet(tmp_cfg)
            end
            dur = round(Dates.value(now() - t_job) / 1000, digits=1)
            n_ok += 1
            push!(durations, Float64(dur))
            last_status = "$tag · OK $(_format_eta_s(dur))"
            println(log_io, "$(job.subject_id),$(job.bids_id),$(job.group),$(job.bids_session),$(job.condition),OK,$(dur),,$(now())")
        catch e
            dur = round(Dates.value(now() - t_job) / 1000, digits=1)
            errmsg = replace(string(e), "," => ";", "\n" => " ")
            n_err += 1
            push!(durations, Float64(dur))
            last_status = "$tag · ERROR $(_format_eta_s(dur))"
            # Error siempre visible (fuera del silencio del pipeline)
            if in_place && dash_lines[] > 0
                println()  # no pisar el dashboard con el aviso
                dash_lines[] = 0
            end
            println("  ✗ ERROR en $tag: $e")
            println(log_io, "$(job.subject_id),$(job.bids_id),$(job.group),$(job.bids_session),$(job.condition),ERROR,$(dur),\"$(errmsg)\",$(now())")
        end
        flush(log_io)

        group_seen[job.group] = get(group_seen, job.group, 0) + 1
        group_ok[job.group]   = get(group_ok, job.group, 0) + 1
        group_pend_live = Dict(g => get(n_by_group, g, 0) - get(group_seen, g, 0)
                               for g in keys(n_by_group))
        print_batch_dashboard(idx, max_n, n_ok, n_skip, n_err, durations,
                              group_ok, group_pend_live;
                              dash_lines=dash_lines, current="—",
                              last_status=last_status, in_place=in_place)
    end

    close(log_io)
    if in_place && dash_lines[] > 0
        println()  # dejar el último dashboard fijo antes del resumen
    end

    t_total_s = Dates.value(now() - t_start) / 1000
    t_total_min = round(t_total_s / 60, digits=1)
    mean_dur = isempty(durations) ? 0.0 : sum(durations) / length(durations)

    # Resumen por grupo / condición desde el CSV log (recontar en memoria)
    # Usamos group_seen y contadores globales; para OK/SKIP/ERROR por grupo
    # re-leemos no es necesario si trackeamos — añadimos contadores dedicados.
    println("\n═══ Resumen batch ═══════════════════════════════════════")
    println("  Pendientes al inicio : $n_pending_at_start")
    println("  Completados (OK)     : $n_ok")
    println("  Saltados (SKIP)      : $n_skip")
    println("  Errores              : $n_err")
    println("  Jobs en esta corrida : $max_n")
    println("  surrogates           : $(surr_on ? "ON" : "OFF")")
    println("  Tiempo total         : $(_format_eta_s(t_total_s)) ($(t_total_min) min)")
    if !isempty(durations)
        println("  Media por job ejec.  : $(_format_eta_s(mean_dur))")
    end
    println("  Log guardado         : $log_path")

    # Tabla por grupo — recontar desde jobs_run + log es costoso; usamos un
    # segundo pase ligero sobre el fichero de log.
    _print_summary_tables(log_path)
    println("═════════════════════════════════════════════════════════")

    if max_n == 0
        println("\n⚠ Ningún job tras filtros — nada que ejecutar")
    elseif n_ok + n_skip == max_n && n_err == 0
        println("\n✅ Fase C completada")
        println("\nPróximos pasos:")
        println("  julia --project=. scripts/run_transversal_analysis.jl")
        println("  julia --project=. scripts/run_longitudinal_analysis.jl")
    elseif n_err > 0
        println("\n⚠ Fase C terminó con $n_err error(es) — revisar log")
    end
end

"""
    _print_summary_tables(log_path)

Lee el CSV de la corrida y muestra tablas por grupo y condición.
"""
function _print_summary_tables(log_path::String)
    # status counts: group => (ok, skip, err), cond => (ok, skip, err)
    gstat = Dict{String, Dict{String,Int}}()
    cstat = Dict{String, Dict{String,Int}}()

    for (i, line) in enumerate(readlines(log_path))
        i == 1 && continue
        isempty(strip(line)) && continue
        parts = split(line, ",")
        length(parts) < 6 && continue
        group = strip(parts[3])
        cond  = strip(parts[5])
        status = strip(parts[6])
        key = status == "OK" ? "ok" : status == "SKIPPED" ? "skip" : "err"

        if !haskey(gstat, group)
            gstat[group] = Dict("ok" => 0, "skip" => 0, "err" => 0)
        end
        if !haskey(cstat, cond)
            cstat[cond] = Dict("ok" => 0, "skip" => 0, "err" => 0)
        end
        gstat[group][key] = get(gstat[group], key, 0) + 1
        cstat[cond][key]  = get(cstat[cond], key, 0) + 1
    end

    println("\n  Grupo  OK  SKIP  ERROR")
    for g in sort(collect(keys(gstat)))
        s = gstat[g]
        println("  $(rpad(g, 5))  $(lpad(string(s["ok"]), 2))  $(lpad(string(s["skip"]), 4))  $(lpad(string(s["err"]), 5))")
    end
    println("\n  Cond   OK  SKIP  ERROR")
    for c in sort(collect(keys(cstat)))
        s = cstat[c]
        println("  $(rpad(c, 5))  $(lpad(string(s["ok"]), 2))  $(lpad(string(s["skip"]), 4))  $(lpad(string(s["err"]), 5))")
    end
end

# ─── Punto de entrada ─────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    cli = parse_cli_args()
    run_batch(cli)
end
