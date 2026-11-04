#!/usr/bin/env julia

# Publica GitHub Pages desde main:/docs con un flujo simple y predecible.
# Uso:
#   julia tools/pages/publish.jl
#   julia tools/pages/publish.jl "Mensaje personalizado"

const DEFAULT_COMMIT_MESSAGE = "Update GitHub Pages (auto build)"

function project_root()::String
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function run_capture(cmd::Cmd; cwd::String = project_root())
    cmd_in_dir = Cmd(cmd; dir = cwd)
    output_buffer = IOBuffer()
    process = run(
        pipeline(ignorestatus(cmd_in_dir), stdout = output_buffer, stderr = output_buffer),
        wait = true,
    )
    output = String(take!(output_buffer))
    return success(process), output
end

function run_or_error(step_name::String, cmd::Cmd; cwd::String = project_root())
    println("→ ", step_name)
    ok, output = run_capture(cmd; cwd = cwd)
    if !isempty(strip(output))
        println(output)
    end
    if !ok
        error("Fallo en paso: $(step_name)")
    end
end

function current_branch()::String
    ok, output = run_capture(`git rev-parse --abbrev-ref HEAD`)
    if !ok
        error("No se pudo detectar la rama actual.")
    end
    return strip(output)
end

function push_with_upstream_fallback(branch::String)
    println("→ Haciendo push")
    ok, output = run_capture(`git push`)
    if !isempty(strip(output))
        println(output)
    end
    if ok
        return
    end

    # Fallback unico: rama sin upstream configurado.
    if occursin("has no upstream branch", output)
        println("⚠️  No hay upstream configurado. Reintentando con --set-upstream.")
        run_or_error("Configurando upstream y haciendo push", `git push --set-upstream origin $branch`)
        return
    end

    error("Fallo en paso: Haciendo push.")
end

function git_status_short()::String
    ok, output = run_capture(`git status --short`)
    if !ok
        error("No se pudo consultar git status --short.")
    end
    return chomp(output)
end

function has_staged_changes()::Bool
    ok, _ = run_capture(`git diff --cached --quiet`)
    if !ok
        # Exit code 1: hay cambios staged. Exit code 0: no hay cambios staged.
        return true
    end
    return false
end

function main()
    root = project_root()
    commit_message = isempty(ARGS) ? DEFAULT_COMMIT_MESSAGE : join(ARGS, " ")
    branch = current_branch()

    println("==========================================")
    println(" NeuroSmart-EEG · Publish GitHub Pages")
    println("==========================================")
    println("Repositorio: ", root)
    println("Rama actual: ", branch)
    if branch != "main"
        println("⚠️  Advertencia: la rama recomendada para publicar es 'main'.")
    end
    println("Mensaje de commit: \"", commit_message, "\"")
    println("")

    # 1) Ejecutar build
    build_script = joinpath(root, "tools", "pages", "build.jl")
    if !isfile(build_script)
        error("No se encontro tools/pages/build.jl")
    end
    run_or_error("Construyendo Pages", `julia $build_script`)

    # 2) Mostrar git status --short y 3) Comprobar si hay cambios
    println("→ git status --short")
    status_output = git_status_short()
    if isempty(strip(status_output))
        println("(sin cambios)")
        println("No hay cambios para publicar. Fin del proceso.")
        return
    end
    println(status_output)
    println("")

    # 4) Si hay cambios: add, commit, push
    run_or_error("Agregando cambios (git add .)", `git add .`)

    if !has_staged_changes()
        println("No hay cambios para publicar. Fin del proceso.")
        return
    end
    run_or_error("Creando commit", `git commit -m $commit_message`)

    # 5) Push con fallback de upstream.
    push_with_upstream_fallback(branch)

    println("✅ Publicacion completada.")
end

try
    main()
catch err
    println("❌ Error durante publish: ", err)
    println("Proceso detenido. Revisa el mensaje anterior para corregir el problema.")
    rethrow(err)
end
