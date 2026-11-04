#!/usr/bin/env julia

# Construye contenido de GitHub Pages desde exportaciones HTML de Pluto.
#
# Flujo:
#   - Fuente cientifica: Pluto/<Modulo>/<Notebook>.jl
#   - Staging HTML manual: exports/Pluto/<Modulo>/<Notebook>.html
#   - Publicacion web: docs/Pluto/<Modulo>/<Notebook>.html
#
# Uso:
#   julia tools/pages/build.jl
#   julia tools/pages/build.jl /ruta/custom/exports/Pluto

using Dates

const MODULES = [
    (
        name = "Setup",
        source = "Pluto/Setup/Project_Setup.jl",
        target = "docs/Pluto/Setup/Project_Setup.html",
        description = "Configuracion general y contexto del proyecto",
    ),
    (
        name = "BIDS",
        source = "Pluto/BIDS/BIDS.jl",
        target = "docs/Pluto/BIDS/BIDS.html",
        description = "Estandarizacion BIDS para datasets EEG",
    ),
    (
        name = "Preprocessing",
        source = "Pluto/Preprocessing/Preprocessing.jl",
        target = "docs/Pluto/Preprocessing/Preprocessing.html",
        description = "Limpieza y preparacion de senales EEG",
    ),
    (
        name = "ICA",
        source = "Pluto/ICA/ICA.jl",
        target = "docs/Pluto/ICA/ICA.html",
        description = "Separacion de componentes independientes para limpieza EEG",
    ),
    (
        name = "Processing",
        source = "Pluto/Processing/Processing.jl",
        target = "docs/Pluto/Processing/Processing.html",
        description = "Extraccion y analisis de caracteristicas",
    ),
    (
        name = "Spectral",
        source = "Pluto/Spectral/Spectral.jl",
        target = "docs/Pluto/Spectral/Spectral.html",
        description = "Analisis espectral y dinamica en frecuencia para EEG",
    ),
    (
        name = "Connectivity",
        source = "Pluto/Connectivity/Connectivity.jl",
        target = "docs/Pluto/Connectivity/Connectivity.html",
        description = "Analisis de conectividad funcional en EEG",
    ),
    (
        name = "Surrogate",
        source = "Pluto/Surrogate/Surrogate.jl",
        target = "docs/Pluto/Surrogate/Surrogate.html",
        description = "Modelo surrogate para significacion estadistica de conectividad",
    ),
]

function project_root()::String
    return normpath(joinpath(@__DIR__, "..", ".."))
end

function write_text(path::String, content::String)::Nothing
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, content)
    end
    return nothing
end

function json_escape(value::String)::String
    escaped = replace(value, "\\" => "\\\\")
    escaped = replace(escaped, "\"" => "\\\"")
    escaped = replace(escaped, "\n" => "\\n")
    escaped = replace(escaped, "\r" => "\\r")
    escaped = replace(escaped, "\t" => "\\t")
    return escaped
end

function html_name_from_source(source_path::String)::String
    return replace(basename(source_path), ".jl" => ".html")
end

function placeholder_html(module_name::String, source_path::String)::String
    generated_at = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
    return """
<!DOCTYPE html>
<html lang="es">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>NeuroSmart-EEG | $(module_name)</title>
  <style>
    :root {
      --bg: #f4f7fb;
      --surface: #ffffff;
      --text: #1f2937;
      --muted: #4b5563;
      --primary: #1d4ed8;
      --primary-hover: #1e40af;
      --border: #dbe4f0;
      --shadow: 0 10px 24px rgba(16, 24, 40, 0.08);
      --radius: 14px;
      --max-width: 860px;
    }

    * { box-sizing: border-box; }
    body {
      margin: 0;
      padding: 0;
      background: var(--bg);
      color: var(--text);
      font-family: "Inter", "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      line-height: 1.6;
    }
    main {
      min-height: 100dvh;
      display: grid;
      place-items: center;
      padding: 2rem 1rem;
    }
    .card {
      width: min(100%, var(--max-width));
      background: var(--surface);
      border: 1px solid var(--border);
      border-radius: var(--radius);
      box-shadow: var(--shadow);
      padding: 1.6rem;
    }
    h1 {
      margin-top: 0;
      margin-bottom: 0.6rem;
      font-size: clamp(1.5rem, 3vw, 2rem);
      color: #0f172a;
    }
    p {
      margin: 0.5rem 0;
      color: var(--muted);
    }
    .badge {
      display: inline-block;
      margin-bottom: 0.7rem;
      padding: 0.28rem 0.7rem;
      border-radius: 999px;
      background: #eff6ff;
      border: 1px solid #cfe0ff;
      color: #1e3a8a;
      font-size: 0.86rem;
      font-weight: 600;
    }
    .source {
      margin-top: 0.85rem;
      padding: 0.75rem;
      background: #f8fafc;
      border: 1px solid var(--border);
      border-radius: 10px;
      font-family: "SFMono-Regular", Menlo, Consolas, monospace;
      font-size: 0.88rem;
      color: #334155;
      overflow-x: auto;
    }
    .actions {
      margin-top: 1.2rem;
      display: flex;
      gap: 0.7rem;
      flex-wrap: wrap;
    }
    .btn {
      display: inline-block;
      padding: 0.55rem 0.9rem;
      border-radius: 10px;
      border: 1px solid transparent;
      background: var(--primary);
      color: #fff;
      text-decoration: none;
      font-weight: 600;
      font-size: 0.92rem;
    }
    .btn:hover {
      background: var(--primary-hover);
      text-decoration: none;
    }
    footer {
      margin-top: 1rem;
      color: #64748b;
      font-size: 0.86rem;
    }
  </style>
</head>
<body>
  <main>
    <article class="card">
      <span class="badge">NeuroSmart-EEG</span>
      <h1>Modulo $(module_name)</h1>
      <p><strong>Aun no hay exportacion HTML disponible para este modulo.</strong></p>
      <p>
        Exporta el notebook desde Pluto.jl y coloca el archivo en:
      </p>
      <div class="source">exports/Pluto/$(module_name)/$(html_name_from_source(source_path))</div>
      <p>Notebook fuente:</p>
      <div class="source">$(source_path)</div>
      <div class="actions">
        <a class="btn" href="../../index.html">Volver a la portada</a>
      </div>
      <footer>Placeholder generado automaticamente el $(generated_at).</footer>
    </article>
  </main>
</body>
</html>
"""
end

function incomplete_export_html(module_name::String, source_path::String, missing_assets::Vector{String})::String
    generated_at = Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS")
    missing_list = join(["<li><code>$(asset)</code></li>" for asset in missing_assets], "\n")
    return """
<!DOCTYPE html>
<html lang="es">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>NeuroSmart-EEG | $(module_name)</title>
  <style>
    body { font-family: "Inter", "Segoe UI", Roboto, Helvetica, Arial, sans-serif; margin: 0; background: #f4f7fb; color: #1f2937; }
    main { min-height: 100dvh; display: grid; place-items: center; padding: 2rem 1rem; }
    .card { width: min(100%, 860px); background: #fff; border: 1px solid #dbe4f0; border-radius: 14px; padding: 1.4rem; }
    h1 { margin: 0 0 .7rem 0; }
    code { background: #f8fafc; padding: 0.1rem 0.35rem; border-radius: 6px; }
    .btn { display: inline-block; margin-top: 1rem; padding: .55rem .9rem; border-radius: 10px; background: #1d4ed8; color: #fff; text-decoration: none; font-weight: 600; }
    ul { margin-top: .5rem; }
    .source { margin-top: .8rem; color: #475569; }
    footer { margin-top: 1rem; color: #64748b; font-size: .86rem; }
  </style>
</head>
<body>
  <main>
    <article class="card">
      <h1>Modulo $(module_name)</h1>
      <p><strong>La exportacion HTML existe, pero esta incompleta.</strong></p>
      <p>Faltan recursos locales requeridos por el HTML exportado del modulo:</p>
      <ul>
        $(missing_list)
      </ul>
      <p class="source">Notebook fuente: <code>$(source_path)</code></p>
      <p class="source">Reexporta desde Pluto y copia la carpeta completa en <code>exports/Pluto/$(module_name)/</code>.</p>
      <a class="btn" href="../../index.html">Volver a la portada</a>
      <footer>Diagnostico generado automaticamente el $(generated_at).</footer>
    </article>
  </main>
</body>
</html>
"""
end

function default_exports_dir(root::String)::String
    return joinpath(root, "exports", "Pluto")
end

function staged_export_path(exports_pluto_dir::String, module_name::String, source_path::String)::String
    return joinpath(exports_pluto_dir, module_name, html_name_from_source(source_path))
end

function staged_export_dir(exports_pluto_dir::String, module_name::String)::String
    return joinpath(exports_pluto_dir, module_name)
end

function relative_to_root(path::String, root::String)::String
    rel = relpath(path, root)
    return replace(rel, "\\" => "/")
end

function build_manifest_entry(name::String, path::String, source::String, published::Bool)::String
    published_text = published ? "true" : "false"
    return string(
        "  {\n",
        "    \"name\": \"", json_escape(name), "\",\n",
        "    \"path\": \"", json_escape(path), "\",\n",
        "    \"source\": \"", json_escape(source), "\",\n",
        "    \"published\": ", published_text, "\n",
        "  }",
    )
end

function copy_file(src::String, dst::String)::Nothing
    mkpath(dirname(dst))
    cp(src, dst; force = true)
    return nothing
end

function copy_dir_replace(src_dir::String, dst_dir::String)::Nothing
    rm(dst_dir; recursive = true, force = true)
    mkpath(dirname(dst_dir))
    cp(src_dir, dst_dir; force = true)
    return nothing
end

function is_local_asset_ref(ref::AbstractString)::Bool
    value = strip(ref)
    isempty(value) && return false

    lower = lowercase(value)
    startswith(lower, "http://") && return false
    startswith(lower, "https://") && return false
    startswith(lower, "data:") && return false
    startswith(lower, "mailto:") && return false
    startswith(lower, "javascript:") && return false
    startswith(value, "//") && return false
    startswith(value, "#") && return false
    startswith(value, "/") && return false

    return true
end

function normalize_local_ref(ref::AbstractString)::String
    value = strip(ref)
    value = split(value, '#'; limit = 2)[1]
    value = split(value, '?'; limit = 2)[1]
    return value
end

function missing_local_assets(index_html::String, base_dir::String)::Vector{String}
    # Validamos solo referencias reales en atributos HTML (src/href),
    # para evitar falsos positivos de texto embebido (ej. importmap JSON).
    attr_matches = eachmatch(r"""(?i)\b(?:src|href)\s*=\s*["']([^"']+)["']""", index_html)
    refs = Set{String}()

    for m in attr_matches
        raw_ref = m.captures[1]
        is_local_asset_ref(raw_ref) || continue
        ref = normalize_local_ref(raw_ref)
        isempty(ref) && continue
        push!(refs, ref)
    end

    missing = String[]
    for ref in sort(collect(refs))
        if !isfile(joinpath(base_dir, ref))
            push!(missing, ref)
        end
    end
    return missing
end

function main()
    root = project_root()
    docs_dir = joinpath(root, "docs")
    exports_pluto_dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : default_exports_dir(root)

    if !isdir(exports_pluto_dir)
        println("Aviso: no existe directorio de staging HTML: ", exports_pluto_dir)
        println("Se generaran placeholders para modulos sin exportacion.")
    end

    println("== NeuroSmart-EEG Pages Build ==")
    println("Raiz del proyecto: ", root)
    println("Staging de exportaciones Pluto: ", exports_pluto_dir)
    println("Destino GitHub Pages: ", joinpath(docs_dir, "Pluto"))

    manifest_entries = String[]
    real_modules = String[]
    placeholder_modules = String[]

    for mod in MODULES
        target_html_abs = joinpath(root, mod.target)
        staged_html_abs = staged_export_path(exports_pluto_dir, mod.name, mod.source)
        staged_dir_abs = staged_export_dir(exports_pluto_dir, mod.name)
        target_dir_abs = dirname(target_html_abs)

        published_real = false
        if isfile(staged_html_abs)
            staged_html_content = read(staged_html_abs, String)
            missing_assets = missing_local_assets(staged_html_content, staged_dir_abs)

            if isempty(missing_assets)
                copy_dir_replace(staged_dir_abs, target_dir_abs)
                # Preferimos el HTML nominal del notebook y limpiamos index legacy.
                legacy_index = joinpath(target_dir_abs, "index.html")
                if isfile(legacy_index)
                    rm(legacy_index; force=true)
                end
                published_real = true
            else
                println("Aviso: exportacion incompleta para ", mod.name, " (faltan assets).")
                println("       Recursos faltantes: ", join(missing_assets, ", "))
                html = incomplete_export_html(mod.name, mod.source, missing_assets)
                write_text(target_html_abs, html)
            end
        else
            html = placeholder_html(mod.name, mod.source)
            write_text(target_html_abs, html)
        end

        rel_published_path = relative_to_root(target_html_abs, docs_dir)
        push!(
            manifest_entries,
            build_manifest_entry(mod.name, rel_published_path, mod.source, true),
        )

        if published_real
            push!(real_modules, mod.name)
        else
            push!(placeholder_modules, mod.name)
        end
    end

    manifest_path = joinpath(docs_dir, "pages_manifest.json")
    manifest_content = "[\n" * join(manifest_entries, ",\n") * "\n]\n"
    write_text(manifest_path, manifest_content)

    println("\nResumen:")
    println(" - Manifest: ", manifest_path)
    println(" - Modulos con HTML real: ", isempty(real_modules) ? "-" : join(real_modules, ", "))
    println(" - Modulos en placeholder: ", isempty(placeholder_modules) ? "-" : join(placeholder_modules, ", "))
end

main()
