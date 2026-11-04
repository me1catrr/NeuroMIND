#!/usr/bin/env julia

# Exporta notebooks Pluto a HTML estatico en exports/Pluto/<Modulo>/.
#
# Uso:
#   julia --project=. tools/pages/export_pluto.jl
#   julia --project=. tools/pages/export_pluto.jl ICA Spectral

using Pkg
import PlutoSliderServer

function project_root()::String
    return normpath(joinpath(@__DIR__, "..", ".."))
end

const MODULE_NOTEBOOK = Dict(
    "Setup" => "Pluto/Setup/Project_Setup.jl",
    "BIDS" => "Pluto/BIDS/BIDS.jl",
    "Preprocessing" => "Pluto/Preprocessing/Preprocessing.jl",
    "ICA" => "Pluto/ICA/ICA.jl",
    "Processing" => "Pluto/Processing/Processing.jl",
    "Spectral" => "Pluto/Spectral/Spectral.jl",
    "Connectivity" => "Pluto/Connectivity/Connectivity.jl",
    "Surrogate" => "Pluto/Surrogate/Surrogate.jl",
)

function normalize_modules(args::Vector{String})::Vector{String}
    if isempty(args)
        return collect(keys(MODULE_NOTEBOOK))
    end

    unknown = filter(m -> !haskey(MODULE_NOTEBOOK, m), args)
    if !isempty(unknown)
        error("Modulos no validos: $(join(unknown, ", ")). Usa: $(join(sort(collect(keys(MODULE_NOTEBOOK))), ", "))")
    end
    return args
end

function export_module(root::String, module_name::String)
    notebook_rel = MODULE_NOTEBOOK[module_name]
    notebook_abs = joinpath(root, notebook_rel)
    output_dir = joinpath(root, "exports", "Pluto", module_name)

    if !isfile(notebook_abs)
        error("No existe notebook para $module_name: $notebook_abs")
    end

    mkpath(output_dir)
    PlutoSliderServer.export_notebook(notebook_abs; Export_output_dir=output_dir)

    # Limpia index legacy para mantener solo <Notebook>.html en staging.
    legacy_index = joinpath(output_dir, "index.html")
    if isfile(legacy_index)
        rm(legacy_index; force=true)
    end

    return output_dir
end

function exported_html_name(module_name::String)::String
    return replace(basename(MODULE_NOTEBOOK[module_name]), ".jl" => ".html")
end

function main()
    root = project_root()
    Pkg.activate(root)
    modules = normalize_modules(ARGS)

    println("== Pluto HTML Export ==")
    println("Proyecto: $root")
    println("Modulos: $(join(modules, ", "))")

    for m in modules
        println("\n→ Exportando $m")
        out = export_module(root, m)
        println("  ✓ OK: $out/$(exported_html_name(m))")
    end

    println("\nListo. Siguiente paso:")
    println("  julia --project=. tools/pages/build.jl")
    println("  julia --project=. tools/pages/publish.jl \"mensaje\"")
end

main()
