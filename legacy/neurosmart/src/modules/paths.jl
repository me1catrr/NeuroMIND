"""
Utilidades de rutas del proyecto.

Objetivo: evitar strings hardcodeadas dispersas y construir rutas de forma robusta
con `joinpath` y una única noción de "raíz del proyecto".
"""

"""
    project_root() -> String

Devuelve la raíz del proyecto (carpeta que contiene `src/`, `data/`, `results/`, etc.).
"""
project_root() = normpath(joinpath(@__DIR__, "..", ".."))

"""`data_root(cfg=nothing)` devuelve la carpeta base de `data/`."""
function data_root(cfg=nothing)
    if cfg !== nothing && hasproperty(cfg, :data_dir)
        return normpath(getproperty(cfg, :data_dir))
    end
    return joinpath(project_root(), "data")
end

"""`results_root(cfg=nothing)` devuelve la carpeta base de `results/`."""
function results_root(cfg=nothing)
    if cfg !== nothing && hasproperty(cfg, :output_dir)
        return normpath(getproperty(cfg, :output_dir))
    end
    return joinpath(project_root(), "results")
end

"""
    stage_dir(stage::Symbol; kind=:data, cfg=nothing) -> String

Ruta a la carpeta de una etapa, bajo `data/` o bajo `results/`.

Ejemplos:
- `stage_dir(:IO)` → `data/Preprocessing/IO`
- `stage_dir(:FFT; kind=:figures)` → `results/Processing/figures/FFT`
"""
const STAGE_PHASE_MAP = Dict{Symbol, Symbol}(
    :build_participants => :BIDS,
    :build_eeg_bids => :BIDS,
    :validate_bids => :BIDS,
    :IO => :Preprocessing,
    :filtering => :Preprocessing,
    :ICA => :Processing,
    :ICA_cleaning => :Processing,
    :segmentation => :Processing,
    :baseline => :Processing,
    :artifact_rejection => :Processing,
    :baseline_2st => :Processing,
    :FFT => :Processing,
    :CSD => :Connectivity,
    :wPLI => :Connectivity,
)

function stage_dir(stage::Symbol; kind::Symbol = :data, cfg=nothing, phase::Union{Nothing,Symbol}=nothing)
    st = String(stage)
    if kind === :data
        ph = isnothing(phase) ? get(STAGE_PHASE_MAP, stage, nothing) : phase
        return isnothing(ph) ? joinpath(data_root(cfg), st) : joinpath(data_root(cfg), String(ph), st)
    elseif kind === :results
        ph = isnothing(phase) ? get(STAGE_PHASE_MAP, stage, nothing) : phase
        return isnothing(ph) ? joinpath(results_root(cfg), st) : joinpath(results_root(cfg), String(ph), st)
    elseif kind === :figures
        ph = isnothing(phase) ? get(STAGE_PHASE_MAP, stage, nothing) : phase
        return isnothing(ph) ? joinpath(results_root(cfg), "figures", st) : joinpath(results_root(cfg), String(ph), "figures", st)
    elseif kind === :tables
        ph = isnothing(phase) ? get(STAGE_PHASE_MAP, stage, nothing) : phase
        return isnothing(ph) ? joinpath(results_root(cfg), "tables", st) : joinpath(results_root(cfg), String(ph), "tables", st)
    elseif kind === :logs
        ph = isnothing(phase) ? get(STAGE_PHASE_MAP, stage, nothing) : phase
        return isnothing(ph) ? joinpath(results_root(cfg), "logs", st) : joinpath(results_root(cfg), String(ph), "logs", st)
    else
        throw(ArgumentError("kind inválido: $kind (usa :data, :results, :figures, :tables, :logs)"))
    end
end

"""
Ruta base de BIDS.

Si `data_root(cfg)` ya termina en `BIDS`, se usa tal cual.
Si no, se asume estructura `data/BIDS`.
"""
function bids_root(cfg=nothing)
    base = normpath(data_root(cfg))
    return basename(base) == "BIDS" ? base : joinpath(base, "BIDS")
end

"""
Ruta de subcarpeta BIDS con compatibilidad legacy.

Prioriza la ruta nueva `data/BIDS/<subdir>`. Si no existe, intenta
`data/<subdir>` para mantener compatibilidad con estructuras anteriores.
"""
function bids_subdir(cfg, subdir::String)
    primary = joinpath(bids_root(cfg), subdir)
    fallback = joinpath(data_root(cfg), subdir)
    return isdir(primary) ? primary : (isdir(fallback) ? fallback : primary)
end

"""Ruta a carpeta raw (nueva: `data/BIDS/raw`)."""
raw_dir(cfg=nothing) = bids_subdir(cfg, "raw")

"""Ruta a carpeta electrodes (nueva: `data/BIDS/electrodes`)."""
electrodes_dir(cfg=nothing) = bids_subdir(cfg, "electrodes")

