"""
Utilidades comunes del pipeline (IO ligero, logging simple, limpieza de carpetas).

Nota: se prioriza no añadir dependencias nuevas. Todo se basa en stdlib y paquetes ya presentes.
"""

using Dates
using Serialization

"""Asegura que `dir` existe y devuelve `dir`."""
ensure_dir(dir::AbstractString) = (isdir(dir) || mkpath(dir); String(dir))

"""
    clean_dir_by_ext(dir, ext; verbose=true)

Elimina archivos en `dir` cuyo nombre termine en `ext` (p. ej. ".png", ".bin", ".csv", ".log").
"""
function clean_dir_by_ext(dir::AbstractString, ext::AbstractString; verbose::Bool = true)
    if !isdir(dir)
        return 0
    end
    files = filter(f -> endswith(f, ext), readdir(dir))
    for f in files
        rm(joinpath(dir, f))
        verbose && println("  ✓ Eliminado: $f")
    end
    return length(files)
end

"""Timestamp compacto para nombres de ejecución/archivos."""
run_id() = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

"""Escribe un mensaje con hora en un `io` de log."""
logmsg(io, s) = println(io, "[$(Dates.format(now(), "HH:MM:SS"))] $s")

"""
    open_log(dir, prefix) -> (path, io)

Crea `dir` si hace falta y abre un archivo de log en modo escritura.
Devuelve la ruta y el `IO` abierto (responsabilidad del caller cerrarlo).
"""
function open_log(dir::AbstractString, prefix::AbstractString)
    ensure_dir(dir)
    path = joinpath(dir, "$(prefix)_$(run_id()).log")
    io = open(path, "w")
    return path, io
end

"""Serializa `x` en `path` asegurando directorio padre."""
function save_bin(path::AbstractString, x)
    ensure_dir(dirname(path))
    Serialization.serialize(path, x)
    return path
end

"""Deserializa desde `path`."""
load_bin(path::AbstractString) = Serialization.deserialize(path)

"""
    pick_first_existing(paths...) -> Union{String,Nothing}

Devuelve el primer path existente (archivo) o `nothing` si ninguno existe.
"""
function pick_first_existing(paths::AbstractString...)
    for p in paths
        isfile(p) && return String(p)
    end
    return nothing
end

