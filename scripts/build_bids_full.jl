# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Construcción del dataset BIDS-EEG
# ═══════════════════════════════════════════════════════════════
#
#  Fase B — Lee data/tables/inventory.csv (audit_full_dataset.jl) y
#  genera en data/bids/ un dataset BIDS-EEG real, validable con
#  bids-validator, a partir de los tripletes BrainVision de
#  data/sourcedata/brainvision/ (que no se modifican nunca).
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/build_bids_full.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 15-10-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/<este-script>.jl …
#  Sin shebang: #!/usr/bin/env julia no activaría --project=.
#
#  Uso
#  ───
#    julia --project=. scripts/build_bids_full.jl
#    julia --project=. scripts/build_bids_full.jl --subjects M05,M04
#    julia --project=. scripts/build_bids_full.jl --force
#
#  Opciones
#  ────────
#    --inventory PATH   inventario (defecto: data/tables/inventory.csv)
#    --subjects  LIST   solo estos sujetos BIDS (M05,M04,…); no limpia obsoletos
#    --force            recopia también los .eeg ya existentes de igual tamaño
#
#  Qué escribe (todo regenerable; nada se borra)
#  ─────────────────────────────────────────────
#    data/bids/
#    ├── dataset_description.json · README.md · participants.json   (config/bids/)
#    ├── participants.tsv                                           (audit_full_dataset.jl)
#    ├── task-{eyesclosed,eyesopen}_eeg.json · _channels.tsv        (herencia BIDS)
#    └── sub-*/ses-*/eeg/
#        ├── sub-*_ses-*_task-*_run-01_eeg.{vhdr,vmrk,eeg}          señal renombrada
#        ├── sub-*_ses-*_task-*_run-01_eeg.json                     solo registros recortados
#        └── sub-*_ses-*_{electrodes.tsv,coordsystem.json}
#
#  Registros recortados (config/raw_selection.csv, acción `segment`)
#  ─────────────────────────────────────────────────────────────────
#  Un raw con ojos abiertos + cerrados se parte ANTES de cualquier
#  filtro en dos grabaciones independientes con intervalos
#  semiabiertos [tmin, tmax): se copian los bytes de esas muestras
#  (IEEE_FLOAT_32, sin reescalar) y el .vmrk recibe un `New Segment`
#  con la hora de adquisición correspondiente.
#
#  Los archivos que dejan de ser válidos se mueven a
#  data/_archive/bids_descartados/ (nunca se borran).
#
#  Siguiente paso: validar con bids-validator (ver README, §Datos).

using Dates

# ─── Rutas ────────────────────────────────────────────────────

const PROJ_ROOT   = dirname(@__DIR__)
const BIDS_DIR    = joinpath(PROJ_ROOT, "data", "bids")
const INVENTORY   = joinpath(PROJ_ROOT, "data", "tables", "inventory.csv")
const TEMPLATES   = joinpath(PROJ_ROOT, "config", "bids")
const DISCARD_DIR = joinpath(PROJ_ROOT, "data", "_archive", "bids_descartados")

const TASKS = Dict(
    "EC" => ("eyesclosed", "Resting-state EEG with eyes closed."),
    "EO" => ("eyesopen",   "Resting-state EEG with eyes open."),
)

# ─── Inventario ───────────────────────────────────────────────

struct InventoryRow
    filename     :: String
    subject_id   :: String
    bids_id      :: String
    group        :: String
    session      :: String
    bids_session :: String
    condition    :: String
    excluded     :: Bool
    filepath     :: String                    # relativa a la raíz del proyecto
    tmin_s       :: Union{Nothing,Float64}    # ventana de corte (M4/M5/M6 T1)
    tmax_s       :: Union{Nothing,Float64}
end

function load_inventory(path::String)::Vector{InventoryRow}
    isfile(path) || error("Inventario no encontrado: $path\n  → Ejecuta primero: julia --project=. scripts/audit_full_dataset.jl")
    rows = InventoryRow[]
    # Cabecera: filename,subject_id,bids_id,group,session,session_sub,
    # bids_session,condition_raw,condition,date_str,initials,excluded,sexo,
    # edad,nivel_edu,filepath,note,tmin_s,tmax_s
    for line in readlines(path)[2:end]
        isempty(strip(line)) && continue
        p = split(line, ",")
        length(p) >= 19 || error("inventory.csv: se esperaban 19 columnas, hay $(length(p)): $line")
        push!(rows, InventoryRow(
            strip(p[1]), strip(p[2]), strip(p[3]), strip(p[4]), strip(p[5]),
            strip(p[7]), strip(p[9]), strip(p[12]) == "true", strip(p[16]),
            tryparse(Float64, strip(p[18])), tryparse(Float64, strip(p[19])),
        ))
    end
    return rows
end

# ─── BrainVision ──────────────────────────────────────────────

struct VhdrInfo
    n_channels :: Int
    fs         :: Float64
    ch_names   :: Vector{String}
    unit       :: String
    binfmt     :: String
    orient     :: String
    data_file  :: String
    marker_file:: String
end

function read_vhdr(path::String)::VhdrInfo
    n_ch, fs, binfmt, orient, dfile, mfile = 0, 500.0, "", "", "", ""
    names, units = String[], String[]
    for raw in readlines(path)
        line = strip(raw)
        if startswith(line, "NumberOfChannels=");  n_ch = parse(Int, split(line, "=")[2])
        elseif startswith(line, "SamplingInterval="); fs = 1e6 / parse(Float64, split(line, "=")[2])
        elseif startswith(line, "BinaryFormat=");     binfmt = String(split(line, "="; limit=2)[2])
        elseif startswith(line, "DataOrientation=");  orient = String(split(line, "="; limit=2)[2])
        elseif startswith(line, "DataFile=");         dfile  = String(split(line, "="; limit=2)[2])
        elseif startswith(line, "MarkerFile=");       mfile  = String(split(line, "="; limit=2)[2])
        elseif (m = match(r"^Ch(\d+)=([^,]*),[^,]*,[^,]*,(.*)$", line)) !== nothing
            push!(names, String(m[2])); push!(units, String(strip(m[3])))
        end
    end
    n_ch > 0 && length(names) >= n_ch || error("Cabecera sin canales válidos: $path")
    return VhdrInfo(n_ch, fs, names[1:n_ch], units[1], binfmt, orient, dfile, mfile)
end

"Sustituye una clave `Clave=valor` de una cabecera BrainVision, respetando saltos de línea."
_set_key(txt::String, key::String, val::String) =
    replace(txt, Regex("(?m)^" * key * "=[^\\r\\n]*") => key * "=" * val)

# Marcadores: Mk<N>=<Tipo>,<Descripción>,<Posición>,<Tamaño>,<Canal>,<Fecha 20 dígitos>
struct Marker
    type :: String
    desc :: String
    pos  :: Int
    size :: Int
    chan :: Int
    date :: String
end

function parse_markers(txt::String)::Vector{Marker}
    out = Marker[]
    for line in split(txt, r"\r?\n")
        m = match(r"^Mk\d+=([^,]*),([^,]*),(\d+),(\d+),(\d+)(?:,(\d*))?", strip(line))
        m === nothing && continue
        push!(out, Marker(m[1], m[2], parse(Int, m[3]), parse(Int, m[4]),
                          parse(Int, m[5]), m[6] === nothing ? "" : String(m[6])))
    end
    return out
end

"Suma `us` microsegundos a una fecha BrainVision `YYYYMMDDhhmmssuuuuuu`."
function shift_bv_date(date::String, us::Float64)::String
    length(date) == 20 || return date
    dt   = DateTime(date[1:4], "yyyy") + Month(parse(Int, date[5:6]) - 1) + Day(parse(Int, date[7:8]) - 1) +
           Hour(parse(Int, date[9:10])) + Minute(parse(Int, date[11:12])) + Second(parse(Int, date[13:14]))
    tot  = parse(Int, date[15:20]) + round(Int, us)
    dt  += Second(fld(tot, 1_000_000))
    return Dates.format(dt, "yyyymmddHHMMSS") * lpad(mod(tot, 1_000_000), 6, '0')
end

function marker_line(i::Int, m::Marker)
    base = "Mk$(i)=$(m.type),$(m.desc),$(m.pos),$(m.size),$(m.chan)"
    return isempty(m.date) ? base : base * "," * m.date
end

# ─── Escritura de un registro ─────────────────────────────────

"""
    write_recording(row, src_vhdr, dst_dir, base; force) -> Symbol

Escribe el tríptico BIDS de un registro. Devuelve :copied, :kept (el .eeg ya
estaba) o :cropped. La señal original no se toca.
"""
function write_recording(row::InventoryRow, src_vhdr::String, dst_dir::String,
                         base::String; force::Bool=false)::Symbol
    info = read_vhdr(src_vhdr)
    info.binfmt == "IEEE_FLOAT_32" || error("Formato no soportado ($(info.binfmt)): $src_vhdr")
    info.orient == "MULTIPLEXED"   || error("Orientación no soportada ($(info.orient)): $src_vhdr")

    srcdir  = dirname(src_vhdr)
    src_eeg = joinpath(srcdir, info.data_file)
    src_vmrk = joinpath(srcdir, info.marker_file)
    isfile(src_eeg)  || error(".eeg no encontrado: $src_eeg")
    isfile(src_vmrk) || error(".vmrk no encontrado: $src_vmrk")

    dst_eeg  = joinpath(dst_dir, base * "_eeg.eeg")
    dst_vhdr = joinpath(dst_dir, base * "_eeg.vhdr")
    dst_vmrk = joinpath(dst_dir, base * "_eeg.vmrk")

    vhdr_txt = read(src_vhdr, String)
    vhdr_txt = _set_key(vhdr_txt, "DataFile",   base * "_eeg.eeg")
    vhdr_txt = _set_key(vhdr_txt, "MarkerFile", base * "_eeg.vmrk")
    vmrk_txt = read(src_vmrk, String)

    bytes_per_frame = info.n_channels * 4
    total_bytes = filesize(src_eeg)
    total_bytes % bytes_per_frame == 0 ||
        error("Tamaño de .eeg no múltiplo de $(bytes_per_frame) B: $src_eeg")
    n_total = total_bytes ÷ bytes_per_frame

    if row.tmin_s === nothing || row.tmax_s === nothing
        # ── Registro completo: copia byte a byte con enlaces reescritos ──
        status = :kept
        if force || !(isfile(dst_eeg) && filesize(dst_eeg) == total_bytes)
            cp(src_eeg, dst_eeg; force=true)
            status = :copied
        end
        write(dst_vhdr, vhdr_txt)
        write(dst_vmrk, _set_key(vmrk_txt, "DataFile", base * "_eeg.eeg"))
        return status
    end

    # ── Recorte [tmin, tmax) en muestras ────────────────────────────────
    i0 = round(Int, row.tmin_s * info.fs)
    i1 = round(Int, row.tmax_s * info.fs)
    0 <= i0 < i1 <= n_total ||
        error("Ventana fuera de rango en $(row.filename): [$(row.tmin_s), $(row.tmax_s)) s " *
              "→ muestras [$i0, $i1) de $n_total")
    buf = open(src_eeg) do io
        seek(io, i0 * bytes_per_frame)
        read(io, (i1 - i0) * bytes_per_frame)
    end
    length(buf) == (i1 - i0) * bytes_per_frame || error("Lectura truncada en $src_eeg")
    write(dst_eeg, buf)
    write(dst_vhdr, vhdr_txt)

    # .vmrk nuevo: cabecera original + marcadores de la ventana (posiciones 1-based)
    markers = parse_markers(vmrk_txt)
    seg = filter(m -> m.type == "New Segment" && m.pos <= i0 + 1, markers)
    t0_date = isempty(seg) ? "" :
        shift_bv_date(last(seg).date, (i0 + 1 - last(seg).pos) / info.fs * 1e6)
    inwin = [Marker(m.type, m.desc, m.pos - i0, m.size, m.chan, m.date)
             for m in markers if m.type != "New Segment" && i0 < m.pos <= i1]
    head = String[]
    for line in split(vmrk_txt, r"\r?\n")
        startswith(strip(line), "Mk") && break
        push!(head, startswith(line, "DataFile=") ? "DataFile=" * base * "_eeg.eeg" : String(line))
    end
    out = vcat(head,
               [marker_line(1, Marker("New Segment", "", 1, 1, 0, t0_date))],
               [marker_line(i + 1, m) for (i, m) in enumerate(inwin)])
    write(dst_vmrk, join(out, "\n") * "\n")
    return :cropped
end

# ─── Ficheros de dataset y de sesión ──────────────────────────

_json_s(s) = "\"" * replace(s, "\\" => "\\\\", "\"" => "\\\"") * "\""

function write_dataset_files(ch_names::Vector{String}, unit::String, fs::Float64)
    mkpath(BIDS_DIR)
    for f in ("dataset_description.json", "README.md", "participants.json")
        cp(joinpath(TEMPLATES, f), joinpath(BIDS_DIR, f); force=true)
    end
    tpl = read(joinpath(TEMPLATES, "task_eeg.template.json"), String)
    for (task, desc) in values(TASKS)
        write(joinpath(BIDS_DIR, "task-$(task)_eeg.json"),
              replace(tpl, "@TASK@" => task, "@DESCRIPTION@" => desc))
        open(joinpath(BIDS_DIR, "task-$(task)_channels.tsv"), "w") do io
            println(io, "name\ttype\tunits\tsampling_frequency\treference")
            for c in ch_names
                println(io, "$c\tEEG\t$unit\t$(Int(fs))\tFCz")
            end
        end
    end
end

function write_session_files(dst_dir::String, sub::String, ses::String)
    pre = "sub-$(sub)_ses-$(ses)"
    cp(joinpath(TEMPLATES, "montage_eci_31ch.tsv"), joinpath(dst_dir, pre * "_electrodes.tsv"); force=true)
    cp(joinpath(TEMPLATES, "coordsystem.json"),     joinpath(dst_dir, pre * "_coordsystem.json"); force=true)
end

function write_crop_sidecar(path::String, row::InventoryRow, fs::Float64)
    i0, i1 = round(Int, row.tmin_s * fs), round(Int, row.tmax_s * fs)
    write(path, """{
    "Description": "The eyes-open and eyes-closed blocks were recorded in one raw file. This recording is the block between SourceWindow.TminS and TmaxS of that file (half-open interval, cut before any filtering).",
    "RecordingDuration": $((i1 - i0) / fs),
    "SourceWindow": {
        "TminS": $(row.tmin_s),
        "TmaxS": $(row.tmax_s),
        "SourceSamples": [$i0, $i1],
        "Interval": "[TminS, TmaxS)"
    }
}
""")
end

# ─── Limpieza de lo que ya no es válido (nada se borra) ───────

function sweep_stale!(expected::Set{String})
    moved = 0
    for sub in filter(d -> startswith(d, "sub-") && isdir(joinpath(BIDS_DIR, d)), readdir(BIDS_DIR))
        for (root, _, files) in walkdir(joinpath(BIDS_DIR, sub)), f in files
            full = joinpath(root, f)
            rel  = relpath(full, BIDS_DIR)
            (rel in expected || startswith(f, ".")) && continue
            dst = joinpath(DISCARD_DIR, rel)
            mkpath(dirname(dst))
            mv(full, dst; force=true)
            open(joinpath(DISCARD_DIR, "README.md"), "a") do io
                println(io, "- $(Dates.format(now(), "yyyy-mm-dd")) `$rel` — ya no figura como válido en inventory.csv (ver config/raw_selection.csv)")
            end
            moved += 1
        end
    end
    # carpetas vacías
    for sub in filter(d -> startswith(d, "sub-"), readdir(BIDS_DIR))
        for (root, _, _) in walkdir(joinpath(BIDS_DIR, sub); topdown=false)
            isempty(readdir(root)) && rm(root)
        end
    end
    return moved
end

# ─── Función principal ────────────────────────────────────────

function main()
    println("=" ^ 60)
    println(" NeuroMIND — Fase B: Construcción BIDS-EEG")
    println(" $(now())")
    println("=" ^ 60)

    inv_path = INVENTORY
    only_subjects = String[]
    force = "--force" in ARGS
    for (i, a) in enumerate(ARGS)
        a == "--inventory" && i < length(ARGS) && (inv_path = ARGS[i+1])
        a == "--subjects"  && i < length(ARGS) && (only_subjects = String.(split(ARGS[i+1], ",")))
    end

    rows  = load_inventory(inv_path)
    valid = filter(r -> !r.excluded && haskey(TASKS, r.condition), rows)
    isempty(only_subjects) || (valid = filter(r -> r.bids_id in only_subjects, valid))
    println("✓ Inventario: $(length(rows)) entradas · $(length(valid)) grabaciones válidas a construir")
    isempty(valid) && error("Nada que construir")

    # Las grabaciones válidas deben compartir canales y orden (cierra MONT-07)
    ref = nothing
    for r in valid
        info = read_vhdr(joinpath(PROJ_ROOT, r.filepath))
        if ref === nothing
            ref = info
        else
            info.ch_names == ref.ch_names ||
                error("Canales/orden distintos en $(r.filename): $(info.ch_names)")
            info.fs == ref.fs || error("fs distinta en $(r.filename): $(info.fs)")
        end
    end
    println("✓ $(length(valid)) cabeceras con los mismos $(length(ref.ch_names)) canales (último: $(last(ref.ch_names))) @ $(ref.fs) Hz")

    write_dataset_files(ref.ch_names, ref.unit, ref.fs)

    counts = Dict(:copied => 0, :kept => 0, :cropped => 0)
    expected = Set{String}()
    sessions = Set{Tuple{String,String}}()

    for r in sort(valid, by = r -> (r.bids_id, r.bids_session, r.condition))
        task, _ = TASKS[r.condition]
        base    = "sub-$(r.bids_id)_ses-$(r.bids_session)_task-$(task)_run-01"
        dst_dir = joinpath(BIDS_DIR, "sub-$(r.bids_id)", "ses-$(r.bids_session)", "eeg")
        mkpath(dst_dir)

        status = write_recording(r, joinpath(PROJ_ROOT, r.filepath), dst_dir, base; force)
        counts[status] += 1
        for ext in ("eeg", "vhdr", "vmrk")
            push!(expected, relpath(joinpath(dst_dir, "$(base)_eeg.$(ext)"), BIDS_DIR))
        end
        if status == :cropped
            side = joinpath(dst_dir, base * "_eeg.json")
            write_crop_sidecar(side, r, ref.fs)
            push!(expected, relpath(side, BIDS_DIR))
        end
        if (r.bids_id, r.bids_session) ∉ sessions
            push!(sessions, (r.bids_id, r.bids_session))
            write_session_files(dst_dir, r.bids_id, r.bids_session)
            pre = "sub-$(r.bids_id)_ses-$(r.bids_session)"
            push!(expected, relpath(joinpath(dst_dir, pre * "_electrodes.tsv"),  BIDS_DIR))
            push!(expected, relpath(joinpath(dst_dir, pre * "_coordsystem.json"), BIDS_DIR))
        end
    end

    n_moved = isempty(only_subjects) ? sweep_stale!(expected) : 0

    println("\n── Resumen ──────────────────────────────────────────────")
    println("  Registros copiados       : $(counts[:copied])")
    println("  Ya existían (sin recopiar): $(counts[:kept])")
    println("  Registros recortados     : $(counts[:cropped])")
    println("  Sujetos / sesiones       : $(length(unique(first.(sessions)))) / $(length(sessions))")
    println("  Apartados a _archive     : $n_moved")
    println("  Dataset BIDS en          : $BIDS_DIR")
    println("\n✅ Fase B completada: $(Dates.format(now(), "HH:MM:SS"))")
    println("\nValidar:  bids-validator-deno $(relpath(BIDS_DIR, pwd()))")
end

# ─── Punto de entrada ─────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
