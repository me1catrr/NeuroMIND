# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Construcción BIDS ligera (dataset completo)
# ═══════════════════════════════════════════════════════════════
#
#  Fase B — Lee el inventario de audit_full_dataset.jl y crea:
#    · metadata JSON por grabación (apunta al .vhdr; no copia datos)
#    · electrodes TSV por sujeto/sesión (plantilla 10-20, 31 ch)
#    · dataset_description.json
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/build_bids_full.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 22-07-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/<este-script>.jl …
#  Sin shebang: #!/usr/bin/env julia no activaría --project=.
#
#  Prerrequisito
#  ─────────────
#    data/full_data/inventory.csv
#
#  Uso
#  ───
#    julia --project=. scripts/build_bids_full.jl
#    julia --project=. scripts/build_bids_full.jl --inventory ruta.csv
#
#  Opciones
#  ────────
#    --inventory  PATH   CSV de inventario (defecto: data/full_data/inventory.csv)
#
#  Salida
#  ──────
#    data/bids/sub-*/ses-*/…  (metadatos BIDS ligeros)

using Dates

# ─── Rutas ────────────────────────────────────────────────────

const PROJ_ROOT    = dirname(@__DIR__)
const BIDS_DIR     = joinpath(PROJ_ROOT, "data", "bids")
const INVENTORY    = joinpath(PROJ_ROOT, "data", "full_data", "inventory.csv")
const ELEC_TMPL    = joinpath(BIDS_DIR, "electrodes",
                              "sub-M05_ses-T2_electrodes.tsv")   # plantilla existente

# Metadata de los canales (obtenidos del vhdr de M11)
const CH_NAMES_31 = [
    "Fz","F3","F7","FT9","FC5","FC1","C3","T7","TP9","CP5","CP1",
    "Pz","P3","P7","O1","Oz","O2","P4","P8","TP10","CP6","CP2",
    "Cz","C4","T8","FT10","FC6","FC2","F4","F8","Fp2"
]

# ─── Lectura del inventario ───────────────────────────────────

struct InventoryRow
    filename     :: String
    subject_id   :: String
    bids_id      :: String
    group        :: String
    session      :: String
    bids_session :: String
    condition    :: String
    excluded     :: Bool
    filepath     :: String
end

function load_inventory(path::String)::Vector{InventoryRow}
    rows = InventoryRow[]
    isfile(path) || error("Inventario no encontrado: $path\n  → Ejecuta primero: julia scripts/audit_full_dataset.jl")
    lines = readlines(path)
    length(lines) < 2 && return rows

    # Cabecera esperada:
    # filename,subject_id,bids_id,group,session,session_sub,bids_session,
    # condition_raw,condition,date_str,initials,excluded,sexo,edad,nivel_edu,filepath,note
    for line in lines[2:end]
        isempty(strip(line)) && continue
        parts = split(line, ",")
        length(parts) < 16 && continue
        excluded = strip(String(parts[12])) == "true"
        push!(rows, InventoryRow(
            strip(String(parts[1])),   # filename
            strip(String(parts[2])),   # subject_id
            strip(String(parts[3])),   # bids_id
            strip(String(parts[4])),   # group
            strip(String(parts[5])),   # session
            strip(String(parts[7])),   # bids_session
            strip(String(parts[9])),   # condition (normalizada)
            excluded,
            strip(String(parts[16])),  # filepath
        ))
    end
    return rows
end

# ─── Lectura de cabecera .vhdr ────────────────────────────────

function read_vhdr_basics(vhdr_path::String)
    n_ch = 31
    fs   = 500.0
    try
        for line in readlines(vhdr_path; keep=false)
            startswith(line, "NumberOfChannels=") &&
                (n_ch = parse(Int, split(line,"=")[2]); continue)
            startswith(line, "SamplingInterval=") &&
                (fs = 1_000_000.0 / parse(Float64, split(line,"=")[2]); continue)
        end
    catch
    end
    return n_ch, fs
end

# ─── Escritura de metadata JSON ───────────────────────────────

function write_metadata_json(path::String, row::InventoryRow,
                              n_ch::Int, fs::Float64)
    task = row.condition == "EC" ? "eyesclosed" : "eyesopen"
    open(path, "w") do io
        print(io, """{
    "data_format": "brainvision",
    "vhdr_path": "$(escape_string(row.filepath))",
    "fs": $(fs),
    "n_channels": $(n_ch),
    "channel_names": $(json_string_array(CH_NAMES_31[1:min(n_ch,31)])),
    "subject": "$(row.bids_id)",
    "session": "$(row.bids_session)",
    "task": "$(task)",
    "run": 1,
    "group": "$(row.group)",
    "PowerLineFrequency": 50,
    "EEGReference": "Cz",
    "EEGPlacementScheme": "10-20",
    "Manufacturer": "Brain Products",
    "ManufacturersModelName": "BrainAmp",
    "SamplingFrequency": $(fs),
    "created_at": "$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))"
}
""")
    end
end

function json_string_array(v::Vector{String})::String
    "[" * join(["\"$s\"" for s in v], ", ") * "]"
end

# ─── Copia de electrodes TSV ──────────────────────────────────

function ensure_electrodes_tsv(elec_dir::String, bids_id::String, sess::String)
    dest = joinpath(elec_dir, "sub-$(bids_id)_ses-$(sess)_electrodes.tsv")
    isfile(dest) && return dest
    if isfile(ELEC_TMPL)
        cp(ELEC_TMPL, dest; force=false)
    else
        # Crear electrodes TSV mínimo si no existe plantilla
        open(dest, "w") do io
            println(io, "name\tx\ty\tz\ttype")
            positions = [
                ("Fz",0.0,0.7,0.7141),("F3",-0.35,0.68,0.6443),
                ("F7",-0.85,0.68,0.0),("FT9",-0.95,0.2,0.2398),
                ("FC5",-0.6,0.68,0.4214),("FC1",-0.6,0.68,0.4214),
                ("C3",-0.35,0.0,0.9367),("T7",-0.95,0.0,0.3122),
                ("TP9",-0.95,-0.2,0.2398),("CP5",-0.6,0.0,0.8),
                ("CP1",-0.6,0.0,0.8),("Pz",0.0,-0.7,0.7141),
                ("P3",-0.35,-0.68,0.6443),("P7",-0.85,-0.68,0.0),
                ("O1",-0.6,-0.9,0.0),("Oz",0.0,-0.9,0.4359),
                ("O2",0.6,-0.9,0.0),("P4",0.35,-0.68,0.6443),
                ("P8",0.85,-0.68,0.0),("TP10",0.95,-0.2,0.2398),
                ("CP6",0.6,0.0,0.8),("CP2",0.6,0.0,0.8),
                ("Cz",0.0,0.0,1.0),("C4",0.35,0.0,0.9367),
                ("T8",0.95,0.0,0.3122),("FT10",0.95,0.2,0.2398),
                ("FC6",0.6,0.68,0.4214),("FC2",0.6,0.68,0.4214),
                ("F4",0.35,0.68,0.6443),("F8",0.85,0.68,0.0),
                ("Fp2",0.35,0.9,0.2598),
            ]
            for (name,x,y,z) in positions
                println(io, "$(name)\t$(x)\t$(y)\t$(z)\tEEG")
            end
            println(io, "FCz\t0.0\t0.35\t0.9367\tREF")
            println(io, "Fpz\t0.0\t0.9\t0.4359\tGND")
        end
    end
    return dest
end

# ─── dataset_description.json ─────────────────────────────────

function write_dataset_description()
    path = joinpath(BIDS_DIR, "dataset_description.json")
    isfile(path) && return
    open(path, "w") do io
        print(io, """{
    "Name": "MINDEM-IMIBIC Resting-State EEG",
    "BIDSVersion": "1.7.0",
    "License": "CC-BY-4.0",
    "Authors": ["Rafael Castro Triguero"],
    "Acknowledgements": "Grupo BRAIN, Hospital Universitario Reina Sofía, Córdoba",
    "HowToAcknowledge": "Please cite the NeuroMIND paper when using this dataset.",
    "DatasetType": "raw",
    "Description": "Resting-state EEG (eyes-closed and eyes-open) recorded in MS patients (n=41) and healthy controls (n=37) at two time points (T1 baseline, T2 follow-up ~12 months). 31-channel 10-20 system, 500 Hz, BrainAmp amplifier.",
    "ReferencesAndLinks": ["https://github.com/me1catrr/NeuroMIND"]
}
""")
    end
    println("✓ dataset_description.json")
end

# ─── Función principal ────────────────────────────────────────

function main()
    println("=" ^ 60)
    println(" NeuroMIND — Fase B: Construcción BIDS Ligera")
    println(" $(now())")
    println("=" ^ 60)

    rows = load_inventory(INVENTORY)
    println("✓ Inventario cargado: $(length(rows)) entradas")

    valid = filter(r -> !r.excluded && r.condition ∈ ("EC","EO"), rows)
    println("  Grabaciones válidas para pipeline: $(length(valid))")

    # Crear directorios raíz
    raw_dir  = joinpath(BIDS_DIR, "raw")
    elec_dir = joinpath(BIDS_DIR, "electrodes")
    mkpath(raw_dir)
    mkpath(elec_dir)

    write_dataset_description()

    n_created  = 0
    n_existing = 0
    n_errors   = 0

    for row in sort(valid, by = r -> (r.bids_id, r.bids_session, r.condition))
        task = row.condition == "EC" ? "eyesclosed" : "eyesopen"
        base = "sub-$(row.bids_id)_ses-$(row.bids_session)_task-$(task)_run-01"

        meta_path = joinpath(raw_dir, "$(base)_eeg_metadata.json")

        # Electrodes TSV
        elec_path = ensure_electrodes_tsv(elec_dir, row.bids_id, row.bids_session)

        if isfile(meta_path)
            n_existing += 1
            continue
        end

        # Verificar que el .vhdr existe
        if !isfile(row.filepath)
            @warn "Archivo fuente no encontrado: $(row.filepath)"
            n_errors += 1
            continue
        end

        n_ch, fs = read_vhdr_basics(row.filepath)
        try
            write_metadata_json(meta_path, row, n_ch, fs)
            n_created += 1
        catch e
            @warn "Error creando metadata para $(row.bids_id)/$(row.bids_session)/$(task): $e"
            n_errors += 1
        end
    end

    println("\n── Resumen ──────────────────────────────────────────────")
    println("  Metadata JSON creados    : $n_created")
    println("  Ya existían              : $n_existing")
    println("  Errores                  : $n_errors")
    println("  Estructura BIDS en       : $BIDS_DIR")

    println("\n✅ Fase B completada: $(Dates.format(now(), "HH:MM:SS"))")
    println("\nPróximo paso:")
    println("  julia --project=. scripts/run_batch_pipeline.jl --condition EC")
end

# ─── Punto de entrada ─────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
