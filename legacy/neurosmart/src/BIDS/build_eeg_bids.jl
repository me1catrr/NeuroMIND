#!/usr/bin/env julia

using ArgParse
using DataFrames
using CSV
using Dates
using JSON
using Glob

# ---------- funciones auxiliares ----------

function parse_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--sourcedata"
            help = "Carpeta con archivos EEG originales"
            arg_type = String
            default = "sourcedata"
        "--bids"
            help = "Carpeta de salida BIDS"
            arg_type = String
            default = "bids"
        "--participants"
            help = "Archivo participants.tsv generado previamente"
            arg_type = String
            default = "bids/participants.tsv"
    end
    return ArgParse.parse_args(s)
end

# Extraer información del nombre de archivo
function parse_filename(filename)
    basename_file = basename(filename)
    
    # Patrón 1: M{ID}_T{sesión}_{iniciales}_{fecha}_{condición}
    pattern1 = r"^M(\d+)_T(\d+)_([^_]+)_(\d+)_(.+)\.(vhdr|vmrk|eeg)$"
    m1 = match(pattern1, basename_file)
    
    # Patrón 2: MC{ID}_{iniciales}_{fecha}_{condición}
    pattern2 = r"^MC(\d+)_([^_]+)_(\d+)_(.+)\.(vhdr|vmrk|eeg)$"
    m2 = match(pattern2, basename_file)
    
    # Patrón 3: M{ID}_T{sesión}_{iniciales}_{fecha}_{condición} (con formato de fecha diferente)
    pattern3 = r"^M(\d+)_T(\d+)_(\d+)_([^_]+)_(\d+)_(.+)\.(vhdr|vmrk|eeg)$"
    m3 = match(pattern3, basename_file)
    
    if m1 !== nothing
        participant_id = "M$(m1.captures[1])"
        session = parse(Int, m1.captures[2])
        initials = m1.captures[3]
        date_str = m1.captures[4]
        condition = m1.captures[5]
        file_type = m1.captures[6]
    elseif m2 !== nothing
        participant_id = "MC$(m2.captures[1])"
        session = 1  # Asumir sesión 1 para archivos MC
        initials = m2.captures[2]
        date_str = m2.captures[3]
        condition = m2.captures[4]
        file_type = m2.captures[5]
    elseif m3 !== nothing
        participant_id = "M$(m3.captures[1])"
        session = parse(Int, m3.captures[2])
        initials = m3.captures[4]
        date_str = m3.captures[5]
        condition = m3.captures[6]
        file_type = m3.captures[7]
    else
        return nothing
    end
    
    # Normalizar condición (maneja separadores y acentos comunes)
    condition_normalized = lowercase(condition)
    condition_normalized = replace(
        condition_normalized,
        "á" => "a",
        "é" => "e",
        "í" => "i",
        "ó" => "o",
        "ú" => "u",
        "_" => "",
        "-" => "",
        " " => ""
    )

    if occursin("ojosabiertos", condition_normalized)
        task = "eyesopen"
    elseif occursin("ojoscerrados", condition_normalized)
        task = "eyesclosed"
    else
        task = "unknown"
    end
    
    return Dict(
        "participant_id" => participant_id,
        "session" => session,
        "initials" => initials,
        "date" => date_str,
        "condition" => condition,
        "task" => task,
        "file_type" => file_type
    )
end

# Leer metadatos del archivo .vhdr
function read_vhdr_metadata(vhdr_path)
    metadata = Dict()
    
    if !isfile(vhdr_path)
        return metadata
    end
    
    lines = readlines(vhdr_path)
    
    for line in lines
        line = strip(line)
        
        # Información común
        if startswith(line, "DataFile=")
            metadata["data_file"] = split(line, "=")[2]
        elseif startswith(line, "MarkerFile=")
            metadata["marker_file"] = split(line, "=")[2]
        elseif startswith(line, "NumberOfChannels=")
            metadata["n_channels"] = parse(Int, split(line, "=")[2])
        elseif startswith(line, "SamplingInterval=")
            sampling_interval = parse(Int, split(line, "=")[2])
            metadata["sampling_rate"] = 1000000 / sampling_interval  # convertir de µs a Hz
        elseif startswith(line, "DataFormat=")
            metadata["data_format"] = split(line, "=")[2]
        elseif startswith(line, "DataOrientation=")
            metadata["data_orientation"] = split(line, "=")[2]
        elseif startswith(line, "BinaryFormat=")
            metadata["binary_format"] = split(line, "=")[2]
        end
        
        # Información de canales
        if startswith(line, "Ch") && contains(line, "=")
            parts = split(line, "=")
            if length(parts) >= 2
                ch_info = split(parts[2], ",")
                if length(ch_info) >= 4
                    ch_num = parse(Int, replace(parts[1], "Ch" => ""))
                    ch_name = ch_info[1]
                    ch_ref = ch_info[2]
                    ch_resolution = parse(Float64, ch_info[3])
                    ch_unit = ch_info[4]
                    
                    if !haskey(metadata, "channels")
                        metadata["channels"] = []
                    end
                    
                    push!(metadata["channels"], Dict(
                        "number" => ch_num,
                        "name" => ch_name,
                        "reference" => ch_ref,
                        "resolution" => ch_resolution,
                        "unit" => ch_unit
                    ))
                end
            end
        end
    end
    
    return metadata
end

# Leer marcadores del archivo .vmrk
function read_vmrk_markers(vmrk_path)
    markers = []
    
    if !isfile(vmrk_path)
        return markers
    end
    
    lines = readlines(vmrk_path)
    in_marker_section = false
    
    for line in lines
        line = strip(line)
        
        if line == "[Marker Infos]"
            in_marker_section = true
            continue
        elseif startswith(line, "[")
            in_marker_section = false
            continue
        end
        
        if in_marker_section && startswith(line, "Mk") && contains(line, "=")
            parts = split(line, "=")
            if length(parts) >= 2
                marker_info = split(parts[2], ",")
                if length(marker_info) >= 3
                    marker_type = marker_info[1]
                    marker_desc = marker_info[2]
                    marker_pos = parse(Int, marker_info[3])
                    
                    push!(markers, Dict(
                        "type" => marker_type,
                        "description" => marker_desc,
                        "position" => marker_pos
                    ))
                end
            end
        end
    end
    
    return markers
end

# Crear estructura de carpetas BIDS
function create_bids_structure(bids_dir, participant_id, session)
    # Crear carpetas principales
    sub_dir = joinpath(bids_dir, "sub-$(participant_id)")
    ses_dir = joinpath(sub_dir, "ses-$(session)")
    eeg_dir = joinpath(ses_dir, "eeg")
    
    mkpath(eeg_dir)
    
    return sub_dir, ses_dir, eeg_dir
end

# Convertir archivo BrainVision a BIDS
function convert_brainvision_to_bids(vhdr_path, bids_eeg_dir, participant_id, session, task)
    # Leer metadatos
    metadata = read_vhdr_metadata(vhdr_path)
    markers = read_vmrk_markers(replace(vhdr_path, ".vhdr" => ".vmrk"))
    
    # Crear nombres de archivos BIDS
    base_name = "sub-$(participant_id)_ses-$(session)_task-$(task)"
    
    # Copiar archivos principales
    eeg_file = replace(vhdr_path, ".vhdr" => ".eeg")
    vmrk_file = replace(vhdr_path, ".vhdr" => ".vmrk")
    
    bids_eeg_file = joinpath(bids_eeg_dir, "$(base_name)_eeg.eeg")
    bids_vhdr_file = joinpath(bids_eeg_dir, "$(base_name)_eeg.vhdr")
    bids_vmrk_file = joinpath(bids_eeg_dir, "$(base_name)_eeg.vmrk")
    
    # Copiar archivos si existen
    if isfile(eeg_file)
        cp(eeg_file, bids_eeg_file, force=true)
    end
    if isfile(vhdr_path)
        cp(vhdr_path, bids_vhdr_file, force=true)
    end
    if isfile(vmrk_file)
        cp(vmrk_file, bids_vmrk_file, force=true)
    end
    
    # Crear archivo _eeg.json
    eeg_json = Dict(
        "TaskName" => task,
        "SamplingFrequency" => get(metadata, "sampling_rate", 500),
        "PowerLineFrequency" => 50,
        "SoftwareFilters" => Dict(
            "HighPassFilter" => Dict(
                "FilterType" => "Butterworth",
                "HighPassCutoffHz" => 0.1
            ),
            "LowPassFilter" => Dict(
                "FilterType" => "Butterworth", 
                "LowPassCutoffHz" => 140
            )
        ),
        "EEGChannelCount" => get(metadata, "n_channels", 31),
        "EOGChannelCount" => 0,
        "ECGChannelCount" => 0,
        "EMGChannelCount" => 0,
        "MiscChannelCount" => 0,
        "TriggerChannelCount" => 0,
        "RecordingType" => "continuous",
        "DeviceSoftwareVersion" => "BrainVision Recorder Professional V. 1.21.0303",
        "EEGPlacementScheme" => "10-20",
        "EEGGround" => "Fpz",
        "EEGReference" => "Cz",
        "Manufacturer" => "Brain Products",
        "ManufacturersModelName" => "BrainAmp",
        "CapManufacturer" => "EasyCap",
        "CapManufacturersModelName" => "actiCAP",
        "HardwareFilters" => Dict(
            "HighPassFilter" => Dict(
                "FilterType" => "DC",
                "HighPassCutoffHz" => 0
            ),
            "LowPassFilter" => Dict(
                "FilterType" => "Butterworth",
                "LowPassCutoffHz" => 140
            ),
            "NotchFilter" => Dict(
                "FilterType" => "Off"
            )
        )
    )
    
    # Agregar información de canales si está disponible
    if haskey(metadata, "channels")
        eeg_json["EEGChannelCount"] = length(metadata["channels"])
    end
    
    # Escribir archivo JSON
    eeg_json_path = joinpath(bids_eeg_dir, "$(base_name)_eeg.json")
    open(eeg_json_path, "w") do f
        JSON.print(f, eeg_json, 4)
    end
    
    # Crear archivo _channels.tsv
    if haskey(metadata, "channels")
        channels_df = DataFrame()
        channels_df.name = [ch["name"] for ch in metadata["channels"]]
        channels_df.type = fill("EEG", length(metadata["channels"]))
        channels_df.units = [ch["unit"] for ch in metadata["channels"]]
        channels_df.low_cutoff = fill("n/a", length(metadata["channels"]))
        channels_df.high_cutoff = fill(140.0, length(metadata["channels"]))
        channels_df.description = fill("n/a", length(metadata["channels"]))
        channels_df.sampling_frequency = fill(get(metadata, "sampling_rate", 500), length(metadata["channels"]))
        channels_df.manufacturer = fill("Brain Products", length(metadata["channels"]))
        channels_df.manufacturer_sn = fill("n/a", length(metadata["channels"]))
        channels_df.group = fill("n/a", length(metadata["channels"]))
        channels_df.status = fill("good", length(metadata["channels"]))
        channels_df.status_description = fill("n/a", length(metadata["channels"]))
        channels_df.bad_channel = fill(false, length(metadata["channels"]))
        
        channels_tsv_path = joinpath(bids_eeg_dir, "$(base_name)_channels.tsv")
        CSV.write(channels_tsv_path, channels_df, delim='\t')
    end
    
    # Crear archivo _events.tsv si hay marcadores
    if !isempty(markers)
        events_df = DataFrame()
        events_df.onset = [Float64(m["position"]) / get(metadata, "sampling_rate", 500) for m in markers]
        events_df.duration = fill(0.0, length(markers))
        events_df.trial_type = [m["type"] for m in markers]
        events_df.value = [isempty(m["description"]) ? "n/a" : m["description"] for m in markers]
        events_df.sample = [m["position"] for m in markers]
        
        events_tsv_path = joinpath(bids_eeg_dir, "$(base_name)_events.tsv")
        CSV.write(events_tsv_path, events_df, delim='\t')
    end
    
    return base_name
end

# ---------- función principal ----------

function main()
    args = parse_args()
    
    sourcedata_dir = args["sourcedata"]
    bids_dir = args["bids"]
    participants_file = args["participants"]
    
    # Leer participantes
    if isfile(participants_file)
        participants_df = CSV.read(participants_file, DataFrame)
        println("✓ Cargados $(nrow(participants_df)) participantes")
    else
        error("No se encontró el archivo de participantes: $participants_file")
    end
    
    # Buscar archivos .vhdr
    vhdr_files = glob("*.vhdr", sourcedata_dir)
    println("✓ Encontrados $(length(vhdr_files)) archivos .vhdr")
    
    processed_files = 0
    skipped_files = 0
    
    for vhdr_file in vhdr_files
        file_info = parse_filename(vhdr_file)
        
        if file_info === nothing
            println("⚠ Saltando archivo con formato no reconocido: $(basename(vhdr_file))")
            skipped_files += 1
            continue
        end
        
        participant_id = file_info["participant_id"]
        session = file_info["session"]
        task = file_info["task"]
        
        # Buscar el participante en la lista (puede estar como M4 o sub-M4)
        bids_participant_id = nothing
        if participant_id in participants_df.participant_id
            bids_participant_id = participant_id
        elseif "sub-$participant_id" in participants_df.participant_id
            bids_participant_id = "sub-$participant_id"
        else
            println("⚠ Participante $participant_id no encontrado en participants.tsv")
            skipped_files += 1
            continue
        end
        
        # Extraer el ID numérico para crear la estructura de carpetas
        if startswith(bids_participant_id, "sub-")
            folder_participant_id = replace(bids_participant_id, "sub-" => "")
        else
            folder_participant_id = participant_id
        end
        
        # Crear estructura de carpetas
        sub_dir, ses_dir, eeg_dir = create_bids_structure(bids_dir, folder_participant_id, session)
        
        # Convertir archivo
        try
            base_name = convert_brainvision_to_bids(vhdr_file, eeg_dir, folder_participant_id, session, task)
            println("✓ Procesado: $base_name")
            processed_files += 1
        catch e
            println("✗ Error procesando $(basename(vhdr_file)): $e")
            skipped_files += 1
        end
    end
    
    println("\n=== Resumen ===")
    println("✓ Archivos procesados: $processed_files")
    println("⚠ Archivos saltados: $skipped_files")
    println("✓ Estructura BIDS creada en: $bids_dir")
end

# Ejecutar si se llama directamente
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end