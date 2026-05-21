# NeuroMIND/src/io/BIDSLoader.jl
# Carga de datos EEG en formato BIDS.
# Formato real del proyecto: TSV con filas=canales, columnas=muestras.
# Primera columna = "Channel" con nombres de canal.

"""
    load_eeg_bids(cfg, subject_id, session_id, condition; run=1) -> EEGRecording

Carga una grabación EEG desde la estructura BIDS.
Condición: "EC" (eyes closed) o "EO" (eyes open) → mapeado a task name.
"""
function load_eeg_bids(
    cfg::PipelineConfig,
    subject_id::String,
    session_id::String,
    condition::String;
    run::Int = 1
)::EEGRecording

    task = condition == "EO" ? "eyesopen" : "eyesclosed"
    bids_dir = bids_root_dir(cfg)
    prefix   = "sub-$(subject_id)_ses-$(session_id)_task-$(task)_run-$(lpad(run,2,'0'))_eeg"
    raw_dir  = joinpath(bids_dir, "raw")

    data_path = joinpath(raw_dir, "$(prefix)_data.tsv")
    meta_path = joinpath(raw_dir, "$(prefix)_metadata.json")
    elec_path = joinpath(bids_dir, "electrodes",
                         "sub-$(subject_id)_ses-$(session_id)_electrodes.tsv")

    isfile(data_path) || error("Archivo EEG no encontrado: $data_path")
    isfile(meta_path) || error("Metadatos no encontrados: $meta_path")

    meta_raw = _parse_bids_json(meta_path)
    fs = Float64(meta_raw["fs"])

    # TSV: filas=canales, primera columna="Channel", resto=muestras temporales
    df       = CSV.read(data_path, DataFrame; delim='\t')
    ch_names = string.(df[:, 1])
    data     = Float64.(Matrix(df[:, 2:end]))   # (channels × samples)

    ch_pos = isfile(elec_path) ? _load_electrode_positions(elec_path) : nothing

    meta = RecordingMeta(
        subject_id, session_id, condition, run,
        fs, length(ch_names), ch_names, ch_pos, data_path
    )
    times = collect(0.0:(1/fs):(size(data, 2) - 1) / fs)

    return EEGRecording(meta, data, times)
end

# ─── Helpers privados ─────────────────────────────────────────

function _load_electrode_positions(tsv_path::String)
    df  = CSV.read(tsv_path, DataFrame; delim='\t')
    pos = Dict{String,Tuple{Float64,Float64}}()
    for r in eachrow(df)
        type_col = hasproperty(r, :type) ? string(r.type) : "EEG"
        if type_col == "EEG"
            name = uppercase(strip(string(r.name)))
            pos[name] = (Float64(r.x), Float64(r.y))
        end
    end
    return pos
end

"""
    _parse_bids_json(path) -> Dict

Parser robusto para los JSON de metadata del proyecto.
Soporta los campos reales: "fs", "channel_names", "n_samples", etc.
"""
function _parse_bids_json(path::String)::Dict{String,Any}
    raw = read(path, String)
    d   = Dict{String,Any}()

    # Sampling frequency: "fs" (formato real) o "SamplingFrequency" (BIDS estándar)
    for pat in [r"\"fs\"\s*:\s*([0-9]+\.?[0-9]*)",
                r"\"SamplingFrequency\"\s*:\s*([0-9]+\.?[0-9]*)"]
        m = match(pat, raw)
        if m !== nothing
            d["fs"] = parse(Float64, m.captures[1])
            break
        end
    end

    # n_samples
    m = match(r"\"n_samples\"\s*:\s*([0-9]+)", raw)
    m !== nothing && (d["n_samples"] = parse(Int, m.captures[1]))

    # duration
    m = match(r"\"duration_s\"\s*:\s*([0-9]+\.?[0-9]*)", raw)
    m !== nothing && (d["duration_s"] = parse(Float64, m.captures[1]))

    # n_channels
    m = match(r"\"n_channels\"\s*:\s*([0-9]+)", raw)
    m !== nothing && (d["n_channels"] = parse(Int, m.captures[1]))

    # subject / session / task
    for (key, pat) in [("subject", r"\"subject\"\s*:\s*\"([^\"]+)\""),
                       ("session", r"\"session\"\s*:\s*\"([^\"]+)\""),
                       ("task",    r"\"task\"\s*:\s*\"([^\"]+)\"")]
        m = match(pat, raw)
        m !== nothing && (d[key] = String(m.captures[1]))
    end

    # channel_names (puede ser multilínea)
    for pat in [r"\"channel_names\"\s*:\s*\[([^\]]+)\]"s,
                r"\"ChannelNames\"\s*:\s*\[([^\]]+)\]"s]
        m = match(pat, raw)
        if m !== nothing
            parts = split(m.captures[1], ',')
            names = [strip(strip(p), ['"', ' ', '\n', '\r', '\t']) for p in parts]
            d["channel_names"] = filter(!isempty, names)
            break
        end
    end

    return d
end

# Alias de compatibilidad hacia atrás
JSON_parse = _parse_bids_json
