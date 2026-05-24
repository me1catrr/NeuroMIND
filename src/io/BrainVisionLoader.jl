# NeuroMIND/src/io/BrainVisionLoader.jl
#
# Carga de datos EEG en formato BrainVision nativo (.vhdr + .eeg).
# Soporta IEEE_FLOAT_32 en orientación MULTIPLEXED (el formato estándar
# del amplificador BrainAmp usado en este proyecto).
#
# IMPORTANTE: BrainVision Recorder guarda los datos en unidades digitales
# (ADC) incluso en formato IEEE_FLOAT_32. La resolución del .vhdr DEBE
# aplicarse siempre: value_µV = raw_value × resolution_from_vhdr.
# Esto aplica tanto para IEEE_FLOAT_32 como para INT_16.

# ─── Lectura de cabecera .vhdr ────────────────────────────────

"""
    read_vhdr_header(vhdr_path) -> Dict{String,Any}

Parsea el archivo de cabecera BrainVision (.vhdr).
Claves devueltas:
  "n_channels"     Int     — número de canales en el .eeg
  "fs"             Float64 — tasa de muestreo en Hz
  "ch_names"       Vector{String}  — nombres de canales (en orden del fichero)
  "resolutions"    Vector{Float64} — factor de escala µV/unidad (aplicar siempre)
  "binary_format"  String  — "IEEE_FLOAT_32" | "INT_16" | ...
  "orientation"    String  — "MULTIPLEXED" | "VECTORIZED"
  "eeg_file"       String  — ruta relativa/absoluta al archivo .eeg
"""
function read_vhdr_header(vhdr_path::String)::Dict{String,Any}
    d = Dict{String,Any}(
        "n_channels"    => 0,
        "fs"            => 500.0,
        "ch_names"      => String[],
        "resolutions"   => Float64[],
        "binary_format" => "IEEE_FLOAT_32",
        "orientation"   => "MULTIPLEXED",
        "eeg_file"      => ""
    )

    lines = readlines(vhdr_path; keep=false)
    for line in lines
        line = strip(line)
        isempty(line) && continue
        startswith(line, ";") && continue   # comentario

        if startswith(line, "DataFile=")
            d["eeg_file"] = String(split(line, "="; limit=2)[2])
        elseif startswith(line, "NumberOfChannels=")
            d["n_channels"] = parse(Int, strip(split(line, "=")[2]))
        elseif startswith(line, "SamplingInterval=")
            us = parse(Float64, strip(split(line, "=")[2]))
            d["fs"] = 1_000_000.0 / us
        elseif startswith(line, "BinaryFormat=")
            d["binary_format"] = strip(String(split(line, "="; limit=2)[2]))
        elseif startswith(line, "DataOrientation=")
            d["orientation"] = strip(String(split(line, "="; limit=2)[2]))
        else
            # Canales: Ch1=Fz,,0.0488281,µV
            m = match(r"^Ch(\d+)=(.+)$", line)
            if m !== nothing
                ch_info = split(m[2], ",")
                push!(d["ch_names"], String(ch_info[1]))
                res = length(ch_info) >= 3 && !isempty(strip(ch_info[3])) ?
                      parse(Float64, strip(ch_info[3])) : 1.0
                push!(d["resolutions"], res)
            end
        end
    end

    # Asegura consistencia de longitud de vectores de canales
    n = d["n_channels"]
    if n > 0
        while length(d["ch_names"]) < n
            push!(d["ch_names"], "CH$(length(d["ch_names"])+1)")
        end
        while length(d["resolutions"]) < n
            push!(d["resolutions"], 1.0)
        end
        # Truncar si hay más (cabecera puede tener canales AUX no incluidos en .eeg)
        d["ch_names"]   = d["ch_names"][1:n]
        d["resolutions"] = d["resolutions"][1:n]
    end

    return d
end

# ─── Carga completa desde BrainVision ─────────────────────────

"""
    load_eeg_brainvision(vhdr_path, subj_id, sess_id, task;
                         run=1, ch_pos=nothing) -> EEGRecording

Lee una grabación EEG completa desde formato BrainVision binario.

# Argumentos
- `vhdr_path`  Ruta al archivo de cabecera (.vhdr)
- `subj_id`    ID del sujeto (ej. "M11")
- `sess_id`    ID de sesión (ej. "T1")
- `task`       Nombre de tarea BIDS: "eyesclosed" o "eyesopen"
- `run`        Número de run (por defecto 1)
- `ch_pos`     Diccionario de posiciones de electrodos (opcional)

# Notas de formato
- Soporta IEEE_FLOAT_32 e INT_16. En ambos casos se aplica la resolución
  del .vhdr (µV/unidad) para convertir unidades ADC a µV.
- La orientación esperada es MULTIPLEXED: [ch₁t₁, ch₂t₁, …, chₙtₙ].
  En Julia (column-major), `reshape(raw, n_ch, n_samples)` produce la
  matriz correcta con ch como primer índice.
"""
function load_eeg_brainvision(
    vhdr_path::String,
    subj_id::String,
    sess_id::String,
    task::String;
    run::Int    = 1,
    ch_pos      = nothing
)::EEGRecording

    isfile(vhdr_path) || error("Archivo .vhdr no encontrado: $vhdr_path")

    hdr  = read_vhdr_header(vhdr_path)
    n_ch = hdr["n_channels"]
    fs   = hdr["fs"]
    ch_names  = hdr["ch_names"]
    ch_res    = hdr["resolutions"]
    binfmt    = hdr["binary_format"]
    orient    = hdr["orientation"]

    n_ch > 0  || error("Cabecera sin canales: $vhdr_path")

    # ── Ruta al archivo .eeg ──────────────────────────────────
    eeg_rel  = hdr["eeg_file"]
    eeg_path = isabspath(eeg_rel) ? eeg_rel :
               joinpath(dirname(vhdr_path), eeg_rel)
    isfile(eeg_path) || error("Archivo .eeg no encontrado: $eeg_path")

    # ── Lectura binaria ───────────────────────────────────────
    n_bytes = filesize(eeg_path)

    data = if binfmt == "IEEE_FLOAT_32"
        bytes_per_sample = 4
        n_total = div(n_bytes, bytes_per_sample)
        n_t     = div(n_total, n_ch)
        raw     = Vector{Float32}(undef, n_ch * n_t)
        read!(eeg_path, raw)
        # MULTIPLEXED → column-major reshape produce (n_ch × n_t) correcto
        d = Float64.(reshape(raw, n_ch, n_t))
        # BrainVision Recorder guarda en unidades ADC (no µV) incluso como
        # float32. Aplicar resolución canal a canal para convertir a µV.
        for i in 1:n_ch
            d[i, :] .*= ch_res[i]
        end
        d

    elseif binfmt == "INT_16"
        bytes_per_sample = 2
        n_total = div(n_bytes, bytes_per_sample)
        n_t     = div(n_total, n_ch)
        raw     = Vector{Int16}(undef, n_ch * n_t)
        read!(eeg_path, raw)
        d = Float64.(reshape(raw, n_ch, n_t))
        # Aplicar resolución canal a canal
        for i in 1:n_ch
            d[i, :] .*= ch_res[i]
        end
        d

    else
        error("Formato binario no soportado: $binfmt  (soportados: IEEE_FLOAT_32, INT_16)")
    end

    # ── Construcción de EEGRecording ──────────────────────────
    n_samples = size(data, 2)
    condition = task == "eyesclosed" ? "EC" : (task == "eyesopen" ? "EO" : task)

    meta = RecordingMeta(
        subj_id, sess_id, condition, run,
        fs, n_ch, ch_names, ch_pos, vhdr_path
    )
    times = collect(0.0:(1.0/fs):(n_samples - 1) / fs)

    @info "BV cargado: sub-$(subj_id) ses-$(sess_id) $(task) — $(n_ch)ch × $(n_samples) muestras @ $(fs) Hz ($(round(n_samples/fs, digits=1)) s)"

    return EEGRecording(meta, data, times)
end

# ─── Helper: lee posiciones de electrodos para BV ─────────────

"""
    bv_electrode_positions(electrodes_tsv_path) -> Dict{String,Tuple{Float64,Float64}}

Carga posiciones de electrodos desde un TSV BIDS.
Alias conveniente para usar desde el loader de BrainVision.
"""
function bv_electrode_positions(tsv_path::String)
    isfile(tsv_path) || return nothing
    _load_electrode_positions(tsv_path)
end
