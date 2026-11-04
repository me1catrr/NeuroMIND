#!/usr/bin/env julia
# Build BIDS participants.tsv/json from two demographics files.
#
# Usage (desde la raíz del repo):
#   julia --project=. scripts/build_participants.jl \
#       --base config/demographics.csv \
#       --clinical config/demographics_v2.csv \
#       --out bids
#
# Requisitos en Project.toml: CSV, DataFrames, JSON3, ArgParse, FilePathsBase

using ArgParse
using CSV, DataFrames
using JSON3
using FilePathsBase: joinpath
import Base.Filesystem: mkpath

# ---------- helpers ----------
function guess_delim(path::AbstractString)
    # Heurística simple: si la primera línea tiene más ';' que ',' usamos ';'
    open(path, "r") do io
        line = readline(io)
        count_semis = count(==(';'), line)
        count_commas = count(==(','), line)
        return count_semis > count_commas ? ';' : ','
    end
end

striplower(s) = ismissing(s) ? missing : lowercase(strip(String(s)))

# Mapear sexo y grupo a BIDS
function map_sex(x)
    s = striplower(x)
    if ismissing(s)                    ; return missing
    elseif s in ("hombre","m","male")  ; return "M"
    elseif s in ("mujer","f","female") ; return "F"
    elseif isempty(s) || s == "na"     ; return missing
    else                               ; return uppercase(s)[1:1]  # fallback
    end
end

function map_group(x)
    s = striplower(x)
    if ismissing(s)                                       ; return missing
    elseif s in ("pac","patient","ms","esclerosis multiple"); return "patient"
    elseif s in ("con","control","healthy");              return "control"
    elseif isempty(s) || s == "na";                       return missing
    else                                                  return s
    end
end

# Normaliza nombres de columnas a identificadores simples (ASCII aproximado)
function normalize_colnames!(df::DataFrame)
    repl = Dict(
        'Á'=>'A','É'=>'E','Í'=>'I','Ó'=>'O','Ú'=>'U','Ü'=>'U','Ñ'=>'N',
        'á'=>'a','é'=>'e','í'=>'i','ó'=>'o','ú'=>'u','ü'=>'u','ñ'=>'n'
    )
    function clean(name::AbstractString)
        s = String(name)
        s = replace(s, repl...)
        s = replace(s, r"[\.()/]" => " ")
        s = replace(s, r"\s+" => "_")
        s = lowercase(strip(s))
        return s
    end
    rename!(df, Dict(n => Symbol(clean(string(n))) for n in names(df)))
    return df
end

# Selección segura de columna (por múltiples alias)
function col(df::DataFrame, aliases::Vector{String})
    for a in aliases
        if a in names(df); return df[!, a] end
    end
    return fill(missing, nrow(df))
end

# ---------- main ----------
function parse_args()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--base"
            help = "Ruta a demographics.csv (base con todos los sujetos)"
            arg_type = String
            required = true
        "--clinical"
            help = "Ruta a demographics_v2.csv (clínico opcional)"
            arg_type = String
            required = false
            default = ""
        "--out"
            help = "Carpeta de salida BIDS (donde escribir participants.*)"
            arg_type = String
            default = "bids"
        "--powerline"
            help = "No usado aquí; solo recordatorio de entorno EEG"
            arg_type = Int
            default = 50
    end
    return ArgParse.parse_args(s)
end

function main()
    args = parse_args()

    # --- leer base ---
    delim_base = guess_delim(args["base"])
    base = CSV.read(args["base"], DataFrame; delim=delim_base, ignorerepeated=true, missingstring=["NA","N/A","na",""])
    normalize_colnames!(base)

    # columnas esperadas (con alias)
    codigo_base = col(base, ["codigo", "id", "subject", "subject_id"])
    sexo_base   = col(base, ["sexo", "sex"])
    edad_base   = col(base, ["edad", "age"])
    tipo_base   = col(base, ["tipo", "group"])
    edu_base    = col(base, ["educacion", "nivel_educativo", "education", "education_level"])

    # armar DF participants con base
    participants = DataFrame(
        participant_id = "sub-" .* string.(codigo_base),
        age            = tryparse.(Float64, string.(edad_base)),
        sex            = map(map_sex, sexo_base),
        group          = map(map_group, tipo_base),
        education_level = string.(edu_base),
    )

    # --- leer clínico y mezclar (opcional) ---
    if !isempty(args["clinical"])
        delim_c = guess_delim(args["clinical"])
        clin = CSV.read(args["clinical"], DataFrame; delim=delim_c, ignorerepeated=true, missingstring=["NA","N/A","na",""])
        normalize_colnames!(clin)

        # claves de unión - buscar columna de identificador
        if "codigo" ∈ names(clin)
            # Ya existe codigo, asegurar que esté en la primera posición
            if names(clin)[1] != "codigo"
                select!(clin, "codigo", Not("codigo"))
            end
        else
            # Buscar otra columna de identificador y renombrarla a codigo
            for alias in ["id", "subject", "subject_id"]
                if alias ∈ names(clin)
                    rename!(clin, alias => "codigo")
                    # Mover codigo a la primera posición
                    select!(clin, "codigo", Not("codigo"))
                    break
                end
            end
            if "codigo" ∉ names(clin)
                error("No se encontró columna de identificador en el archivo clínico")
            end
        end

        # seleccionar y renombrar variables clínicas a ASCII/snake_case
        select!(clin, intersect(names(clin), [
            "codigo",
            "sexo", "edad", "nivel_educativo", "tipo",
            "anos_evol", "t1_edss_", "t1_tratam", "t1_diagn", "t2_edss", "t2_tratam",
            "discapacidad_t_1", "discapacidad_t_2", "mejora_discapacidad"
        ]), renamecols=false)

        # si vienen con otros nombres (por ejemplo "años_evol"), crear alias
        if "años_evol" in names(clin); clin.anos_evol = clin[!, "años_evol"]; end

        # renombrados finales BIDS-friendly
        ren = Dict(
            "anos_evol"=>"disease_duration_years",
            "t1_edss_"=>"edss_t1",
            "t1_tratam"=>"treatment_t1",
            "t1_diagn"=>"diagnosis_t1",
            "t2_edss"=>"edss_t2",
            "t2_tratam"=>"treatment_t2",
            "discapacidad_t_1"=>"disability_t1",
            "discapacidad_t_2"=>"disability_t2",
            "mejora_discapacidad"=>"disability_improvement"
        )
        rename!(clin, intersect(keys(ren), names(clin)) .=> getindex.(Ref(ren), intersect(keys(ren), names(clin))))

        # nos quedamos solo con columnas de interés para el join
        keep = ["codigo", "disease_duration_years", "edss_t1", "treatment_t1", "diagnosis_t1",
                "edss_t2", "treatment_t2", "disability_t1", "disability_t2", "disability_improvement"]
        keep = filter(in(names(clin)), keep)
        clin2 = select(clin, keep, renamecols=false)

        # añadir clave participant_id para el join
        clin2.participant_id = "sub-" .* string.(clin2.codigo)
        select!(clin2, Not(:codigo))

        # left join (base tiene todos los sujetos)
        participants = leftjoin(participants, clin2, on=:participant_id)
    end

    # --- limpiar tipos y faltantes ---
    # age a Float64, education_level a String con "n/a" si vacío
    if :age ∈ names(participants)
        participants.age = tryparse.(Float64, string.(participants.age))
    end
    if :education_level ∈ names(participants)
        participants.education_level = replace.(string.(participants.education_level), "" => "n/a")
    end

    # --- escribir salida ---
    outdir = args["out"]
    mkpath(outdir)

    # participants.tsv
    tsv_path = joinpath(outdir, "participants.tsv")
    CSV.write(tsv_path, participants; delim='\t', missingstring="n/a", quotechar='"')
    println("✔ Wrote $(tsv_path)  (rows=$(nrow(participants)), cols=$(ncol(participants)))")

    # participants.json (diccionario de columnas)
    dict = Dict(
        "participant_id" => Dict("Description"=>"Unique participant identifier (BIDS)", "Pattern"=>"^sub-[A-Za-z0-9]+\$"),
        "age"            => Dict("Description"=>"Age at first session", "Units"=>"years"),
        "sex"            => Dict("Description"=>"Biological sex", "Levels"=>Dict("M"=>"male","F"=>"female")),
        "group"          => Dict("Description"=>"Study group", "Levels"=>Dict("patient"=>"Multiple Sclerosis","control"=>"Healthy control")),
        "education_level"=> Dict("Description"=>"Education level (study-specific coding)", "Levels"=>Dict(
                                   "1"=>"primary","2"=>"secondary","3"=>"higher","n/a"=>"not available"))
    )
    # añadir descripciones clínicas si existen
    haskey_col = (c)->(c ∈ names(participants))
    if haskey_col(:disease_duration_years)
        dict["disease_duration_years"] = Dict("Description"=>"Disease duration since diagnosis", "Units"=>"years")
    end
    if haskey_col(:edss_t1)
        dict["edss_t1"] = Dict("Description"=>"EDSS at timepoint T1", "Units"=>"score")
    end
    if haskey_col(:treatment_t1)
        dict["treatment_t1"] = Dict("Description"=>"Treatment at timepoint T1")
    end
    if haskey_col(:diagnosis_t1)
        dict["diagnosis_t1"] = Dict("Description"=>"Diagnosis at timepoint T1")
    end
    if haskey_col(:edss_t2)
        dict["edss_t2"] = Dict("Description"=>"EDSS at timepoint T2", "Units"=>"score")
    end
    if haskey_col(:treatment_t2)
        dict["treatment_t2"] = Dict("Description"=>"Treatment at timepoint T2")
    end
    if haskey_col(:disability_t1)
        dict["disability_t1"] = Dict("Description"=>"Disability at timepoint T1")
    end
    if haskey_col(:disability_t2)
        dict["disability_t2"] = Dict("Description"=>"Disability at timepoint T2")
    end
    if haskey_col(:disability_improvement)
        dict["disability_improvement"] = Dict("Description"=>"Improvement in disability between T1 and T2")
    end

    json_path = joinpath(outdir, "participants.json")
    open(json_path, "w") do io
        JSON3.pretty(io, dict)
    end
    println("✔ Wrote $(json_path)")
    println("Done.")
end

# tryparse helper that returns missing on failure
function tryparse(::Type{T}, s) where {T}
    x = try Base.tryparse(T, s) catch; nothing end
    x === nothing ? missing : x
end

main()
main()