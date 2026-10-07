# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Auditoría del dataset MINDEM-IMIBIC
# ═══════════════════════════════════════════════════════════════
#
#  Fase A — Escanea los .vhdr brutos, normaliza condiciones
#  (>15 variantes ortográficas), detecta pares T1/T2 y genera
#  los metadatos necesarios para el pipeline.
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/audit_full_dataset.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 22-07-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/<este-script>.jl …
#  Sin shebang: #!/usr/bin/env julia no activaría --project=.
#
#  Uso
#  ───
#    julia --project=. scripts/audit_full_dataset.jl
#    julia --project=. scripts/audit_full_dataset.jl --data ruta/raw
#
#  Opciones
#  ────────
#    --data  PATH   carpeta de .vhdr brutos (defecto interno del script)
#
#  Salida
#  ──────
#    data/full_data/inventory.csv
#    data/bids/participants.tsv
#    data/bids/groups.csv
#    data/bids/longitudinal_pairs.csv

using Dates

# ─── Rutas ────────────────────────────────────────────────────

const RAW_DATA_DIR = joinpath(
    dirname(@__DIR__),
    "data", "full_data",
    "Pacientes MINDEM_IMIBIC_27 03 25"
)
const DEMOG_CSV = joinpath(
    dirname(@__DIR__),
    "data", "full_data",
    "demográficos Rafael MIND EM_pac vs cont t1.csv"
)
const BIDS_DIR      = joinpath(dirname(@__DIR__), "data", "bids")
const FULL_DATA_DIR = joinpath(dirname(@__DIR__), "data", "full_data")

# ─── Normalización de nombres de condición ────────────────────

"""
    normalize_condition(raw_cond) -> "EC" | "EO" | "ODDBALL" | "UNKNOWN"

Mapea >30 variantes ortográficas a una de 4 etiquetas canónicas.
La normalización es insensible a mayúsculas, acentos, espacios y guiones.
"""
function normalize_condition(raw::String)::String
    # Quitar prefijo de paradigma cognitivo si lo hay
    s = raw
    s = replace(s, r"(?i)cog-?10-?1-?" => "")
    s = replace(s, r"(?i)Cog-?10-?1[-\s]*" => "")

    # Normalización: minúsculas, quitar acentos y caracteres no alfanuméricos
    s = lowercase(s)
    s = replace(s, "á" => "a", "é" => "e", "í" => "i", "ó" => "o", "ú" => "u")
    s = replace(s, r"[\s\-_]" => "")

    # Patrones EC (ojos cerrados)
    ec_patterns = [
        "ojoscerrados", "ojoscerrados",   # correcto
        "ojoscerradps",   # typo: "ps" en lugar de "dos"
        "ojoscerradpsb",  # typo + "b" extra
        "ojoscerreados",  # typo: "reados"
        "ojoscerrads",    # typo: falta "o" final
        "ojoscerrado",    # typo: falta "s"
        "ojoscerredos",   # typo
        "ojoscerradoss",  # typo: doble "s"
        "ojoscerradso",   # typo: orden
        "ojocerrados",    # typo: falta "s" en "ojos"
        "ojoscerraados",  # typo: doble "a"
        "ojosCerrados",   # camelCase → ya normalizado
        "ojoscerr",       # abreviado
    ]

    # Patrones EO (ojos abiertos)
    eo_patterns = [
        "ojosabiertos",   # correcto
        "ojosaciertos",   # typo: "aciertos"
        "ojosabioertos",  # typo: "io" invertido
        "ojosabiertes",   # typo: "es" en lugar de "os"
        "ojosabierto",    # typo: falta "s"
        "ojosabierttos",  # typo: doble "t"
        "ojoscabiertos",  # typo: "c" extra
        "ojosabiertos",   # correcto (repetido por seguridad)
        "ojosbiertos",    # typo: falta "a"
        "ojosabieto",     # typo: falta "r"
        "ojosabietos",    # typo: falta "r" (M9_T2 real)
        "ojoabiert",      # abreviado
    ]

    if any(startswith(s, p) || s == p for p in ec_patterns) ||
       occursin("cerrado", s)
        return "EC"
    elseif any(startswith(s, p) || s == p for p in eo_patterns) ||
           occursin("abierto", s) || occursin("acierto", s)
        return "EO"
    elseif isempty(s) || s == "cog101" || occursin("cog", s)
        # Solo paradigma cognitivo sin condición explícita
        return "ODDBALL"
    else
        return "UNKNOWN"
    end
end

# ─── Parseo de nombres de archivo ─────────────────────────────

"""
    parse_vhdr_filename(filename) -> NamedTuple o nothing

Extrae subject_id, session, session_sub (para T2_1/T2_2), initials,
date, condition_raw, group ("MS"|"HC") de un nombre de archivo .vhdr.
"""
function parse_vhdr_filename(fname::String)
    b = basename(fname)
    endswith(b, ".vhdr") || return nothing

    name_noext = b[1:end-5]

    # ── Patrón MS: M{n}_T{sess}[_{sub}]_{init}_{date}_{cond}
    #    M37_T2_1_JLP_270722_Ojosabiertos
    m = match(r"^(M\d+)_T(\d+)(?:_(\d+))?_([^_]+)_(\d+)_(.+)$", name_noext)
    if m !== nothing
        subj_raw  = String(m[1])
        sess_n    = String(m[2])
        sess_sub  = m[3] === nothing ? "" : String(m[3])
        cond_raw  = String(m[6])
        return (
            subject_id   = subj_raw,
            group        = "MS",
            session      = "T" * sess_n,
            session_sub  = sess_sub,       # "" | "1" | "2" ...
            initials     = String(m[4]),
            date_str     = String(m[5]),
            condition_raw = cond_raw,
            condition    = normalize_condition(cond_raw),
            filename     = b,
            filepath     = fname,
        )
    end

    # ── Patrón early MS (M4-M7 formato COG): M{n}_T{sess}_{init}_{date}_{task}
    #    donde task puede ser "COG-10-1" sin condición de ojos
    m2 = match(r"^(M\d+)_T(\d+)_([^_]+)_(\d+)_(COG-10-1)$"i, name_noext)
    if m2 !== nothing
        return (
            subject_id   = String(m2[1]),
            group        = "MS",
            session      = "T" * String(m2[2]),
            session_sub  = "",
            initials     = String(m2[3]),
            date_str     = String(m2[4]),
            condition_raw = String(m2[5]),
            condition    = "ODDBALL",
            filename     = b,
            filepath     = fname,
        )
    end

    # ── Patrón Control: MC{n}_{init}_{date}_{cond}
    m3 = match(r"^(MC\d+)_([^_]+)_(\d+)_(.+)$", name_noext)
    if m3 !== nothing
        cond_raw = String(m3[4])
        return (
            subject_id   = String(m3[1]),
            group        = "HC",
            session      = "T1",       # controles: sesión única → T1
            session_sub  = "",
            initials     = String(m3[2]),
            date_str     = String(m3[3]),
            condition_raw = cond_raw,
            condition    = normalize_condition(cond_raw),
            filename     = b,
            filepath     = fname,
        )
    end

    return nothing
end

# ─── Sujetos excluidos ────────────────────────────────────────

const EXCLUDED_SUBJECTS = Set(["MC5", "MC6", "MC11"])
const EXCLUDED_SUBJ_SESS = Set([
    ("M3",  "T2"),   # en carpeta EXCLUIDOS
    ("M4",  "T1"),   # solo Oddball
    ("M5",  "T1"),   # solo Oddball
    ("M6",  "T1"),   # solo Oddball
    ("M6",  "T2"),   # no tiene T2
])

function is_excluded(subj::String, sess::String, cond::String)::Bool
    subj in EXCLUDED_SUBJECTS && return true
    (subj, sess) in EXCLUDED_SUBJ_SESS && return true
    # M4 y M6 nunca tienen resting state
    subj == "M4" && return true
    subj == "M6" && return true
    cond == "ODDBALL" && return true
    cond == "UNKNOWN" && return true
    return false
end

# ─── Leer demografías ─────────────────────────────────────────

struct SubjectDemog
    codigo::String
    tipo::String        # PAC / CON
    sexo::String        # Hombre / Mujer
    edad::Int
    nivel_edu::Int      # 1-3
end

function load_demographics()::Dict{String,SubjectDemog}
    d = Dict{String,SubjectDemog}()
    isfile(DEMOG_CSV) || (@warn "Demografías no encontradas: $DEMOG_CSV"; return d)
    for line in readlines(DEMOG_CSV)[2:end]   # saltar cabecera
        isempty(strip(line)) && continue
        parts = split(line, ";")
        length(parts) < 5 && continue
        codigo = strip(String(parts[1]))
        tipo   = strip(String(parts[2]))
        sexo   = strip(String(parts[3]))
        edad   = tryparse(Int, strip(String(parts[4])))
        edu    = tryparse(Int, strip(String(parts[5])))
        isnothing(edad) && continue
        d[codigo] = SubjectDemog(
            codigo, tipo, sexo,
            something(edad, 0),
            something(edu, 0)
        )
    end
    return d
end

# ─── Función principal ────────────────────────────────────────

function main()
    println("=" ^ 60)
    println(" NeuroMIND — Fase A: Auditoría del Dataset")
    println(" $(now())")
    println("=" ^ 60)

    isdir(RAW_DATA_DIR) || error("Directorio de datos no encontrado: $RAW_DATA_DIR")

    demog = load_demographics()
    println("✓ Demografías cargadas: $(length(demog)) sujetos")

    # ── Escanear .vhdr ────────────────────────────────────────
    vhdr_files = filter(f -> endswith(f, ".vhdr"),
                        readdir(RAW_DATA_DIR; join=true))
    println("✓ Archivos .vhdr encontrados: $(length(vhdr_files))")

    # ── Parsear y clasificar ──────────────────────────────────
    rows = []
    n_skipped = 0

    for fpath in sort(vhdr_files)
        info = parse_vhdr_filename(fpath)
        if info === nothing
            @warn "No se pudo parsear: $(basename(fpath))"
            n_skipped += 1
            continue
        end

        subj  = info.subject_id
        sess  = info.session
        cond  = info.condition
        excl  = is_excluded(subj, sess, cond)
        note  = ""

        # Notas específicas
        if !isempty(info.session_sub)
            note = "T2_$(info.session_sub) duplicado: usar T2_1 como primario"
            # Si es T2_2, marcar como duplicado
            info.session_sub == "2" && (excl = true; note = "DUPLICADO T2: usar T2_1")
        end
        if info.condition == "ODDBALL"
            note = "Paradigma Oddball — no es resting state"
        end

        # ID de sujeto BIDS (zero-padded)
        bids_id = _bids_subject_id(subj)
        # Session BIDS
        bids_sess = info.session_sub == "1" ? sess : sess   # T2_1 → ses-T2

        # Demografías
        dem = get(demog, subj, nothing)
        sexo  = dem === nothing ? "n/a" : dem.sexo
        edad  = dem === nothing ? 0     : dem.edad
        edu   = dem === nothing ? 0     : dem.nivel_edu
        grupo = dem === nothing ? info.group : (dem.tipo == "PAC" ? "MS" : "HC")

        push!(rows, (
            filename      = info.filename,
            subject_id    = subj,
            bids_id       = bids_id,
            group         = grupo,
            session       = sess,
            session_sub   = info.session_sub,
            bids_session  = bids_sess,
            condition_raw = info.condition_raw,
            condition     = cond,
            date_str      = info.date_str,
            initials      = info.initials,
            excluded      = excl,
            sexo          = sexo,
            edad          = edad,
            nivel_edu     = edu,
            filepath      = fpath,
            note          = note,
        ))
    end

    # ── Estadísticas ──────────────────────────────────────────
    valid  = filter(r -> !r.excluded, rows)
    ms_t1  = filter(r -> !r.excluded && r.group == "MS" && r.session == "T1", rows)
    ms_t2  = filter(r -> !r.excluded && r.group == "MS" && r.session == "T2", rows)
    hc     = filter(r -> !r.excluded && r.group == "HC", rows)

    println("\n── Resumen ──────────────────────────────────────────────")
    println("  Total archivos parseados : $(length(rows))")
    println("  Total excluidos          : $(count(r -> r.excluded, rows))")
    println("  Válidos para pipeline    : $(length(valid))")
    println("  MS T1 (EC+EO)           : $(length(ms_t1))")
    println("  MS T2 (EC+EO)           : $(length(ms_t2))")
    println("  HC sesión única (EC+EO) : $(length(hc))")

    # Condiciones normalizadas con problemas
    unknown = filter(r -> r.condition == "UNKNOWN", rows)
    if !isempty(unknown)
        println("\n⚠ Condiciones no reconocidas ($(length(unknown))):")
        for r in unknown
            println("  ⚠ $(r.filename)  → '$(r.condition_raw)'")
        end
    end

    # ── Guardar inventory.csv ─────────────────────────────────
    mkpath(FULL_DATA_DIR)
    inv_path = joinpath(FULL_DATA_DIR, "inventory.csv")
    open(inv_path, "w") do io
        println(io, "filename,subject_id,bids_id,group,session,session_sub," *
                    "bids_session,condition_raw,condition,date_str,initials," *
                    "excluded,sexo,edad,nivel_edu,filepath,note")
        for r in sort(rows, by = r -> (r.subject_id, r.session, r.condition))
            ex_str = r.excluded ? "true" : "false"
            line = join([
                r.filename, r.subject_id, r.bids_id, r.group,
                r.session, r.session_sub, r.bids_session,
                r.condition_raw, r.condition, r.date_str, r.initials,
                ex_str, r.sexo, string(r.edad), string(r.nivel_edu),
                r.filepath, r.note
            ], ",")
            println(io, line)
        end
    end
    println("\n✓ Inventario guardado: $inv_path")

    # ── Generar participants.tsv ──────────────────────────────
    mkpath(BIDS_DIR)
    seen_subj = Set{String}()
    parts_path = joinpath(BIDS_DIR, "participants.tsv")
    open(parts_path, "w") do io
        println(io, "participant_id\tgroup\tsex\tage\teducation_level\thas_t1\thas_t2")
        for r in sort(rows, by = r -> r.bids_id)
            r.excluded && continue
            r.bids_id in seen_subj && continue
            push!(seen_subj, r.bids_id)
            has_t1 = any(x -> x.bids_id == r.bids_id && x.session == "T1" && !x.excluded, rows)
            has_t2 = any(x -> x.bids_id == r.bids_id && x.session == "T2" && !x.excluded, rows)
            sex_bids = r.sexo == "Hombre" ? "M" : (r.sexo == "Mujer" ? "F" : "n/a")
            println(io, "sub-$(r.bids_id)\t$(r.group)\t$(sex_bids)\t$(r.edad)\t$(r.nivel_edu)\t$(has_t1)\t$(has_t2)")
        end
    end
    println("✓ participants.tsv guardado: $parts_path")

    # ── Generar groups.csv ────────────────────────────────────
    groups_path = joinpath(BIDS_DIR, "groups.csv")
    open(groups_path, "w") do io
        println(io, "subject_id,group,session,condition,bids_id")
        for r in sort(rows, by = r -> (r.bids_id, r.session, r.condition))
            r.excluded && continue
            println(io, "$(r.subject_id),$(r.group),$(r.session),$(r.condition),$(r.bids_id)")
        end
    end
    println("✓ groups.csv guardado: $groups_path")

    # ── Generar longitudinal_pairs.csv ────────────────────────
    # Solo pacientes MS con EC y EO en ambos T1 y T2
    ms_subjects = unique(r.subject_id for r in rows if r.group == "MS" && !r.excluded)
    long_pairs_path = joinpath(BIDS_DIR, "longitudinal_pairs.csv")
    open(long_pairs_path, "w") do io
        println(io, "subject_id,bids_id,has_t1_ec,has_t1_eo,has_t2_ec,has_t2_eo,include_longitudinal")
        for subj in sort(ms_subjects)
            subj_rows = filter(r -> r.subject_id == subj && !r.excluded, rows)
            bids_id = isempty(subj_rows) ? _bids_subject_id(subj) : first(subj_rows).bids_id
            t1_ec = any(r -> r.session == "T1" && r.condition == "EC", subj_rows)
            t1_eo = any(r -> r.session == "T1" && r.condition == "EO", subj_rows)
            t2_ec = any(r -> r.session == "T2" && r.condition == "EC", subj_rows)
            t2_eo = any(r -> r.session == "T2" && r.condition == "EO", subj_rows)
            include = t1_ec && t1_eo && t2_ec && t2_eo
            println(io, "$(subj),$(bids_id),$(t1_ec),$(t1_eo),$(t2_ec),$(t2_eo),$(include)")
        end
    end
    println("✓ longitudinal_pairs.csv guardado: $long_pairs_path")

    n_long = count(begin
        subj_rows = filter(r -> r.subject_id == subj && !r.excluded, rows)
        t1_ec = any(r -> r.session == "T1" && r.condition == "EC", subj_rows)
        t1_eo = any(r -> r.session == "T1" && r.condition == "EO", subj_rows)
        t2_ec = any(r -> r.session == "T2" && r.condition == "EC", subj_rows)
        t2_eo = any(r -> r.session == "T2" && r.condition == "EO", subj_rows)
        t1_ec && t1_eo && t2_ec && t2_eo
    end for subj in ms_subjects)
    println("  Pares longitudinales completos (T1+T2 EC+EO): $n_long")

    println("\n✅ Fase A completada: $(Dates.format(now(), "HH:MM:SS"))")
    return inv_path
end

# ─── Helpers ──────────────────────────────────────────────────

"""
    _bids_subject_id(raw_id) -> String

Convierte "M5" → "M05", "MC9" → "MC09", "M44" → "M44", etc.
Zero-padding hasta 2 dígitos del número.
"""
function _bids_subject_id(raw::String)::String
    m = match(r"^(M?C?)(\d+)$", raw)
    if m !== nothing
        prefix = String(m[1])
        n      = String(m[2])
        padded = length(n) < 2 ? lpad(n, 2, '0') : n
        return prefix * padded
    end
    return raw   # fallback: sin cambios
end

# ─── Punto de entrada ─────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
