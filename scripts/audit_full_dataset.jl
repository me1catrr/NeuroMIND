# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Auditoría del dataset MINDEM-IMIBIC
# ═══════════════════════════════════════════════════════════════
#
#  Fase A — Escanea los .vhdr brutos, normaliza condiciones
#  (>15 variantes ortográficas), aplica las decisiones de
#  config/raw_selection.csv, fusiona los datos clínicos y genera
#  el inventario y participants.tsv.
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/audit_full_dataset.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 15-10-2026
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
#  Entradas (data/sourcedata/, solo lectura)
#  ─────────────────────────────────────────
#    brainvision/                         tripletes .vhdr/.vmrk/.eeg originales
#    clinical/anteriores/BASE DE DATOS MIND-EM*.xlsx EDSS, D-FIS, BDI-II, fechas, tratamiento
#                                         (si hay varias versiones, la más reciente)
#    clinical/anteriores/Javier_NeuroMIND.xlsx       categorías: EDSS leve/grave, D-FIS alta/baja,
#                                         respuesta (T1 → transversal, T1+T2 → longitudinal)
#    clinical/anteriores/demográficos … .csv         nivel educativo y demografía de controles
#
#  Salida
#  ──────
#    data/tables/inventory.csv            una fila por grabación (incluye excluidas y motivo)
#    data/bids/participants.tsv           una fila por participante con grabación válida
#
#  Precedencia de fuentes en participants.tsv: Excel MIND-EM para pacientes
#  (sexo, edad, fechas, escalas); CSV demográfico para controles y para el
#  nivel educativo; Javier para las categorías. Las discrepancias entre
#  fuentes se imprimen al final.

using Dates, XLSX

# ─── Rutas ────────────────────────────────────────────────────

const PROJ_ROOT    = dirname(@__DIR__)
const SOURCE_DIR   = joinpath(PROJ_ROOT, "data", "sourcedata")
const RAW_DATA_DIR = joinpath(SOURCE_DIR, "brainvision")
const CLINICAL_DIR = joinpath(SOURCE_DIR, "clinical", "anteriores")
const DEMOG_CSV    = joinpath(CLINICAL_DIR, "demográficos Rafael MIND EM_pac vs cont t1.csv")
const JAVIER_XLSX  = joinpath(CLINICAL_DIR, "Javier_NeuroMIND.xlsx")
const BIDS_DIR     = joinpath(PROJ_ROOT, "data", "bids")
const TABLES_DIR   = joinpath(PROJ_ROOT, "data", "tables")

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
    ("M6",  "T2"),   # no tiene T2
])

function is_excluded(subj::String, sess::String, cond::String)::Bool
    subj in EXCLUDED_SUBJECTS && return true
    (subj, sess) in EXCLUDED_SUBJ_SESS && return true
    cond == "ODDBALL" && return true
    cond == "UNKNOWN" && return true
    return false
end

# ─── Selección de raws (decisiones de Javier) ─────────────────
#
#  config/raw_selection.csv es la fuente única de las decisiones sobre
#  qué archivo raw usar. Se aplica por nombre de archivo (sin .vhdr):
#    segment  → el raw único contiene abiertos y cerrados: genera una fila
#               por condición con la ventana tmin_s–tmax_s (no excluida)
#    use      → fuerza incluir el archivo
#    discard  → duplicado descartado (excluded = true)
#    exclude  → grabación no válida (excluded = true)
#  (`use` sobre un archivo sin sufijo de ojos, clasificado ODDBALL por nombre, lo reetiqueta
#  con la condición indicada: es el caso de M7 T1, cuyo archivo sin sufijo es el EO.)
#  Los archivos con decisión sustituyen a las reglas por defecto de arriba
#  (p. ej. el ODDBALL por nombre en M4/M5/M6, o T2_2 = duplicado en M37).
#
#  Privacidad: los nombres originales llevan las iniciales del participante
#  y la fecha de grabación (M4_T1_<iniciales>_<ddmmaa>_COG-10-1). El CSV está
#  versionado en un repositorio público, así que ese tramo se sustituye por
#  «*» (M4_T1_*_COG-10-1) y la comparación se hace con `_anon_stem`.

const RAW_SELECTION_CSV = joinpath(dirname(@__DIR__), "config", "raw_selection.csv")

"""
    _anon_stem(stem) -> String

Sustituye `<iniciales>_<fecha>` por `*` en el nombre de un raw:
`M12_T2_2_ABC_010122_Ojoscerrados` → `M12_T2_2_*_Ojoscerrados` (ejemplo ficticio).
"""
_anon_stem(stem::AbstractString) =
    replace(String(stem), r"^(MC?\d+_T\d(?:_\d)?)_[^_]+_\d{6,8}_" => s"\1_*_")

struct SelectionEntry
    file_stem :: String
    subject_id:: String
    session   :: String
    condition :: String
    action    :: String
    tmin_s    :: Union{Nothing,Float64}
    tmax_s    :: Union{Nothing,Float64}
    reason    :: String
end

function load_raw_selection(path::String = RAW_SELECTION_CSV)::Dict{String,Vector{SelectionEntry}}
    sel = Dict{String,Vector{SelectionEntry}}()
    isfile(path) || (@warn "raw_selection.csv no encontrado: $path"; return sel)
    for (i, line) in enumerate(readlines(path)[2:end])
        isempty(strip(line)) && continue
        p = split(line, ",")
        length(p) >= 9 || error("raw_selection.csv línea $(i + 1): se esperaban 9 campos, hay $(length(p))")
        action = strip(String(p[5]))
        action in ("segment", "use", "discard", "exclude") ||
            error("raw_selection.csv línea $(i + 1): acción desconocida '$action'")
        tmin = tryparse(Float64, strip(String(p[6])))
        tmax = tryparse(Float64, strip(String(p[7])))
        if action == "segment"
            (tmin !== nothing && tmax !== nothing && tmax > tmin >= 0) ||
                error("raw_selection.csv línea $(i + 1): 'segment' requiere 0 ≤ tmin_s < tmax_s")
        end
        stem = strip(String(p[1]))
        push!(get!(sel, stem, SelectionEntry[]),
              SelectionEntry(stem, strip(String(p[2])), strip(String(p[3])),
                             strip(String(p[4])), action, tmin, tmax,
                             strip(String(p[8]))))
    end
    return sel
end

_fmt_s(x::Union{Nothing,Float64}) = x === nothing ? "" : string(x)

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

# ─── Datos clínicos: Excel MIND-EM + categorías de Javier ─────

_norm(s) = lowercase(strip(string(coalesce(s, ""))))

function _num(x)::Union{Nothing,Float64}
    (x === nothing || ismissing(x)) && return nothing
    x isa Number && return Float64(x)
    return tryparse(Float64, replace(strip(string(x)), "," => "."))
end

function _date(x)::Union{Nothing,Date}
    (x === nothing || ismissing(x)) && return nothing
    x isa DateTime && return Date(x)
    x isa Date && return x
    v = _num(x)
    return v === nothing ? nothing : Date(1899, 12, 30) + Day(round(Int, v))
end

# "M1-T2" → "M1"; las hojas mezclan "M42" y "M42-T1"
_sid(x) = (x === nothing || ismissing(x)) ? "" : replace(uppercase(strip(string(x))), r"-T[12]$" => "")

# Hoja → (cabeceras normalizadas, filas)
function _table(M)
    hdr  = [_norm(M[1, j]) for j in 1:size(M, 2)]
    body = [[M[i, j] for j in 1:size(M, 2)] for i in 2:size(M, 1)]
    return hdr, body
end
_col(hdr, pat::Regex) = findfirst(h -> occursin(pat, h), hdr)

function _neuro_kind(name::AbstractString)
    u = strip(uppercase(name))
    occursin("NEUROPSICOLOGICA", u) || return nothing
    endswith(u, "MC") && return :mc
    endswith(u, "T1") && return :t1
    (endswith(u, "T2") || endswith(u, " T")) && return :t2   # la hoja de T2 se llama «… NEUROPSICOLOGICA T»
    return nothing
end

const TREATMENT_FIX = Dict("natallizumab" => "natalizumab")
function _treatment(x)
    (x === nothing || ismissing(x)) && return nothing
    s = lowercase(strip(string(x)))
    isempty(s) && return nothing
    return uppercasefirst(get(TREATMENT_FIX, s, s))
end

"""
    load_clinical() -> Dict{String,Dict{String,Any}}

Lee la base clínica MIND-EM más reciente de data/sourcedata/clinical/anteriores/ y
devuelve, por código de sujeto original (M5, MC9…), sexo, edad, fechas,
EDSS, tratamiento, diagnóstico, D-FIS y BDI-II de T1 y T2.
"""
function load_clinical()::Dict{String,Dict{String,Any}}
    clin = Dict{String,Dict{String,Any}}()
    files = isdir(CLINICAL_DIR) ?
        filter(f -> occursin(r"^BASE DE DATOS MIND-EM.*\.xlsx$", f), readdir(CLINICAL_DIR)) : String[]
    if isempty(files)
        @warn "Sin base clínica MIND-EM en $CLINICAL_DIR; participants.tsv saldrá sin variables clínicas"
        return clin
    end
    sort!(files, by = f -> mtime(joinpath(CLINICAL_DIR, f)))
    length(files) > 1 && @warn "Varias versiones de la base clínica; se usa la más reciente: $(last(files))"
    path = joinpath(CLINICAL_DIR, last(files))
    println("✓ Base clínica: $(last(files))")

    rec(id) = get!(clin, id, Dict{String,Any}())
    XLSX.openxlsx(path) do xf
        names = XLSX.sheetnames(xf)
        sheet(pred) = (i = findfirst(pred, names);
                       i === nothing ? error("Hoja no encontrada en $(basename(path))") : xf[names[i]])

        for (suf, tag) in (("T1", "t1"), ("T2", "t2"))
            ws = sheet(n -> occursin("ANTROPOM", uppercase(n)) && endswith(strip(uppercase(n)), suf))
            hdr, body = _table(XLSX.getdata(ws))
            ci, cinc, cbd = _col(hdr, r"paciente"), _col(hdr, r"inclus"), _col(hdr, r"nacim")
            csex, cage, cedss = _col(hdr, r"^sexo$"), _col(hdr, r"^edad$"), _col(hdr, r"^edss$")
            ctr, cdx, cint = _col(hdr, r"^trat"), _col(hdr, r"diagn"), _col(hdr, r"intervenci")
            for r in body
                id = _sid(r[ci])
                isempty(id) && continue
                (_date(r[cinc]) === nothing && _num(r[cedss]) === nothing) && continue   # fila vacía (p. ej. sin T2)
                d = rec(id)
                d["inclusion_date_$tag"] = _date(r[cinc])
                d["edss_$tag"]           = _num(r[cedss])
                d["treatment_$tag"]      = ctr === nothing ? nothing : _treatment(r[ctr])
                d["birth_date"]  = something(get(d, "birth_date", nothing), cbd === nothing ? nothing : _date(r[cbd]), Some(nothing))
                if tag == "t1"
                    d["sex"]       = r[csex] === missing ? nothing : strip(string(r[csex]))
                    d["age"]       = _num(r[cage])
                    d["diagnosis"] = cdx === nothing || r[cdx] === missing ? nothing : strip(string(r[cdx]))
                elseif cint !== nothing
                    d["months_t1_t2"] = _num(r[cint])
                end
            end
        end

        for n in names
            kind = _neuro_kind(n)
            kind === nothing && continue
            tag = kind == :t2 ? "t2" : "t1"          # controles: única sesión, se guarda como t1
            hdr, body = _table(XLSX.getdata(xf[n]))
            ci, cb, cf = _col(hdr, r"participante"), _col(hdr, r"bdi"), _col(hdr, r"d-?fis")
            for r in body
                id = _sid(r[ci])
                isempty(id) && continue
                (_num(r[cb]) === nothing && _num(r[cf]) === nothing) && continue
                d = rec(id)
                d["bdi_$tag"]  = _num(r[cb])
                d["dfis_$tag"] = _num(r[cf])
            end
        end
    end
    return clin
end

"""
    load_javier() -> Dict{String,Dict{String,String}}

Categorías de Javier (Hoja1 de Javier_NeuroMIND.xlsx): EDSS leve/grave,
D-FIS alta/baja con la mediana de T1 (transversal) o de T1+T2 (longitudinal)
y respuesta mayor/menor. Se traducen a las etiquetas de participants.json.
"""
function load_javier()::Dict{String,Dict{String,String}}
    out = Dict{String,Dict{String,String}}()
    if !isfile(JAVIER_XLSX)
        @warn "No encontrado: $JAVIER_XLSX; participants.tsv saldrá sin categorías"
        return out
    end
    lvl = Dict("LEVE" => "mild", "GRAVE" => "severe", "ALTA" => "high", "BAJA" => "low",
               "MAYOR" => "greater", "MENOR" => "lesser")
    XLSX.openxlsx(JAVIER_XLSX) do xf
        hdr, body = _table(XLSX.getdata(xf["Hoja1"]))
        cs, ce = _col(hdr, r"^sujeto"), _col(hdr, r"^edss$")
        cd1, cd12, cr = _col(hdr, r"^dfis t1$"), _col(hdr, r"^dfis t1\+t2$"), _col(hdr, r"^respuesta")
        for r in body
            id = _sid(r[cs])
            isempty(id) && continue
            tr(c) = (c === nothing || ismissing(r[c])) ? nothing : get(lvl, uppercase(strip(string(r[c]))), nothing)
            d = Dict{String,String}()
            for (k, c) in (("edss_group", ce), ("dfis_group_t1", cd1),
                           ("dfis_group_t1t2", cd12), ("fatigue_response", cr))
                v = tr(c)
                v === nothing || (d[k] = v)
            end
            out[id] = d
        end
    end
    return out
end

_fmt(x) = x === nothing || ismissing(x) ? "n/a" :
          x isa AbstractFloat ? (isinteger(x) ? string(Int(x)) : string(round(x; digits=2))) : string(x)

const PARTICIPANT_COLS = [
    "participant_id", "group", "sex", "age", "birth_date", "education_level", "has_t1", "has_t2",
    "inclusion_date_t1", "inclusion_date_t2", "months_t1_t2", "diagnosis", "treatment_t1", "treatment_t2",
    "edss_t1", "edss_t2", "dfis_t1", "dfis_t2", "bdi_t1", "bdi_t2",
    "edss_group", "dfis_group_t1", "dfis_group_t1t2", "fatigue_response",
]

"""
    write_participants(rows, demog)

Escribe data/bids/participants.tsv (una fila por participante con alguna
grabación válida) fusionando inventario, CSV demográfico, base clínica y
categorías de Javier, y lista las discrepancias entre fuentes.
"""
function write_participants(rows, demog)
    clin = load_clinical()
    jav  = load_javier()
    mkpath(BIDS_DIR)
    seen  = Set{String}()
    notes = String[]
    n_clin = Dict("MS" => 0, "HC" => 0); n_tot = Dict("MS" => 0, "HC" => 0)

    path = joinpath(BIDS_DIR, "participants.tsv")
    open(path, "w") do io
        println(io, join(PARTICIPANT_COLS, "\t"))
        for r in sort(rows, by = r -> r.bids_id)
            (r.excluded || r.bids_id in seen) && continue
            push!(seen, r.bids_id)
            sid = r.subject_id
            c = get(clin, sid, Dict{String,Any}())
            j = get(jav, sid, Dict{String,String}())
            g(k) = get(c, k, nothing)

            # sexo y edad: Excel MIND-EM (pacientes) > CSV demográfico
            sex_x = g("sex")
            sex   = sex_x !== nothing ? sex_x : r.sexo
            sex_b = sex == "Hombre" ? "M" : (sex == "Mujer" ? "F" : "n/a")
            age_x = g("age")
            age   = age_x !== nothing ? age_x : (r.edad > 0 ? Float64(r.edad) : nothing)

            sex_x !== nothing && r.sexo != "n/a" && sex_x != r.sexo &&
                push!(notes, "$sid: sexo Excel=$(sex_x) vs CSV=$(r.sexo)")
            age_x !== nothing && r.edad > 0 && abs(age_x - r.edad) >= 1 &&
                push!(notes, "$sid: edad Excel=$(_fmt(age_x)) vs CSV=$(r.edad)")
            if age_x !== nothing && g("birth_date") !== nothing && g("inclusion_date_t1") !== nothing
                calc = floor(Int, Dates.value(g("inclusion_date_t1") - g("birth_date")) / 365.2425)
                abs(calc - age_x) >= 1 &&
                    push!(notes, "$sid: edad EDAD=$(_fmt(age_x)) vs fechas nacimiento/inclusión=$(calc)")
            end
            dm = get(demog, sid, nothing)
            dm === nothing && push!(notes, "$sid: sin fila en el CSV demográfico (grupo/educación/edad sin fuente)")
            isempty(c) && r.group == "MS" && push!(notes, "$sid: paciente sin datos en la base clínica MIND-EM")
            isempty(j) && r.group == "MS" && push!(notes, "$sid: paciente sin categorías de Javier")

            has_t1 = any(x -> x.bids_id == r.bids_id && x.session == "T1" && !x.excluded, rows)
            has_t2 = any(x -> x.bids_id == r.bids_id && x.session == "T2" && !x.excluded, rows)
            haskey(n_tot, r.group) && (n_tot[r.group] += 1; isempty(c) || (n_clin[r.group] += 1))

            vals = Dict{String,Any}(
                "participant_id" => "sub-$(r.bids_id)", "group" => r.group, "sex" => sex_b,
                "age" => age, "birth_date" => g("birth_date"),
                "education_level" => r.nivel_edu > 0 ? r.nivel_edu : nothing,
                "has_t1" => has_t1, "has_t2" => has_t2,
            )
            for k in ("inclusion_date_t1", "inclusion_date_t2", "months_t1_t2", "diagnosis", "treatment_t1",
                      "treatment_t2", "edss_t1", "edss_t2", "dfis_t1", "dfis_t2", "bdi_t1", "bdi_t2")
                vals[k] = g(k)
            end
            for k in ("edss_group", "dfis_group_t1", "dfis_group_t1t2", "fatigue_response")
                vals[k] = get(j, k, nothing)
            end
            println(io, join((_fmt(vals[k]) for k in PARTICIPANT_COLS), "\t"))
        end
    end
    println("✓ participants.tsv guardado: $path ($(length(seen)) participantes)")
    for grp in ("MS", "HC")
        println("  $(grp): $(n_tot[grp]) participantes · $(n_clin[grp]) con datos en la base clínica")
    end
    if !isempty(notes)
        println("\n⚠ Discrepancias / huecos entre fuentes ($(length(notes))):")
        foreach(n -> println("  · ", n), first(notes, 40))
        length(notes) > 40 && println("  … y $(length(notes) - 40) más")
    end
    return path
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
    selection = load_raw_selection()
    used_stems = Set{String}()
    println("✓ Decisiones de raw_selection.csv: $(sum(length, values(selection); init=0)) filas / $(length(selection)) archivos")

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

        # Notas específicas (reglas por defecto; raw_selection.csv las sustituye)
        if !isempty(info.session_sub)
            note = "T2_$(info.session_sub) duplicado: sin decision en raw_selection.csv; por defecto T2_1 primario"
            # Si es T2_2, marcar como duplicado
            info.session_sub == "2" && (excl = true; note = "DUPLICADO T2: sin decision en raw_selection.csv; por defecto usar T2_1")
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

        # Decisiones de raw_selection.csv: 0 entradas = reglas por defecto,
        # 1 = use/discard/exclude, 2 = segment (una fila por condición).
        stem = _anon_stem(splitext(info.filename)[1])
        push!(used_stems, stem)
        entries = get(selection, stem, SelectionEntry[])
        variants = isempty(entries) ? [nothing] : entries

        for e in variants
            cond_r, excl_r, note_r = cond, excl, note
            tmin_r = nothing
            tmax_r = nothing
            if e !== nothing
                if e.subject_id != subj || e.session != sess
                    error("raw_selection.csv: $(stem) declarado como $(e.subject_id)/$(e.session) pero el archivo es $(subj)/$(sess)")
                end
                if e.action == "segment"
                    cond_r, excl_r = e.condition, false
                    tmin_r, tmax_r = e.tmin_s, e.tmax_s
                else
                    if e.action == "use" && cond == "ODDBALL" && e.condition in ("EC", "EO")
                        # archivo sin sufijo de ojos que el equipo clínico confirma como reposo
                        # (p. ej. M7 T1: es el EO); la regla por nombre lo habría marcado oddball
                        cond_r = e.condition
                    else
                        e.condition == cond || error("raw_selection.csv: $(stem) declarado como $(e.condition) pero el archivo es $(cond)")
                    end
                    excl_r = e.action in ("discard", "exclude")
                end
                note_r = "[$(e.action)] $(e.reason)"
            end

            push!(rows, (
                filename      = info.filename,
                subject_id    = subj,
                bids_id       = bids_id,
                group         = grupo,
                session       = sess,
                session_sub   = info.session_sub,
                bids_session  = bids_sess,
                condition_raw = info.condition_raw,
                condition     = cond_r,
                date_str      = info.date_str,
                initials      = info.initials,
                excluded      = excl_r,
                sexo          = sexo,
                edad          = edad,
                nivel_edu     = edu,
                filepath      = relpath(fpath, PROJ_ROOT),   # relativa: portable entre Macs
                note          = note_r,
                tmin_s        = tmin_r,
                tmax_s        = tmax_r,
            ))
        end
    end

    # Toda decisión de raw_selection.csv debe corresponder a un archivo real
    stale = setdiff(keys(selection), used_stems)
    isempty(stale) || error("raw_selection.csv: archivos sin .vhdr correspondiente: $(join(sort(collect(stale)), ", "))")

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
    mkpath(TABLES_DIR)
    inv_path = joinpath(TABLES_DIR, "inventory.csv")
    open(inv_path, "w") do io
        println(io, "filename,subject_id,bids_id,group,session,session_sub," *
                    "bids_session,condition_raw,condition,date_str,initials," *
                    "excluded,sexo,edad,nivel_edu,filepath,note,tmin_s,tmax_s")
        for r in sort(rows, by = r -> (r.subject_id, r.session, r.condition))
            ex_str = r.excluded ? "true" : "false"
            # tmin_s/tmax_s van al final (cols 18-19) para no desplazar las
            # columnas que ya leen run_batch_pipeline.jl y build_bids_full.jl.
            # `note` no debe contener comas.
            line = join([
                r.filename, r.subject_id, r.bids_id, r.group,
                r.session, r.session_sub, r.bids_session,
                r.condition_raw, r.condition, r.date_str, r.initials,
                ex_str, r.sexo, string(r.edad), string(r.nivel_edu),
                r.filepath, replace(r.note, "," => ";"),
                _fmt_s(r.tmin_s), _fmt_s(r.tmax_s)
            ], ",")
            println(io, line)
        end
    end
    println("\n✓ Inventario guardado: $inv_path")

    # ── Generar participants.tsv (fusión de fuentes clínicas) ──
    # groups.csv y longitudinal_pairs.csv ya no se generan: los análisis de
    # grupo calculan cohorte y pares T1/T2 desde este inventario.
    write_participants(rows, demog)

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
