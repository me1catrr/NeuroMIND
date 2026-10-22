# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Documento maestro clínico + selección EEG
# ═══════════════════════════════════════════════════════════════
#
#  Une en un único documento (para compartir con Javier y Alejandro)
#  lo que hoy está repartido entre la base clínica MIND-EM (EDSS,
#  D-FIS, BDI-II numéricos), las categorías e indicaciones finales de
#  Javier y la selección de grabaciones EEG. No lee las fuentes
#  originales: parte de las salidas ya fusionadas de la Fase A.
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/build_clinical_master.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 15-10-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/build_clinical_master.jl [--check]
#
#  --check   tras escribir, cruza BIDS en disco, inventario, raw_selection.csv,
#            participants.tsv, el CSV/xlsx recién escritos, la base clínica original
#            y las cifras del capítulo 4 del informe; imprime las filas de los
#            Cuadros 4.1 y 4.2 y sale con error si algo no coincide.
#
#  Entradas (regenerar antes con audit_full_dataset.jl si cambian las fuentes)
#  ───────────────────────────────────────────────────────────────────────────
#    data/bids/participants.tsv                 clínica fusionada (una fila por sujeto)
#    data/tables/inventory.csv                  grabaciones EEG y decisiones de selección
#    data/sourcedata/clinical/anteriores/Javier_NeuroMIND.xlsx   solo observaciones de Javier (col. H/I)
#
#  Salidas (data/ está en .gitignore)
#  ──────────────────────────────────
#    data/sourcedata/clinical/datos_clinicos_unificados.csv    hoja «Participantes» (UTF-8, coma, punto decimal)
#    data/sourcedata/clinical/datos_clinicos_unificados.xlsx   Léeme · Participantes · Grabaciones EEG ·
#                                                 Control de calidad · Diccionario
#
#  La fecha de nacimiento no se incluye (dato identificativo; basta la edad).
#  El .xlsx se escribe sin dependencias nuevas: XML SpreadsheetML + `zip` del sistema.

using CSV, Dates, Statistics, XLSX

# ─── Rutas ────────────────────────────────────────────────────

const PROJ_ROOT    = dirname(@__DIR__)
const PARTIC_TSV   = joinpath(PROJ_ROOT, "data", "bids", "participants.tsv")
const INVENTORY    = joinpath(PROJ_ROOT, "data", "tables", "inventory.csv")
const CLINICAL_DIR = joinpath(PROJ_ROOT, "data", "sourcedata", "clinical", "anteriores")
const JAVIER_XLSX  = joinpath(CLINICAL_DIR, "Javier_NeuroMIND.xlsx")
const OUT_BASE     = joinpath(PROJ_ROOT, "data", "sourcedata", "clinical", "datos_clinicos_unificados")

# ─── Lectura ──────────────────────────────────────────────────

"participants.tsv → filas Dict (n/a → nothing)."
function read_participants()
    L = readlines(PARTIC_TSV)
    hdr = split(L[1], '\t')
    return [Dict(String(h) => (v == "n/a" ? nothing : String(v)) for (h, v) in zip(hdr, split(l, '\t')))
            for l in L[2:end]]
end

"sub-M05 → M5 · sub-MC09 → MC9 (código original del proyecto)."
function original_code(bids_id::AbstractString)
    m = match(r"^sub-(MC?)0*(\d+)$", bids_id)
    return m === nothing ? String(bids_id) : m[1] * m[2]
end

"Inventario de grabaciones EEG: Vector de NamedTuples (una por .vhdr)."
function read_inventory()
    s(x) = ismissing(x) ? nothing : String(x)
    f(x) = ismissing(x) ? nothing : parse(Float64, String(x))
    return [(code = String(r.subject_id), bids = String(r.bids_id), bses = String(r.bids_session), group = String(r.group), session = String(r.session),
             sub = s(r.session_sub), cond = String(r.condition), cond_raw = s(r.condition_raw),
             excluded = String(r.excluded) == "true", note = s(r.note),
             tmin = f(r.tmin_s), tmax = f(r.tmax_s), sexo = s(r.sexo),
             edad = ismissing(r.edad) ? 0 : parse(Int, String(r.edad)))
            for r in CSV.Rows(INVENTORY; types = String)]
end

"Observaciones de texto libre de Javier (Hoja1, columnas H e I), por código."
function read_javier_notes()
    notes = Dict{String,String}()
    XLSX.openxlsx(JAVIER_XLSX) do xf
        M = XLSX.getdata(xf["Hoja1"])
        for i in 2:size(M, 1)
            ismissing(M[i, 1]) && continue
            parts = String[]
            for j in (8, 9)                        # H = cerrados T1, I = T2
                v = M[i, j]
                (ismissing(v) || !(v isa AbstractString)) && continue
                t = strip(v)
                (isempty(t) || t == "FALTA" || occursin(r"^[\d.]+ s\s*-", t)) && continue
                push!(parts, "Javier (" * (j == 8 ? "T1 cerrados" : "T2") * "): " * t)
            end
            isempty(parts) || (notes[uppercase(strip(string(M[i, 1])))] = join(parts, "; "))
        end
    end
    return notes
end

# ─── Utilidades ───────────────────────────────────────────────

num(x) = x === nothing ? nothing : parse(Float64, x)
dat(x) = x === nothing ? nothing : Date(x)
delta(a, b) = (a === nothing || b === nothing) ? nothing : b - a
const CAT = Dict("mild" => "Leve", "severe" => "Grave", "high" => "Alta", "low" => "Baja",
                 "greater" => "Mayor", "lesser" => "Menor")
cat_es(x) = x === nothing ? nothing : CAT[x]
fmtn(x) = isinteger(x) ? string(Int(x)) : string(round(x; digits = 2))

# ─── Tabla de participantes ───────────────────────────────────

const COLUMNS = [   # (clave, anchura, descripción, fuente)
    ("codigo", 9, "Código original del proyecto (M = paciente EM, MC = control).", "Proyecto"),
    ("id_bids", 11, "Identificador BIDS (M5 → sub-M05); el mismo que en data/bids.", "Proyecto"),
    ("grupo", 10, "Paciente (EM) o Control.", "Base clínica MIND-EM / CSV demográfico"),
    ("sexo", 8, "Hombre / Mujer.", "Base clínica (pacientes); CSV demográfico (controles)"),
    ("edad_t1", 8, "Edad en años en la primera evaluación (T1).", "Base clínica (pacientes); CSV demográfico (controles)"),
    ("nivel_educativo", 10, "Nivel educativo, codificación ordinal 1–3. El significado de cada categoría no está documentado en las fuentes.", "CSV demográfico"),
    ("fecha_inclusion_t1", 13, "Fecha de inclusión / primera evaluación. Solo pacientes.", "Base clínica MIND-EM"),
    ("fecha_inclusion_t2", 13, "Fecha de la evaluación de seguimiento. Solo pacientes con T2.", "Base clínica MIND-EM"),
    ("intervalo_meses_t1_t2", 11, "Meses entre T1 y T2 (columna T.INTERVENCIÓN de la base).", "Base clínica MIND-EM"),
    ("diagnostico", 11, "EMRR = remitente-recurrente; EMPP = primaria progresiva.", "Base clínica MIND-EM"),
    ("tratamiento_t1", 14, "Tratamiento modificador de la enfermedad en T1 (nombres normalizados).", "Base clínica MIND-EM"),
    ("tratamiento_t2", 14, "Tratamiento modificador de la enfermedad en T2 (nombres normalizados).", "Base clínica MIND-EM"),
    ("edss_t1", 8, "EDSS (0–10) en T1. Valor numérico.", "Base clínica MIND-EM"),
    ("edss_t2", 8, "EDSS (0–10) en T2. Valor numérico.", "Base clínica MIND-EM"),
    ("dfis_t1", 8, "D-FIS (escala de impacto de la fatiga) en T1. También disponible para controles.", "Base clínica MIND-EM, hoja de valoración neuropsicológica"),
    ("dfis_t2", 8, "D-FIS en T2.", "Base clínica MIND-EM, hoja de valoración neuropsicológica"),
    ("bdi2_t1", 8, "BDI-II (depresión) en T1. También disponible para controles.", "Base clínica MIND-EM, hoja de valoración neuropsicológica"),
    ("bdi2_t2", 8, "BDI-II en T2.", "Base clínica MIND-EM, hoja de valoración neuropsicológica"),
    ("delta_edss", 9, "CALCULADO: edss_t2 − edss_t1.", "Calculado"),
    ("delta_dfis", 9, "CALCULADO: dfis_t2 − dfis_t1 (negativo = menos fatiga).", "Calculado"),
    ("edss_categoria", 11, "Leve / Grave (categoría de Javier, análisis transversal).", "Javier (Hoja1, col. C)"),
    ("dfis_cat_transversal", 12, "Alta / Baja según mediana de D-FIS T1 (categoría de Javier; usar en transversal).", "Javier (Hoja1, col. D)"),
    ("dfis_cat_longitudinal", 12, "Alta / Baja según mediana recalculada solo con sujetos T1+T2 (usar en longitudinal).", "Javier (Hoja1, col. E)"),
    ("respuesta_fatiga", 11, "Mayor / Menor: respuesta de la fatiga entre T1 y T2 (solo sujetos con T2).", "Javier (Hoja1, col. F)"),
    ("eeg_t1_ec", 18, "¿Hay EEG utilizable en T1, ojos cerrados? «Sí (segmento a–b s)» = recorte de un raw único.", "Inventario EEG / raw_selection.csv"),
    ("eeg_t1_eo", 18, "Idem, T1 ojos abiertos.", "Inventario EEG / raw_selection.csv"),
    ("eeg_t2_ec", 18, "Idem, T2 ojos cerrados.", "Inventario EEG / raw_selection.csv"),
    ("eeg_t2_eo", 18, "Idem, T2 ojos abiertos.", "Inventario EEG / raw_selection.csv"),
    ("observaciones", 80, "Indicaciones de Javier, decisiones sobre los raw y avisos de control de calidad a confirmar.", "Javier / Control de calidad"),
]
const KEYS = first.(COLUMNS)

function eeg_status(recs)
    ok = filter(r -> !r.excluded, recs)
    isempty(ok) && return "No"
    length(ok) > 1 && return "Sí (¡$(length(ok)) archivos sin resolver!)"
    r = only(ok)
    return r.tmin === nothing ? "Sí" : "Sí (segmento $(fmtn(r.tmin))–$(fmtn(r.tmax)) s)"
end

"Aplica las comprobaciones de coherencia; devuelve (filas del informe, avisos por código)."
function quality_checks(P, inv)
    rows = Vector{NamedTuple}()
    flags = Dict{String,Vector{String}}()
    add!(chk, res, who, det) = (push!(rows, (comprobacion = chk, resultado = res, participante = who, detalle = det));
                                res == "Revisar" && who != "" && push!(get!(flags, who, String[]), det))
    ms = filter(p -> p["group"] == "MS" && p["code"] != "", P)

    # 1) EDSS: la categoría debe separar por un corte único
    mild   = [p["edss_t1"] for p in ms if p["edss_group"] == "mild"   && p["edss_t1"] !== nothing]
    severe = [p["edss_t1"] for p in ms if p["edss_group"] == "severe" && p["edss_t1"] !== nothing]
    mx, mn = maximum(num.(mild)), minimum(num.(severe))
    if mx < mn
        add!("EDSS T1 vs categoría (Leve/Grave)", "OK", "", "Corte coherente: Leve ≤ $(fmtn(mx)), Grave ≥ $(fmtn(mn)) (N=$(length(mild) + length(severe)))")
    else
        for p in ms
            e = num(p["edss_t1"]); e === nothing && continue
            ((p["edss_group"] == "mild" && e >= mn) || (p["edss_group"] == "severe" && e <= mx)) &&
                add!("EDSS T1 vs categoría (Leve/Grave)", "Revisar", p["code"], "EDSS T1=$(fmtn(e)) pero categoría $(cat_es(p["edss_group"]))")
        end
    end

    # 2) D-FIS transversal: Alta ⇔ D-FIS T1 ≥ mediana de todos los pacientes con T1
    d1 = [num(p["dfis_t1"]) for p in ms if p["dfis_t1"] !== nothing]
    med1 = median(d1); bad = 0
    for p in ms
        v = num(p["dfis_t1"]); (v === nothing || p["dfis_group_t1"] === nothing) && continue
        if (v >= med1) != (p["dfis_group_t1"] == "high")
            bad += 1
            add!("D-FIS T1 vs categoría transversal", "Revisar", p["code"], "D-FIS T1=$(fmtn(v)) pero categoría $(cat_es(p["dfis_group_t1"])) (mediana=$(fmtn(med1)))")
        end
    end
    bad == 0 && add!("D-FIS T1 vs categoría transversal", "OK", "", "Coherente: Alta ⇔ D-FIS T1 ≥ mediana = $(fmtn(med1)) (N=$(length(d1)); los empates con la mediana cuentan como Alta)")

    # 3) D-FIS longitudinal: mediana de T1 solo en sujetos con categoría T1+T2
    lg = [p for p in ms if p["dfis_group_t1t2"] !== nothing && p["dfis_t1"] !== nothing]
    med12 = median(num.(getindex.(lg, "dfis_t1"))); bad = 0
    for p in lg
        v = num(p["dfis_t1"])
        if (v >= med12) != (p["dfis_group_t1t2"] == "high")
            bad += 1
            add!("D-FIS T1 vs categoría longitudinal", "Revisar", p["code"], "D-FIS T1=$(fmtn(v)) pero categoría longitudinal $(cat_es(p["dfis_group_t1t2"])) (mediana del subgrupo=$(fmtn(med12)))")
        end
    end
    bad == 0 && add!("D-FIS T1 vs categoría longitudinal", "OK", "", "Coherente: Alta ⇔ D-FIS T1 ≥ mediana del subgrupo = $(fmtn(med12)) (N=$(length(lg)))")

    # 4) Respuesta: debe ser monótona en ΔD-FIS; se busca el corte que menos discrepancias deja
    rs = [(p, delta(num(p["dfis_t1"]), num(p["dfis_t2"]))) for p in ms if p["fatigue_response"] !== nothing]
    rs = [(p, d) for (p, d) in rs if d !== nothing]
    cuts = sort(unique(last.(rs)))
    miss(t) = count(((p, d),) -> (d <= t) != (p["fatigue_response"] == "greater"), rs)
    t = cuts[argmin(miss.(cuts))]
    viol = [(p, d) for (p, d) in rs if (d <= t) != (p["fatigue_response"] == "greater")]
    for (p, d) in viol
        add!("Respuesta de fatiga vs ΔD-FIS", "Revisar", p["code"],
             "ΔD-FIS=$(fmtn(d)) (T1=$(p["dfis_t1"]), T2=$(p["dfis_t2"])) pero respuesta $(cat_es(p["fatigue_response"])); la regla que mejor encaja es Mayor ⇔ ΔD-FIS ≤ $(fmtn(t))")
    end
    isempty(viol) && add!("Respuesta de fatiga vs ΔD-FIS", "OK", "", "Coherente: Mayor ⇔ ΔD-FIS ≤ $(fmtn(t)) (N=$(length(rs)))")

    # 5) Misma pareja (T1,T2) de D-FIS con categorías distintas = contradicción directa
    for (nombre, key) in (("D-FIS longitudinal", "dfis_group_t1t2"), ("Respuesta de fatiga", "fatigue_response"))
        g = Dict{Tuple,Vector{Any}}()
        for p in ms
            (p[key] === nothing || p["dfis_t1"] === nothing || p["dfis_t2"] === nothing) && continue
            push!(get!(g, (p["dfis_t1"], p["dfis_t2"]), Any[]), p)
        end
        for ((a, b), ps) in g
            length(unique(p[key] for p in ps)) > 1 &&
                add!("Mismos valores D-FIS, distinta categoría", "Revisar", join(sort(getindex.(ps, "code")), " y "),
                     "$nombre: D-FIS T1=$a / T2=$b en todos ellos, pero categorías " * join(("$(p["code"])=$(cat_es(p[key]))" for p in ps), ", "))
        end
    end

    # 6) Edad y sexo entre fuentes
    byinv = Dict(r.code => r for r in inv)
    for p in P
        r = get(byinv, p["code"], nothing)
        a = num(p["age"])
        if r !== nothing && a !== nothing && r.edad > 0 && abs(a - r.edad) >= 1
            add!("Edad: base clínica vs CSV demográfico", "Revisar", p["code"], "Edad base=$(fmtn(a)) vs CSV demográfico=$(r.edad)")
        end
        b, i = dat(p["birth_date"]), dat(p["inclusion_date_t1"])
        if a !== nothing && b !== nothing && i !== nothing
            calc = floor(Int, Dates.value(i - b) / 365.2425)
            abs(calc - a) >= 1 && add!("Edad: columna EDAD vs fechas", "Revisar", p["code"], "Edad=$(fmtn(a)) vs $(calc) calculada con fechas de nacimiento e inclusión")
        end
        if r !== nothing && r.sexo !== nothing && p["sex"] !== nothing &&
           (p["sex"] == "M") != (r.sexo == "Hombre")
            add!("Sexo: base clínica vs CSV demográfico", "Revisar", p["code"], "Sexo base=$(p["sex"]) vs CSV=$(r.sexo)")
        end
    end

    # 7) Datos clínicos de T2 sin EEG T2 utilizable
    for p in P
        has_clin = p["inclusion_date_t2"] !== nothing || p["edss_t2"] !== nothing || p["dfis_t2"] !== nothing
        has_clin && p["has_t2"] == "false" &&
            add!("Datos clínicos T2 sin EEG T2", "Info", p["code"], "Hay EDSS/D-FIS T2 en la base, pero el EEG T2 no se analiza (decisión de Javier)")
    end

    # 8) Valores clínicos que faltan
    for p in ms, k in ("edss_t1", "dfis_t1", "bdi_t1", "edss_group", "dfis_group_t1")
        p[k] === nothing && add!("Valor clínico ausente en paciente", "Revisar", p["code"], "Falta $k")
    end
    return rows, flags
end

function build_rows()
    P = read_participants()
    foreach(p -> p["code"] = original_code(p["participant_id"]), P)
    inv  = read_inventory()
    jnot = read_javier_notes()
    qa, flags = quality_checks(P, inv)

    recs = Dict{Tuple{String,String,String},Vector}()
    foreach(r -> push!(get!(recs, (r.code, r.session, r.cond), []), r), inv)

    rows = Vector{Vector{Any}}()
    for p in P
        c = p["code"]
        e1, e2, d1, d2 = num(p["edss_t1"]), num(p["edss_t2"]), num(p["dfis_t1"]), num(p["dfis_t2"])
        obs = String[]
        haskey(jnot, c) && push!(obs, jnot[c])
        any(r -> r.code == c && r.session == "T1" && r.tmin !== nothing, inv) &&
            push!(obs, "T1: raw único dividido en abiertos/cerrados por segundos (indicación de Javier)")
        any(r -> r.code == c && !r.excluded && r.note !== nothing && occursin("sin sufijo de ojos", r.note), inv) &&
            push!(obs, "El archivo sin sufijo de ojos es su EO de T1 (confirmado por Alejandro, 28-09-2026)")
        haskey(flags, c) && push!(obs, "REVISAR: " * join(flags[c], " · "))
        push!(rows, Any[
            c, p["participant_id"], p["group"] == "MS" ? "Paciente" : "Control",
            p["sex"] == "M" ? "Hombre" : p["sex"] == "F" ? "Mujer" : nothing,
            num(p["age"]), num(p["education_level"]),
            dat(p["inclusion_date_t1"]), dat(p["inclusion_date_t2"]), num(p["months_t1_t2"]),
            p["diagnosis"], p["treatment_t1"], p["treatment_t2"],
            e1, e2, d1, d2, num(p["bdi_t1"]), num(p["bdi_t2"]),
            delta(e1, e2), delta(d1, d2),
            cat_es(p["edss_group"]), cat_es(p["dfis_group_t1"]), cat_es(p["dfis_group_t1t2"]), cat_es(p["fatigue_response"]),
            eeg_status(get(recs, (c, "T1", "EC"), [])), eeg_status(get(recs, (c, "T1", "EO"), [])),
            eeg_status(get(recs, (c, "T2", "EC"), [])), eeg_status(get(recs, (c, "T2", "EO"), [])),
            isempty(obs) ? nothing : join(obs, " | "),
        ])
    end
    return rows, qa, inv
end

"Hoja de grabaciones: una fila por .vhdr con su estado final (sin nombres de archivo originales)."
function recording_rows(inv)
    ord(r) = (r.group == "MS" ? 0 : 1, parse(Int, match(r"\d+", r.code).match), r.session, something(r.sub, ""), r.cond)
    out = Vector{Vector{Any}}()
    for r in sort(inv; by = ord)
        m = r.note === nothing ? nothing : match(r"^\[(\w+)\]\s*(.*)$", r.note)
        decision = m === nothing ? nothing : m[1]
        motivo   = m === nothing ? r.note : m[2]
        estado = r.excluded ? "Excluida" : (r.tmin !== nothing ? "Usada (segmento)" : "Usada")
        push!(out, Any[r.code, r.group == "MS" ? "Paciente" : "Control", r.session, r.sub, r.cond, r.cond_raw,
                       estado, decision, r.tmin, r.tmax, motivo])
    end
    return out
end
const REC_HEADER = ["codigo", "grupo", "sesion", "sub_sesion", "condicion", "condicion_original",
                    "estado", "decision", "tmin_s", "tmax_s", "motivo"]
const REC_WIDTHS = [9, 10, 8, 10, 10, 26, 18, 10, 9, 9, 70]

# ─── Escritura CSV ────────────────────────────────────────────

csvcell(x) = x === nothing ? "" : x isa Date ? string(x) : x isa Number ? fmtn(x) :
             occursin(r"[\",\n]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : String(x)

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        foreach(r -> println(io, join(csvcell.(r), ",")), rows)
    end
end

# ─── Escritura XLSX (SpreadsheetML mínimo) ────────────────────

xesc(s) = replace(replace(String(s), r"[\x00-\x08\x0b\x0c\x0e-\x1f]" => ""),
                  "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
function colname(n)
    s = ""
    while n > 0; n, r = divrem(n - 1, 26); s = string(Char('A' + r)) * s; end
    return s
end

# Estilos: 0 normal · 1 cabecera · 2 fecha · 3 texto con ajuste · 4 negrita con ajuste
const STYLES_XML = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><numFmts count="1"><numFmt numFmtId="164" formatCode="yyyy\\-mm\\-dd"/></numFmts><fonts count="3"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts><fills count="3"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill><fill><patternFill patternType="solid"><fgColor rgb="FF1F4E78"/><bgColor indexed="64"/></patternFill></fill></fills><borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders><cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs><cellXfs count="5"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1" applyAlignment="1"><alignment vertical="center" wrapText="1"/></xf><xf numFmtId="164" fontId="0" fillId="0" borderId="0" xfId="0" applyNumberFormat="1"/><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf><xf numFmtId="0" fontId="2" fillId="0" borderId="0" xfId="0" applyFont="1" applyAlignment="1"><alignment vertical="top" wrapText="1"/></xf></cellXfs><cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>"""

"""
Hoja: (nombre, filas, anchuras; header=true → 1.ª fila es cabecera con filtro y panel fijo;
`freezecols` columnas inmovilizadas; `colstyle[j]` estilo de las celdas de texto de la columna j).
"""
function sheet_xml(rows, widths; header = true, freezecols = 0, colstyle = Dict{Int,Int}(), grid = true)
    io = IOBuffer()
    nc = maximum(length, rows)
    print(io, """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">""")
    print(io, "<dimension ref=\"A1:$(colname(nc))$(length(rows))\"/><sheetViews><sheetView workbookViewId=\"0\"", grid ? "" : " showGridLines=\"0\"", ">")
    if header
        if freezecols > 0
            print(io, "<pane xSplit=\"$freezecols\" ySplit=\"1\" topLeftCell=\"$(colname(freezecols + 1))2\" activePane=\"bottomRight\" state=\"frozen\"/><selection pane=\"bottomRight\"/>")
        else
            print(io, "<pane ySplit=\"1\" topLeftCell=\"A2\" activePane=\"bottomLeft\" state=\"frozen\"/><selection pane=\"bottomLeft\"/>")
        end
    end
    print(io, "</sheetView></sheetViews><sheetFormatPr defaultRowHeight=\"15\"/><cols>")
    for (j, w) in enumerate(widths)
        print(io, "<col min=\"$j\" max=\"$j\" width=\"$w\" customWidth=\"1\"/>")
    end
    print(io, "</cols><sheetData>")
    for (i, row) in enumerate(rows)
        print(io, "<row r=\"$i\">")
        for (j, v) in enumerate(row)
            (v === nothing || (v isa AbstractString && isempty(v))) && continue
            ref = colname(j) * string(i)
            if v isa Date
                print(io, "<c r=\"$ref\" s=\"2\"><v>$(Dates.value(v - Date(1899, 12, 30)))</v></c>")
            elseif v isa Number
                print(io, "<c r=\"$ref\"><v>$(Float64(v))</v></c>")
            else
                st = (header && i == 1) ? 1 : get(colstyle, j, 0)
                print(io, "<c r=\"$ref\" t=\"inlineStr\"", st == 0 ? "" : " s=\"$st\"", "><is><t xml:space=\"preserve\">", xesc(v), "</t></is></c>")
            end
        end
        print(io, "</row>")
    end
    print(io, "</sheetData>")
    header && print(io, "<autoFilter ref=\"A1:$(colname(nc))$(length(rows))\"/>")
    print(io, "</worksheet>")
    return String(take!(io))
end

function write_xlsx(path, sheets)   # sheets: Vector{Pair{String,String}} (nombre => xml de hoja)
    mktempdir() do dir
        mkpath.(joinpath.(dir, ["_rels", "xl/_rels", "xl/worksheets"]))
        n = length(sheets)
        write(joinpath(dir, "[Content_Types].xml"), """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/><Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>""" *
            join(("<Override PartName=\"/xl/worksheets/sheet$i.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>" for i in 1:n)) * "</Types>")
        write(joinpath(dir, "_rels/.rels"), """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>""")
        write(joinpath(dir, "xl/workbook.xml"), """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>""" *
            join(("<sheet name=\"$(xesc(first(s)))\" sheetId=\"$i\" r:id=\"rId$i\"/>" for (i, s) in enumerate(sheets))) * "</sheets></workbook>")
        write(joinpath(dir, "xl/_rels/workbook.xml.rels"), """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">""" *
            join(("<Relationship Id=\"rId$i\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet$i.xml\"/>" for i in 1:n)) *
            "<Relationship Id=\"rId$(n + 1)\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles\" Target=\"styles.xml\"/></Relationships>")
        write(joinpath(dir, "xl/styles.xml"), STYLES_XML)
        for (i, s) in enumerate(sheets)
            write(joinpath(dir, "xl/worksheets/sheet$i.xml"), last(s))
        end
        rm(path; force = true)
        cd(dir) do
            run(pipeline(`zip -q -X -D -r $path .`; stdout = devnull))
        end
    end
    return path
end

# ─── Hojas de texto: Léeme y Diccionario ──────────────────────

function readme_rows(rows, qa, nrec)
    idx(k) = findfirst(==(k), KEYS)
    col(k) = [r[idx(k)] for r in rows]
    isyes(x) = x isa AbstractString && startswith(x, "Sí")
    pac = findall(==("Paciente"), col("grupo")); ctl = findall(==("Control"), col("grupo"))
    pair(c) = count(i -> isyes(col("eeg_t1_" * c)[i]) && isyes(col("eeg_t2_" * c)[i]), eachindex(rows))
    nrev = length(unique(q.participante for q in qa if q.resultado == "Revisar" && q.participante != ""))
    srcmod(f) = Dates.format(unix2datetime(mtime(f)), "yyyy-mm-dd")
    base = filter(f -> occursin(r"^BASE DE DATOS MIND-EM.*\.xlsx$", f), readdir(CLINICAL_DIR))
    return Vector{Any}[
        ["NeuroMIND — Datos clínicos y selección EEG unificados", nothing],
        ["Generado", string(today()) * " (scripts/build_clinical_master.jl)"],
        ["Para qué sirve", "Un único documento con los valores numéricos (EDSS, D-FIS, BDI-II), las categorías e indicaciones finales de Javier y la selección de grabaciones EEG, actualizado a la fecha de generación. El mismo contenido de la hoja «Participantes» está en el .csv."],
        ["Hojas", "Participantes: una fila por sujeto.\nGrabaciones EEG: una fila por archivo raw, con su estado final (usada / segmento / excluida) y el motivo.\nControl de calidad: comprobaciones de coherencia entre categorías y valores numéricos, y discrepancias entre fuentes.\nDiccionario: significado, unidades y fuente de cada columna."],
        ["Muestra", "$(length(rows)) sujetos: $(length(pac)) pacientes EM y $(length(ctl)) controles. $nrec grabaciones EEG registradas.\nParejas T1+T2 con EEG utilizable: $(pair("ec")) en ojos cerrados y $(pair("eo")) en ojos abiertos.\nSujetos con algún aviso por revisar: $nrev (ver Control de calidad)."],
        ["Fuentes", "• Base clínica MIND-EM ($(isempty(base) ? "no encontrada" : last(base)), modificada $(isempty(base) ? "—" : srcmod(joinpath(CLINICAL_DIR, last(base))))): fechas, sexo, edad, EDSS, tratamiento, diagnóstico, D-FIS y BDI-II de T1 y T2.\n• Javier_NeuroMIND.xlsx (= segmentos.xlsx, modificado $(srcmod(JAVIER_XLSX))): categorías EDSS / D-FIS / respuesta, segmentos de los raw únicos de T1 y observaciones.\n• CSV demográfico: nivel educativo y demografía de los controles.\n• Selección de raws (config/raw_selection.csv): decisiones de Javier aplicadas al inventario de grabaciones."],
        ["Cómo leerlo", "• Celda vacía = dato no disponible (p. ej. controles sin EDSS ni fechas; pacientes sin T2).\n• Los valores numéricos y las categorías son los de las fuentes, sin modificar. Las únicas columnas calculadas son delta_edss y delta_dfis.\n• Usar dfis_cat_transversal en el análisis transversal y dfis_cat_longitudinal en el longitudinal (la mediana cambia, así que un mismo sujeto puede ser Baja en uno y Alta en otro).\n• M4, M5 y M6 tienen un único raw en T1: se recortan por segundos en ojos abiertos y cerrados (columnas eeg_t1_*).\n• M31: hay D-FIS/EDSS de T2, pero su raw de T2 no vale y solo se analiza T1.\n• M37: de los 4 archivos de T2 se usan los T2_2; M16: en T1 cerrados se usa el que acaba en _ojoscerrados.\n• M7: el archivo T1 sin sufijo de ojos es su ojos abiertos (confirmado por Alejandro el 28-09-2026).\n• Los controles MC5, MC6, MC11 y MC12 no tienen datos."],
        ["Privacidad", "No se incluye la fecha de nacimiento; se mantienen la edad y las fechas de inclusión necesarias para el intervalo T1–T2. Los nombres de archivo originales (con iniciales) tampoco aparecen."],
    ]
end

# ─── Comprobaciones de consistencia (--check) ─────────────────

const BIDS_DIR      = joinpath(PROJ_ROOT, "data", "bids")
const RAW_SELECTION = joinpath(PROJ_ROOT, "config", "raw_selection.csv")
const REPORT_TEX    = [joinpath(PROJ_ROOT, "Report_Pre", "chapters", "04_diseno_adquisicion.tex"),
                       joinpath(PROJ_ROOT, "Report_Pre", "assets", "04_diseno_adquisicion", "report", "tikz", "bids_folder_tree.tikz"),
                       joinpath(PROJ_ROOT, "Report_Pre", "assets", "04_diseno_adquisicion", "report", "tikz", "bids_brainvision.tikz")]
const TASK_OF       = Dict("EC" => "eyesclosed", "EO" => "eyesopen")

"Valores numéricos de la base MIND-EM leídos celda a celda: código → (edss_t1, dfis_t1, bdi2_t1, …)."
function read_base_values()
    files = isdir(CLINICAL_DIR) ? filter(f -> occursin(r"^BASE DE DATOS MIND-EM.*\.xlsx$", f), readdir(CLINICAL_DIR)) : String[]
    isempty(files) && return nothing
    path = joinpath(CLINICAL_DIR, last(sort(files, by = f -> mtime(joinpath(CLINICAL_DIR, f)))))
    out = Dict{String,Dict{String,Float64}}()
    sid(x) = replace(uppercase(strip(string(x))), r"-T[12]$" => "")
    val(x) = ismissing(x) ? nothing : x isa Number ? Float64(x) : tryparse(Float64, replace(strip(string(x)), "," => "."))
    XLSX.openxlsx(path) do xf
        for n in XLSX.sheetnames(xf)
            u = strip(uppercase(n)); d = XLSX.getdata(xf[n])
            hdr = [lowercase(strip(string(coalesce(x, "")))) for x in d[1, :]]
            if occursin("ANTROPOM", u) && (endswith(u, "T1") || endswith(u, "T2"))
                tag = endswith(u, "T1") ? "t1" : "t2"
                ci = findfirst(h -> occursin("paciente", h), hdr); ce = findfirst(==("edss"), hdr)
                for i in 2:size(d, 1)
                    ismissing(d[i, ci]) && continue
                    e = val(d[i, ce]); e === nothing && continue
                    get!(out, sid(d[i, ci]), Dict{String,Float64}())["edss_$tag"] = e
                end
            elseif occursin("NEUROPSICOLOGICA", u)
                tag = (endswith(u, "T1") || endswith(u, "MC")) ? "t1" : "t2"
                cb = findfirst(h -> occursin("bdi", h), hdr); cf = findfirst(h -> occursin(r"d-?fis", h), hdr)
                for i in 2:size(d, 1)
                    ismissing(d[i, 1]) && continue
                    b, f = val(d[i, cb]), val(d[i, cf]); (b === nothing && f === nothing) && continue
                    r = get!(out, sid(d[i, 1]), Dict{String,Float64}())
                    b === nothing || (r["bdi2_$tag"] = b); f === nothing || (r["dfis_$tag"] = f)
                end
            end
        end
    end
    return out
end

squash(s) = replace(String(s), r"\s+" => "")
q25_75(v) = (quantile(v, 0.25), median(v), quantile(v, 0.75))   # cuantiles inclusivos, como el informe

"Filas del Cuadro 4.1 y 4.2 del capítulo 4, con el formato LaTeX del informe."
function report_tables(rows)
    idx = Dict(k => i for (i, k) in enumerate(KEYS))
    yes(x) = x isa AbstractString && startswith(x, "Sí")
    function t41(label, S; clin = true, t2 = true)
        c(f) = count(f, S)
        g(k) = [r[idx[k]] for r in S]
        ages = Float64.(g("edad_t1"))
        fmt(v) = lpad(string(round(v; digits = 1)), 0)
        edad = "$(fmt(mean(ages))) \$\\pm\$ $(fmt(std(ages))) ($(Int(minimum(ages)))--$(Int(maximum(ages))))"
        edu = join((c(r -> r[idx["nivel_educativo"]] == k) for k in 1:3), " / ")
        pair(a, b) = "$(c(r -> yes(r[idx[a]]) && (b === nothing || yes(r[idx[b]])))) / $(c(r -> yes(r[idx[replace(a, "_ec" => "_eo")]]) && (b === nothing || yes(r[idx[replace(b, "_ec" => "_eo")]]))))"
        dx   = clin ? "$(c(r -> r[idx["diagnostico"]] == "EMRR")) / $(c(r -> r[idx["diagnostico"]] == "EMPP"))" : "--"
        edss = clin ? "$(c(r -> r[idx["edss_categoria"]] == "Leve")) / $(c(r -> r[idx["edss_categoria"]] == "Grave"))" : "--"
        dfis = clin ? "$(c(r -> r[idx["dfis_cat_transversal"]] == "Alta")) / $(c(r -> r[idx["dfis_cat_transversal"]] == "Baja"))" : "--"
        clin2 = t2 ? string(c(r -> r[idx["edss_t2"]] !== nothing || r[idx["dfis_t2"]] !== nothing)) : "--"
        par = t2 ? pair("eeg_t1_ec", "eeg_t2_ec") : "--"
        "$label & $(length(S)) & $edad & $(c(r -> r[idx["sexo"]] == "Mujer")) / $(c(r -> r[idx["sexo"]] == "Hombre")) & $edu & $dx & $edss & $dfis & $clin2 & $(pair("eeg_t1_ec", nothing)) & $(pair("eeg_t2_ec", nothing)) & $par \\\\"
    end
    em = filter(r -> r[idx["grupo"]] == "Paciente", rows); hc = filter(r -> r[idx["grupo"]] == "Control", rows)
    l41 = [t41("EM / MS", em), t41("Control / HC", hc; clin = false, t2 = false),
           t41("Total", rows; clin = false)]
    function t42(label, k, hcflag)
        s(S, t) = (v = Float64[x for x in (r[idx["$(k)_$t"]] for r in S) if x !== nothing];
                   isempty(v) ? "--" : (q = q25_75(v); "$(fmtn(q[2])) [$(fmtn(q[1]))--$(fmtn(q[3]))] ($(length(v)))"))
        "$label & $(s(em, "t1")) & $(s(em, "t2")) & $(hcflag ? s(hc, "t1") : "--") \\\\"
    end
    l42 = [t42("EDSS", "edss", false), t42("D-FIS", "dfis", true), t42("BDI-II", "bdi2", true)]
    return l41, l42
end

"""
    run_checks(rows, inv) -> Int

Cruza BIDS en disco, inventario, `raw_selection.csv`, `participants.tsv`, el CSV y el
libro recién escritos, la base clínica original y las cifras del capítulo 4 del informe.
Devuelve el número de fallos.
"""
function run_checks(rows, inv)
    res = Tuple{Bool,String}[]
    ck(ok, msg) = push!(res, (ok, msg))
    idx = Dict(k => i for (i, k) in enumerate(KEYS))
    P = read_participants()
    valid = filter(r -> !r.excluded && haskey(TASK_OF, r.cond), inv)
    base(r) = "sub-$(r.bids)/ses-$(r.bses)/eeg/sub-$(r.bids)_ses-$(r.bses)_task-$(TASK_OF[r.cond])_run-01_eeg"

    # 1. BIDS en disco frente al inventario
    expected = Set(base(r) * "." * e for r in valid for e in ("eeg", "vhdr", "vmrk"))
    found = Set{String}(); found_json = Set{String}(); dirs = Set{Tuple{String,String}}()
    for sub in filter(d -> startswith(d, "sub-"), readdir(BIDS_DIR))
        for ses in filter(d -> startswith(d, "ses-"), readdir(joinpath(BIDS_DIR, sub)))
            push!(dirs, (sub[5:end], ses[5:end]))
            for f in readdir(joinpath(BIDS_DIR, sub, ses, "eeg"))
                p = "$sub/$ses/eeg/$f"
                occursin(r"_task-.+_eeg\.(eeg|vhdr|vmrk)$", f) && push!(found, p)
                occursin(r"_task-.+_eeg\.json$", f) && push!(found_json, p)
            end
        end
    end
    ck(found == expected, "BIDS: $(length(found) ÷ 3) registros en disco frente a $(length(expected) ÷ 3) utilizables del inventario (sobran $(length(setdiff(found, expected))), faltan $(length(setdiff(expected, found))))")
    crop = Set(base(r) * ".json" for r in valid if r.tmin !== nothing)
    ck(found_json == crop, "BIDS: $(length(found_json)) sidecars de recorte en disco frente a $(length(crop)) filas con ventana temporal")
    ck(dirs == Set((r.bids, r.bses) for r in valid), "BIDS: $(length(dirs)) carpetas ses-* frente a $(length(Set((r.bids, r.bses) for r in valid))) sesiones del inventario")

    # 2. participants.tsv, inventario y selección de grabaciones
    bad = [p["participant_id"] for p in P if (p["has_t1"] == "true") != ((p["participant_id"][5:end], "T1") in dirs) ||
                                              (p["has_t2"] == "true") != ((p["participant_id"][5:end], "T2") in dirs)]
    ck(isempty(bad), "participants.tsv: has_t1/has_t2 coherentes con las carpetas ses-* (discrepan: $(join(bad, ", ")))")
    ck(length(P) == length(rows) == length(Set(first.(dirs))), "$(length(P)) sujetos en participants.tsv, en el libro y en el BIDS ($(length(Set(first.(dirs)))))")
    sel = filter(l -> !isempty(strip(l)), readlines(RAW_SELECTION)[2:end])
    nnote = count(r -> r.note !== nothing && startswith(r.note, "["), inv)
    ck(length(sel) == nnote, "raw_selection.csv: $(length(sel)) decisiones frente a $nnote filas del inventario con decisión")

    # 3. Libro y CSV escritos frente al BIDS y entre sí
    miss = String[]
    for r in rows, (s, c) in (("T1", "EC"), ("T1", "EO"), ("T2", "EC"), ("T2", "EO"))
        v = r[idx["eeg_$(lowercase(s))_$(lowercase(c))"]]; id = r[idx["id_bids"]]
        isfile(joinpath(BIDS_DIR, id, "ses-$s", "eeg", "$(id)_ses-$(s)_task-$(TASK_OF[c])_run-01_eeg.eeg")) ==
            startswith(v, "Sí") || push!(miss, "$(r[1]) $s $c")
    end
    ck(isempty(miss), "libro: $(4 * length(rows)) celdas eeg_* coinciden con el BIDS (discrepan: $(join(miss, ", ")))")
    ck(readlines(OUT_BASE * ".csv") == vcat([join(KEYS, ",")], [join(csvcell.(r), ",") for r in rows]),
       "CSV en disco idéntico a la tabla en memoria")
    xl = XLSX.openxlsx(xf -> XLSX.getdata(xf["Participantes"]), OUT_BASE * ".xlsx")
    same(a, b) = a === nothing ? ismissing(b) : a isa Number ? (b isa Number && Float64(a) == Float64(b)) :
                 a isa Date ? (b isa Date && a == b) : (b isa AbstractString && a == b)
    ndiff = count(!same(rows[i][j], j <= size(xl, 2) ? xl[i + 1, j] : missing) for i in eachindex(rows), j in eachindex(KEYS))
    ck(size(xl, 1) == length(rows) + 1 && ndiff == 0, "xlsx: hoja Participantes idéntica a la tabla en memoria ($(size(xl, 1) - 1) filas, $ndiff celdas distintas)")

    # 4. Valores clínicos frente a la base original
    bv = read_base_values()
    if bv === nothing
        ck(false, "base clínica MIND-EM no encontrada en $CLINICAL_DIR")
    else
        n = 0; nbad = String[]
        for r in rows, k in ("edss_t1", "edss_t2", "dfis_t1", "dfis_t2", "bdi2_t1", "bdi2_t2")
            a = get(get(bv, r[idx["codigo"]], Dict{String,Float64}()), k, nothing); b = r[idx[k]]; n += 1
            (a === nothing ? b === nothing : (b !== nothing && Float64(a) == Float64(b))) || push!(nbad, "$(r[1]) $k")
        end
        ck(isempty(nbad), "valores EDSS/D-FIS/BDI-II: $n comparados con la base original (discrepan: $(join(first(nbad, 8), ", ")))")
    end

    # 5. Cifras del capítulo 4 del informe
    tex = join((isfile(f) ? squash(read(f, String)) : "" for f in REPORT_TEX[1:1]))
    l41, l42 = report_tables(rows)
    for l in vcat(l41, l42)
        ck(occursin(squash(l), tex), "informe, Cuadro 4.1/4.2: fila «" * first(split(l, " & ")) * "» coincide con el capítulo 4")
        occursin(squash(l), tex) || println("   esperada: ", l)
    end
    nv = length(valid); nec = count(r -> r.cond == "EC", valid); neo = count(r -> r.cond == "EO", valid)
    nses = length(dirs); nt2 = count(p -> p[2] == "T2", dirs)
    nsel_files = length(Set(first(split(l, ",")) for l in sel))
    phrases = [(1, "quedan $nv grabaciones válidas de $(length(P)) participantes"),
               (1, "$nses sesiones (T1 de los $(length(P)) participantes y T2 de $nt2 pacientes) y $nv grabaciones, $nec en ojos cerrados y $neo en ojos abiertos"),
               (1, "($(length(inv)) filas: $nv utilizables y $(length(inv) - nv) excluidas)"),
               (1, "$(length(sel)) filas sobre $nsel_files archivos")]
    for (i, ph) in phrases
        ck(occursin(squash(ph), tex), "informe, capítulo 4: «" * ph * "»")
    end
    for (f, ph) in ((REPORT_TEX[2], "(T1: $(length(P)) · T2: $nt2 pacientes)"), (REPORT_TEX[3], "$(length(P)) sujetos, $nv grabaciones"))
        ck(isfile(f) && occursin(squash(ph), squash(read(f, String))), "figura $(basename(f)): «$ph»")
    end

    println("\nComprobaciones de consistencia:")
    for (ok, msg) in res
        println("  ", ok ? "✓" : "✗", " ", msg)
    end
    nfail = count(!first, res)
    println(nfail == 0 ? "\n✅ $(length(res)) comprobaciones, sin fallos." : "\n❌ $nfail de $(length(res)) comprobaciones fallan.")
    return nfail
end

# ─── Principal ────────────────────────────────────────────────

function main()
    rows, qa, inv = build_rows()
    recs = recording_rows(inv)

    write_csv(OUT_BASE * ".csv", KEYS, rows)

    qa_sorted = sort(qa; by = q -> (q.resultado == "Revisar" ? 0 : q.resultado == "Info" ? 1 : 2, q.comprobacion, q.participante))
    qa_rows = vcat([["comprobacion", "resultado", "participante", "detalle"]],
                   [Any[q.comprobacion, q.resultado, q.participante, q.detalle] for q in qa_sorted])
    dic_rows = vcat([["columna", "descripcion", "fuente"]],
                    [Any[k, d, s] for (k, _, d, s) in COLUMNS],
                    [Any["(Grabaciones EEG) estado", "Usada / Usada (segmento) / Excluida.", "Inventario EEG"],
                     Any["(Grabaciones EEG) decision", "Acción de raw_selection.csv: use, segment, discard, exclude (vacío = regla por defecto).", "raw_selection.csv"],
                     Any["(Grabaciones EEG) tmin_s, tmax_s", "Segundos del raw usados cuando es un recorte (intervalo [tmin, tmax)).", "Javier (Hoja1, col. G y H)"],
                     Any["(Grabaciones EEG) condicion", "EC = ojos cerrados; EO = ojos abiertos; ODDBALL = paradigma distinto, no se analiza.", "Inventario EEG"]])

    sheets = [
        "Léeme" => sheet_xml(readme_rows(rows, qa, length(inv)), [26, 120]; header = false, grid = false,
                             colstyle = Dict(1 => 4, 2 => 3)),
        "Participantes" => sheet_xml(vcat([collect(KEYS)], rows), [w for (_, w, _, _) in COLUMNS];
                                     freezecols = 2, colstyle = Dict(length(KEYS) => 3)),
        "Grabaciones EEG" => sheet_xml(vcat([REC_HEADER], recs), REC_WIDTHS; freezecols = 1, colstyle = Dict(11 => 3)),
        "Control de calidad" => sheet_xml(qa_rows, [38, 11, 16, 120]; colstyle = Dict(4 => 3)),
        "Diccionario" => sheet_xml(dic_rows, [30, 100, 50]; colstyle = Dict(2 => 3, 3 => 3)),
    ]
    write_xlsx(OUT_BASE * ".xlsx", sheets)

    println("✓ ", OUT_BASE, ".csv  ($(length(rows)) sujetos × $(length(KEYS)) columnas)")
    println("✓ ", OUT_BASE, ".xlsx ($(length(recs)) grabaciones, $(length(qa)) filas de control de calidad)")
    println("\nControl de calidad:")
    for q in qa_sorted
        q.resultado == "OK" || println("  [", q.resultado, "] ", q.comprobacion, q.participante == "" ? "" : " · " * q.participante, " — ", q.detalle)
    end
    foreach(q -> q.resultado == "OK" && println("  [OK] ", q.comprobacion, " — ", q.detalle), qa_sorted)

    if "--check" in ARGS
        l41, l42 = report_tables(rows)
        println("\nFilas para los Cuadros 4.1 y 4.2 del informe (formato LaTeX):")
        foreach(l -> println("  ", l), vcat(l41, l42))
        run_checks(rows, inv) == 0 || exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
