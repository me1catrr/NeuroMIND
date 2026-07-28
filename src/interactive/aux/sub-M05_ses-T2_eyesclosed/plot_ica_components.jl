#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Panel de revisión ICA (componentes)
# ═══════════════════════════════════════════════════════════════
#
#  Panel de revisión: componente ACTIVO (topo/PSD/features) vs
#  componentes COMPARADOS (activaciones) vs DECISIÓN
#  (Conservar / Rechazar / Revisar).
#
#  Entrada:
#    ../../tables/ica/ica_components.csv
#    ../../tables/ica/ica_component_features.csv
#    ../../tables/ica/ica_activations.csv
#    ../../json/ica_summary.json
#    ../../figures/ica/ica_topomap_NNN.png
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_ica_components.jl
#    # → http://127.0.0.1:8770/  (Ctrl+C para salir)
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_ica_components.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      24-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using CSV, DataFrames, CairoMakie, Sockets, Dates, FFTW, Statistics

const HERE    = @__DIR__
const ROOT    = normpath(joinpath(HERE, "..", "..", "..", "..",
                "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))   # …/eyesclosed
const TABLES  = joinpath(ROOT, "tables", "ica")
const FIGS    = joinpath(ROOT, "figures", "ica")
const JSON_P  = joinpath(ROOT, "json", "ica_summary.json")
const HOST    = "127.0.0.1"
const PORT    = 8770

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK    = "eyesclosed"

const TYPE_COLOR = Dict(
    "brain"      => "#2563eb",
    "eye/blink"  => "#7c3aed",
    "ocular"     => "#7c3aed",
    "muscle"     => "#ea580c",
    "line_noise" => "#dc2626",
    "jump"       => "#ca8a04",
    "unknown"    => "#64748b",
)

# Bandas NeuroMIND (sombreado PSD)
const BANDS = [
    ("δ", 0.5, 4.0, "#94a3b8"),
    ("θ", 4.0, 8.0, "#64748b"),
    ("α", 7.8, 11.7, "#2563eb"),
    ("β", 12.0, 30.0, "#0f766e"),
    ("γ", 30.0, 50.0, "#7c3aed"),
]

# ── Store ──────────────────────────────────────────────────────

mutable struct ICACompStore
    components::DataFrame
    features::DataFrame
    t_act::Vector{Float64}
    act::Dict{Int,Vector{Float64}}
    summary::Dict{String,Any}
    n_comp::Int
    t_min::Float64
    t_max::Float64
    fs::Float64
end

function _json_num(s::String, key::String, default)
    m = match(Regex("\"" * key * "\"\\s*:\\s*([-\\d.eE]+)"), s)
    m === nothing && return default
    v = tryparse(Float64, m.captures[1])
    return v === nothing ? default : v
end
function _json_int(s::String, key::String, default::Int)
    Int(round(_json_num(s, key, Float64(default))))
end
function _json_bool(s::String, key::String, default::Bool)
    m = match(Regex("\"" * key * "\"\\s*:\\s*(true|false)"), s)
    m === nothing && return default
    return m.captures[1] == "true"
end
function _json_str(s::String, key::String, default::String)
    m = match(Regex("\"" * key * "\"\\s*:\\s*\"([^\"]*)\""), s)
    m === nothing && return default
    return String(m.captures[1])
end

function load_summary(path::String)::Dict{String,Any}
    d = Dict{String,Any}(
        "n_components" => 0, "n_rejected" => 0, "n_accepted" => 0,
        "variance_retained" => NaN, "rejection_source" => "?",
        "artifact_threshold" => 1.5, "converged" => false,
    )
    isfile(path) || return d
    s = read(path, String)
    d["n_components"] = _json_int(s, "n_components", 0)
    d["n_rejected"] = _json_int(s, "n_rejected", 0)
    d["n_accepted"] = _json_int(s, "n_accepted", 0)
    d["variance_retained"] = _json_num(s, "variance_retained", NaN)
    d["rejection_source"] = _json_str(s, "rejection_source", "?")
    d["artifact_threshold"] = _json_num(s, "artifact_threshold", 1.5)
    d["converged"] = _json_bool(s, "converged", false)
    d["n_iter"] = _json_int(s, "n_iter", 0)
    return d
end

function load_store()::ICACompStore
    comp_p = joinpath(TABLES, "ica_components.csv")
    feat_p = joinpath(TABLES, "ica_component_features.csv")
    act_p  = joinpath(TABLES, "ica_activations.csv")
    isfile(comp_p) || error("No encontrado: $comp_p")
    isfile(act_p)  || error("No encontrado: $act_p")

    components = CSV.read(comp_p, DataFrame)
    features = isfile(feat_p) ? CSV.read(feat_p, DataFrame) : DataFrame()

    adf = CSV.read(act_p, DataFrame)
    t_act = Float64.(adf.t_s)
    act = Dict{Int,Vector{Float64}}()
    for c in 1:nrow(components)
        lab = "IC$c"
        hasproperty(adf, Symbol(lab)) || continue
        act[c] = Float64.(adf[!, Symbol(lab)])
    end

    summary = load_summary(JSON_P)
    n_comp = nrow(components)
    t_min = isempty(t_act) ? 0.0 : minimum(t_act)
    t_max = isempty(t_act) ? 10.0 : maximum(t_act)
    fs = length(t_act) ≥ 2 ? 1.0 / (t_act[2] - t_act[1]) : 500.0

    return ICACompStore(components, features, t_act, act, summary, n_comp, t_min, t_max, fs)
end

function topomap_path(comp::Int)::String
    joinpath(FIGS, "ica_topomap_" * lpad(comp, 3, '0') * ".png")
end

# ── PSD periodograma (misma fórmula que _bandpower) ────────────

function compute_ic_psd(store::ICACompStore, comp::Int; fmax::Float64 = 100.0)
    haskey(store.act, comp) || error("IC$comp sin activación")
    x = store.act[comp] .- mean(store.act[comp])
    N = length(x)
    N < 8 && error("Serie demasiado corta para PSD")
    X = rfft(x)
    freqs = collect((0:length(X)-1) .* (store.fs / N))
    psd = abs.(X) .^ 2 ./ (store.fs * N)
    keep = freqs .<= fmax
    f = freqs[keep]
    p = psd[keep]
    # dB (relativo a 1; floor para evitar -Inf)
    psd_db = 10 .* log10.(max.(p, 1e-30))
    # pico en 1–45 Hz (evita DC)
    band = (f .>= 1.0) .& (f .<= 45.0)
    peak_hz = if any(band)
        fb = f[band]; pb = p[band]
        fb[argmax(pb)]
    else
        NaN
    end
    return (f = f, psd = p, psd_db = psd_db, peak_hz = peak_hz, df = store.fs / N)
end

# ── PNG (CairoMakie): activaciones ─────────────────────────────

function save_activation_png(
    store::ICACompStore,
    comps::Vector{Int},
    t0::Float64,
    t1::Float64,
)::String
    isempty(comps) && error("Selecciona al menos un IC")
    mask = (store.t_act .>= t0) .& (store.t_act .<= t1)
    count(mask) < 2 && error("Ventana vacía")
    t = store.t_act[mask]

    fig = Figure(size = (1100, 420), fontsize = 13)
    ax = Axis(fig[1, 1];
        title = "Activaciones ICA — $SUBJECT / $SESSION / $TASK",
        xlabel = "Tiempo (s)", ylabel = "a.u.",
    )
    ax.xgridvisible = true; ax.ygridvisible = true

    for c in comps
        haskey(store.act, c) || continue
        row = store.components[store.components.component .== c, :]
        typ = nrow(row) > 0 ? String(row.artifact_type[1]) : "unknown"
        rej = nrow(row) > 0 && Bool(row.rejected[1])
        y = store.act[c][mask]
        lab = "IC$c ($typ)" * (rej ? " ✕" : "")
        lines!(ax, t, y; linewidth = rej ? 1.6 : 1.2, label = lab,
               color = typ == "brain" ? :steelblue :
                       typ == "line_noise" ? :red :
                       typ == "jump" ? :goldenrod :
                       (typ == "ocular" || typ == "eye/blink") ? :purple :
                       typ == "muscle" ? :orange : :gray)
    end
    axislegend(ax; position = :rt, framevisible = true, labelsize = 11)
    xlims!(ax, t0, t1)

    tag = join(["IC$c" for c in comps], "_")
    out = joinpath(HERE, "ica_act_" * tag * "_$(round(Int,t0))-$(round(Int,t1))s.png")
    save(out, fig; px_per_unit = 2)
    return out
end

# ── HTTP ───────────────────────────────────────────────────────

function _read_headers(sock)::Dict{String,String}
    headers = Dict{String,String}()
    while true
        line = readline(sock)
        (line == "" || line == "\r") && break
        m = match(r"^([^:]+):\s*(.+?)(?:\r)?$", line)
        m === nothing && continue
        headers[lowercase(m.captures[1])] = strip(m.captures[2])
    end
    return headers
end

function _read_body(sock, headers)::String
    n = parse(Int, get(headers, "content-length", "0"))
    n <= 0 && return ""
    buf = Vector{UInt8}(undef, n)
    filled = 0
    while filled < n
        nb = readbytes!(sock, view(buf, filled+1:n))
        nb == 0 && break
        filled += nb
    end
    return String(view(buf, 1:filled))
end

function _send(sock, status::Int, body::String; content_type::String = "text/plain; charset=utf-8")
    write(sock, "HTTP/1.1 $status\r\nContent-Type: $content_type\r\n" *
                "Content-Length: $(sizeof(body))\r\nConnection: close\r\n" *
                "Access-Control-Allow-Origin: *\r\n\r\n" * body)
end

function _send_bytes(sock, status::Int, body::Vector{UInt8}; content_type::String)
    write(sock, "HTTP/1.1 $status\r\nContent-Type: $content_type\r\n" *
                "Content-Length: $(length(body))\r\nConnection: close\r\n" *
                "Access-Control-Allow-Origin: *\r\n\r\n")
    write(sock, body)
end

function _parse_qs(path::AbstractString)::Dict{String,String}
    qs = Dict{String,String}()
    occursin('?', path) || return qs
    for pair in split(split(path, '?', limit=2)[2], '&')
        kv = split(pair, '=', limit=2)
        length(kv) == 2 || continue
        qs[String(kv[1])] = String(replace(replace(kv[2], "%2C" => ","), "%20" => " "))
    end
    return qs
end

function _parse_save_json(s::String)
    comps = Int[]
    m = match(Regex("\"comps\"\\s*:\\s*\\[(.*?)\\]", "s"), s)
    if m !== nothing
        for cap in eachmatch(r"(\d+)", m.captures[1])
            push!(comps, parse(Int, cap.captures[1]))
        end
    end
    _num = (key, default) -> begin
        m2 = match(Regex("\"" * key * "\"\\s*:\\s*([-\\d.eE]+)"), s)
        m2 === nothing && return Float64(default)
        something(tryparse(Float64, m2.captures[1]), Float64(default))
    end
    return (comps = comps, t0 = _num("t0", 0.0), t1 = _num("t1", 10.0))
end

function _components_json(store::ICACompStore)::String
    parts = String[]
    for r in eachrow(store.components)
        c = Int(r.component)
        typ = String(r.artifact_type)
        decision = Bool(r.rejected) ? "reject" : "keep"
        push!(parts,
            "{\"comp\":$c,\"label\":\"$(r.label)\",\"variance_pct\":$(r.variance_pct)," *
            "\"rejected\":$(r.rejected ? "true" : "false"),\"artifact_type\":\"$typ\"," *
            "\"decision\":\"$decision\"," *
            "\"has_topomap\":$(isfile(topomap_path(c)) ? "true" : "false")}")
    end
    return "[" * join(parts, ",") * "]"
end

function _feature_json(store::ICACompStore, comp::Int)::String
    nrow(store.features) == 0 && return "{}"
    rows = store.features[store.features.component .== comp, :]
    isempty(rows) && return "{}"
    r = rows[1, :]
    keys_keep = [
        "frontal_ratio","temporal_ratio","blink_ratio","emg_ratio","line_ratio",
        "kurtosis","extreme_frac","ocular_score","muscle_score","line_score",
        "jump_score","artifact_score","artifact_type",
    ]
    parts = String[]
    for k in keys_keep
        sym = Symbol(k)
        hasproperty(r, sym) || continue
        v = r[sym]
        if v isa AbstractString
            push!(parts, "\"$k\":\"$v\"")
        elseif v isa Bool
            push!(parts, "\"$k\":$(v ? "true" : "false")")
        elseif v isa Number && isfinite(Float64(v))
            push!(parts, "\"$k\":$(round(Float64(v); digits=4))")
        else
            push!(parts, "\"$k\":null")
        end
    end
    # Frase de decisión determinista
    expl = _decision_explanation(r, Float64(store.summary["artifact_threshold"]))
    push!(parts, "\"explanation\":$(repr(expl))")
    push!(parts, "\"artifact_threshold\":$(Float64(store.summary["artifact_threshold"]))")
    return "{" * join(parts, ",") * "}"
end

function _decision_explanation(r, thresh::Float64)::String
    typ = hasproperty(r, :artifact_type) ? String(r.artifact_type) : "unknown"
    score = hasproperty(r, :artifact_score) ? Float64(r.artifact_score) : NaN
    scores = Dict{String,Float64}()
    for (k, lab) in (("ocular_score","ocular"), ("muscle_score","muscle"),
                     ("line_score","line"), ("jump_score","jump"))
        hasproperty(r, Symbol(k)) || continue
        scores[lab] = Float64(r[Symbol(k)])
    end
    if typ == "brain" || (!isfinite(score) || score ≤ thresh)
        return "Clasificado como brain: artifact_score " *
               (isfinite(score) ? string(round(score; digits=3)) : "?") *
               " ≤ umbral $thresh (z-score relativo entre ICs; mayor ⇒ más artefacto)."
    end
    isempty(scores) && return "Etiquetado $typ (artifact_score=$(round(score; digits=3)) > $thresh)."
    winner = argmax(scores)
    return "Etiquetado $typ porque $(winner)_score ($(round(scores[winner]; digits=3))) " *
           "supera el umbral $thresh y es el máximo entre scores de artefacto."
end

function handle_request(sock, store::ICACompStore)
    req_line = try readline(sock) catch; return end
    isempty(strip(req_line, ['\r','\n',' '])) && return
    parts = split(strip(req_line, ['\r']), ' ')
    length(parts) < 2 && return
    method, path = parts[1], parts[2]
    path_only = split(path, '?')[1]
    headers = _read_headers(sock)

    method == "OPTIONS" && (_send(sock, 204, ""); return)

    try
        if method == "GET" && path_only in ("/", "/index.html")
            _send(sock, 200, html_page(store); content_type = "text/html; charset=utf-8")

        elseif method == "GET" && path_only == "/api/meta"
            s = store.summary
            var_ret = Float64(s["variance_retained"])
            var_rej = isfinite(var_ret) ? round(100.0 - var_ret; digits=2) : NaN
            body = "{" *
                "\"t_min\":$(store.t_min),\"t_max\":$(store.t_max),\"fs\":$(store.fs)," *
                "\"n_comp\":$(store.n_comp)," *
                "\"n_rejected\":$(s["n_rejected"]),\"n_accepted\":$(s["n_accepted"])," *
                "\"variance_retained\":$(var_ret),\"variance_rejected\":$var_rej," *
                "\"rejection_source\":\"$(s["rejection_source"])\"," *
                "\"artifact_threshold\":$(Float64(s["artifact_threshold"]))," *
                "\"components\":" * _components_json(store) *
                "}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/activation"
            qs = _parse_qs(path)
            comps = Int[parse(Int, c) for c in split(get(qs, "comps", ""), ',') if !isempty(c)]
            t0 = something(tryparse(Float64, get(qs, "t0", "0")), store.t_min)
            t1 = something(tryparse(Float64, get(qs, "t1", string(store.t_max))), store.t_max)
            max_pts = something(tryparse(Int, get(qs, "max_points", "4000")), 4000)
            mask = findall(i -> store.t_act[i] >= t0 && store.t_act[i] <= t1, eachindex(store.t_act))
            if length(mask) > max_pts
                step = ceil(Int, length(mask) / max_pts)
                mask = mask[1:step:end]
            end
            t_out = store.t_act[mask]
            parts_j = String["\"t\":[" * join(round.(t_out; digits=4), ",") * "]"]
            for c in comps
                haskey(store.act, c) || continue
                push!(parts_j, "\"IC$c\":[" * join(round.(store.act[c][mask]; digits=4), ",") * "]")
            end
            _send(sock, 200, "{" * join(parts_j, ",") * "}"; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/psd"
            qs = _parse_qs(path)
            comp = something(tryparse(Int, get(qs, "comp", "1")), 1)
            ps = compute_ic_psd(store, comp)
            # Downsample for JSON if needed
            max_pts = 800
            idx = collect(eachindex(ps.f))
            if length(idx) > max_pts
                step = ceil(Int, length(idx) / max_pts)
                idx = idx[1:step:end]
            end
            body = "{\"f\":[" * join(round.(ps.f[idx]; digits=3), ",") * "]," *
                   "\"psd_db\":[" * join(round.(ps.psd_db[idx]; digits=2), ",") * "]," *
                   "\"peak_hz\":$(isfinite(ps.peak_hz) ? round(ps.peak_hz; digits=2) : "null")," *
                   "\"df\":$(round(ps.df; digits=4))}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/features"
            qs = _parse_qs(path)
            comp = something(tryparse(Int, get(qs, "comp", "1")), 1)
            _send(sock, 200, _feature_json(store, comp); content_type = "application/json")

        elseif method == "GET" && path_only == "/api/topomap"
            qs = _parse_qs(path)
            comp = something(tryparse(Int, get(qs, "comp", "1")), 1)
            p = topomap_path(comp)
            isfile(p) || error("Topomap no encontrado: $p")
            bytes = read(p)
            _send_bytes(sock, 200, bytes; content_type = "image/png")

        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_save_json(raw)
            out = save_activation_png(store, req.comps, req.t0, req.t1)
            body = "{\"ok\":true,\"file\":$(repr(basename(out)))}"
            _send(sock, 200, body; content_type = "application/json")
            println("[$(Dates.format(now(), "HH:MM:SS"))] PNG → $out")
        else
            _send(sock, 404, "{\"error\":\"not found\"}"; content_type = "application/json")
        end
    catch e
        msg = sprint(showerror, e)
        @warn "Request error" exception=e
        _send(sock, 400, "{\"ok\":false,\"error\":$(repr(msg))}"; content_type = "application/json")
    end
end

# ── HTML ───────────────────────────────────────────────────────

function html_page(store::ICACompStore)::String
    s = store.summary
    var_ret = Float64(s["variance_retained"])
    var_rej = isfinite(var_ret) ? round(100.0 - var_ret; digits=1) : NaN
    thresh = Float64(s["artifact_threshold"])
    rows = String[]
    for r in eachrow(store.components)
        c = Int(r.component)
        typ = String(r.artifact_type)
        rej = Bool(r.rejected)
        vp = Float64(r.variance_pct)
        dec = rej ? "reject" : "keep"
        push!(rows,
            """<tr class="ic-row$(rej ? " rej" : "")" data-comp="$c" data-type="$typ" data-decision="$dec" onclick="setActive($c,event)">
                 <td><input type="checkbox" class="ic-check" value="$c" title="Comparar" onclick="event.stopPropagation(); onCompareChange()"/></td>
                 <td>IC$c</td>
                 <td title="% varianza explicada por el componente (variance_explained del ICA)"><div class="vbar"><span style="width:$(min(vp*4,100))%"></span></div> $(round(vp;digits=1))%</td>
                 <td><span class="badge" data-type="$typ">$typ</span></td>
                 <td onclick="event.stopPropagation()">
                   <select class="dec-sel" data-comp="$c" onchange="onDecisionChange($c,this.value)">
                     <option value="keep"$(dec=="keep" ? " selected" : "")>Conservar</option>
                     <option value="reject"$(dec=="reject" ? " selected" : "")>Rechazar</option>
                     <option value="review">Revisar</option>
                   </select>
                 </td>
               </tr>""")
    end
    rows_html = join(rows, "\n")
    type_colors_js = join(["'$k':'$v'" for (k,v) in TYPE_COLOR], ",")
    bands_js = join(["{name:'$(b[1])',lo:$(b[2]),hi:$(b[3]),color:'$(b[4])'}" for b in BANDS], ",")

    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — ICA revisión</title>
<style>
  :root {
    --bg:#f4f6f8; --card:#fff; --border:#d8dee6; --text:#1e293b;
    --muted:#64748b; --accent:#7c3aed; --accent2:#0f766e;
    --keep:#047857; --reject:#b91c1c; --review:#a16207;
  }
  * { box-sizing: border-box; }
  body { margin:0; font-family:"IBM Plex Sans","Segoe UI",sans-serif; background:var(--bg); color:var(--text); }
  header { padding:14px 20px; background:#1e1b4b; color:#e2e8f0; display:flex; flex-wrap:wrap; gap:12px 18px; align-items:baseline; }
  header h1 { margin:0; font-size:16px; font-weight:600; }
  header span { font-size:12px; color:#a5b4fc; }
  .wrap { max-width:1360px; margin:16px auto; padding:0 16px; }
  .grid { display:grid; grid-template-columns: 360px 1fr 320px; gap:14px; }
  @media (max-width: 1100px) { .grid { grid-template-columns: 1fr 1fr; } }
  @media (max-width: 720px) { .grid { grid-template-columns: 1fr; } }
  .card { background:var(--card); border:1px solid var(--border); border-radius:10px; padding:12px 14px; margin-bottom:14px; }
  .toolbar { display:flex; flex-wrap:wrap; gap:10px 14px; align-items:end; }
  .field { display:flex; flex-direction:column; gap:4px; }
  .field label { font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:.04em; }
  select, input[type=number] {
    height:34px; padding:0 8px; border:1px solid var(--border);
    border-radius:6px; background:#fff; font-size:13px; min-width:80px;
  }
  .dec-sel { height:28px; font-size:11px; min-width:96px; padding:0 4px; }
  .btn { height:34px; padding:0 14px; border:none; border-radius:6px; font-size:13px; font-weight:600; cursor:pointer; }
  .btn-primary { background:var(--accent); color:#fff; }
  .btn-save { background:var(--accent2); color:#fff; }
  .btn:disabled { opacity:.5; cursor:wait; }
  table.ics { width:100%; border-collapse:collapse; font-size:12px; }
  table.ics th { text-align:left; color:var(--muted); font-weight:600; padding:4px 6px; border-bottom:1px solid var(--border); }
  table.ics td { padding:5px 6px; border-bottom:1px solid #f1f5f9; vertical-align:middle; }
  tr.ic-row { cursor:pointer; }
  tr.ic-row:hover { background:#f8fafc; }
  tr.ic-row.active { background:#ede9fe; outline:2px solid var(--accent); outline-offset:-2px; }
  tr.ic-row.compared td:first-child { box-shadow: inset 3px 0 0 var(--accent2); }
  tr.ic-row.rej .badge { opacity:.9; }
  .vbar { display:inline-block; width:40px; height:6px; background:#e2e8f0; border-radius:3px; vertical-align:middle; margin-right:4px; overflow:hidden; }
  .vbar span { display:block; height:100%; background:var(--accent); }
  .badge { font-size:10px; padding:2px 6px; border-radius:4px; color:#fff; background:#64748b; white-space:nowrap; }
  .badge[data-type="brain"] { background:#2563eb; }
  .badge[data-type="line_noise"] { background:#dc2626; }
  .badge[data-type="jump"] { background:#ca8a04; }
  .badge[data-type="ocular"], .badge[data-type="eye/blink"] { background:#7c3aed; }
  .badge[data-type="muscle"] { background:#ea580c; }
  .scroll { max-height:520px; overflow:auto; }
  canvas { width:100%; display:block; background:#fff; border:1px solid var(--border); border-radius:8px; }
  #topo { width:100%; max-width:300px; border:1px solid var(--border); border-radius:8px; background:#fff; }
  .status { font-size:12px; color:var(--muted); margin-top:8px; min-height:1.2em; }
  .status.ok { color:#047857; } .status.err { color:#b91c1c; }
  .active-banner {
    font-size:13px; font-weight:600; padding:10px 12px; background:#ede9fe;
    border:1px solid #c4b5fd; border-radius:8px; margin-bottom:14px;
  }
  .active-banner .sub { font-weight:400; color:var(--muted); font-size:12px; margin-top:4px; }
  .feat-block { margin-bottom:12px; }
  .feat-block h4 { margin:0 0 6px; font-size:12px; color:var(--muted); text-transform:uppercase; letter-spacing:.04em; }
  .feat { display:grid; grid-template-columns:repeat(auto-fill,minmax(100px,1fr)); gap:8px; }
  .feat .box { background:#f8fafc; border:1px solid var(--border); border-radius:8px; padding:8px; }
  .feat .k { font-size:10px; color:var(--muted); text-transform:uppercase; }
  .feat .v { font-size:14px; font-weight:600; margin-top:2px; }
  .feat .dir { font-size:10px; color:var(--muted); }
  .hint { font-size:12px; color:var(--muted); }
  .expl { font-size:12px; line-height:1.45; background:#f8fafc; border-left:3px solid var(--accent); padding:8px 10px; margin-top:8px; }
  .score-bar { position:relative; height:10px; background:#e2e8f0; border-radius:5px; margin:8px 0 4px; }
  .score-bar .fill { position:absolute; left:0; top:0; bottom:0; background:var(--accent); border-radius:5px; max-width:100%; }
  .score-bar .thresh { position:absolute; top:-3px; bottom:-3px; width:2px; background:#b91c1c; }
  .quick { display:flex; flex-wrap:wrap; gap:6px; margin-top:8px; }
  .chip { border:1px solid var(--border); background:#fff; border-radius:999px; padding:3px 10px; font-size:12px; cursor:pointer; color:var(--muted); }
  .chip:hover { border-color:var(--accent); color:var(--accent); }
  .dec-keep { color:var(--keep); } .dec-reject { color:var(--reject); } .dec-review { color:var(--review); }
</style>
</head>
<body>
<header>
  <h1>ICA · Revisión</h1>
  <span>$SUBJECT / $SESSION / $TASK</span>
  <span>$(s["n_components"]) ICs · $(s["n_rejected"]) rechazadas · $(var_rej)% var. eliminada · $(var_ret)% retenida · fuente $(s["rejection_source"]) · umbral $thresh</span>
</header>
<div class="wrap">
  <div class="card">
    <div class="toolbar">
      <div class="field"><label>Inicio (s)</label>
        <input type="number" id="t0" value="0" min="$(store.t_min)" max="$(store.t_max)" step="0.5"/></div>
      <div class="field"><label>Ventana</label>
        <select id="win">
          <option value="2">2 s</option>
          <option value="5">5 s</option>
          <option value="10" selected>10 s</option>
          <option value="full">Completo ($(round(store.t_max; digits=1)) s)</option>
        </select></div>
      <button class="btn btn-primary" onclick="refreshAll()">Actualizar</button>
      <button class="btn btn-save" id="btn-save" onclick="savePng()">Guardar PNG activaciones</button>
    </div>
    <p class="hint" style="margin:8px 0 0">Clic en fila = <strong>activo</strong> (topo · PSD · features). Checkbox = <strong>comparar</strong> activaciones. Decisión no reescribe el pipeline.</p>
  </div>

  <div class="active-banner" id="active-banner">Activo: —</div>

  <div class="grid">
    <div>
      <div class="card">
        <div style="display:flex;justify-content:space-between;align-items:baseline">
          <strong style="font-size:13px">Componentes</strong>
          <span class="hint" id="cmp-count">0 en comparación</span>
        </div>
        <div class="quick">
          <button class="chip" type="button" onclick="presetFilter('rejected')">Rechazados</button>
          <button class="chip" type="button" onclick="presetFilter('brain')">Brain</button>
          <button class="chip" type="button" onclick="presetFilter('clear')">Desmarcar comparación</button>
          <button class="chip" type="button" onclick="presetTopVar()">Top 5 var.</button>
        </div>
        <div class="scroll" style="margin-top:8px">
          <table class="ics">
            <thead><tr><th title="Comparar">Cmp</th><th>IC</th><th title="% varianza explicada">Var.</th><th>Tipo auto</th><th>Decisión</th></tr></thead>
            <tbody>$rows_html</tbody>
          </table>
        </div>
      </div>
    </div>

    <div>
      <div class="card">
        <strong style="font-size:13px">Activación temporal</strong>
        <span class="hint" id="act-hint"> · activo</span>
        <canvas id="cv-act" height="260"></canvas>
      </div>
      <div class="card">
        <strong style="font-size:13px">PSD · IC activo</strong>
        <canvas id="cv-psd" height="220"></canvas>
        <p class="hint" id="psd-hint">PSD sobre los primeros ~10 s exportados.</p>
      </div>
      <div class="status" id="status">Clic en una fila para activar un IC.</div>
    </div>

    <div>
      <div class="card">
        <strong style="font-size:13px">Topomap</strong>
        <div style="margin-top:8px;text-align:center">
          <img id="topo" alt="topomap" src=""/>
          <div class="hint" id="topo-label">—</div>
          <p class="hint" style="margin-top:6px">Escala por componente (pipeline). Pesos de mezcla ICA (a.u.).</p>
          <p class="hint">La polaridad del topomap es arbitraria; interprete la geometría espacial, no el signo rojo/azul.</p>
        </div>
      </div>
      <div class="card">
        <strong style="font-size:13px">Clasificación · IC activo</strong>
        <div id="feat-panel"><span class="hint">—</span></div>
      </div>
    </div>
  </div>
</div>

<script>
const T_MIN = $(store.t_min), T_MAX = $(store.t_max);
const TYPE_COLORS = {$type_colors_js};
const BANDS = [$bands_js];
const THRESH = $thresh;
const COMPS = $(_components_json(store));
const DECISIONS = {};
COMPS.forEach(c => { DECISIONS[c.comp] = c.decision; });
let activeComp = COMPS.length ? COMPS[0].comp : 1;

function setStatus(msg, kind) {
  const el = document.getElementById('status');
  el.textContent = msg; el.className = 'status' + (kind ? ' ' + kind : '');
}
function comparedComps() {
  return [...document.querySelectorAll('.ic-check:checked')].map(el => parseInt(el.value, 10));
}
function decisionLabel(d) {
  return d === 'reject' ? 'Rechazar' : (d === 'review' ? 'Revisar' : 'Conservar');
}
function decisionClass(d) {
  return d === 'reject' ? 'dec-reject' : (d === 'review' ? 'dec-review' : 'dec-keep');
}
function updateActiveBanner() {
  const meta = COMPS.find(x => x.comp === activeComp) || {};
  const d = DECISIONS[activeComp] || meta.decision || 'keep';
  const cmp = comparedComps();
  document.getElementById('active-banner').innerHTML =
    'Activo: <strong>IC' + activeComp + '</strong> · tipo auto: <strong>' + (meta.artifact_type||'?') +
    '</strong> · decisión: <span class="' + decisionClass(d) + '">' + decisionLabel(d) +
    '</span> · var ' + (meta.variance_pct != null ? meta.variance_pct : '?') + ' %' +
    '<div class="sub">' + (cmp.length ? ('Comparando activaciones: IC' + cmp.join(', IC')) : 'Sin comparación (se muestra solo el activo)') + '</div>';
  document.getElementById('cmp-count').textContent = cmp.length + ' en comparación';
  document.getElementById('act-hint').textContent = cmp.length ? (' · comparando ' + cmp.length) : ' · solo activo';
}
function highlightRows() {
  const cmp = new Set(comparedComps());
  document.querySelectorAll('.ic-row').forEach(tr => {
    const c = parseInt(tr.dataset.comp, 10);
    tr.classList.toggle('active', c === activeComp);
    tr.classList.toggle('compared', cmp.has(c));
  });
}
function onCompareChange() {
  updateActiveBanner();
  highlightRows();
  refreshAct();
}
function onDecisionChange(c, val) {
  DECISIONS[c] = val;
  const tr = document.querySelector('tr.ic-row[data-comp="'+c+'"]');
  if (tr) tr.dataset.decision = val;
  if (c === activeComp) updateActiveBanner();
  setStatus('Decisión IC' + c + ' → ' + decisionLabel(val) + ' (solo sesión; no altera el pipeline)', 'ok');
}
function setActive(c, ev) {
  activeComp = c;
  highlightRows();
  updateActiveBanner();
  loadTopomap(c);
  loadFeatures(c);
  loadPsd(c);
  refreshAct();
}

function presetFilter(mode) {
  document.querySelectorAll('.ic-check').forEach(el => {
    const row = COMPS.find(x => x.comp === parseInt(el.value,10));
    if (!row) { el.checked = false; return; }
    if (mode === 'rejected') el.checked = row.rejected;
    else if (mode === 'brain') el.checked = row.artifact_type === 'brain' && !row.rejected;
    else el.checked = false; // clear
  });
  onCompareChange();
  if (mode === 'rejected' || mode === 'brain') {
    const sel = comparedComps();
    if (sel.length) setActive(sel[0]);
  }
}
function presetTopVar() {
  const top = [...COMPS].sort((a,b) => b.variance_pct - a.variance_pct).slice(0,5);
  const set = new Set(top.map(x => x.comp));
  document.querySelectorAll('.ic-check').forEach(el => {
    el.checked = set.has(parseInt(el.value,10));
  });
  setActive(top[0].comp);
  onCompareChange();
}

function timeRange() {
  const win = document.getElementById('win').value;
  let t0 = parseFloat(document.getElementById('t0').value);
  if (Number.isNaN(t0)) t0 = T_MIN;
  t0 = Math.max(T_MIN, Math.min(t0, T_MAX));
  if (win === 'full') return [T_MIN, T_MAX];
  const w = parseFloat(win);
  return [t0, Math.min(T_MAX, t0 + w)];
}

function setupCanvas(id, cssH) {
  const canvas = document.getElementById(id);
  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || 700;
  canvas.style.height = cssH + 'px';
  canvas.width = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, W: cssW, H: cssH };
}

function drawSeries(canvasId, cssH, t, series, t0, t1) {
  const { ctx, W, H } = setupCanvas(canvasId, cssH);
  ctx.fillStyle = '#fff'; ctx.fillRect(0,0,W,H);
  const pad = {l:52, r:14, t:18, b:34};
  const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
  if (!t || t.length < 2 || !series.length) {
    ctx.fillStyle = '#94a3b8'; ctx.font = '13px sans-serif'; ctx.textAlign = 'center';
    ctx.fillText('Sin datos', W/2, H/2); return;
  }
  let ymin = Infinity, ymax = -Infinity;
  series.forEach(s => { for (const v of s.y) { if (v < ymin) ymin = v; if (v > ymax) ymax = v; } });
  const padY = 0.08 * Math.max(ymax - ymin, 1e-6);
  ymin -= padY; ymax += padY;
  const xAt = u => pad.l + ((u - t0) / (t1 - t0 || 1)) * pw;
  const yAt = v => pad.t + (1 - (v - ymin) / (ymax - ymin || 1)) * ph;

  ctx.strokeStyle = '#e2e8f0'; ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i++) {
    const y = pad.t + (i/4)*ph;
    ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(pad.l+pw,y); ctx.stroke();
    const v = ymax - (i/4)*(ymax-ymin);
    ctx.fillStyle = '#64748b'; ctx.font = '10px sans-serif'; ctx.textAlign = 'right';
    ctx.fillText(v.toFixed(2), pad.l - 4, y + 3);
  }
  for (let i = 0; i <= 4; i++) {
    const u = t0 + (i/4)*(t1-t0);
    const x = xAt(u);
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t+ph); ctx.stroke();
    ctx.fillStyle = '#64748b'; ctx.textAlign = 'center';
    ctx.fillText(u.toFixed(1), x, H - 12);
  }
  // zero line
  if (ymin < 0 && ymax > 0) {
    ctx.strokeStyle = '#cbd5e1'; ctx.setLineDash([4,3]);
    ctx.beginPath(); ctx.moveTo(pad.l, yAt(0)); ctx.lineTo(pad.l+pw, yAt(0)); ctx.stroke();
    ctx.setLineDash([]);
  }
  ctx.strokeStyle = '#94a3b8'; ctx.strokeRect(pad.l, pad.t, pw, ph);

  series.forEach(s => {
    ctx.strokeStyle = s.color; ctx.lineWidth = s.lw || 1.3;
    ctx.beginPath();
    for (let k = 0; k < t.length; k++) {
      const x = xAt(t[k]), y = yAt(s.y[k]);
      if (k === 0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
    }
    ctx.stroke();
  });

  let lx = pad.l + 8, ly = pad.t + 12;
  ctx.font = '11px sans-serif';
  series.forEach(s => {
    ctx.strokeStyle = s.color; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(lx, ly); ctx.lineTo(lx+14, ly); ctx.stroke();
    ctx.fillStyle = '#1e293b'; ctx.textAlign = 'left';
    ctx.fillText(s.label, lx + 18, ly + 3);
    lx += ctx.measureText(s.label).width + 40;
  });
  ctx.fillStyle = '#475569'; ctx.textAlign = 'center';
  ctx.fillText('Tiempo (s)', pad.l + pw/2, H - 2);
}

function drawPsd(f, psdDb, peakHz) {
  const { ctx, W, H } = setupCanvas('cv-psd', 220);
  ctx.fillStyle = '#fff'; ctx.fillRect(0,0,W,H);
  const pad = {l:52, r:14, t:14, b:34};
  const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
  if (!f || f.length < 2) {
    ctx.fillStyle = '#94a3b8'; ctx.font = '13px sans-serif'; ctx.textAlign = 'center';
    ctx.fillText('Sin PSD', W/2, H/2); return;
  }
  const f0 = 0, f1 = Math.min(100, f[f.length-1]);
  let ymin = Infinity, ymax = -Infinity;
  for (const v of psdDb) { if (v < ymin) ymin = v; if (v > ymax) ymax = v; }
  const padY = 0.08 * Math.max(ymax - ymin, 1);
  ymin -= padY; ymax += padY;
  const xAt = u => pad.l + ((u - f0) / (f1 - f0 || 1)) * pw;
  const yAt = v => pad.t + (1 - (v - ymin) / (ymax - ymin || 1)) * ph;

  // band shades
  BANDS.forEach(b => {
    const x1 = xAt(Math.max(f0, b.lo)), x2 = xAt(Math.min(f1, b.hi));
    if (x2 <= x1) return;
    ctx.fillStyle = b.color + '22';
    ctx.fillRect(x1, pad.t, x2-x1, ph);
  });

  // grid
  ctx.strokeStyle = '#e2e8f0'; ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i++) {
    const y = pad.t + (i/4)*ph;
    ctx.beginPath(); ctx.moveTo(pad.l,y); ctx.lineTo(pad.l+pw,y); ctx.stroke();
    const v = ymax - (i/4)*(ymax-ymin);
    ctx.fillStyle = '#64748b'; ctx.font = '10px sans-serif'; ctx.textAlign = 'right';
    ctx.fillText(v.toFixed(0), pad.l - 4, y + 3);
  }
  for (const u of [0,10,20,30,40,50,60,80,100]) {
    if (u > f1) continue;
    const x = xAt(u);
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t+ph); ctx.stroke();
    ctx.fillStyle = '#64748b'; ctx.textAlign = 'center';
    ctx.fillText(String(u), x, H - 12);
  }
  // 50 / 100 Hz markers
  [50,100].forEach(hz => {
    if (hz > f1) return;
    ctx.strokeStyle = '#dc2626'; ctx.setLineDash([3,3]);
    ctx.beginPath(); ctx.moveTo(xAt(hz), pad.t); ctx.lineTo(xAt(hz), pad.t+ph); ctx.stroke();
    ctx.setLineDash([]);
  });
  ctx.strokeStyle = '#94a3b8'; ctx.strokeRect(pad.l, pad.t, pw, ph);

  ctx.strokeStyle = '#7c3aed'; ctx.lineWidth = 1.4;
  ctx.beginPath();
  for (let k = 0; k < f.length; k++) {
    if (f[k] > f1) break;
    const x = xAt(f[k]), y = yAt(psdDb[k]);
    if (k === 0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
  }
  ctx.stroke();

  if (peakHz != null && !Number.isNaN(peakHz)) {
    ctx.fillStyle = '#1e293b'; ctx.font = '11px sans-serif'; ctx.textAlign = 'left';
    ctx.fillText('Pico ≈ ' + peakHz.toFixed(1) + ' Hz', pad.l + 8, pad.t + 14);
  }
  ctx.fillStyle = '#475569'; ctx.textAlign = 'center';
  ctx.fillText('Frecuencia (Hz)', pad.l + pw/2, H - 2);
  ctx.save(); ctx.translate(12, pad.t + ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText('dB', 0, 0); ctx.restore();
}

async function refreshAct() {
  const cmp = comparedComps();
  const comps = cmp.length ? cmp : [activeComp];
  const [t0, t1] = timeRange();
  const url = '/api/activation?comps=' + comps.join(',') + '&t0=' + t0 + '&t1=' + t1;
  const res = await fetch(url); const data = await res.json();
  if (!res.ok) { setStatus(data.error || 'Error', 'err'); return; }
  const series = comps.map(c => {
    const meta = COMPS.find(x => x.comp === c) || {artifact_type:'unknown'};
    const isActive = c === activeComp;
    return {
      y: data['IC'+c] || [],
      color: TYPE_COLORS[meta.artifact_type] || '#64748b',
      label: 'IC' + c + (isActive ? ' ★' : ''),
      lw: isActive ? 1.8 : 1.2,
    };
  });
  drawSeries('cv-act', 260, data.t, series, t0, t1);
}

async function loadPsd(c) {
  const res = await fetch('/api/psd?comp=' + c);
  const data = await res.json();
  if (!res.ok) { setStatus(data.error || 'Error PSD', 'err'); return; }
  drawPsd(data.f, data.psd_db, data.peak_hz);
  document.getElementById('psd-hint').textContent =
    'PSD IC' + c + ' sobre ~10 s exportados · Δf ≈ ' + (data.df != null ? data.df.toFixed(3) : '?') + ' Hz · líneas en 50/100 Hz';
}

function loadTopomap(c) {
  const img = document.getElementById('topo');
  img.src = '/api/topomap?comp=' + c + '&_=' + Date.now();
  const meta = COMPS.find(x => x.comp === c);
  const d = DECISIONS[c] || (meta && meta.decision) || 'keep';
  document.getElementById('topo-label').textContent =
    meta ? ('IC' + c + ' · ' + meta.artifact_type + ' · ' + meta.variance_pct + '% · ' + decisionLabel(d)) : ('IC' + c);
}

async function loadFeatures(c) {
  const res = await fetch('/api/features?comp=' + c);
  const f = await res.json();
  const panel = document.getElementById('feat-panel');
  if (!f || !Object.keys(f).length) { panel.innerHTML = '<span class="hint">Sin features</span>'; return; }

  const score = f.artifact_score;
  const thr = (f.artifact_threshold != null) ? f.artifact_threshold : THRESH;
  // Map score to bar: assume typical [-1, 4] → 0–100%
  const barMin = -1, barMax = 4;
  const pct = Math.max(0, Math.min(100, ((score - barMin) / (barMax - barMin)) * 100));
  const thrPct = Math.max(0, Math.min(100, ((thr - barMin) / (barMax - barMin)) * 100));

  let html = '<div class="feat-block"><h4>Clasificación</h4>';
  html += '<div class="feat"><div class="box"><div class="k">Tipo auto</div><div class="v">' + (f.artifact_type||'?') + '</div></div>';
  html += '<div class="box"><div class="k">artifact_score</div><div class="v">' + (typeof score==='number'?score.toFixed(3):'?') + '</div><div class="dir">↑ artefacto</div></div>';
  html += '<div class="box"><div class="k">Umbral</div><div class="v">' + thr + '</div></div></div>';
  html += '<div class="score-bar"><div class="fill" style="width:'+pct+'%"></div><div class="thresh" style="left:'+thrPct+'%"></div></div>';
  html += '<p class="hint">z-score relativo entre ICs del sujeto; mayor ⇒ más artefacto. Línea roja = umbral.</p></div>';

  const scores = [
    ['ocular_score','Ocular','↑'],
    ['muscle_score','Músculo','↑'],
    ['line_score','Línea 50 Hz','↑'],
    ['jump_score','Jump (temporal)','↑'],
  ];
  html += '<div class="feat-block"><h4>Scores (↑ artefacto)</h4><div class="feat">';
  scores.forEach(([k,lab,dir]) => {
    if (f[k] === undefined) return;
    html += '<div class="box"><div class="k">'+lab+'</div><div class="v">'+Number(f[k]).toFixed(3)+'</div><div class="dir">'+dir+'</div></div>';
  });
  html += '</div></div>';

  const evid = [
    ['frontal_ratio','Frontalidad','|A| front / no-front'],
    ['temporal_ratio','Temporalidad','|A| temp / no-temp'],
    ['blink_ratio','Blink (δ/mid)','P(0.5–4)/P(4–40)'],
    ['emg_ratio','EMG (HF/LF)','P(30–80)/P(1–30)'],
    ['line_ratio','Línea','P(48–52)/P(0.5–100)'],
    ['kurtosis','Curtosis','Pearson (~3 Gauss)'],
    ['extreme_frac','Extremos','fracción >5σ'],
  ];
  html += '<div class="feat-block"><h4>Evidencias</h4><div class="feat">';
  evid.forEach(([k,lab,hint]) => {
    if (f[k] === undefined) return;
    const v = f[k];
    const txt = (typeof v === 'number') ? (k==='extreme_frac' ? (v*100).toFixed(2)+'%' : v.toFixed(3)) : String(v);
    html += '<div class="box"><div class="k">'+lab+'</div><div class="v">'+txt+'</div><div class="dir">'+hint+'</div></div>';
  });
  html += '</div></div>';

  if (f.explanation) {
    html += '<div class="expl">' + f.explanation + '</div>';
  }
  panel.innerHTML = html;
}

async function refreshAll() {
  await refreshAct();
  loadTopomap(activeComp);
  loadFeatures(activeComp);
  loadPsd(activeComp);
  updateActiveBanner();
  setStatus('Activo IC' + activeComp + (comparedComps().length ? (' · comparando ' + comparedComps().join(',')) : ''));
}

async function savePng() {
  let comps = comparedComps();
  if (!comps.length) comps = [activeComp];
  const [t0, t1] = timeRange();
  const btn = document.getElementById('btn-save');
  btn.disabled = true;
  setStatus('Generando PNG…');
  try {
    const res = await fetch('/api/save', {
      method:'POST', headers:{'Content-Type':'application/json'},
      body: JSON.stringify({comps, t0, t1})
    });
    const data = await res.json();
    if (!res.ok || !data.ok) throw new Error(data.error || 'Error');
    setStatus('Guardado: ' + data.file, 'ok');
  } catch (e) {
    setStatus(String(e.message || e), 'err');
  } finally { btn.disabled = false; }
}

// Arranque: activo = mayor varianza; sin comparación previa
(function(){
  const top = [...COMPS].sort((a,b) => b.variance_pct - a.variance_pct);
  activeComp = top.length ? top[0].comp : 1;
  setActive(activeComp);
})();
</script>
</body>
</html>
"""
end

# ── Main ───────────────────────────────────────────────────────

function open_browser(url::String)
    try
        if Sys.isapple(); run(`open $url`)
        elseif Sys.islinux(); run(`xdg-open $url`)
        elseif Sys.iswindows(); run(`cmd /c start $url`)
        end
    catch
        @warn "No se pudo abrir el navegador automáticamente"
    end
end

function main()
    println("Cargando panel revisión ICA desde $TABLES …")
    store = load_store()
    s = store.summary
    var_ret = Float64(s["variance_retained"])
    var_rej = isfinite(var_ret) ? round(100.0 - var_ret; digits=1) : NaN
    println("  $(store.n_comp) ICs · rechazados $(s["n_rejected"]) · retenida $(var_ret)% · eliminada $(var_rej)%")
    println("  Activaciones: $(length(store.t_act)) muestras · fs≈$(round(store.fs; digits=1)) Hz")
    n_topo = count(c -> isfile(topomap_path(c)), 1:store.n_comp)
    println("  Topomaps: $n_topo PNG en $FIGS")

    server = listen(IPv4(HOST), PORT)
    url = "http://$HOST:$PORT/"
    println()
    println("UI ICA revisión → $url")
    println("Ctrl+C para detener.")
    println()
    open_browser(url)

    try
        while true
            sock = accept(server)
            @async begin
                try
                    handle_request(sock, store)
                finally
                    close(sock)
                end
            end
        end
    catch e
        isa(e, InterruptException) || rethrow(e)
        println("\nServidor detenido.")
    finally
        close(server)
    end
end

main()
