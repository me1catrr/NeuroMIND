#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Raw vs filtrado (superposición por etapa)
# ═══════════════════════════════════════════════════════════════
#
#  Superpone la señal cruda y cada etapa de filtrado
#  (notch → bandreject → highpass → lowpass) leyendo los CSV
#  generados por el pipeline en tables/.
#
#  Entrada:
#    ../../tables/raw_signal.csv
#    ../../tables/filtered_signal_notch.csv
#    ../../tables/filtered_signal_bandreject.csv
#    ../../tables/filtered_signal_highpass.csv
#    ../../tables/filtered_signal_lowpass.csv
#
#  Salida:  filt_vs_raw_<canales>[_etapas][_t0-t1].png
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_filtered_vs_raw.jl
#    # → http://127.0.0.1:8769/  (Ctrl+C para salir)
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_filtered_vs_raw.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      23-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using CSV, DataFrames, CairoMakie, Sockets, Dates

const HERE         = @__DIR__
const RESULTS_UNIT = normpath(joinpath(HERE, "..", "..", "..", "..",
                     "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))
const TABLES = joinpath(RESULTS_UNIT, "tables")
const HOST   = "127.0.0.1"
const PORT   = 8769   # raw=8765, butterfly=8766, PSD=8767, hist=8768

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK    = "eyesclosed"

# Orden de la cadena eeg_julia (+ raw al inicio)
const STAGE_ORDER = ["raw", "notch", "bandreject", "highpass", "lowpass"]
const STAGE_FILE = Dict(
    "raw"         => "raw_signal.csv",
    "notch"       => "filtered_signal_notch.csv",
    "bandreject"  => "filtered_signal_bandreject.csv",
    "highpass"    => "filtered_signal_highpass.csv",
    "lowpass"     => "filtered_signal_lowpass.csv",
)
const STAGE_LABEL = Dict(
    "raw"         => "Raw",
    "notch"       => "Notch",
    "bandreject"  => "Bandreject",
    "highpass"    => "High-pass",
    "lowpass"     => "Low-pass",
)
# Colores fijos por etapa (canvas + CairoMakie)
const STAGE_HEX = Dict(
    "raw"         => "#64748b",
    "notch"       => "#dc2626",
    "bandreject"  => "#ea580c",
    "highpass"    => "#2563eb",
    "lowpass"     => "#16a34a",
)
const STAGE_RGB = Dict(
    "raw"         => RGBf(0.39, 0.45, 0.55),
    "notch"       => RGBf(0.86, 0.15, 0.15),
    "bandreject"  => RGBf(0.92, 0.35, 0.05),
    "highpass"    => RGBf(0.15, 0.39, 0.92),
    "lowpass"     => RGBf(0.09, 0.64, 0.29),
)

# ── Datos ──────────────────────────────────────────────────────

mutable struct MultiStore
    t::Vector{Float64}
    channels::Vector{String}
    stages::Vector{String}                              # disponibles
    data::Dict{String,Dict{String,Vector{Float64}}}     # stage → ch → y
    t_min::Float64
    t_max::Float64
end

function _load_stage_csv(path::String)
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    hasproperty(df, :t_s) || error("Columna 't_s' ausente en $path")
    t = Float64.(df.t_s)
    channels = String[string(n) for n in names(df) if string(n) != "t_s"]
    data = Dict{String,Vector{Float64}}(ch => Float64.(df[!, Symbol(ch)]) for ch in channels)
    return (t=t, channels=channels, data=data)
end

function load_store(tables_dir::String)::MultiStore
    raw_path = joinpath(tables_dir, STAGE_FILE["raw"])
    raw = _load_stage_csv(raw_path)
    raw === nothing && error("No encontrado: $raw_path")

    data = Dict{String,Dict{String,Vector{Float64}}}()
    stages = String[]
    t = raw.t
    channels = raw.channels

    for key in STAGE_ORDER
        path = joinpath(tables_dir, STAGE_FILE[key])
        st = _load_stage_csv(path)
        if st === nothing
            @warn "Etapa ausente, se omite" stage=key path=path
            continue
        end
        if length(st.t) != length(t) || !all(isapprox.(st.t, t; atol=1e-6))
            @warn "Grid temporal distinto; se alinea por índice mínimo" stage=key
        end
        # Canales: intersección con raw
        ch_ok = intersect(channels, st.channels)
        isempty(ch_ok) && continue
        data[key] = Dict{String,Vector{Float64}}(ch => st.data[ch] for ch in ch_ok)
        push!(stages, key)
    end
    isempty(stages) && error("Ninguna etapa cargada desde $tables_dir")
    # Canales comunes a todas las etapas disponibles
    common = copy(channels)
    for key in stages
        common = intersect(common, collect(keys(data[key])))
    end
    isempty(common) && error("Sin canales comunes entre etapas")
    sort!(common; by = ch -> findfirst(==(ch), channels))

    return MultiStore(t, common, stages, data, minimum(t), maximum(t))
end

# ── PNG ────────────────────────────────────────────────────────

function save_png(
    store::MultiStore,
    selected_ch::Vector{String},
    selected_st::Vector{String},
    t0::Float64,
    t1::Float64,
    scale::String,
)::String
    isempty(selected_ch) && error("Selecciona al menos un canal")
    isempty(selected_st) && error("Selecciona al menos una etapa")
    t0 >= t1 && error("t0 debe ser < t1")
    for ch in selected_ch
        ch in store.channels || error("Canal desconocido: $ch")
    end
    for st in selected_st
        st in store.stages || error("Etapa desconocida: $st")
    end

    mask = (store.t .>= t0) .& (store.t .<= t1)
    count(mask) < 2 && error("Ventana temporal vacía")
    t = store.t[mask]

    fig = Figure(size = (1100, 440), fontsize = 13)
    ax = Axis(fig[1, 1];
        title  = "Raw vs filtrado — $SUBJECT / $SESSION / $TASK",
        xlabel = "Tiempo (s)",
        ylabel = "Amplitud (µV)",
    )
    ax.xgridvisible = true
    ax.ygridvisible = true

    ymin, ymax = Inf, -Inf
    n_ch = length(selected_ch)
    for (ci, ch) in enumerate(selected_ch)
        # Estilo de trazo por canal si hay varios
        ls = n_ch == 1 ? :solid :
             (ci == 1 ? :solid : (ci == 2 ? :dash : :dot))
        for st in selected_st
            y = store.data[st][ch][mask]
            col = STAGE_RGB[st]
            lw = st == "raw" ? 1.6 : 1.2
            lab = n_ch == 1 ? STAGE_LABEL[st] : "$ch · $(STAGE_LABEL[st])"
            lines!(ax, t, y; color = col, linewidth = lw, linestyle = ls, label = lab)
            ymin = min(ymin, minimum(y))
            ymax = max(ymax, maximum(y))
        end
    end
    axislegend(ax; position = :rt, framevisible = true, labelsize = 11, nbanks = 1)

    if scale == "auto"
        pad = 0.05 * max(ymax - ymin, 1.0)
        ylims!(ax, ymin - pad, ymax + pad)
    else
        lim = parse(Float64, scale)
        ylims!(ax, -lim, lim)
    end
    xlims!(ax, t0, t1)

    tag_t = (isapprox(t0, store.t_min; atol=1e-3) && isapprox(t1, store.t_max; atol=1e-3)) ?
            "" : "_$(round(Int, t0))-$(round(Int, t1))s"
    ch_tag = length(selected_ch) <= 3 ? join(selected_ch, "_") : "$(length(selected_ch))ch"
    st_tag = length(selected_st) == length(store.stages) ? "" :
             "_" * join(selected_st, "-")
    out_name = "filt_vs_raw_" * ch_tag * st_tag * tag_t * ".png"
    out_path = joinpath(HERE, out_name)
    save(out_path, fig; px_per_unit = 2)
    return out_path
end

# ── HTTP ───────────────────────────────────────────────────────

function _read_headers(sock)::Dict{String,String}
    headers = Dict{String,String}()
    while true
        line = readline(sock)
        line == "" && break
        line == "\r" && break
        m = match(r"^([^:]+):\s*(.+?)(?:\r)?$", line)
        m === nothing && continue
        headers[lowercase(m.captures[1])] = strip(m.captures[2])
    end
    return headers
end

function _read_body(sock, headers::Dict{String,String})::String
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
    msg = "HTTP/1.1 $status\r\n" *
          "Content-Type: $content_type\r\n" *
          "Content-Length: $(sizeof(body))\r\n" *
          "Connection: close\r\n" *
          "Access-Control-Allow-Origin: *\r\n" *
          "\r\n" * body
    write(sock, msg)
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
    _list = function (key)
        out = String[]
        m = match(Regex("\"" * key * "\"\\s*:\\s*\\[(.*?)\\]", "s"), s)
        m === nothing && return out
        for cap in eachmatch(r"\"([^\"]+)\"", m.captures[1])
            push!(out, cap.captures[1])
        end
        return out
    end
    _num = function (key, default)
        m2 = match(Regex("\"" * key * "\"\\s*:\\s*([-\\d.eE]+)"), s)
        m2 === nothing && return Float64(default)
        v = tryparse(Float64, m2.captures[1])
        return v === nothing ? Float64(default) : v
    end
    _str = function (key, default)
        m2 = match(Regex("\"" * key * "\"\\s*:\\s*\"([^\"]*)\""), s)
        return m2 === nothing ? String(default) : String(m2.captures[1])
    end
    return (
        channels = _list("channels"),
        stages   = _list("stages"),
        t0 = _num("t0", 0.0),
        t1 = _num("t1", 10.0),
        scale = _str("scale", "auto"),
    )
end

function handle_request(sock, store::MultiStore)
    req_line = try readline(sock) catch; return end
    isempty(strip(req_line, ['\r','\n',' '])) && return
    parts = split(strip(req_line, ['\r']), ' ')
    length(parts) < 2 && return
    method, path = parts[1], parts[2]
    path_only = split(path, '?')[1]
    headers = _read_headers(sock)

    if method == "OPTIONS"
        _send(sock, 204, "")
        return
    end

    try
        if method == "GET" && path_only in ("/", "/index.html")
            _send(sock, 200, html_page(store); content_type = "text/html; charset=utf-8")

        elseif method == "GET" && path_only == "/api/meta"
            body = "{\"t_min\":$(store.t_min),\"t_max\":$(store.t_max)," *
                   "\"n_samples\":$(length(store.t))," *
                   "\"channels\":[" * join(["\"$c\"" for c in store.channels], ",") * "]," *
                   "\"stages\":[" * join(["\"$s\"" for s in store.stages], ",") * "]}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/signal"
            qs = _parse_qs(path)
            chs = String[String(c) for c in split(get(qs, "chs", ""), ',') if !isempty(c)]
            sts = String[String(c) for c in split(get(qs, "stages", ""), ',') if !isempty(c)]
            isempty(sts) && (sts = copy(store.stages))
            t0 = something(tryparse(Float64, get(qs, "t0", "0")), store.t_min)
            t1 = something(tryparse(Float64, get(qs, "t1", string(store.t_max))), store.t_max)
            max_pts = something(tryparse(Int, get(qs, "max_points", "5000")), 5000)

            mask = findall(i -> store.t[i] >= t0 && store.t[i] <= t1, eachindex(store.t))
            if length(mask) > max_pts
                step = ceil(Int, length(mask) / max_pts)
                mask = mask[1:step:end]
            end
            t_out = store.t[mask]
            parts_j = String["\"t\":[" * join(round.(t_out; digits=4), ",") * "]"]
            for st in sts
                st in store.stages || continue
                for ch in chs
                    haskey(store.data[st], ch) || continue
                    ys = store.data[st][ch][mask]
                    push!(parts_j, "\"$(st)__$ch\":[" * join(round.(ys; digits=3), ",") * "]")
                end
            end
            _send(sock, 200, "{" * join(parts_j, ",") * "}"; content_type = "application/json")

        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_save_json(raw)
            isempty(req.channels) && error("channels vacío")
            isempty(req.stages) && error("stages vacío")
            out = save_png(store, req.channels, req.stages, req.t0, req.t1, req.scale)
            body = "{\"ok\":true,\"path\":$(repr(out)),\"file\":$(repr(basename(out)))}"
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

function html_page(store::MultiStore)::String
    ch_opts = join(["<label class=\"ch\"><input type=\"checkbox\" value=\"$c\">$c</label>"
                    for c in store.channels], "\n")
    st_opts = join([
        let hex = STAGE_HEX[s], lab = STAGE_LABEL[s], checked = s in ("raw", "lowpass") ? " checked" : ""
            """<label class="st" style="--st:$hex">
                 <input type="checkbox" value="$s"$checked>
                 <span class="dot"></span>$lab
               </label>"""
        end
        for s in store.stages], "\n")
    stage_colors_js = join(["'$s':'$(STAGE_HEX[s])'" for s in store.stages], ",")
    stage_labels_js = join(["'$s':'$(STAGE_LABEL[s])'" for s in store.stages], ",")
    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — Raw vs filtrado</title>
<style>
  :root {
    --bg:#f4f6f8; --card:#fff; --border:#d8dee6; --text:#1e293b;
    --muted:#64748b; --accent:#2563eb; --accent2:#0f766e;
  }
  * { box-sizing: border-box; }
  body { margin:0; font-family:"IBM Plex Sans","Segoe UI",sans-serif; background:var(--bg); color:var(--text); }
  header { padding:14px 20px; background:#0f172a; color:#e2e8f0; display:flex; align-items:baseline; gap:14px; flex-wrap:wrap; }
  header h1 { margin:0; font-size:16px; font-weight:600; }
  header span { font-size:12px; color:#94a3b8; }
  .wrap { max-width:1200px; margin:16px auto; padding:0 16px; }
  .card { background:var(--card); border:1px solid var(--border); border-radius:10px; padding:14px 16px; margin-bottom:14px; }
  .toolbar { display:flex; flex-wrap:wrap; gap:12px 18px; align-items:end; }
  .field { display:flex; flex-direction:column; gap:4px; }
  .field label { font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:.04em; }
  select, input[type=number] {
    height:34px; padding:0 10px; border:1px solid var(--border);
    border-radius:6px; background:#fff; font-size:13px; min-width:100px;
  }
  .btn { height:34px; padding:0 16px; border:none; border-radius:6px; font-size:13px; font-weight:600; cursor:pointer; }
  .btn-primary { background:var(--accent); color:#fff; }
  .btn-primary:hover { background:#1d4ed8; }
  .btn-save { background:var(--accent2); color:#fff; }
  .btn-save:hover { background:#0d9488; }
  .btn:disabled { opacity:.5; cursor:wait; }
  .channels {
    display:grid; grid-template-columns:repeat(auto-fill,minmax(88px,1fr));
    gap:6px 8px; max-height:140px; overflow:auto; padding:4px 0;
  }
  .stages { display:flex; flex-wrap:wrap; gap:8px 14px; padding:4px 0; }
  .ch, .st {
    display:flex; align-items:center; gap:6px; font-size:13px;
    padding:4px 8px; border-radius:5px; cursor:pointer; user-select:none;
  }
  .ch:hover, .st:hover { background:#f1f5f9; }
  .ch input, .st input { accent-color:var(--accent); }
  .st .dot {
    width:10px; height:10px; border-radius:50%; background:var(--st);
    display:inline-block; flex-shrink:0;
  }
  .quick { display:flex; flex-wrap:wrap; gap:6px; margin-top:8px; }
  .chip {
    border:1px solid var(--border); background:#fff; border-radius:999px;
    padding:3px 10px; font-size:12px; cursor:pointer; color:var(--muted);
  }
  .chip:hover { border-color:var(--accent); color:var(--accent); }
  canvas {
    width:100%; height:400px; background:#fff;
    border:1px solid var(--border); border-radius:8px;
  }
  .status { font-size:12px; color:var(--muted); margin-top:8px; min-height:1.2em; }
  .status.ok { color:#047857; }
  .status.err { color:#b91c1c; }
  .hint { font-size:12px; color:var(--muted); margin-top:4px; }
</style>
</head>
<body>
<header>
  <h1>Raw vs filtrado</h1>
  <span>$SUBJECT / $SESSION / $TASK · etapas: $(join(store.stages, " → "))</span>
</header>
<div class="wrap">
  <div class="card">
    <div class="toolbar">
      <div class="field">
        <label>Inicio (s)</label>
        <input type="number" id="t0" value="0" min="$(store.t_min)" max="$(store.t_max)" step="0.5"/>
      </div>
      <div class="field">
        <label>Ventana</label>
        <select id="win">
          <option value="5">5 s</option>
          <option value="10" selected>10 s</option>
          <option value="30">30 s</option>
          <option value="50">50 s</option>
          <option value="100">100 s</option>
          <option value="full">Completo</option>
        </select>
      </div>
      <div class="field">
        <label>Escala</label>
        <select id="scale">
          <option value="auto" selected>Auto</option>
          <option value="50">±50 µV</option>
          <option value="100">±100 µV</option>
          <option value="200">±200 µV</option>
          <option value="500">±500 µV</option>
        </select>
      </div>
      <button class="btn btn-primary" onclick="refreshPlot()">Actualizar</button>
      <button class="btn btn-save" id="btn-save" onclick="savePng()">Guardar PNG</button>
    </div>
  </div>

  <div class="card">
    <div style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:6px">
      <strong style="font-size:13px">Etapas a superponer</strong>
      <span style="font-size:12px;color:var(--muted)" id="st-count">—</span>
    </div>
    <div class="stages" id="stages">$st_opts</div>
    <div class="quick">
      <button class="chip" type="button" onclick="presetStages(['raw','lowpass'])">Raw + Low-pass</button>
      <button class="chip" type="button" onclick="presetStages(['raw','notch','lowpass'])">Raw + Notch + LP</button>
      <button class="chip" type="button" onclick="presetStages($(replace(repr(store.stages), "\"" => "'")))">Todas</button>
      <button class="chip" type="button" onclick="presetStages([])">Ninguna</button>
    </div>
    <p class="hint">Orden pipeline <code>eeg_julia</code>: Notch → Bandreject → High-pass → Low-pass. Cada CSV es la señal acumulada tras ese paso.</p>
  </div>

  <div class="card">
    <div style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:6px">
      <strong style="font-size:13px">Canales</strong>
      <span style="font-size:12px;color:var(--muted)" id="ch-count">0 seleccionados</span>
    </div>
    <div class="channels" id="channels">$ch_opts</div>
    <div class="quick">
      <button class="chip" type="button" onclick="presetCh(['Cz'])">Cz</button>
      <button class="chip" type="button" onclick="presetCh(['Fp2','Cz'])">Fp2 + Cz</button>
      <button class="chip" type="button" onclick="presetCh(['O1','Oz','O2'])">Occipitales</button>
      <button class="chip" type="button" onclick="presetCh([])">Ninguno</button>
    </div>
  </div>

  <div class="card">
    <canvas id="cv" width="1100" height="400"></canvas>
    <div class="status" id="status">Elige canal(es) y etapas, luego Actualizar.</div>
  </div>
</div>

<script>
const T_MIN = $(store.t_min);
const T_MAX = $(store.t_max);
const STAGE_COLORS = {$stage_colors_js};
const STAGE_LABELS = {$stage_labels_js};
const STAGE_ORDER = $(replace(repr(store.stages), "\"" => "'"));

function selectedChannels() {
  return [...document.querySelectorAll('#channels input:checked')].map(el => el.value);
}
function selectedStages() {
  const checked = [...document.querySelectorAll('#stages input:checked')].map(el => el.value);
  return STAGE_ORDER.filter(s => checked.includes(s));
}
function updateCounts() {
  const n = selectedChannels().length;
  const m = selectedStages().length;
  document.getElementById('ch-count').textContent = n + ' seleccionado' + (n === 1 ? '' : 's');
  document.getElementById('st-count').textContent = m + ' etapa' + (m === 1 ? '' : 's');
}
document.querySelectorAll('#channels input, #stages input').forEach(el =>
  el.addEventListener('change', updateCounts));

function presetCh(chs) {
  document.querySelectorAll('#channels input').forEach(el => { el.checked = chs.includes(el.value); });
  updateCounts();
  if (chs.length) refreshPlot();
}
function presetStages(sts) {
  document.querySelectorAll('#stages input').forEach(el => { el.checked = sts.includes(el.value); });
  updateCounts();
  if (selectedChannels().length && sts.length) refreshPlot();
}

function timeRange() {
  const win = document.getElementById('win').value;
  let t0 = parseFloat(document.getElementById('t0').value);
  if (Number.isNaN(t0)) t0 = T_MIN;
  t0 = Math.max(T_MIN, Math.min(t0, T_MAX));
  if (win === 'full') return [T_MIN, T_MAX];
  const w = parseFloat(win);
  let t1 = t0 + w;
  if (t1 > T_MAX) { t1 = T_MAX; t0 = Math.max(T_MIN, t1 - w); }
  return [t0, t1];
}

function setStatus(msg, kind) {
  const el = document.getElementById('status');
  el.textContent = msg;
  el.className = 'status' + (kind ? ' ' + kind : '');
}

async function refreshPlot() {
  const chs = selectedChannels();
  const sts = selectedStages();
  if (!chs.length) { setStatus('Selecciona al menos un canal.', 'err'); return; }
  if (!sts.length) { setStatus('Selecciona al menos una etapa.', 'err'); return; }
  const [t0, t1] = timeRange();
  setStatus('Cargando…');
  const url = '/api/signal?chs=' + encodeURIComponent(chs.join(',')) +
              '&stages=' + encodeURIComponent(sts.join(',')) +
              '&t0=' + t0 + '&t1=' + t1 + '&max_points=5000';
  const res = await fetch(url);
  const data = await res.json();
  if (!res.ok) { setStatus(data.error || 'Error', 'err'); return; }
  draw(data, chs, sts, t0, t1);
  setStatus(chs.join(', ') + ' · ' + sts.map(s => STAGE_LABELS[s] || s).join(' + ') +
            ' · ' + t0.toFixed(1) + '–' + t1.toFixed(1) + ' s');
}

function draw(data, chs, sts, t0, t1) {
  const canvas = document.getElementById('cv');
  const ctx = canvas.getContext('2d');
  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || 1100;
  const cssH = 400;
  canvas.style.height = cssH + 'px';
  canvas.width  = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  const W = cssW, H = cssH;

  const pad = {l: 58, r: 18, t: 28, b: 42};
  const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);

  const t = data.t || [];
  if (t.length < 2) return;

  let ymin = Infinity, ymax = -Infinity;
  const scale = document.getElementById('scale').value;
  const series = [];
  chs.forEach((ch, ci) => {
    sts.forEach(st => {
      const key = st + '__' + ch;
      const y = data[key] || [];
      series.push({ch, st, y, ci});
      if (scale === 'auto') {
        for (const v of y) { if (v < ymin) ymin = v; if (v > ymax) ymax = v; }
      }
    });
  });
  if (scale === 'auto') {
    const padY = 0.05 * Math.max(ymax - ymin, 1);
    ymin -= padY; ymax += padY;
  } else {
    const lim = parseFloat(scale);
    ymin = -lim; ymax = lim;
  }

  const xAt = u => pad.l + ((u - t0) / (t1 - t0)) * pw;
  const yAt = v => pad.t + (1 - (v - ymin) / (ymax - ymin)) * ph;

  ctx.strokeStyle = '#e2e8f0'; ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i++) {
    const v = ymin + (i / 4) * (ymax - ymin);
    const y = yAt(v);
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + pw, y); ctx.stroke();
    ctx.fillStyle = '#64748b'; ctx.font = '11px sans-serif';
    ctx.textAlign = 'right'; ctx.fillText(v.toFixed(0), pad.l - 6, y + 3);
  }
  for (let i = 0; i <= 5; i++) {
    const u = t0 + (i / 5) * (t1 - t0);
    const x = xAt(u);
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
    ctx.fillStyle = '#64748b'; ctx.textAlign = 'center';
    ctx.fillText(u.toFixed(1), x, H - 16);
  }
  ctx.strokeStyle = '#94a3b8';
  ctx.strokeRect(pad.l, pad.t, pw, ph);

  // series: color = etapa; dash = canal (si >1)
  series.forEach(s => {
    ctx.strokeStyle = STAGE_COLORS[s.st] || '#334155';
    ctx.lineWidth = s.st === 'raw' ? 1.8 : 1.3;
    if (chs.length === 1) ctx.setLineDash([]);
    else if (s.ci === 0) ctx.setLineDash([]);
    else if (s.ci === 1) ctx.setLineDash([6, 4]);
    else ctx.setLineDash([2, 3]);
    ctx.beginPath();
    for (let k = 0; k < t.length; k++) {
      const x = xAt(t[k]), yy = yAt(s.y[k]);
      if (k === 0) ctx.moveTo(x, yy); else ctx.lineTo(x, yy);
    }
    ctx.stroke();
  });
  ctx.setLineDash([]);

  // leyenda
  let lx = pad.l + 8, ly = pad.t + 14;
  ctx.font = '12px "IBM Plex Sans", sans-serif';
  series.forEach(s => {
    const lab = chs.length === 1 ? (STAGE_LABELS[s.st] || s.st)
                                 : (s.ch + ' · ' + (STAGE_LABELS[s.st] || s.st));
    ctx.strokeStyle = STAGE_COLORS[s.st] || '#334155';
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(lx, ly - 3); ctx.lineTo(lx + 14, ly - 3); ctx.stroke();
    ctx.fillStyle = '#1e293b'; ctx.textAlign = 'left';
    ctx.fillText(lab, lx + 18, ly);
    lx += ctx.measureText(lab).width + 40;
    if (lx > pad.l + pw - 80) { lx = pad.l + 8; ly += 16; }
  });

  ctx.fillStyle = '#475569'; ctx.font = '12px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('Tiempo (s)', pad.l + pw / 2, H - 4);
  ctx.save();
  ctx.translate(14, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('Amplitud (µV)', 0, 0);
  ctx.restore();
}

async function savePng() {
  const chs = selectedChannels();
  const sts = selectedStages();
  if (!chs.length) { setStatus('Selecciona al menos un canal.', 'err'); return; }
  if (!sts.length) { setStatus('Selecciona al menos una etapa.', 'err'); return; }
  const [t0, t1] = timeRange();
  const scale = document.getElementById('scale').value;
  const btn = document.getElementById('btn-save');
  btn.disabled = true;
  setStatus('Generando PNG…');
  try {
    const res = await fetch('/api/save', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({channels: chs, stages: sts, t0, t1, scale})
    });
    const data = await res.json();
    if (!res.ok || !data.ok) throw new Error(data.error || 'Error al guardar');
    setStatus('Guardado: ' + data.file, 'ok');
  } catch (e) {
    setStatus(String(e.message || e), 'err');
  } finally {
    btn.disabled = false;
  }
}

updateCounts();
presetCh(['Cz']);
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
    println("Cargando tablas desde $TABLES …")
    store = load_store(TABLES)
    println("  $(length(store.channels)) canales · $(length(store.t)) muestras · ",
            "$(round(store.t_min; digits=1))–$(round(store.t_max; digits=1)) s")
    println("  Etapas: " * join(store.stages, " → "))

    server = listen(IPv4(HOST), PORT)
    url = "http://$HOST:$PORT/"
    println()
    println("UI raw vs filtrado → $url")
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
