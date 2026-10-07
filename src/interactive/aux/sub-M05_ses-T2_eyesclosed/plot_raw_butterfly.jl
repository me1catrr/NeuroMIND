#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor butterfly / montaje apilado (GUI navegador)
# ═══════════════════════════════════════════════════════════════
#
#  Vista multicanal con offset vertical (escala de amplitud común).
#  Misma interfaz interactiva que plot_raw.jl: canales, ventana,
#  escala y «Guardar PNG» (CairoMakie).
#
#  Entrada:  results/subjects/sub-M05/ses-T2/eyesclosed/tables/raw_signal.csv
#  Salida:   raw_butterfly[_t0-t1][_Nch].png (en este directorio)
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_raw_butterfly.jl
#    # → http://127.0.0.1:8766/  (Ctrl+C para salir)
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_raw_butterfly.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      23-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics

const HERE         = @__DIR__
const RESULTS_UNIT = normpath(joinpath(HERE, "..", "..", "..", "..",
                     "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))
const CSV_PATH     = joinpath(RESULTS_UNIT, "tables", "raw_signal.csv")
const HOST         = "127.0.0.1"
const PORT         = 8766   # distinto de plot_raw.jl (8765)

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK_ID = "EC"          # etiqueta estilo BIDS del informe
const FS_HZ   = 500.0

# Paleta cualitativa (trazas distinguibles en montaje)
const LINE_COLORS = [
    RGBf(0.89, 0.10, 0.11), RGBf(0.21, 0.49, 0.72), RGBf(0.30, 0.69, 0.29),
    RGBf(0.60, 0.31, 0.64), RGBf(1.00, 0.50, 0.00), RGBf(0.65, 0.34, 0.16),
    RGBf(0.97, 0.51, 0.75), RGBf(0.40, 0.40, 0.40), RGBf(0.70, 0.70, 0.15),
    RGBf(0.00, 0.69, 0.73), RGBf(0.80, 0.20, 0.40), RGBf(0.15, 0.55, 0.35),
    RGBf(0.45, 0.25, 0.70), RGBf(0.90, 0.60, 0.10), RGBf(0.20, 0.40, 0.70),
    RGBf(0.55, 0.15, 0.15), RGBf(0.10, 0.60, 0.50), RGBf(0.70, 0.40, 0.10),
    RGBf(0.35, 0.35, 0.70), RGBf(0.60, 0.60, 0.20), RGBf(0.80, 0.30, 0.55),
    RGBf(0.25, 0.55, 0.75), RGBf(0.75, 0.25, 0.25), RGBf(0.40, 0.65, 0.30),
    RGBf(0.55, 0.30, 0.55), RGBf(0.85, 0.45, 0.15), RGBf(0.20, 0.50, 0.55),
    RGBf(0.65, 0.20, 0.45), RGBf(0.30, 0.45, 0.20), RGBf(0.50, 0.50, 0.75),
    RGBf(0.90, 0.35, 0.35),
]

mutable struct RawStore
    t::Vector{Float64}
    channels::Vector{String}
    data::Dict{String,Vector{Float64}}
    t_min::Float64
    t_max::Float64
    fs::Float64
end

function load_store(path::String)::RawStore
    isfile(path) || error("No encontrado: $path")
    df = CSV.read(path, DataFrame)
    hasproperty(df, :t_s) || error("Columna 't_s' ausente")
    t = Float64.(df.t_s)
    channels = String[string(n) for n in names(df) if string(n) != "t_s"]
    data = Dict{String,Vector{Float64}}(ch => Float64.(df[!, Symbol(ch)]) for ch in channels)
    fs = length(t) >= 2 ? 1.0 / (t[2] - t[1]) : FS_HZ
    return RawStore(t, channels, data, minimum(t), maximum(t), fs)
end

"""Espaciado vertical (µV entre baselines) según escala común."""
function _spacing_uv(store::RawStore, selected::AbstractVector{<:AbstractString}, mask, scale::AbstractString)::Float64
    if scale != "auto"
        return parse(Float64, scale)
    end
    stds = Float64[]
    for ch in selected
        y = store.data[String(ch)][mask]
        s = std(y)
        isfinite(s) && s > 0 && push!(stds, s)
    end
    isempty(stds) && return 50.0
    # 3.5·mediana(σ) deja margen visual sin solape habitual
    return max(3.5 * median(stds), 25.0)
end

function save_butterfly_png(
    store::RawStore,
    selected::AbstractVector{<:AbstractString},
    t0::Float64,
    t1::Float64,
    scale::AbstractString,
)::String
    selected = String[String(ch) for ch in selected]
    isempty(selected) && error("Selecciona al menos un canal")
    t0 >= t1 && error("t0 debe ser < t1")
    for ch in selected
        haskey(store.data, ch) || error("Canal desconocido: $ch")
    end

    mask = (store.t .>= t0) .& (store.t .<= t1)
    count(mask) < 2 && error("Ventana temporal vacía")
    t = store.t[mask]
    n = length(selected)
    spacing = _spacing_uv(store, selected, mask, scale)

    # Offset: canal 1 arriba
    offsets = Float64[(n - i) * spacing for i in 1:n]

    # Tipografía y tamaño adaptados al nº de canales (PNG nítido)
    row_px = n <= 16 ? 34 : (n <= 24 ? 30 : 26)
    fig_h  = max(520, min(row_px * n + 110, 1800))
    fig_w  = 1200
    tick_fs = clamp(Int(round(0.42 * row_px)), 9, 14)
    title_fs = 15
    fig = Figure(size = (fig_w, fig_h), fontsize = tick_fs)
    n_lab = "$n chans"
    ax = Axis(fig[1, 1];
        title  = "DATOS CRUDOS — $(SUBJECT)_$(SESSION)_task-$(TASK_ID) — fs=$(round(Int, store.fs)) Hz, $n_lab",
        xlabel = "Tiempo (s)",
        ylabel = "Canal (offset)",
        yticks = (offsets, selected),
        titlesize = title_fs,
        xlabelsize = tick_fs + 2,
        ylabelsize = tick_fs + 2,
        xticklabelsize = tick_fs,
        yticklabelsize = tick_fs,
        yticklabelspace = 58,
    )
    ax.xgridvisible = true
    ax.ygridvisible = true

    ymin = Inf
    ymax = -Inf
    for (i, ch) in enumerate(selected)
        y = store.data[ch][mask] .+ offsets[i]
        col = LINE_COLORS[mod1(i, length(LINE_COLORS))]
        lines!(ax, t, y; color = col, linewidth = 1.0)
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))
    end

    # Incluir la amplitud real (p. ej. Fp2 con artefactos oculares),
    # no solo ±0.6·spacing — evita trazas “flotando” bajo el eje.
    pad = max(0.05 * (ymax - ymin), 0.25 * spacing)
    ylims!(ax, ymin - pad, ymax + pad)
    xlims!(ax, t0, t1)

    tag_t = (isapprox(t0, store.t_min; atol=1e-3) && isapprox(t1, store.t_max; atol=1e-3)) ?
            "" : "_$(round(Int, t0))-$(round(Int, t1))s"
    out_name = "raw_butterfly$(tag_t)_$(n)ch.png"
    out_path = joinpath(HERE, out_name)
    save(out_path, fig; px_per_unit = 3)
    return out_path
end

# ── HTTP mínimo ────────────────────────────────────────────────

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
    write(sock,
        "HTTP/1.1 $status\r\nContent-Type: $content_type\r\n" *
        "Content-Length: $(sizeof(body))\r\nConnection: close\r\n" *
        "Access-Control-Allow-Origin: *\r\n\r\n" * body)
end

function _parse_save_json(s::String)
    chs = String[]
    m = match(r"\"channels\"\s*:\s*\[(.*?)\]"s, s)
    if m !== nothing
        for cap in eachmatch(r"\"([^\"]+)\"", m.captures[1])
            push!(chs, cap.captures[1])
        end
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
    return (channels = chs, t0 = _num("t0", 0.0), t1 = _num("t1", 10.0), scale = _str("scale", "auto"))
end

function handle_request(sock, store::RawStore)
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
            body = "{\"t_min\":$(store.t_min),\"t_max\":$(store.t_max)," *
                   "\"fs\":$(store.fs),\"n_samples\":$(length(store.t))," *
                   "\"channels\":[" * join(["\"$c\"" for c in store.channels], ",") * "]}"
            _send(sock, 200, body; content_type = "application/json")
        elseif method == "GET" && path_only == "/api/signal"
            qs = Dict{String,String}()
            if occursin('?', path)
                for pair in split(split(path, '?', limit=2)[2], '&')
                    kv = split(pair, '=', limit=2)
                    length(kv) == 2 && (qs[kv[1]] = kv[2])
                end
            end
            chs = String[String(c) for c in split(replace(get(qs, "chs", ""), "%2C" => ","), ',') if !isempty(c)]
            isempty(chs) && (chs = copy(store.channels))
            t0 = something(tryparse(Float64, get(qs, "t0", "0")), store.t_min)
            t1 = something(tryparse(Float64, get(qs, "t1", "10")), min(10.0, store.t_max))
            max_pts = something(tryparse(Int, get(qs, "max_points", "2500")), 2500)
            scale = String(get(qs, "scale", "auto"))

            mask_idx = findall(i -> store.t[i] >= t0 && store.t[i] <= t1, eachindex(store.t))
            if length(mask_idx) > max_pts
                step = ceil(Int, length(mask_idx) / max_pts)
                mask_idx = mask_idx[1:step:end]
            end
            t_out = store.t[mask_idx]
            # spacing para preview JS
            full_mask = (store.t .>= t0) .& (store.t .<= t1)
            spacing = _spacing_uv(store, chs, full_mask, scale)

            parts = String[
                "\"t\":[" * join(round.(t_out; digits=4), ",") * "]",
                "\"spacing\":$(round(spacing; digits=2))",
            ]
            for ch in chs
                haskey(store.data, ch) || continue
                ys = store.data[ch][mask_idx]
                push!(parts, "\"$ch\":[" * join(round.(ys; digits=2), ",") * "]")
            end
            _send(sock, 200, "{" * join(parts, ",") * "}"; content_type = "application/json")
        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_save_json(raw)
            isempty(req.channels) && error("channels vacío")
            out = save_butterfly_png(store, req.channels, req.t0, req.t1, req.scale)
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

function html_page(store::RawStore)::String
    ch_opts = join([
        "<label class=\"ch\"><input type=\"checkbox\" value=\"$c\" checked>$c</label>"
        for c in store.channels], "\n")
    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — Butterfly / montaje</title>
<style>
  :root {
    --bg:#f4f6f8; --card:#fff; --border:#d8dee6; --text:#1e293b;
    --muted:#64748b; --accent:#2563eb; --accent2:#0f766e;
  }
  * { box-sizing: border-box; }
  body { margin:0; font-family:"IBM Plex Sans","Segoe UI",sans-serif; background:var(--bg); color:var(--text); }
  header { padding:14px 20px; background:#0f172a; color:#e2e8f0; display:flex; align-items:baseline; gap:14px; }
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
    gap:6px 8px; max-height:160px; overflow:auto; padding:4px 0;
  }
  .ch { display:flex; align-items:center; gap:6px; font-size:13px; padding:4px 6px; border-radius:5px; cursor:pointer; user-select:none; }
  .ch:hover { background:#f1f5f9; }
  .ch input { accent-color:var(--accent); }
  .quick { display:flex; flex-wrap:wrap; gap:6px; margin-top:8px; }
  .chip { border:1px solid var(--border); background:#fff; border-radius:999px; padding:3px 10px; font-size:12px; cursor:pointer; color:var(--muted); }
  .chip:hover { border-color:var(--accent); color:var(--accent); }
  canvas {
    width: 100%;
    display: block;
    background: #fff;
    border: 1px solid var(--border);
    border-radius: 8px;
    /* height se fija en JS según nº de canales + DPR */
  }
  .status { font-size:12px; color:var(--muted); margin-top:8px; min-height:1.2em; }
  .status.ok { color:#047857; }
  .status.err { color:#b91c1c; }
</style>
</head>
<body>
<header>
  <h1>Butterfly / montaje apilado</h1>
  <span>$SUBJECT / $SESSION / task-$TASK_ID · fs=$(round(Int, store.fs)) Hz</span>
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
        <label>Escala (µV / canal)</label>
        <select id="scale">
          <option value="auto" selected>Auto</option>
          <option value="25">25 µV</option>
          <option value="50">50 µV</option>
          <option value="100">100 µV</option>
          <option value="200">200 µV</option>
        </select>
      </div>
      <button class="btn btn-primary" id="btn-plot" onclick="refreshPlot()">Actualizar</button>
      <button class="btn btn-save" id="btn-save" onclick="savePng()">Guardar PNG</button>
    </div>
  </div>

  <div class="card">
    <div style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:6px">
      <strong style="font-size:13px">Canales (orden de montaje)</strong>
      <span style="font-size:12px;color:var(--muted)" id="ch-count">—</span>
    </div>
    <div class="channels" id="channels">$ch_opts</div>
    <div class="quick">
      <button class="chip" type="button" onclick="selectAll(true)">Todos</button>
      <button class="chip" type="button" onclick="selectAll(false)">Ninguno</button>
      <button class="chip" type="button" onclick="preset(['Fz','Cz','Pz','Oz','Fp2'])">Línea media</button>
      <button class="chip" type="button" onclick="preset(['O1','Oz','O2','P7','P8'])">Occipital / parietal</button>
    </div>
  </div>

  <div class="card">
    <canvas id="cv"></canvas>
    <div class="status" id="status">Cargando montaje completo (0–10 s)…</div>
  </div>
</div>

<script>
const T_MIN = $(store.t_min);
const T_MAX = $(store.t_max);
const ALL_CHS = $(replace(repr(store.channels), "\"" => "'"));
const COLORS = [
  '#e31a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#666666',
  '#b2b21a','#00b0ba','#cc3366','#279158','#7340b3','#e6991a','#3366b3','#8c2626',
  '#1a9980','#b3661a','#5959b3','#999933','#cc4d8c','#408cbc','#bf4040','#66a64d',
  '#8c4d8c','#d97326','#33808c','#a63373','#4d7333','#8080bf','#e65959'
];

let _lastPlot = null;  // para redibujar al redimensionar

function selectedChannels() {
  return [...document.querySelectorAll('#channels input:checked')].map(el => el.value);
}
function updateCount() {
  const n = selectedChannels().length;
  document.getElementById('ch-count').textContent = n + ' / ' + ALL_CHS.length;
}
document.querySelectorAll('#channels input').forEach(el => el.addEventListener('change', updateCount));

function selectAll(on) {
  document.querySelectorAll('#channels input').forEach(el => el.checked = on);
  updateCount();
  if (on) refreshPlot();
}
function preset(chs) {
  document.querySelectorAll('#channels input').forEach(el => { el.checked = chs.includes(el.value); });
  updateCount();
  refreshPlot();
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

/** Canvas HiDPI: resolución interna = CSS × devicePixelRatio (texto nítido). */
function setupHiDPICanvas(canvas, cssH) {
  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || canvas.parentElement.clientWidth || 1100;
  canvas.style.height = cssH + 'px';
  canvas.width  = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);  // dibujar en coords CSS
  return { ctx, W: cssW, H: cssH, dpr };
}

async function refreshPlot() {
  const chs = selectedChannels();
  if (!chs.length) { setStatus('Selecciona al menos un canal.', 'err'); return; }
  const [t0, t1] = timeRange();
  const scale = document.getElementById('scale').value;
  setStatus('Cargando…');
  const url = '/api/signal?chs=' + encodeURIComponent(chs.join(',')) +
              '&t0=' + t0 + '&t1=' + t1 + '&scale=' + encodeURIComponent(scale) +
              '&max_points=2500';
  const res = await fetch(url);
  const data = await res.json();
  if (!res.ok) { setStatus(data.error || 'Error', 'err'); return; }
  _lastPlot = { data, chs, t0, t1 };
  drawButterfly(data, chs, t0, t1);
  const sp = data.spacing || 0;
  setStatus(chs.length + ' canales · ' + t0.toFixed(1) + '–' + t1.toFixed(1) +
            ' s · spacing ≈ ' + sp.toFixed(1) + ' µV');
}

function drawButterfly(data, chs, t0, t1) {
  const canvas = document.getElementById('cv');
  const n = chs.length;
  const rowPx = n <= 16 ? 34 : (n <= 24 ? 30 : 26);
  const cssH = Math.max(480, Math.min(rowPx * n + 90, 1600));
  const { ctx, W, H } = setupHiDPICanvas(canvas, cssH);

  // Tipografía proporcional al alto por canal
  const tickFs  = Math.max(10, Math.min(14, Math.round(0.45 * rowPx)));
  const titleFs = Math.max(13, Math.min(16, tickFs + 2));
  const labelFs = tickFs + 1;
  const pad = {
    l: Math.max(64, tickFs * 5.2),
    r: 18,
    t: titleFs + 22,
    b: labelFs + 28,
  };
  const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;

  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);

  const t = data.t || [];
  if (t.length < 2) return;
  const spacing = data.spacing || 50;
  const offsets = chs.map((_, i) => (n - 1 - i) * spacing);

  let yMin = Infinity, yMax = -Infinity;
  chs.forEach((ch, i) => {
    const y = data[ch] || [];
    const off = offsets[i];
    for (let k = 0; k < y.length; k++) {
      const v = y[k] + off;
      if (v < yMin) yMin = v;
      if (v > yMax) yMax = v;
    }
  });
  const yPad = Math.max(0.05 * (yMax - yMin), 0.25 * spacing);
  yMin -= yPad;
  yMax += yPad;

  const xAt = u => pad.l + ((u - t0) / (t1 - t0)) * pw;
  const yAt = v => pad.t + (1 - (v - yMin) / (yMax - yMin)) * ph;

  // grid + ticks X
  ctx.strokeStyle = '#e2e8f0'; ctx.lineWidth = 1;
  const xTicks = 4;
  for (let i = 0; i <= xTicks; i++) {
    const u = t0 + (i / xTicks) * (t1 - t0);
    const x = xAt(u);
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
    ctx.fillStyle = '#475569';
    ctx.font = tickFs + 'px "IBM Plex Sans", "Segoe UI", sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(u.toFixed(1), x, H - 12);
  }
  // baselines + labels Y
  chs.forEach((ch, i) => {
    const y = yAt(offsets[i]);
    ctx.strokeStyle = '#eef2f7';
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + pw, y); ctx.stroke();
    ctx.fillStyle = '#1e293b';
    ctx.font = tickFs + 'px "IBM Plex Sans", "Segoe UI", sans-serif';
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText(ch, pad.l - 10, y);
  });

  ctx.strokeStyle = '#64748b';
  ctx.lineWidth = 1.25;
  ctx.strokeRect(pad.l, pad.t, pw, ph);

  ctx.save();
  ctx.beginPath();
  ctx.rect(pad.l, pad.t, pw, ph);
  ctx.clip();

  chs.forEach((ch, i) => {
    const y = data[ch] || [];
    const off = offsets[i];
    ctx.strokeStyle = COLORS[i % COLORS.length];
    ctx.lineWidth = 1.1;
    ctx.beginPath();
    for (let k = 0; k < t.length; k++) {
      const x = xAt(t[k]), yy = yAt(y[k] + off);
      if (k === 0) ctx.moveTo(x, yy); else ctx.lineTo(x, yy);
    }
    ctx.stroke();
  });
  ctx.restore();

  // título y ejes
  ctx.fillStyle = '#0f172a';
  ctx.font = '600 ' + titleFs + 'px "IBM Plex Sans", "Segoe UI", sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'alphabetic';
  ctx.fillText('DATOS CRUDOS — sub-M05_ses-T2_task-EC — fs=$(round(Int, store.fs)) Hz, ' + n + ' chans',
               pad.l + pw / 2, titleFs + 4);

  ctx.fillStyle = '#334155';
  ctx.font = labelFs + 'px "IBM Plex Sans", "Segoe UI", sans-serif';
  ctx.fillText('Tiempo (s)', pad.l + pw / 2, H - 2);

  ctx.save();
  ctx.translate(16, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.fillText('Canal (offset)', 0, 0);
  ctx.restore();
}

async function savePng() {
  const chs = selectedChannels();
  if (!chs.length) { setStatus('Selecciona al menos un canal.', 'err'); return; }
  const [t0, t1] = timeRange();
  const scale = document.getElementById('scale').value;
  const btn = document.getElementById('btn-save');
  btn.disabled = true;
  setStatus('Generando PNG…');
  try {
    const res = await fetch('/api/save', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({channels: chs, t0, t1, scale})
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

let _resizeTimer = null;
window.addEventListener('resize', () => {
  clearTimeout(_resizeTimer);
  _resizeTimer = setTimeout(() => {
    if (_lastPlot) drawButterfly(_lastPlot.data, _lastPlot.chs, _lastPlot.t0, _lastPlot.t1);
  }, 120);
});

updateCount();
refreshPlot();
</script>
</body>
</html>
"""
end

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
    println("Cargando $(CSV_PATH)…")
    store = load_store(CSV_PATH)
    println("  $(length(store.channels)) canales · $(length(store.t)) muestras · ",
            "fs=$(round(store.fs; digits=1)) Hz · ",
            "$(round(store.t_min; digits=1))–$(round(store.t_max; digits=1)) s")

    server = listen(IPv4(HOST), PORT)
    url = "http://$HOST:$PORT/"
    println()
    println("UI butterfly → $url")
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
