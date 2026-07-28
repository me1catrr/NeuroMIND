#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor interactivo de señal cruda (GUI en navegador)
# ═══════════════════════════════════════════════════════════════
#
#  Abre una interfaz con botones/selectores para elegir canales,
#  rango temporal y escala. El botón «Guardar PNG» genera la
#  figura con CairoMakie en este mismo directorio.
#
#  Entrada:  results/subjects/sub-M05/ses-T2/eyesclosed/tables/raw_signal.csv
#  Salida:   raw_<canales>[_t0-t1].png (en este directorio)
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_raw.jl
#    # → http://127.0.0.1:8765/  (Ctrl+C para salir)
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_raw.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      23-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using CSV, DataFrames, CairoMakie, Sockets, Dates

const HERE         = @__DIR__
const RESULTS_UNIT = normpath(joinpath(HERE, "..", "..", "..", "..",
                     "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))
const CSV_PATH     = joinpath(RESULTS_UNIT, "tables", "raw_signal.csv")
const HOST         = "127.0.0.1"
const PORT         = 8765

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK    = "eyesclosed"

const LINE_COLORS = [
    :red, :blue, :green, :orange, :purple,
    :brown, :teal, :magenta, :olive, :navy,
    :coral, :darkcyan, :maroon, :darkgreen, :indigo,
]

# ── Datos en memoria ───────────────────────────────────────────

mutable struct RawStore
    t::Vector{Float64}
    channels::Vector{String}
    data::Dict{String,Vector{Float64}}  # canal → amplitud
    t_min::Float64
    t_max::Float64
end

function load_store(path::String)::RawStore
    isfile(path) || error("No encontrado: $path")
    df = CSV.read(path, DataFrame)
    hasproperty(df, :t_s) || error("Columna 't_s' ausente")
    t = Float64.(df.t_s)
    channels = String[string(n) for n in names(df) if string(n) != "t_s"]
    data = Dict{String,Vector{Float64}}(ch => Float64.(df[!, Symbol(ch)]) for ch in channels)
    return RawStore(t, channels, data, minimum(t), maximum(t))
end

# ── PNG (CairoMakie) ───────────────────────────────────────────

function save_png(
    store::RawStore,
    selected::Vector{String},
    t0::Float64,
    t1::Float64,
    scale::String,   # "auto" | "50" | "100" | "200" | "500"
)::String
    isempty(selected) && error("Selecciona al menos un canal")
    t0 >= t1 && error("t0 debe ser < t1")
    for ch in selected
        haskey(store.data, ch) || error("Canal desconocido: $ch")
    end

    mask = (store.t .>= t0) .& (store.t .<= t1)
    count(mask) < 2 && error("Ventana temporal vacía")
    t = store.t[mask]

    fig = Figure(size = (1100, 420), fontsize = 14)
    ax  = Axis(fig[1, 1];
        title  = "Señal cruda — $SUBJECT / $SESSION / $TASK",
        xlabel = "Tiempo (s)",
        ylabel = "Amplitud (µV)",
    )
    ax.xgridvisible = true
    ax.ygridvisible = true

    ymin = Inf
    ymax = -Inf
    for (i, ch) in enumerate(selected)
        y = store.data[ch][mask]
        col = LINE_COLORS[mod1(i, length(LINE_COLORS))]
        lines!(ax, t, y; color = col, linewidth = 1.1, label = ch)
        ymin = min(ymin, minimum(y))
        ymax = max(ymax, maximum(y))
    end
    axislegend(ax; position = :rt, framevisible = true, labelsize = 12)

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
    out_name = "raw_" * join(selected, "_") * tag_t * ".png"
    out_path = joinpath(HERE, out_name)
    save(out_path, fig; px_per_unit = 2)
    return out_path
end

# ── HTTP mínimo (stdlib Sockets) ───────────────────────────────

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

function _send_bytes(sock, status::Int, body::Vector{UInt8}; content_type::String)
    msg = "HTTP/1.1 $status\r\n" *
          "Content-Type: $content_type\r\n" *
          "Content-Length: $(length(body))\r\n" *
          "Connection: close\r\n" *
          "Access-Control-Allow-Origin: *\r\n" *
          "\r\n"
    write(sock, msg)
    write(sock, body)
end

"""Parseo mínimo de JSON: {\"channels\":[...],\"t0\":0,\"t1\":10,\"scale\":\"auto\"}."""
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
    return (
        channels = chs,
        t0 = _num("t0", 0.0),
        t1 = _num("t1", 10.0),
        scale = _str("scale", "auto"),
    )
end

function handle_request(sock, store::RawStore)
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
                   "\"channels\":[" * join(["\"$c\"" for c in store.channels], ",") * "]}"
            _send(sock, 200, body; content_type = "application/json")
        elseif method == "GET" && path_only == "/api/signal"
            # ?chs=Cz,Fp2&t0=0&t1=10&max_points=4000
            qs = Dict{String,String}()
            if occursin('?', path)
                for pair in split(split(path, '?', limit=2)[2], '&')
                    kv = split(pair, '=', limit=2)
                    length(kv) == 2 && (qs[kv[1]] = kv[2])
                end
            end
            chs = filter(!isempty, split(replace(get(qs, "chs", ""), "%2C" => ","), ','))
            t0 = something(tryparse(Float64, get(qs, "t0", "0")), store.t_min)
            t1 = something(tryparse(Float64, get(qs, "t1", string(store.t_max))), store.t_max)
            max_pts = something(tryparse(Int, get(qs, "max_points", "5000")), 5000)

            mask = findall(i -> store.t[i] >= t0 && store.t[i] <= t1, eachindex(store.t))
            if length(mask) > max_pts
                step = ceil(Int, length(mask) / max_pts)
                mask = mask[1:step:end]
            end
            t_out = store.t[mask]
            parts = String["\"t\":[" * join(round.(t_out; digits=4), ",") * "]"]
            for ch in chs
                haskey(store.data, ch) || continue
                ys = store.data[ch][mask]
                push!(parts, "\"$ch\":[" * join(round.(ys; digits=3), ",") * "]")
            end
            _send(sock, 200, "{" * join(parts, ",") * "}"; content_type = "application/json")
        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_save_json(raw)
            isempty(req.channels) && error("channels vacío")
            out = save_png(store, req.channels, req.t0, req.t1, req.scale)
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

# ── HTML UI ────────────────────────────────────────────────────

function html_page(store::RawStore)::String
    ch_opts = join(["<label class=\"ch\"><input type=\"checkbox\" value=\"$c\">$c</label>"
                    for c in store.channels], "\n")
    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — Señal cruda</title>
<style>
  :root {
    --bg:#f4f6f8; --card:#fff; --border:#d8dee6; --text:#1e293b;
    --muted:#64748b; --accent:#2563eb; --accent2:#0f766e;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; font-family: "IBM Plex Sans", "Segoe UI", sans-serif;
    background: var(--bg); color: var(--text);
  }
  header {
    padding: 14px 20px; background: #0f172a; color: #e2e8f0;
    display: flex; align-items: baseline; gap: 14px;
  }
  header h1 { margin: 0; font-size: 16px; font-weight: 600; }
  header span { font-size: 12px; color: #94a3b8; }
  .wrap { max-width: 1200px; margin: 16px auto; padding: 0 16px; }
  .card {
    background: var(--card); border: 1px solid var(--border);
    border-radius: 10px; padding: 14px 16px; margin-bottom: 14px;
  }
  .toolbar {
    display: flex; flex-wrap: wrap; gap: 12px 18px; align-items: end;
  }
  .field { display: flex; flex-direction: column; gap: 4px; }
  .field label { font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: .04em; }
  select, input[type=number] {
    height: 34px; padding: 0 10px; border: 1px solid var(--border);
    border-radius: 6px; background: #fff; font-size: 13px; min-width: 100px;
  }
  .btn {
    height: 34px; padding: 0 16px; border: none; border-radius: 6px;
    font-size: 13px; font-weight: 600; cursor: pointer;
  }
  .btn-primary { background: var(--accent); color: #fff; }
  .btn-primary:hover { background: #1d4ed8; }
  .btn-save { background: var(--accent2); color: #fff; }
  .btn-save:hover { background: #0d9488; }
  .btn:disabled { opacity: .5; cursor: wait; }
  .channels {
    display: grid; grid-template-columns: repeat(auto-fill, minmax(88px, 1fr));
    gap: 6px 8px; max-height: 160px; overflow: auto; padding: 4px 0;
  }
  .ch {
    display: flex; align-items: center; gap: 6px; font-size: 13px;
    padding: 4px 6px; border-radius: 5px; cursor: pointer; user-select: none;
  }
  .ch:hover { background: #f1f5f9; }
  .ch input { accent-color: var(--accent); }
  .quick { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 8px; }
  .chip {
    border: 1px solid var(--border); background: #fff; border-radius: 999px;
    padding: 3px 10px; font-size: 12px; cursor: pointer; color: var(--muted);
  }
  .chip:hover { border-color: var(--accent); color: var(--accent); }
  canvas {
    width: 100%; height: 380px; background: #fff;
    border: 1px solid var(--border); border-radius: 8px;
  }
  .status { font-size: 12px; color: var(--muted); margin-top: 8px; min-height: 1.2em; }
  .status.ok { color: #047857; }
  .status.err { color: #b91c1c; }
</style>
</head>
<body>
<header>
  <h1>Señal cruda</h1>
  <span>$SUBJECT / $SESSION / $TASK</span>
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
      <button class="btn btn-primary" id="btn-plot" onclick="refreshPlot()">Actualizar</button>
      <button class="btn btn-save" id="btn-save" onclick="savePng()">Guardar PNG</button>
    </div>
  </div>

  <div class="card">
    <div style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:6px">
      <strong style="font-size:13px">Canales</strong>
      <span style="font-size:12px;color:var(--muted)" id="ch-count">0 seleccionados</span>
    </div>
    <div class="channels" id="channels">
      $ch_opts
    </div>
    <div class="quick">
      <button class="chip" type="button" onclick="preset(['Fp2','Cz'])">Fp2 + Cz</button>
      <button class="chip" type="button" onclick="preset(['Cz','F3'])">Cz + F3</button>
      <button class="chip" type="button" onclick="preset(['O1','Oz','O2'])">Occipitales</button>
      <button class="chip" type="button" onclick="clearCh()">Ninguno</button>
    </div>
  </div>

  <div class="card">
    <canvas id="cv" width="1100" height="380"></canvas>
    <div class="status" id="status">Elige canales y pulsa Actualizar.</div>
  </div>
</div>

<script>
const T_MIN = $(store.t_min);
const T_MAX = $(store.t_max);
const COLORS = ['#dc2626','#2563eb','#16a34a','#ea580c','#7c3aed',
                '#92400e','#0d9488','#db2777','#65a30d','#1e3a8a'];

function selectedChannels() {
  return [...document.querySelectorAll('#channels input:checked')].map(el => el.value);
}
function updateCount() {
  const n = selectedChannels().length;
  document.getElementById('ch-count').textContent = n + ' seleccionado' + (n === 1 ? '' : 's');
}
document.querySelectorAll('#channels input').forEach(el => el.addEventListener('change', updateCount));

function preset(chs) {
  document.querySelectorAll('#channels input').forEach(el => {
    el.checked = chs.includes(el.value);
  });
  updateCount();
  refreshPlot();
}
function clearCh() {
  document.querySelectorAll('#channels input').forEach(el => el.checked = false);
  updateCount();
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
  if (!chs.length) { setStatus('Selecciona al menos un canal.', 'err'); return; }
  const [t0, t1] = timeRange();
  setStatus('Cargando…');
  const url = '/api/signal?chs=' + encodeURIComponent(chs.join(',')) +
              '&t0=' + t0 + '&t1=' + t1 + '&max_points=5000';
  const res = await fetch(url);
  const data = await res.json();
  if (!res.ok) { setStatus(data.error || 'Error', 'err'); return; }
  draw(data, chs, t0, t1);
  setStatus('Vista: ' + chs.join(', ') + ' · ' + t0.toFixed(1) + '–' + t1.toFixed(1) + ' s');
}

function draw(data, chs, t0, t1) {
  const canvas = document.getElementById('cv');
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  const pad = {l: 58, r: 18, t: 28, b: 42};
  const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);

  const t = data.t || [];
  if (t.length < 2) return;

  let ymin = Infinity, ymax = -Infinity;
  const scale = document.getElementById('scale').value;
  if (scale === 'auto') {
    for (const ch of chs) {
      const y = data[ch] || [];
      for (const v of y) { if (v < ymin) ymin = v; if (v > ymax) ymax = v; }
    }
    const padY = 0.05 * Math.max(ymax - ymin, 1);
    ymin -= padY; ymax += padY;
  } else {
    const lim = parseFloat(scale);
    ymin = -lim; ymax = lim;
  }

  const xAt = u => pad.l + ((u - t0) / (t1 - t0)) * pw;
  const yAt = v => pad.t + (1 - (v - ymin) / (ymax - ymin)) * ph;

  // grid
  ctx.strokeStyle = '#e2e8f0'; ctx.lineWidth = 1;
  const yTicks = 4;
  for (let i = 0; i <= yTicks; i++) {
    const v = ymin + (i / yTicks) * (ymax - ymin);
    const y = yAt(v);
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + pw, y); ctx.stroke();
    ctx.fillStyle = '#64748b'; ctx.font = '11px sans-serif';
    ctx.textAlign = 'right'; ctx.fillText(v.toFixed(0), pad.l - 6, y + 3);
  }
  const xTicks = 5;
  for (let i = 0; i <= xTicks; i++) {
    const u = t0 + (i / xTicks) * (t1 - t0);
    const x = xAt(u);
    ctx.strokeStyle = '#e2e8f0';
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
    ctx.fillStyle = '#64748b'; ctx.textAlign = 'center';
    ctx.fillText(u.toFixed(1), x, H - 16);
  }

  // frame
  ctx.strokeStyle = '#94a3b8';
  ctx.strokeRect(pad.l, pad.t, pw, ph);

  // series
  chs.forEach((ch, i) => {
    const y = data[ch] || [];
    ctx.strokeStyle = COLORS[i % COLORS.length];
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    for (let k = 0; k < t.length; k++) {
      const x = xAt(t[k]), yy = yAt(y[k]);
      if (k === 0) ctx.moveTo(x, yy); else ctx.lineTo(x, yy);
    }
    ctx.stroke();
  });

  // legend
  let lx = pad.l + 8, ly = pad.t + 14;
  chs.forEach((ch, i) => {
    ctx.fillStyle = COLORS[i % COLORS.length];
    ctx.fillRect(lx, ly - 7, 12, 3);
    ctx.fillStyle = '#1e293b'; ctx.font = '12px sans-serif'; ctx.textAlign = 'left';
    ctx.fillText(ch, lx + 16, ly);
    lx += ctx.measureText(ch).width + 36;
  });

  // axis labels
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

// defaults: Cz checked for single-channel preview feel; user multi-selects
preset(['Cz']);
</script>
</body>
</html>
"""
end

# ── Main ───────────────────────────────────────────────────────

function open_browser(url::String)
    try
        if Sys.isapple()
            run(`open $url`)
        elseif Sys.islinux()
            run(`xdg-open $url`)
        elseif Sys.iswindows()
            run(`cmd /c start $url`)
        end
    catch
        @warn "No se pudo abrir el navegador automáticamente"
    end
end

function main()
    println("Cargando $(CSV_PATH)…")
    store = load_store(CSV_PATH)
    println("  $(length(store.channels)) canales · $(length(store.t)) muestras · ",
            "$(round(store.t_min; digits=1))–$(round(store.t_max; digits=1)) s")

    server = listen(IPv4(HOST), PORT)
    url = "http://$HOST:$PORT/"
    println()
    println("UI interactiva → $url")
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
