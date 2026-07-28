#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor ICA: before / after por canal EEG
# ═══════════════════════════════════════════════════════════════
#
#  Compara la señal filtrada pre-ICA y la señal limpia post-ICA
#  en un canal seleccionado (primeros ~10 s exportados).
#
#  Entrada:
#    ../../tables/ica/ica_signal_before.csv
#    ../../tables/ica/ica_signal_after.csv
#    ../../json/ica_summary.json
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_ica_before_after.jl
#    # → http://127.0.0.1:8771/  (Ctrl+C para salir)
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_ica_before_after.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      24-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using CSV, DataFrames, CairoMakie, Sockets, Dates

const HERE    = @__DIR__
const ROOT    = normpath(joinpath(HERE, "..", "..", "..", "..",
                "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))   # …/eyesclosed
const TABLES  = joinpath(ROOT, "tables", "ica")
const JSON_P  = joinpath(ROOT, "json", "ica_summary.json")
const HOST    = "127.0.0.1"
const PORT    = 8771   # components=8770

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK    = "eyesclosed"

# ── Store ──────────────────────────────────────────────────────

mutable struct ICASignalStore
    t_sig::Vector{Float64}
    channels::Vector{String}
    before::Dict{String,Vector{Float64}}
    after::Dict{String,Vector{Float64}}
    summary::Dict{String,Any}
    t_min::Float64
    t_max::Float64
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
        "artifact_threshold" => NaN, "converged" => false,
    )
    isfile(path) || return d
    s = read(path, String)
    d["n_components"] = _json_int(s, "n_components", 0)
    d["n_rejected"] = _json_int(s, "n_rejected", 0)
    d["n_accepted"] = _json_int(s, "n_accepted", 0)
    d["variance_retained"] = _json_num(s, "variance_retained", NaN)
    d["rejection_source"] = _json_str(s, "rejection_source", "?")
    d["artifact_threshold"] = _json_num(s, "artifact_threshold", NaN)
    d["converged"] = _json_bool(s, "converged", false)
    d["n_iter"] = _json_int(s, "n_iter", 0)
    return d
end

function load_store()::ICASignalStore
    bef_p = joinpath(TABLES, "ica_signal_before.csv")
    aft_p = joinpath(TABLES, "ica_signal_after.csv")
    isfile(bef_p) || error("No encontrado: $bef_p")
    isfile(aft_p) || error("No encontrado: $aft_p")

    bdf = CSV.read(bef_p, DataFrame)
    adf = CSV.read(aft_p, DataFrame)
    t_sig = Float64.(bdf.t_s)
    channels = String[string(n) for n in names(bdf) if string(n) != "t_s"]
    before = Dict{String,Vector{Float64}}()
    after  = Dict{String,Vector{Float64}}()
    for ch in channels
        before[ch] = Float64.(bdf[!, Symbol(ch)])
        hasproperty(adf, Symbol(ch)) && (after[ch] = Float64.(adf[!, Symbol(ch)]))
    end

    summary = load_summary(JSON_P)
    t_min = isempty(t_sig) ? 0.0 : minimum(t_sig)
    t_max = isempty(t_sig) ? 10.0 : maximum(t_sig)

    return ICASignalStore(t_sig, channels, before, after, summary, t_min, t_max)
end

# ── PNG (CairoMakie) ───────────────────────────────────────────

function save_signal_png(
    store::ICASignalStore,
    ch::String,
    t0::Float64,
    t1::Float64,
)::String
    (!haskey(store.before, ch) || !haskey(store.after, ch)) && error("Canal no disponible: $ch")
    mask = (store.t_sig .>= t0) .& (store.t_sig .<= t1)
    count(mask) < 2 && error("Ventana vacía")
    t = store.t_sig[mask]

    fig = Figure(size = (1100, 420), fontsize = 13)
    ax = Axis(fig[1, 1];
        title = "Before / after ICA — $ch — $SUBJECT / $SESSION / $TASK",
        xlabel = "Tiempo (s)", ylabel = "Amplitud (µV)",
    )
    ax.xgridvisible = true; ax.ygridvisible = true
    lines!(ax, t, store.before[ch][mask]; color = :gray, linewidth = 1.4, label = "$ch before")
    lines!(ax, t, store.after[ch][mask]; color = :purple, linewidth = 1.4, label = "$ch after")
    axislegend(ax; position = :rt, framevisible = true, labelsize = 11)
    xlims!(ax, t0, t1)

    out = joinpath(HERE, "ica_ba_$(ch)_$(round(Int,t0))-$(round(Int,t1))s.png")
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
    _str = (key, default) -> begin
        m2 = match(Regex("\"" * key * "\"\\s*:\\s*\"([^\"]*)\""), s)
        m2 === nothing ? String(default) : String(m2.captures[1])
    end
    _num = (key, default) -> begin
        m2 = match(Regex("\"" * key * "\"\\s*:\\s*([-\\d.eE]+)"), s)
        m2 === nothing && return Float64(default)
        something(tryparse(Float64, m2.captures[1]), Float64(default))
    end
    return (ch = _str("ch", ""), t0 = _num("t0", 0.0), t1 = _num("t1", 10.0))
end

function handle_request(sock, store::ICASignalStore)
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
            body = "{" *
                "\"t_min\":$(store.t_min),\"t_max\":$(store.t_max)," *
                "\"n_rejected\":$(s["n_rejected"]),\"n_accepted\":$(s["n_accepted"])," *
                "\"n_components\":$(s["n_components"])," *
                "\"variance_retained\":$(s["variance_retained"])," *
                "\"rejection_source\":\"$(s["rejection_source"])\"," *
                "\"channels\":[" * join(["\"$c\"" for c in store.channels], ",") * "]" *
                "}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/signal"
            qs = _parse_qs(path)
            ch = get(qs, "ch", isempty(store.channels) ? "" : store.channels[1])
            t0 = something(tryparse(Float64, get(qs, "t0", "0")), store.t_min)
            t1 = something(tryparse(Float64, get(qs, "t1", string(store.t_max))), store.t_max)
            max_pts = something(tryparse(Int, get(qs, "max_points", "4000")), 4000)
            (!haskey(store.before, ch) || !haskey(store.after, ch)) && error("Canal no disponible: $ch")
            mask = findall(i -> store.t_sig[i] >= t0 && store.t_sig[i] <= t1, eachindex(store.t_sig))
            if length(mask) > max_pts
                step = ceil(Int, length(mask) / max_pts)
                mask = mask[1:step:end]
            end
            body = "{\"t\":[" * join(round.(store.t_sig[mask]; digits=4), ",") * "]," *
                   "\"before\":[" * join(round.(store.before[ch][mask]; digits=3), ",") * "]," *
                   "\"after\":[" * join(round.(store.after[ch][mask]; digits=3), ",") * "]}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_save_json(raw)
            ch = isempty(req.ch) ? (isempty(store.channels) ? "" : store.channels[1]) : req.ch
            out = save_signal_png(store, ch, req.t0, req.t1)
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

function html_page(store::ICASignalStore)::String
    s = store.summary
    ch_opts = join(["<option value=\"$c\">$c</option>" for c in store.channels], "\n")

    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — ICA before/after</title>
<style>
  :root {
    --bg:#f4f6f8; --card:#fff; --border:#d8dee6; --text:#1e293b;
    --muted:#64748b; --accent:#7c3aed; --accent2:#0f766e;
  }
  * { box-sizing: border-box; }
  body { margin:0; font-family:"IBM Plex Sans","Segoe UI",sans-serif; background:var(--bg); color:var(--text); }
  header { padding:14px 20px; background:#1e1b4b; color:#e2e8f0; display:flex; flex-wrap:wrap; gap:12px 18px; align-items:baseline; }
  header h1 { margin:0; font-size:16px; font-weight:600; }
  header span { font-size:12px; color:#a5b4fc; }
  .wrap { max-width:1100px; margin:16px auto; padding:0 16px; }
  .card { background:var(--card); border:1px solid var(--border); border-radius:10px; padding:12px 14px; margin-bottom:14px; }
  .toolbar { display:flex; flex-wrap:wrap; gap:10px 14px; align-items:end; }
  .field { display:flex; flex-direction:column; gap:4px; }
  .field label { font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:.04em; }
  select, input[type=number] {
    height:34px; padding:0 10px; border:1px solid var(--border);
    border-radius:6px; background:#fff; font-size:13px; min-width:90px;
  }
  .btn { height:34px; padding:0 14px; border:none; border-radius:6px; font-size:13px; font-weight:600; cursor:pointer; }
  .btn-primary { background:var(--accent); color:#fff; }
  .btn-save { background:var(--accent2); color:#fff; }
  .btn:disabled { opacity:.5; cursor:wait; }
  canvas { width:100%; display:block; background:#fff; border:1px solid var(--border); border-radius:8px; }
  .status { font-size:12px; color:var(--muted); margin-top:8px; min-height:1.2em; }
  .status.ok { color:#047857; } .status.err { color:#b91c1c; }
  .hint { font-size:12px; color:var(--muted); }
</style>
</head>
<body>
<header>
  <h1>ICA · Before / after</h1>
  <span>$SUBJECT / $SESSION / $TASK</span>
  <span>$(s["n_components"]) ICs · rechazados $(s["n_rejected"]) · retenida $(s["variance_retained"])% · fuente $(s["rejection_source"])</span>
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
      <div class="field"><label>Canal</label>
        <select id="ch">$ch_opts</select></div>
      <button class="btn btn-primary" onclick="refreshSig()">Actualizar</button>
      <button class="btn btn-save" id="btn-save" onclick="savePng()">Guardar PNG</button>
    </div>
    <p class="hint" style="margin:8px 0 0">Primeros ~10 s exportados por el pipeline (<code>ica_signal_before/after.csv</code>). Gris = before · violeta = after.</p>
  </div>

  <div class="card">
    <strong style="font-size:13px">Before / after ICA</strong>
    <canvas id="cv-sig" height="360"></canvas>
    <div class="status" id="status">Selecciona un canal y pulsa Actualizar.</div>
  </div>
</div>

<script>
const T_MIN = $(store.t_min), T_MAX = $(store.t_max);

function setStatus(msg, kind) {
  const el = document.getElementById('status');
  el.textContent = msg; el.className = 'status' + (kind ? ' ' + kind : '');
}

function timeRange() {
  const win = document.getElementById('win').value;
  let t0 = parseFloat(document.getElementById('t0').value);
  if (Number.isNaN(t0)) t0 = T_MIN;
  t0 = Math.max(T_MIN, Math.min(t0, T_MAX));
  if (win === 'full') return [T_MIN, T_MAX];
  const w = parseFloat(win);
  let t1 = Math.min(T_MAX, t0 + w);
  return [t0, t1];
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

async function refreshSig() {
  const ch = document.getElementById('ch').value;
  if (!ch) return;
  const [t0, t1] = timeRange();
  const res = await fetch('/api/signal?ch=' + encodeURIComponent(ch) + '&t0=' + t0 + '&t1=' + t1);
  const data = await res.json();
  if (!res.ok) { setStatus(data.error || 'Error signal', 'err'); return; }
  drawSeries('cv-sig', 360, data.t, [
    { y: data.before, color: '#94a3b8', label: ch + ' before', lw: 1.4 },
    { y: data.after,  color: '#7c3aed', label: ch + ' after',  lw: 1.4 },
  ], t0, t1);
  setStatus('Canal ' + ch + ' · ventana ' + t0.toFixed(1) + '–' + t1.toFixed(1) + ' s');
}

async function savePng() {
  const ch = document.getElementById('ch').value;
  const [t0, t1] = timeRange();
  const btn = document.getElementById('btn-save');
  btn.disabled = true;
  setStatus('Generando PNG…');
  try {
    const res = await fetch('/api/save', {
      method:'POST', headers:{'Content-Type':'application/json'},
      body: JSON.stringify({ch, t0, t1})
    });
    const data = await res.json();
    if (!res.ok || !data.ok) throw new Error(data.error || 'Error');
    setStatus('Guardado: ' + data.file, 'ok');
  } catch (e) {
    setStatus(String(e.message || e), 'err');
  } finally { btn.disabled = false; }
}

document.getElementById('ch').addEventListener('change', refreshSig);
(function(){
  const sel = document.getElementById('ch');
  for (const o of sel.options) if (o.value === 'Cz') { sel.value = 'Cz'; break; }
})();
refreshSig();
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
    println("Cargando before/after ICA desde $TABLES …")
    store = load_store()
    s = store.summary
    println("  $(length(store.channels)) canales · $(length(store.t_sig)) muestras")
    println("  ICs: $(s["n_components"]) · rechazados $(s["n_rejected"]) · retenida $(s["variance_retained"])%")

    server = listen(IPv4(HOST), PORT)
    url = "http://$HOST:$PORT/"
    println()
    println("UI ICA before/after → $url")
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
