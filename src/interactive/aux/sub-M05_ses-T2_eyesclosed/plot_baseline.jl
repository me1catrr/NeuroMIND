#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Baseline before / after (por canal y época)
# ═══════════════════════════════════════════════════════════════
#
#  Compara un canal en una época:
#    Antes  = post-ICA segmentada (sin baseline)
#    Después = apply_baseline first_window_mean [0, baseline_end_s]
#
#  Esa 1ª pasada es la que alimenta el AR ±70 µV. Con n_passes=2
#  el mismo método se reaplica tras el rechazo sobre épocas válidas
#  (efecto visual pequeño si la ventana ya tiene media ≈ 0).
#
#  Por defecto: canal C4, época con mayor |offset| medio en la
#  ventana de baseline (como baseline_before_after_C4_ep089.png).
#
#  Entrada:
#    ../../cache/ica_result.jls
#    ../../ica_labels_auto.csv
#    ../../tables/segments_table.csv   (status valid/rejected)
#    ../../config_snapshot.toml        (baseline_end_s, opcional)
#
#  Salida:
#    baseline_before_after_<ch>_epXXX.png
#    baseline_offset_by_channel.csv
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_baseline.jl
#    # → http://127.0.0.1:8773/
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_baseline.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      24-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using NeuroMIND
using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics, Serialization

const HERE   = @__DIR__
const ROOT   = normpath(joinpath(HERE, "..", "..", "..", "..",
               "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))
const TABLES = joinpath(ROOT, "tables")
const HOST   = "127.0.0.1"
const PORT   = 8773

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK    = "eyesclosed"
const TASK_ID = "EC"

const DEFAULT_CHANNEL = "C4"
const EPOCH_LEN_S     = 1.0
const DEFAULT_BL_END  = 0.10

# ── Store ──────────────────────────────────────────────────────

mutable struct BaselineStore
    channels::Vector{String}
    clean::Matrix{Float64}       # ch × samples, post-ICA
    fs::Float64
    duration_s::Float64
    n_epochs::Int
    status::Vector{String}       # por época (1-indexed)
    baseline_end_s::Float64
    default_epoch::Int
    default_channel::String
    rejected_ics::Vector{Int}
    offset_table::DataFrame      # stats |μ| antes/después por canal
    n_valid::Int
end

function _load_rejected_ics(root::String)::Vector{Int}
    auto_p = joinpath(root, "ica_labels_auto.csv")
    if isfile(auto_p)
        df = CSV.read(auto_p, DataFrame)
        hasproperty(df, :label) || return Int[]
        return Int.(df.component[String.(df.label) .== "artifact"])
    end
    return Int[]
end

function _baseline_end_s(root::String)::Float64
    snap = joinpath(root, "config_snapshot.toml")
    isfile(snap) || return DEFAULT_BL_END
    s = read(snap, String)
    m = match(r"baseline_end_s\s*=\s*([0-9.]+)", s)
    m === nothing && return DEFAULT_BL_END
    v = tryparse(Float64, m.captures[1])
    return v === nothing ? DEFAULT_BL_END : v
end

function _reconstruct_clean(root::String)
    cache_p = joinpath(root, "cache", "ica_result.jls")
    isfile(cache_p) || error("No encontrado: $cache_p")
    ica = deserialize(cache_p)
    rej = _load_rejected_ics(root)
    n_comp = size(ica.activations, 1)
    clean = if isempty(rej)
        ica.mixing_matrix * ica.activations
    else
        keep = setdiff(1:n_comp, rej)
        ica.mixing_matrix[:, keep] * ica.activations[keep, :]
    end
    return clean, String.(ica.meta.channel_names), Float64(ica.meta.fs), rej
end

"""Recorta época [start, end) SIN baseline."""
function _epoch_raw(clean::Matrix{Float64}, fs::Float64, epoch::Int;
                    epoch_len_s::Float64 = EPOCH_LEN_S)::Matrix{Float64}
    start_s = (epoch - 1) * epoch_len_s
    end_s   = start_s + epoch_len_s
    n_samp  = size(clean, 2)
    i0 = clamp(Int(round(start_s * fs)) + 1, 1, n_samp)
    i1 = clamp(Int(round(end_s * fs)), i0, n_samp)
    return copy(clean[:, i0:i1])
end

"""Aplica first_window_mean sobre (ch × samples)."""
function _apply_bl(ep::Matrix{Float64}, fs::Float64, baseline_end_s::Float64)::Matrix{Float64}
    out = copy(ep)
    n_bl = clamp(Int(round(baseline_end_s * fs)), 1, size(out, 2))
    out .-= mean(out[:, 1:n_bl]; dims = 2)
    return out
end

"""Índice de época (1-based) con mayor |media| en ventana baseline del canal."""
function _epoch_max_offset(
    clean::Matrix{Float64},
    fs::Float64,
    ci::Int,
    n_epochs::Int,
    baseline_end_s::Float64,
)::Int
    n_bl = clamp(Int(round(baseline_end_s * fs)), 1, Int(round(EPOCH_LEN_S * fs)))
    best_ep, best_abs = 1, -Inf
    for ep in 1:n_epochs
        raw = _epoch_raw(clean, fs, ep)
        size(raw, 2) < n_bl && continue
        μ = mean(raw[ci, 1:n_bl])
        a = abs(μ)
        if a > best_abs
            best_abs = a
            best_ep = ep
        end
    end
    return best_ep
end

"""
Tabla de |offset| medio/máximo por canal en la ventana de baseline,
antes y después de `apply_baseline`, solo épocas válidas.
"""
function compute_offset_table(store_like)::DataFrame
    clean = store_like.clean
    fs = store_like.fs
    channels = store_like.channels
    status = store_like.status
    bl_end = store_like.baseline_end_s
    n_epochs = store_like.n_epochs
    n_ch = length(channels)
    n_bl = clamp(Int(round(bl_end * fs)), 1, Int(round(EPOCH_LEN_S * fs)))
    valid = findall(s -> s == "valid", status)
    isempty(valid) && (valid = collect(1:n_epochs))

    mean_b = zeros(n_ch)
    max_b  = zeros(n_ch)
    mean_a = zeros(n_ch)
    max_a  = zeros(n_ch)

    for (ci, _) in enumerate(channels)
        abs_b = Float64[]
        abs_a = Float64[]
        sizehint!(abs_b, length(valid))
        sizehint!(abs_a, length(valid))
        for ep in valid
            raw = _epoch_raw(clean, fs, ep)
            size(raw, 2) < n_bl && continue
            bl = _apply_bl(raw, fs, bl_end)
            push!(abs_b, abs(mean(raw[ci, 1:n_bl])))
            push!(abs_a, abs(mean(bl[ci, 1:n_bl])))
        end
        mean_b[ci] = mean(abs_b)
        max_b[ci]  = maximum(abs_b)
        mean_a[ci] = mean(abs_a)
        max_a[ci]  = maximum(abs_a)
    end

    df = DataFrame(
        channel = channels,
        offset_mean_before_uv = round.(mean_b; digits = 4),
        offset_max_before_uv  = round.(max_b; digits = 4),
        offset_mean_after_uv  = mean_a,   # ~1e-16; no redondear a 4
        offset_max_after_uv   = max_a,
    )
    sort!(df, :offset_mean_before_uv; rev = true)
    return df
end

function save_offset_table_csv(store::BaselineStore)::String
    out = joinpath(HERE, "baseline_offset_by_channel.csv")
    CSV.write(out, store.offset_table)
    return out
end

function load_store()::BaselineStore
    println("  Reconstruyendo señal post-ICA…")
    clean, channels, fs, rej_ics = _reconstruct_clean(ROOT)
    bl_end = _baseline_end_s(ROOT)
    n_epochs = Int(floor(size(clean, 2) / (EPOCH_LEN_S * fs)))
    n_epochs >= 1 || error("No hay épocas completas en la señal")

    status = fill("valid", n_epochs)
    seg_p = joinpath(TABLES, "segments_table.csv")
    if isfile(seg_p)
        seg = CSV.read(seg_p, DataFrame)
        for r in eachrow(seg)
            e = Int(r.epoch)
            1 <= e <= n_epochs || continue
            status[e] = String(r.status)
        end
    end
    n_valid = count(==("valid"), status)

    ci = findfirst(==(DEFAULT_CHANNEL), channels)
    ci === nothing && error("Canal por defecto no disponible: $DEFAULT_CHANNEL")
    def_ep = _epoch_max_offset(clean, fs, ci, n_epochs, bl_end)

    # Placeholder store para compute_offset_table
    tmp = (
        clean = clean, fs = fs, channels = channels, status = status,
        baseline_end_s = bl_end, n_epochs = n_epochs,
    )
    println("  Calculando tabla de offsets (|μ| ventana baseline)…")
    offset_table = compute_offset_table(tmp)

    return BaselineStore(
        channels, clean, fs, size(clean, 2) / fs, n_epochs, status,
        bl_end, def_ep, DEFAULT_CHANNEL, rej_ics, offset_table, n_valid,
    )
end

# ── Plot ───────────────────────────────────────────────────────

function plot_baseline_before_after(
    store::BaselineStore;
    channel::String = store.default_channel,
    epoch::Int = store.default_epoch,
)::Figure
    ch = String(channel)
    ci = findfirst(==(ch), store.channels)
    ci === nothing && error("Canal no disponible: $ch")
    (1 <= epoch <= store.n_epochs) || error("Época fuera de rango: $epoch")

    raw = _epoch_raw(store.clean, store.fs, epoch)
    bl  = _apply_bl(raw, store.fs, store.baseline_end_s)
    n_samp = size(raw, 2)
    t = collect(range(0.0; step = 1 / store.fs, length = n_samp))
    y0 = raw[ci, :]
    y1 = bl[ci, :]
    st = store.status[epoch]
    bl_end = store.baseline_end_s

    n_bl = clamp(Int(round(bl_end * store.fs)), 1, n_samp)
    μ0 = mean(y0[1:n_bl])
    μ1 = mean(y1[1:n_bl])

    fig = Figure(size = (820, 420), fontsize = 13)
    ax = Axis(fig[1, 1];
        title = "Baseline — canal $ch, época $epoch ($st) · " *
                "offset $(round(μ0; digits=2)) → $(round(μ1; digits=2)) µV",
        xlabel = "Tiempo dentro de la época (s)",
        ylabel = "Amplitud (µV)",
    )
    ax.xgridvisible = true
    ax.ygridvisible = true
    vspan!(ax, 0.0, bl_end; color = (:gray, 0.18))
    hlines!(ax, [0.0]; color = (:gray, 0.65), linestyle = :dash, linewidth = 1.2)
    lines!(ax, t, y0; color = RGBf(0.85, 0.55, 0.10), linewidth = 1.6, label = "antes")
    lines!(ax, t, y1; color = RGBf(0.10, 0.55, 0.55), linewidth = 1.8, label = "después")
    axislegend(ax; position = :rt, framevisible = true, labelsize = 11)
    xlims!(ax, 0.0, maximum(t))

    Label(fig[0, 1],
        "$SUBJECT / $SESSION / task-$TASK_ID · first_window_mean [0, $(bl_end)] s",
        fontsize = 11, color = :gray, tellwidth = false)

    return fig
end

function save_baseline_png(
    store::BaselineStore;
    channel::String = store.default_channel,
    epoch::Int = store.default_epoch,
)::String
    fig = plot_baseline_before_after(store; channel = channel, epoch = epoch)
    out = joinpath(HERE, "baseline_before_after_$(channel)_ep$(lpad(string(epoch), 3, '0')).png")
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
        nb = readbytes!(sock, view(buf, filled + 1:n))
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
    for pair in split(split(path, '?', limit = 2)[2], '&')
        kv = split(pair, '=', limit = 2)
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
    return (ch = _str("ch", DEFAULT_CHANNEL), epoch = Int(round(_num("epoch", 89))))
end

function _listen_available(host::String, port::Int; n_try::Int = 10)
    last_err = nothing
    for p in port:(port + n_try - 1)
        try
            return listen(IPv4(host), p), p
        catch e
            last_err = e
            if e isa Base.IOError && occursin("EADDRINUSE", sprint(showerror, e))
                @warn "Puerto $p ocupado · probando $(p + 1)…"
                continue
            end
            rethrow(e)
        end
    end
    error("No hay puerto libre en $(port)-$(port + n_try - 1). Último error: $last_err")
end

function handle_request(sock, store::BaselineStore)
    req_line = try readline(sock) catch; return end
    isempty(strip(req_line, ['\r', '\n', ' '])) && return
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
            body = "{" *
                "\"n_epochs\":$(store.n_epochs)," *
                "\"n_valid\":$(store.n_valid)," *
                "\"fs\":$(store.fs)," *
                "\"baseline_end_s\":$(store.baseline_end_s)," *
                "\"default_epoch\":$(store.default_epoch)," *
                "\"default_channel\":\"$(store.default_channel)\"," *
                "\"channels\":[" * join(["\"$c\"" for c in store.channels], ",") * "]" *
                "}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/offsets"
            df = store.offset_table
            rows = String[]
            for r in eachrow(df)
                push!(rows,
                    "{\"channel\":\"$(r.channel)\"," *
                    "\"mean_before\":$(r.offset_mean_before_uv)," *
                    "\"max_before\":$(r.offset_max_before_uv)," *
                    "\"mean_after\":$(r.offset_mean_after_uv)," *
                    "\"max_after\":$(r.offset_max_after_uv)}")
            end
            g_mean_b = mean(df.offset_mean_before_uv)
            g_mean_a = mean(df.offset_mean_after_uv)
            body = "{" *
                "\"n_valid\":$(store.n_valid)," *
                "\"baseline_end_s\":$(store.baseline_end_s)," *
                "\"global_mean_before\":$(round(g_mean_b; digits=4))," *
                "\"global_mean_after\":$g_mean_a," *
                "\"rows\":[" * join(rows, ",") * "]" *
                "}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/epoch"
            qs = _parse_qs(path)
            ch = get(qs, "ch", store.default_channel)
            epoch = something(tryparse(Int, get(qs, "epoch", string(store.default_epoch))),
                              store.default_epoch)
            ci = findfirst(==(ch), store.channels)
            ci === nothing && error("Canal no disponible: $ch")
            (1 <= epoch <= store.n_epochs) || error("Época fuera de rango")

            raw = _epoch_raw(store.clean, store.fs, epoch)
            bl  = _apply_bl(raw, store.fs, store.baseline_end_s)
            n_samp = size(raw, 2)
            max_pts = 500
            step = n_samp > max_pts ? ceil(Int, n_samp / max_pts) : 1
            idx = collect(1:step:n_samp)
            t = collect(range(0.0; step = 1 / store.fs, length = n_samp))[idx]
            y0 = round.(raw[ci, idx]; digits = 3)
            y1 = round.(bl[ci, idx]; digits = 3)
            n_bl = clamp(Int(round(store.baseline_end_s * store.fs)), 1, n_samp)
            μ0 = round(mean(raw[ci, 1:n_bl]); digits = 3)
            μ1 = round(mean(bl[ci, 1:n_bl]); digits = 3)
            body = "{" *
                "\"ch\":\"$ch\",\"epoch\":$epoch," *
                "\"status\":\"$(store.status[epoch])\"," *
                "\"baseline_end_s\":$(store.baseline_end_s)," *
                "\"offset_before\":$μ0,\"offset_after\":$μ1," *
                "\"t\":[" * join(round.(t; digits = 4), ",") * "]," *
                "\"before\":[" * join(y0, ",") * "]," *
                "\"after\":[" * join(y1, ",") * "]" *
                "}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_save_json(raw)
            out = save_baseline_png(store; channel = req.ch, epoch = req.epoch)
            body = "{\"ok\":true,\"file\":$(repr(basename(out)))}"
            _send(sock, 200, body; content_type = "application/json")
            println("[$(Dates.format(now(), "HH:MM:SS"))] PNG → $out")
        else
            _send(sock, 404, "{\"error\":\"not found\"}"; content_type = "application/json")
        end
    catch e
        msg = sprint(showerror, e)
        @warn "Request error" exception = e
        _send(sock, 400, "{\"ok\":false,\"error\":$(repr(msg))}"; content_type = "application/json")
    end
end

function html_page(store::BaselineStore)::String
    ch_opts = join(
        ["<option value=\"$c\"$(c == store.default_channel ? " selected" : "")>$c</option>"
         for c in store.channels], "\n")
    def_ep = store.default_epoch
    bl_end = store.baseline_end_s

    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — Baseline before/after</title>
<style>
  :root {
    --bg:#f4f6f8; --card:#fff; --border:#d8dee6; --text:#1e293b;
    --muted:#64748b; --accent:#0f766e; --accent2:#b45309;
  }
  * { box-sizing: border-box; }
  body { margin:0; font-family:"IBM Plex Sans","Segoe UI",sans-serif; background:var(--bg); color:var(--text); }
  header { padding:14px 20px; background:#134e4a; color:#ecfdf5; display:flex; flex-wrap:wrap; gap:12px 18px; align-items:baseline; }
  header h1 { margin:0; font-size:16px; font-weight:600; }
  header span { font-size:12px; color:#99f6e4; }
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
  canvas { width:100%; display:block; background:#fff; border:1px solid var(--border); border-radius:8px; }
  .status { font-size:12px; color:var(--muted); margin-top:8px; min-height:1.2em; }
  .status.ok { color:#047857; } .status.err { color:#b91c1c; }
  .hint { font-size:12px; color:var(--muted); margin:8px 0 0; }
  table.offsets { width:100%; border-collapse:collapse; font-size:13px; }
  table.offsets th, table.offsets td {
    padding:6px 10px; border-bottom:1px solid var(--border); text-align:right;
  }
  table.offsets th:first-child, table.offsets td:first-child { text-align:left; }
  table.offsets th { font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:.03em; font-weight:600; }
  table.offsets tr.highlight td { background:#ecfdf5; font-weight:600; color:#134e4a; }
  table.offsets tr.summary td { background:#f8fafc; font-weight:600; border-top:2px solid var(--border); }
  .tbl-wrap { max-height:360px; overflow:auto; }
</style>
</head>
<body>
<header>
  <h1>Baseline · before / after</h1>
  <span>$SUBJECT / $SESSION / $TASK</span>
  <span>first_window_mean 0–$(bl_end)s · default ep $def_ep (máx |offset|) · $(store.n_valid) válidas</span>
</header>
<div class="wrap">
  <div class="card">
    <div class="toolbar">
      <div class="field"><label>Canal</label>
        <select id="ch">$ch_opts</select></div>
      <div class="field"><label>Época</label>
        <input type="number" id="epoch" value="$def_ep" min="1" max="$(store.n_epochs)" step="1"/></div>
      <button class="btn btn-primary" onclick="refresh()">Actualizar</button>
      <button class="btn btn-save" onclick="savePng()">Guardar PNG</button>
    </div>
    <p class="hint">Antes = post-ICA sin baseline · Después = media de [0, $(bl_end)] s restada (misma corrección que alimenta el AR; n_passes=2 reaplica el método tras el rechazo).</p>
  </div>
  <div class="card">
    <strong id="title" style="font-size:13px">Before / after</strong>
    <canvas id="cv"></canvas>
    <div class="status" id="status"></div>
  </div>
  <div class="card">
    <strong style="font-size:13px">|offset| por canal (ventana baseline, épocas válidas)</strong>
    <p class="hint" id="tbl-caption" style="margin:6px 0 10px"></p>
    <div class="tbl-wrap">
      <table class="offsets">
        <thead>
          <tr>
            <th>Canal</th>
            <th>|offset| medio antes (µV)</th>
            <th>|offset| máx. antes (µV)</th>
            <th>|offset| medio después (µV)</th>
            <th>|offset| máx. después (µV)</th>
          </tr>
        </thead>
        <tbody id="offset-body"></tbody>
      </table>
    </div>
  </div>
</div>
<script>
const BL_END = $bl_end;
const DEFAULT_CH = '$(store.default_channel)';
let data = null;

function setStatus(msg, cls) {
  const el = document.getElementById('status');
  el.textContent = msg || '';
  el.className = 'status' + (cls ? ' ' + cls : '');
}

function setupHiDPICanvas(canvas, cssH) {
  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || canvas.parentElement.clientWidth || 1000;
  canvas.style.height = cssH + 'px';
  canvas.width = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, W: cssW, H: cssH };
}

function draw() {
  if (!data) return;
  const { ctx, W, H } = setupHiDPICanvas(document.getElementById('cv'), 420);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);

  const t = data.t, y0 = data.before, y1 = data.after;
  let ymin = Infinity, ymax = -Infinity;
  for (const v of y0.concat(y1)) { ymin = Math.min(ymin, v); ymax = Math.max(ymax, v); }
  const padY = 0.08 * (ymax - ymin || 1);
  ymin -= padY; ymax += padY;

  const pad = { l: 58, r: 20, t: 36, b: 48 };
  const left = pad.l, top = pad.t;
  const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
  const xAt = (ti) => left + ((ti - t[0]) / (t[t.length - 1] - t[0] || 1)) * pw;
  const yAt = (v) => top + ph - ((v - ymin) / (ymax - ymin || 1)) * ph;

  // shade baseline window
  const xBl = xAt(Math.min(BL_END, t[t.length - 1]));
  ctx.fillStyle = 'rgba(100,116,139,0.15)';
  ctx.fillRect(left, top, Math.max(0, xBl - left), ph);

  // grid
  ctx.strokeStyle = '#eef2f7';
  for (let g = 0; g <= 4; g++) {
    const yy = top + (g / 4) * ph;
    ctx.beginPath(); ctx.moveTo(left, yy); ctx.lineTo(left + pw, yy); ctx.stroke();
  }

  // zero line
  ctx.strokeStyle = 'rgba(100,116,139,0.65)';
  ctx.setLineDash([5, 4]);
  ctx.beginPath(); ctx.moveTo(left, yAt(0)); ctx.lineTo(left + pw, yAt(0)); ctx.stroke();
  ctx.setLineDash([]);

  function strokeTrace(y, color, lw) {
    ctx.strokeStyle = color;
    ctx.lineWidth = lw;
    ctx.beginPath();
    y.forEach((v, i) => i ? ctx.lineTo(xAt(t[i]), yAt(v)) : ctx.moveTo(xAt(t[i]), yAt(v)));
    ctx.stroke();
  }
  strokeTrace(y0, '#d97706', 1.6);
  strokeTrace(y1, '#0f766e', 1.8);

  ctx.strokeStyle = '#94a3b8';
  ctx.lineWidth = 1;
  ctx.strokeRect(left, top, pw, ph);

  ctx.fillStyle = '#0f172a';
  ctx.font = '600 13px "IBM Plex Sans", sans-serif';
  ctx.textAlign = 'left';
  ctx.fillText('Baseline — ' + data.ch + ', época ' + data.epoch +
    ' (' + data.status + ') · offset ' + data.offset_before + ' → ' +
    data.offset_after + ' µV', left, 20);

  // legend
  ctx.font = '12px sans-serif';
  const lx = left + pw - 110, ly = top + 14;
  ctx.fillStyle = '#d97706'; ctx.fillRect(lx, ly - 4, 16, 3);
  ctx.fillStyle = '#334155'; ctx.fillText('antes', lx + 22, ly);
  ctx.fillStyle = '#0f766e'; ctx.fillRect(lx, ly + 14, 16, 3);
  ctx.fillStyle = '#334155'; ctx.fillText('después', lx + 22, ly + 18);

  ctx.fillStyle = '#64748b';
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'right';
  ctx.fillText(ymin.toFixed(0), left - 6, top + ph + 3);
  ctx.fillText(ymax.toFixed(0), left - 6, top + 4);
  ctx.fillText('µV', left - 6, top + ph / 2);
  ctx.textAlign = 'center';
  [0, 0.5, 1].forEach(v => {
    const tt = t[0] + v * (t[t.length - 1] - t[0]);
    ctx.fillText(tt.toFixed(1), xAt(tt), top + ph + 16);
  });
  ctx.fillText('Tiempo dentro de la época (s)', left + pw / 2, H - 10);

  document.getElementById('title').textContent =
    data.ch + ' · ep ' + data.epoch + ' · offset ' +
    data.offset_before + ' → ' + data.offset_after + ' µV · ' + data.status;
}

async function refresh() {
  setStatus('Cargando…');
  try {
    const ch = document.getElementById('ch').value;
    const ep = document.getElementById('epoch').value;
    const r = await fetch('/api/epoch?ch=' + encodeURIComponent(ch) + '&epoch=' + ep);
    const d = await r.json();
    if (!r.ok) throw new Error(d.error || 'error');
    data = d;
    draw();
    setStatus('OK · ' + d.status + ' · Δoffset = ' +
      (d.offset_before - d.offset_after).toFixed(2) + ' µV', 'ok');
  } catch (e) {
    setStatus(String(e.message || e), 'err');
  }
}

async function savePng() {
  setStatus('Generando PNG…');
  try {
    const body = JSON.stringify({
      ch: document.getElementById('ch').value,
      epoch: +document.getElementById('epoch').value,
    });
    const r = await fetch('/api/save', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body,
    });
    const d = await r.json();
    if (!r.ok || !d.ok) throw new Error(d.error || 'error');
    setStatus('Guardado: ' + d.file, 'ok');
  } catch (e) {
    setStatus(String(e.message || e), 'err');
  }
}

document.getElementById('ch').addEventListener('change', refresh);
document.getElementById('epoch').addEventListener('change', refresh);
window.addEventListener('resize', () => { if (data) draw(); });

function fmtAfter(v) {
  if (!isFinite(v)) return '—';
  if (Math.abs(v) < 1e-9) return v.toExponential(2);
  return v.toFixed(4);
}

function renderOffsetTable(d) {
  document.getElementById('tbl-caption').textContent =
    'Épocas válidas: ' + d.n_valid + ' · ventana [0, ' + d.baseline_end_s +
    '] s · media global antes ' + d.global_mean_before.toFixed(2) +
    ' µV → después ' + fmtAfter(d.global_mean_after) + ' µV';
  const body = document.getElementById('offset-body');
  body.innerHTML = '';
  d.rows.forEach(r => {
    const tr = document.createElement('tr');
    if (r.channel === DEFAULT_CH) tr.className = 'highlight';
    tr.innerHTML =
      '<td>' + r.channel + '</td>' +
      '<td>' + r.mean_before.toFixed(2) + '</td>' +
      '<td>' + r.max_before.toFixed(2) + '</td>' +
      '<td>' + fmtAfter(r.mean_after) + '</td>' +
      '<td>' + fmtAfter(r.max_after) + '</td>';
    body.appendChild(tr);
  });
  const sum = document.createElement('tr');
  sum.className = 'summary';
  sum.innerHTML =
    '<td>Los ' + d.rows.length + ' canales (media global)</td>' +
    '<td>' + d.global_mean_before.toFixed(2) + '</td>' +
    '<td>—</td>' +
    '<td>' + fmtAfter(d.global_mean_after) + '</td>' +
    '<td>—</td>';
  body.appendChild(sum);
}

async function loadOffsetTable() {
  try {
    const r = await fetch('/api/offsets');
    const d = await r.json();
    if (!r.ok) throw new Error(d.error || 'error');
    renderOffsetTable(d);
  } catch (e) {
    document.getElementById('tbl-caption').textContent = 'Error tabla: ' + (e.message || e);
  }
}

refresh();
loadOffsetTable();
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
    println("Cargando baseline desde $ROOT …")
    store = load_store()
    println("  $(length(store.channels)) ch · $(store.n_epochs) épocas ($(store.n_valid) válidas) · bl_end=$(store.baseline_end_s)s")
    println("  Default: $(store.default_channel) · época $(store.default_epoch) (máx |offset|)")

    out = save_baseline_png(store)
    println("  → $(basename(out))")
    csv_out = save_offset_table_csv(store)
    println("  → $(basename(csv_out))")
    top = first(store.offset_table)
    println("  Top |offset| medio antes: $(top.channel) = $(top.offset_mean_before_uv) µV")

    server, port = _listen_available(HOST, PORT)
    url = "http://$HOST:$port/"
    println()
    println("UI baseline → $url")
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
