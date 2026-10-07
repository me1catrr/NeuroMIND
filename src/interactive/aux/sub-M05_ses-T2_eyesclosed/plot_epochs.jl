#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Figuras de épocas (quality / overlay / stacked)
# ═══════════════════════════════════════════════════════════════
#
#  Regenera las tres figuras de segmentación del informe:
#    1) Histograma de quality (umbral 0.5)
#    2) Overlay de épocas válidas (canal Oz + media)
#    3) Época rechazada apilada (worst_channel en rojo)
#
#  Señal: reconstrucción post-ICA en memoria (sin CSV de 100 s):
#    cache/ica_result.jls + ica_labels_auto.csv
#      → A_keep * S_keep  (+ baseline first_window_mean 0.10 s)
#
#  Entrada:
#    ../../tables/segments_table.csv
#    ../../tables/rejected_segments.csv
#    ../../cache/ica_result.jls
#    ../../ica_labels_auto.csv
#
#  Salida (este directorio):
#    epoch_quality_histogram.png
#    epochs_overlay_Oz.png
#    epoch_012_rejected_stacked.png
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_epochs.jl
#    # → http://127.0.0.1:8772/  (Ctrl+C para salir)
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_epochs.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      24-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using NeuroMIND
using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics, Serialization

const HERE   = @__DIR__
const ROOT   = normpath(joinpath(HERE, "..", "..", "..", "..",
               "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))   # …/eyesclosed
const TABLES = joinpath(ROOT, "tables")
const HOST   = "127.0.0.1"
const PORT   = 8772   # after before/after = 8771

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK    = "eyesclosed"
const TASK_ID = "EC"

const QUALITY_THRESHOLD = 0.5
const BASELINE_END_S    = 0.10   # pipeline [baseline] first_window_mean
const DEFAULT_OVERLAY_N = 40
const DEFAULT_CHANNEL   = "Oz"

# ── Store ──────────────────────────────────────────────────────

mutable struct EpochStore
    segments::DataFrame
    rejected::DataFrame
    channels::Vector{String}
    clean::Matrix{Float64}          # channels × samples (post-ICA, continuous)
    fs::Float64
    duration_s::Float64
    rejected_ics::Vector{Int}
    subject::String
    session::String
    task::String
end

function _load_rejected_ics(root::String)::Vector{Int}
    auto_p = joinpath(root, "ica_labels_auto.csv")
    if isfile(auto_p)
        df = CSV.read(auto_p, DataFrame)
        hasproperty(df, :label) || return Int[]
        mask = String.(df.label) .== "artifact"
        return Int.(df.component[mask])
    end
    comp_p = joinpath(root, "tables", "ica", "ica_components.csv")
    if isfile(comp_p)
        df = CSV.read(comp_p, DataFrame)
        hasproperty(df, :rejected) || return Int[]
        return Int.(df.component[Bool.(df.rejected)])
    end
    return Int[]
end

"""
Reconstruye la señal ICA-limpia completa en memoria (paso [4/8] sin CSV).
No re-ejecuta FastICA: usa `cache/ica_result.jls` + etiquetas de rechazo.
"""
function _reconstruct_clean(root::String)
    cache_p = joinpath(root, "cache", "ica_result.jls")
    isfile(cache_p) || error("No encontrado: $cache_p — ejecuta el pipeline [4/8] primero")

    ica = deserialize(cache_p)
    rej = _load_rejected_ics(root)
    n_comp = size(ica.activations, 1)
    channels = String.(ica.meta.channel_names)
    fs = Float64(ica.meta.fs)

    clean = if isempty(rej)
        # Sin rechazo → reconstrucción completa A * S ≈ filtrada
        ica.mixing_matrix * ica.activations
    else
        keep = setdiff(1:n_comp, rej)
        ica.mixing_matrix[:, keep] * ica.activations[keep, :]
    end

    # Alinear filas con canales del meta
    if size(clean, 1) != length(channels)
        error("Dimensión ICA ($(size(clean, 1))) ≠ n canales ($(length(channels)))")
    end

    return clean, channels, fs, rej
end

"""Recorta [start_s, end_s) y aplica baseline first_window_mean."""
function _epoch_window(
    clean::Matrix{Float64},
    fs::Float64,
    start_s::Float64,
    end_s::Float64;
    baseline_end_s::Float64 = BASELINE_END_S,
)::Matrix{Float64}
    n_samp = size(clean, 2)
    i0 = clamp(Int(round(start_s * fs)) + 1, 1, n_samp)
    i1 = clamp(Int(round(end_s * fs)), i0, n_samp)
    ep = copy(clean[:, i0:i1])
    n_bl = clamp(Int(round(baseline_end_s * fs)), 1, size(ep, 2))
    ep .-= mean(ep[:, 1:n_bl]; dims = 2)
    return ep
end

function load_store()::EpochStore
    seg_p = joinpath(TABLES, "segments_table.csv")
    rej_p = joinpath(TABLES, "rejected_segments.csv")
    isfile(seg_p) || error("No encontrado: $seg_p")
    isfile(rej_p) || error("No encontrado: $rej_p")

    segments = CSV.read(seg_p, DataFrame)
    rejected = CSV.read(rej_p, DataFrame)
    println("  Reconstruyendo señal post-ICA desde caché…")
    clean, channels, fs, rej_ics = _reconstruct_clean(ROOT)
    dur = size(clean, 2) / fs

    return EpochStore(
        segments, rejected, channels, clean, fs, dur, rej_ics,
        SUBJECT, SESSION, TASK,
    )
end

# ── Plot helpers ───────────────────────────────────────────────

function _task_tag(store::EpochStore)
    return "$(store.subject)_$(store.session)_task-$(TASK_ID)"
end

function _valid_rows(store::EpochStore)
    return store.segments[String.(store.segments.status) .== "valid", :]
end

function _default_rejected_epoch(store::EpochStore)::Int
    isempty(store.rejected) && return 1
    return Int(store.rejected.epoch[1])
end

function _worst_channel(store::EpochStore, epoch::Int)::String
    for df in (store.rejected, store.segments)
        rows = df[Int.(df.epoch) .== epoch, :]
        isempty(rows) && continue
        hasproperty(rows, :worst_channel) || continue
        w = strip(String(rows.worst_channel[1]))
        isempty(w) || return w
    end
    return store.channels[1]
end

"""Histograma de quality ∈ [0,1] con umbral operativo."""
function plot_epoch_quality_histogram(
    store::EpochStore;
    threshold::Float64 = QUALITY_THRESHOLD,
    nbins::Int = 20,
)::Figure
    q = Float64.(store.segments.quality)
    n = length(q)
    nbins = clamp(nbins, 5, 80)
    fig = Figure(size = (720, 480), fontsize = 13)
    ax = Axis(fig[1, 1];
        title = "Distribución de calidad por época ($n épocas)",
        xlabel = "quality ∈ [0, 1]",
        ylabel = "Frecuencia",
    )
    ax.xgridvisible = true
    ax.ygridvisible = true
    hist!(ax, q; bins = nbins, color = RGBf(0.15, 0.55, 0.55), strokewidth = 0.5,
          strokecolor = :white)
    vlines!(ax, [threshold]; color = :red, linestyle = :dash, linewidth = 2)
    xlims!(ax, 0.0, 1.0)
    return fig
end

"""Épocas válidas superpuestas en un canal + media."""
function plot_epochs_overlay(
    store::EpochStore;
    channel::String = DEFAULT_CHANNEL,
    n_max::Int = DEFAULT_OVERLAY_N,
)::Figure
    ch = String(channel)
    ci = findfirst(==(ch), store.channels)
    ci === nothing && error("Canal no disponible: $ch")

    valid = _valid_rows(store)
    n_take = min(n_max, nrow(valid))
    n_take < 1 && error("No hay épocas válidas")
    rows = valid[1:n_take, :]

    fig = Figure(size = (900, 480), fontsize = 13)
    ax = Axis(fig[1, 1];
        title = "Épocas superpuestas — canal $ch ($n_take/$(nrow(valid)) épocas), $(_task_tag(store))",
        xlabel = "Tiempo dentro de la época (s)",
        ylabel = "Amplitud (µV)",
    )
    ax.xgridvisible = true
    ax.ygridvisible = true

    traces = Vector{Vector{Float64}}()
    t_rel = Float64[]
    for r in eachrow(rows)
        ep = _epoch_window(store.clean, store.fs, Float64(r.start_s), Float64(r.end_s))
        y = ep[ci, :]
        push!(traces, y)
        if isempty(t_rel)
            t_rel = collect(range(0.0; step = 1 / store.fs, length = length(y)))
        end
        lines!(ax, t_rel[1:length(y)], y;
               color = (RGBf(0.45, 0.65, 0.85), 0.35), linewidth = 0.8)
    end

    n_min = minimum(length.(traces))
    M = reduce(hcat, [tr[1:n_min] for tr in traces])
    μ = vec(mean(M; dims = 2))
    lines!(ax, t_rel[1:n_min], μ; color = :red, linewidth = 2.2, label = "media")
    axislegend(ax; position = :rt, framevisible = true, labelsize = 11)
    xlims!(ax, 0.0, maximum(t_rel))
    return fig
end

"""Una época, canales apilados; worst_channel en rojo."""
function plot_epochs_stacked(
    store::EpochStore;
    epoch::Int = _default_rejected_epoch(store),
    channels::Union{Nothing,Vector{String}} = nothing,
)::Figure
    rows = store.segments[Int.(store.segments.epoch) .== epoch, :]
    isempty(rows) && error("Época $epoch no está en segments_table")
    r = rows[1, :]
    start_s = Float64(r.start_s)
    end_s   = Float64(r.end_s)
    worst   = _worst_channel(store, epoch)

    ep = _epoch_window(store.clean, store.fs, start_s, end_s)
    ch_list = channels === nothing ? store.channels :
        String[c for c in channels if c in store.channels]
    isempty(ch_list) && error("Ningún canal seleccionado")
    # Mantener orden del montaje
    ch_list = [c for c in store.channels if c in Set(ch_list)]

    n_ch = length(ch_list)
    n_samp = size(ep, 2)
    t_rel = collect(range(0.0; step = 1 / store.fs, length = n_samp))

    scales = Float64[]
    for ch in ch_list
        ci = findfirst(==(ch), store.channels)
        push!(scales, std(ep[ci, :]))
    end
    pos_scales = filter(>(0), scales)
    step = max(3.0 * (isempty(pos_scales) ? 10.0 : mean(pos_scales)), 15.0)

    fig = Figure(size = (900, max(420, 18 * n_ch + 80)), fontsize = 12)
    ax = Axis(fig[1, 1];
        title = "Época $epoch — $(n_ch) canales (worst=$worst), $(_task_tag(store))",
        xlabel = "Tiempo dentro de la época (s)",
        ylabel = "Canal (offset)",
    )
    ax.xgridvisible = true
    ax.ygridvisible = false

    yticks_pos = Float64[]
    yticks_lab = String[]
    for (i, ch) in enumerate(ch_list)
        ci = findfirst(==(ch), store.channels)
        offset = (n_ch - i) * step
        y = ep[ci, :] .+ offset
        is_worst = ch == worst
        lines!(ax, t_rel, y;
               color = is_worst ? :red : RGBf(0.25, 0.45, 0.65),
               linewidth = is_worst ? 1.8 : 0.9)
        push!(yticks_pos, offset)
        push!(yticks_lab, ch)
    end
    ax.yticks = (yticks_pos, yticks_lab)
    xlims!(ax, 0.0, maximum(t_rel))
    return fig
end

function save_view_png(
    store::EpochStore,
    view::String;
    channel::String = DEFAULT_CHANNEL,
    n_max::Int = DEFAULT_OVERLAY_N,
    epoch::Int = _default_rejected_epoch(store),
    nbins::Int = 20,
    channels::Union{Nothing,Vector{String}} = nothing,
)::String
    view = lowercase(strip(view))
    if view in ("hist", "histogram", "quality")
        fig = plot_epoch_quality_histogram(store; nbins = nbins)
        out = joinpath(HERE, "epoch_quality_histogram.png")
    elseif view in ("overlay", "epochs_overlay")
        fig = plot_epochs_overlay(store; channel = channel, n_max = n_max)
        out = joinpath(HERE, "epochs_overlay_$(channel).png")
    elseif view in ("stacked", "epoch_stacked")
        fig = plot_epochs_stacked(store; epoch = epoch, channels = channels)
        out = joinpath(HERE, "epoch_$(lpad(string(epoch), 3, '0'))_rejected_stacked.png")
    else
        error("Vista desconocida: $view (hist|overlay|stacked)")
    end
    save(out, fig; px_per_unit = 2)
    return out
end

function save_all_pngs(store::EpochStore)::Vector{String}
    outs = String[]
    push!(outs, save_view_png(store, "hist"))
    push!(outs, save_view_png(store, "overlay"; channel = DEFAULT_CHANNEL, n_max = DEFAULT_OVERLAY_N))
    push!(outs, save_view_png(store, "stacked"; epoch = _default_rejected_epoch(store)))
    return outs
end

# ── HTTP helpers ───────────────────────────────────────────────

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

function _send_bytes(sock, status::Int, body::Vector{UInt8}; content_type::String)
    write(sock, "HTTP/1.1 $status\r\nContent-Type: $content_type\r\n" *
                "Content-Length: $(length(body))\r\nConnection: close\r\n" *
                "Access-Control-Allow-Origin: *\r\n\r\n")
    write(sock, body)
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
    _int = (key, default) -> Int(round(_num(key, default)))
    chs = String[]
    mchs = match(r"\"channels\"\s*:\s*\[([^\]]*)\]", s)
    if mchs !== nothing
        for m in eachmatch(r"\"([^\"]+)\"" , mchs.captures[1])
            push!(chs, String(m.captures[1]))
        end
    end
    return (
        view = _str("view", "hist"),
        ch = _str("ch", DEFAULT_CHANNEL),
        n_max = _int("n_max", DEFAULT_OVERLAY_N),
        epoch = _int("epoch", 12),
        nbins = _int("nbins", 20),
        channels = chs,
        all = occursin(r"\"all\"\s*:\s*true", s),
    )
end

function handle_request(sock, store::EpochStore)
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
            n_valid = count(String.(store.segments.status) .== "valid")
            n_rej = count(String.(store.segments.status) .== "rejected")
            def_ep = _default_rejected_epoch(store)
            body = "{" *
                "\"n_epochs\":$(nrow(store.segments))," *
                "\"n_valid\":$n_valid,\"n_rejected\":$n_rej," *
                "\"duration_s\":$(round(store.duration_s; digits=2))," *
                "\"fs\":$(store.fs)," *
                "\"rejected_ics\":[" * join(store.rejected_ics, ",") * "]," *
                "\"default_epoch\":$def_ep," *
                "\"worst_channel\":\"$(_worst_channel(store, def_ep))\"," *
                "\"threshold\":$QUALITY_THRESHOLD," *
                "\"channels\":[" * join(["\"$c\"" for c in store.channels], ",") * "]" *
                "}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/quality"
            q = round.(Float64.(store.segments.quality); digits = 4)
            st = String.(store.segments.status)
            ep = Int.(store.segments.epoch)
            body = "{\"quality\":[" * join(q, ",") * "]," *
                   "\"status\":[" * join(["\"$s\"" for s in st], ",") * "]," *
                   "\"epoch\":[" * join(ep, ",") * "]," *
                   "\"threshold\":$QUALITY_THRESHOLD}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/overlay"
            qs = _parse_qs(path)
            ch = get(qs, "ch", DEFAULT_CHANNEL)
            n_max = something(tryparse(Int, get(qs, "n_max", string(DEFAULT_OVERLAY_N))), DEFAULT_OVERLAY_N)
            max_pts = something(tryparse(Int, get(qs, "max_points", "500")), 500)
            ci = findfirst(==(ch), store.channels)
            ci === nothing && error("Canal no disponible: $ch")
            valid = _valid_rows(store)
            n_take = min(max(n_max, 1), nrow(valid))
            n_take < 1 && error("No hay épocas válidas")

            # Longitud común = mínima entre épocas tomadas
            ys = Vector{Vector{Float64}}()
            for r in eachrow(valid[1:n_take, :])
                ep = _epoch_window(store.clean, store.fs, Float64(r.start_s), Float64(r.end_s))
                push!(ys, ep[ci, :])
            end
            n_min = minimum(length.(ys))
            step = n_min > max_pts ? ceil(Int, n_min / max_pts) : 1
            idx = collect(1:step:n_min)
            t_rel = collect(range(0.0; step = 1 / store.fs, length = n_min))[idx]
            traces = ["[" * join(round.(y[idx]; digits = 3), ",") * "]" for y in ys]
            body = "{\"ch\":\"$ch\",\"n\":$n_take,\"n_valid\":$(nrow(valid))," *
                   "\"t\":[" * join(round.(t_rel; digits = 4), ",") * "]," *
                   "\"traces\":[" * join(traces, ",") * "]}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/stacked"
            qs = _parse_qs(path)
            epoch = something(tryparse(Int, get(qs, "epoch", string(_default_rejected_epoch(store)))),
                              _default_rejected_epoch(store))
            max_pts = something(tryparse(Int, get(qs, "max_points", "500")), 500)
            rows = store.segments[Int.(store.segments.epoch) .== epoch, :]
            isempty(rows) && error("Época $epoch no encontrada")
            r = rows[1, :]
            worst = _worst_channel(store, epoch)
            ep = _epoch_window(store.clean, store.fs, Float64(r.start_s), Float64(r.end_s))
            n_samp = size(ep, 2)
            step = n_samp > max_pts ? ceil(Int, n_samp / max_pts) : 1
            idx = collect(1:step:n_samp)
            t_rel = collect(range(0.0; step = 1 / store.fs, length = n_samp))[idx]
            ch_json = String[]
            for (i, ch) in enumerate(store.channels)
                y = round.(ep[i, idx]; digits = 3)
                push!(ch_json, "{\"name\":\"$ch\",\"y\":[" * join(y, ",") * "]}")
            end
            body = "{\"epoch\":$epoch,\"start_s\":$(r.start_s),\"end_s\":$(r.end_s)," *
                   "\"worst\":\"$worst\",\"status\":\"$(r.status)\"," *
                   "\"quality\":$(r.quality)," *
                   "\"t\":[" * join(round.(t_rel; digits = 4), ",") * "]," *
                   "\"channels\":[" * join(ch_json, ",") * "]}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_save_json(raw)
            if req.all
                outs = save_all_pngs(store)
                files = join(["\"$(basename(o))\"" for o in outs], ",")
                body = "{\"ok\":true,\"files\":[$files]}"
                println("[$(Dates.format(now(), "HH:MM:SS"))] PNG ×$(length(outs)) → $HERE")
            else
                chs = isempty(req.channels) ? nothing : req.channels
                out = save_view_png(store, req.view;
                    channel = req.ch, n_max = req.n_max, epoch = req.epoch,
                    nbins = req.nbins, channels = chs)
                body = "{\"ok\":true,\"file\":$(repr(basename(out)))}"
                println("[$(Dates.format(now(), "HH:MM:SS"))] PNG → $out")
            end
            _send(sock, 200, body; content_type = "application/json")
        else
            _send(sock, 404, "{\"error\":\"not found\"}"; content_type = "application/json")
        end
    catch e
        msg = sprint(showerror, e)
        @warn "Request error" exception = e
        _send(sock, 400, "{\"ok\":false,\"error\":$(repr(msg))}"; content_type = "application/json")
    end
end

# ── HTML ───────────────────────────────────────────────────────

function html_page(store::EpochStore)::String
    ch_opts = join(
        ["<option value=\"$c\"$(c == DEFAULT_CHANNEL ? " selected" : "")>$c</option>"
         for c in store.channels], "\n")
    ch_checks = join(
        ["<label class=\"ch\"><input type=\"checkbox\" class=\"chbox\" value=\"$c\" checked>$c</label>"
         for c in store.channels], "\n")
    def_ep = _default_rejected_epoch(store)
    n_valid = count(String.(store.segments.status) .== "valid")
    n_rej = count(String.(store.segments.status) .== "rejected")
    thresh = QUALITY_THRESHOLD

    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — Épocas (interactivo)</title>
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
  .wrap { max-width:1200px; margin:16px auto; padding:0 16px; }
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
  .btn-ghost { background:#fff; color:var(--muted); border:1px solid var(--border); }
  .btn:disabled { opacity:.5; cursor:wait; }
  canvas { width:100%; display:block; background:#fff; border:1px solid var(--border); border-radius:8px; }
  .status { font-size:12px; color:var(--muted); margin-top:8px; min-height:1.2em; }
  .status.ok { color:#047857; } .status.err { color:#b91c1c; }
  .hint { font-size:12px; color:var(--muted); margin:8px 0 0; }
  .tabs { display:flex; gap:6px; margin-bottom:10px; }
  .tab { padding:6px 12px; border:1px solid var(--border); border-radius:6px; background:#fff; cursor:pointer; font-size:12px; font-weight:600; }
  .tab.active { background:#134e4a; color:#fff; border-color:#134e4a; }
  .channels {
    display:grid; grid-template-columns:repeat(auto-fill,minmax(72px,1fr));
    gap:4px 6px; max-height:120px; overflow:auto; margin-top:8px;
  }
  .ch { display:flex; align-items:center; gap:4px; font-size:12px; padding:2px 4px; border-radius:4px; cursor:pointer; }
  .ch:hover { background:#f1f5f9; }
  .ch input { accent-color:var(--accent); }
  .quick { display:flex; flex-wrap:wrap; gap:6px; margin-top:6px; }
  .chip { border:1px solid var(--border); background:#fff; border-radius:999px; padding:3px 10px; font-size:12px; cursor:pointer; color:var(--muted); }
  .chip:hover { border-color:var(--accent); color:var(--accent); }
  .panel { display:none; }
  .panel.active { display:block; }
</style>
</head>
<body>
<header>
  <h1>Épocas · interactivo</h1>
  <span>$SUBJECT / $SESSION / $TASK</span>
  <span>$(nrow(store.segments)) épocas · $n_valid válidas · $n_rej rechazada(s) · ICs rech. [$(join(store.rejected_ics, ","))]</span>
</header>
<div class="wrap">
  <div class="card">
    <div class="tabs">
      <button class="tab active" data-view="hist" onclick="setView('hist')">Histograma quality</button>
      <button class="tab" data-view="overlay" onclick="setView('overlay')">Overlay</button>
      <button class="tab" data-view="stacked" onclick="setView('stacked')">Época apilada</button>
    </div>

    <div id="panel-hist" class="panel active">
      <div class="toolbar">
        <div class="field"><label>Bins</label>
          <input type="number" id="nbins" value="20" min="5" max="80" step="1"/></div>
        <div class="field"><label>Filtro</label>
          <select id="qfilter">
            <option value="all">Todas</option>
            <option value="valid">Solo válidas</option>
            <option value="rejected">Solo rechazadas</option>
          </select></div>
        <div class="field"><label>Umbral</label>
          <input type="number" id="thresh" value="$thresh" min="0" max="1" step="0.05"/></div>
        <button class="btn btn-primary" onclick="refresh()">Actualizar</button>
        <button class="btn btn-save" onclick="savePng(false)">Guardar PNG</button>
      </div>
    </div>

    <div id="panel-overlay" class="panel">
      <div class="toolbar">
        <div class="field"><label>Canal</label>
          <select id="ch">$ch_opts</select></div>
        <div class="field"><label>N épocas</label>
          <input type="number" id="n_max" value="$DEFAULT_OVERLAY_N" min="1" max="$n_valid" step="1"/></div>
        <button class="btn btn-primary" onclick="refresh()">Actualizar</button>
        <button class="btn btn-save" onclick="savePng(false)">Guardar PNG</button>
      </div>
    </div>

    <div id="panel-stacked" class="panel">
      <div class="toolbar">
        <div class="field"><label>Época</label>
          <input type="number" id="epoch" value="$def_ep" min="1" max="$(nrow(store.segments))" step="1"/></div>
        <button class="btn btn-primary" onclick="refresh()">Actualizar</button>
        <button class="btn btn-save" onclick="savePng(false)">Guardar PNG</button>
        <button class="btn btn-save" onclick="savePng(true)">Guardar las 3</button>
      </div>
      <div class="quick">
        <button class="chip" onclick="presetCh('all')">Todos</button>
        <button class="chip" onclick="presetCh('none')">Ninguno</button>
        <button class="chip" onclick="presetCh('worst')">Solo worst</button>
        <button class="chip" onclick="presetCh('occ')">Occipitales</button>
        <button class="chip" onclick="presetCh('central')">Centrales</button>
      </div>
      <div class="channels" id="chgrid">$ch_checks</div>
    </div>

    <p class="hint">Canvas interactivo (datos vía API). «Guardar PNG» exporta con CairoMakie según los controles actuales.</p>
  </div>

  <div class="card">
    <strong id="title" style="font-size:13px">Histograma</strong>
    <canvas id="cv"></canvas>
    <div class="status" id="status"></div>
  </div>
</div>
<script>
let view = 'hist';
let qualityData = null;
let overlayData = null;
let stackedData = null;
const DEFAULT_CH = '$DEFAULT_CHANNEL';

function setStatus(msg, cls) {
  const el = document.getElementById('status');
  el.textContent = msg || '';
  el.className = 'status' + (cls ? ' ' + cls : '');
}

function setView(v) {
  view = v;
  document.querySelectorAll('.tab').forEach(t => t.classList.toggle('active', t.dataset.view === v));
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  document.getElementById('panel-' + v).classList.add('active');
  refresh();
}

function selectedStackedChannels() {
  return Array.from(document.querySelectorAll('.chbox:checked')).map(el => el.value);
}

function presetCh(mode) {
  const boxes = Array.from(document.querySelectorAll('.chbox'));
  const worst = (stackedData && stackedData.worst) || 'C4';
  const occ = new Set(['O1','Oz','O2','P7','P3','Pz','P4','P8']);
  const cen = new Set(['C3','Cz','C4','FC1','FC2','CP1','CP2']);
  boxes.forEach(b => {
    if (mode === 'all') b.checked = true;
    else if (mode === 'none') b.checked = false;
    else if (mode === 'worst') b.checked = (b.value === worst);
    else if (mode === 'occ') b.checked = occ.has(b.value);
    else if (mode === 'central') b.checked = cen.has(b.value);
  });
  if (view === 'stacked') drawStacked();
}

function setupHiDPICanvas(canvas, cssH) {
  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || canvas.parentElement.clientWidth || 1100;
  canvas.style.height = cssH + 'px';
  canvas.width  = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, W: cssW, H: cssH };
}

function drawAxes(ctx, pad, W, H, xlabel, ylabel) {
  const left = pad.l, top = pad.t, pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
  ctx.strokeStyle = '#94a3b8';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(left, top);
  ctx.lineTo(left, top + ph);
  ctx.lineTo(left + pw, top + ph);
  ctx.stroke();
  ctx.fillStyle = '#64748b';
  ctx.font = '12px "IBM Plex Sans", sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText(xlabel, left + pw / 2, H - 10);
  ctx.save();
  ctx.translate(14, top + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(ylabel, 0, 0);
  ctx.restore();
  return { left, top, pw, ph };
}

function drawHist() {
  if (!qualityData) return;
  const filter = document.getElementById('qfilter').value;
  const nbins = Math.max(5, Math.min(80, +document.getElementById('nbins').value || 20));
  const thresh = +document.getElementById('thresh').value;
  let q = qualityData.quality.slice();
  let st = qualityData.status.slice();
  if (filter !== 'all') {
    const qq = [], ss = [];
    for (let i = 0; i < q.length; i++) {
      if (st[i] === filter) { qq.push(q[i]); ss.push(st[i]); }
    }
    q = qq; st = ss;
  }
  const { ctx, W, H } = setupHiDPICanvas(document.getElementById('cv'), 440);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);
  const pad = { l: 58, r: 24, t: 36, b: 48 };
  const { left, top, pw, ph } = drawAxes(ctx, pad, W, H, 'quality ∈ [0, 1]', 'Frecuencia');

  ctx.fillStyle = '#0f172a';
  ctx.font = '600 14px "IBM Plex Sans", sans-serif';
  ctx.textAlign = 'left';
  ctx.fillText('Distribución de calidad por época (' + q.length + ' épocas · ' + nbins + ' bins)', pad.l, 22);

  const lo = 0, hi = 1, bw = (hi - lo) / nbins;
  const counts = new Array(nbins).fill(0);
  for (const v of q) {
    const k = Math.min(nbins - 1, Math.max(0, Math.floor((v - lo) / bw)));
    counts[k]++;
  }
  const maxC = Math.max(1, ...counts);
  // grid
  ctx.strokeStyle = '#eef2f7';
  for (let g = 0; g <= 4; g++) {
    const y = top + (g / 4) * ph;
    ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(left + pw, y); ctx.stroke();
  }
  const barW = pw / nbins;
  ctx.fillStyle = '#268686';
  for (let i = 0; i < nbins; i++) {
    const h = (counts[i] / maxC) * ph;
    ctx.fillRect(left + i * barW + 1, top + ph - h, Math.max(1, barW - 2), h);
  }
  // threshold
  const tx = left + ((thresh - lo) / (hi - lo)) * pw;
  ctx.strokeStyle = '#dc2626';
  ctx.setLineDash([6, 4]);
  ctx.beginPath(); ctx.moveTo(tx, top); ctx.lineTo(tx, top + ph); ctx.stroke();
  ctx.setLineDash([]);
  // ticks
  ctx.fillStyle = '#64748b';
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'center';
  [0, 0.5, 1].forEach(v => {
    ctx.fillText(String(v), left + v * pw, top + ph + 16);
  });
  ctx.textAlign = 'right';
  ctx.fillText('0', left - 6, top + ph + 3);
  ctx.fillText(String(maxC), left - 6, top + 4);
  document.getElementById('title').textContent =
    'Histograma quality · ' + q.length + ' épocas · bins=' + nbins;
}

function drawOverlay() {
  if (!overlayData) return;
  const t = overlayData.t;
  const traces = overlayData.traces;
  const { ctx, W, H } = setupHiDPICanvas(document.getElementById('cv'), 460);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);
  const pad = { l: 58, r: 24, t: 36, b: 48 };
  const { left, top, pw, ph } = drawAxes(ctx, pad, W, H, 'Tiempo dentro de la época (s)', 'Amplitud (µV)');

  ctx.fillStyle = '#0f172a';
  ctx.font = '600 14px "IBM Plex Sans", sans-serif';
  ctx.textAlign = 'left';
  ctx.fillText('Épocas superpuestas — ' + overlayData.ch +
    ' (' + overlayData.n + '/' + overlayData.n_valid + ')', pad.l, 22);

  let ymin = Infinity, ymax = -Infinity;
  for (const tr of traces) for (const v of tr) {
    if (v < ymin) ymin = v;
    if (v > ymax) ymax = v;
  }
  if (!isFinite(ymin)) { ymin = -1; ymax = 1; }
  const padY = 0.05 * (ymax - ymin || 1);
  ymin -= padY; ymax += padY;
  const xAt = (ti) => left + ((ti - t[0]) / (t[t.length - 1] - t[0] || 1)) * pw;
  const yAt = (v) => top + ph - ((v - ymin) / (ymax - ymin || 1)) * ph;

  ctx.strokeStyle = '#eef2f7';
  for (let g = 0; g <= 4; g++) {
    const y = top + (g / 4) * ph;
    ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(left + pw, y); ctx.stroke();
  }

  ctx.strokeStyle = 'rgba(70,130,180,0.35)';
  ctx.lineWidth = 1;
  for (const tr of traces) {
    ctx.beginPath();
    for (let i = 0; i < tr.length; i++) {
      const xx = xAt(t[i]);
      const yy = yAt(tr[i]);
      i ? ctx.lineTo(xx, yy) : ctx.moveTo(xx, yy);
    }
    ctx.stroke();
  }
  // mean
  const mean = t.map((_, i) => {
    let s = 0, n = 0;
    for (const tr of traces) { if (i < tr.length) { s += tr[i]; n++; } }
    return n ? s / n : 0;
  });
  ctx.strokeStyle = '#dc2626';
  ctx.lineWidth = 2.2;
  ctx.beginPath();
  mean.forEach((v, i) => i ? ctx.lineTo(xAt(t[i]), yAt(v)) : ctx.moveTo(xAt(t[i]), yAt(v)));
  ctx.stroke();

  // legend
  ctx.fillStyle = '#dc2626';
  ctx.fillRect(left + pw - 70, top + 8, 14, 3);
  ctx.fillStyle = '#334155';
  ctx.font = '12px sans-serif';
  ctx.textAlign = 'left';
  ctx.fillText('media', left + pw - 52, top + 14);

  ctx.fillStyle = '#64748b';
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'center';
  [0, 0.5, 1].forEach(v => {
    const tt = t[0] + v * (t[t.length - 1] - t[0]);
    ctx.fillText(tt.toFixed(1), xAt(tt), top + ph + 16);
  });
  ctx.textAlign = 'right';
  ctx.fillText(ymin.toFixed(0), left - 6, top + ph + 3);
  ctx.fillText(ymax.toFixed(0), left - 6, top + 4);

  document.getElementById('title').textContent =
    'Overlay · ' + overlayData.ch + ' · ' + overlayData.n + ' épocas';
}

function drawStacked() {
  if (!stackedData) return;
  const sel = new Set(selectedStackedChannels());
  const chs = stackedData.channels.filter(c => sel.has(c.name));
  if (!chs.length) {
    setStatus('Selecciona al menos un canal.', 'err');
    return;
  }
  const t = stackedData.t;
  const worst = stackedData.worst;
  const n = chs.length;
  const cssH = Math.max(420, Math.min(900, 22 * n + 80));
  const { ctx, W, H } = setupHiDPICanvas(document.getElementById('cv'), cssH);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);
  const pad = { l: 52, r: 20, t: 36, b: 40 };
  const left = pad.l, top = pad.t, pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;

  ctx.fillStyle = '#0f172a';
  ctx.font = '600 14px "IBM Plex Sans", sans-serif';
  ctx.textAlign = 'left';
  ctx.fillText('Época ' + stackedData.epoch + ' (' + stackedData.start_s + '-' +
    stackedData.end_s + 's) · ' + stackedData.status + ' · worst=' + worst, pad.l, 22);

  let amp = 0;
  for (const c of chs) for (const v of c.y) amp = Math.max(amp, Math.abs(v));
  amp = Math.max(amp, 1);
  // Espaciado que cabe en el canvas
  const step = ph / Math.max(n, 1);
  const wave = step * 0.38;
  const xAt = (ti) => left + ((ti - t[0]) / (t[t.length - 1] - t[0] || 1)) * pw;

  ctx.strokeStyle = '#eef2f7';
  for (let g = 0; g <= 4; g++) {
    const x = left + (g / 4) * pw;
    ctx.beginPath(); ctx.moveTo(x, top); ctx.lineTo(x, top + ph); ctx.stroke();
  }

  ctx.font = '10px sans-serif';
  chs.forEach((c, i) => {
    const offset = top + (i + 0.5) * step;
    const isW = c.name === worst;
    ctx.strokeStyle = isW ? '#dc2626' : '#3b6ea5';
    ctx.lineWidth = isW ? 1.8 : 0.9;
    ctx.beginPath();
    c.y.forEach((v, j) => {
      const yy = offset - (v / amp) * wave;
      j ? ctx.lineTo(xAt(t[j]), yy) : ctx.moveTo(xAt(t[j]), yy);
    });
    ctx.stroke();
    ctx.fillStyle = isW ? '#dc2626' : '#64748b';
    ctx.textAlign = 'right';
    ctx.fillText(c.name, left - 6, offset + 3);
  });

  ctx.fillStyle = '#64748b';
  ctx.font = '12px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('Tiempo dentro de la época (s)', left + pw / 2, H - 10);
  [0, 0.5, 1].forEach(v => {
    const tt = t[0] + v * (t[t.length - 1] - t[0]);
    ctx.fillText(tt.toFixed(1), xAt(tt), top + ph + 14);
  });

  document.getElementById('title').textContent =
    'Época ' + stackedData.epoch + ' · ' + n + ' canales · worst=' + worst;
}

async function refresh() {
  setStatus('Cargando…');
  try {
    if (view === 'hist') {
      if (!qualityData) {
        const r = await fetch('/api/quality');
        const data = await r.json();
        if (!r.ok) throw new Error(data.error || 'error');
        qualityData = data;
      }
      drawHist();
      setStatus('OK · bins interactivos', 'ok');
    } else if (view === 'overlay') {
      const ch = document.getElementById('ch').value;
      const n = document.getElementById('n_max').value;
      const r = await fetch('/api/overlay?ch=' + encodeURIComponent(ch) + '&n_max=' + n);
      const data = await r.json();
      if (!r.ok) throw new Error(data.error || 'error');
      overlayData = data;
      drawOverlay();
      setStatus('OK · ' + data.n + ' épocas · ' + data.ch, 'ok');
    } else {
      const ep = document.getElementById('epoch').value;
      const r = await fetch('/api/stacked?epoch=' + ep);
      const data = await r.json();
      if (!r.ok) throw new Error(data.error || 'error');
      stackedData = data;
      drawStacked();
      setStatus('OK · época ' + data.epoch + ' · worst=' + data.worst, 'ok');
    }
  } catch (e) {
    setStatus(String(e.message || e), 'err');
  }
}

async function savePng(all) {
  setStatus(all ? 'Generando 3 PNG…' : 'Generando PNG…');
  try {
    const body = all
      ? JSON.stringify({ all: true })
      : JSON.stringify({
          view,
          ch: document.getElementById('ch').value,
          n_max: +document.getElementById('n_max').value,
          epoch: +document.getElementById('epoch').value,
          nbins: +document.getElementById('nbins').value,
          channels: selectedStackedChannels(),
        });
    const r = await fetch('/api/save', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body,
    });
    const data = await r.json();
    if (!r.ok || !data.ok) throw new Error(data.error || 'error');
    if (data.files) setStatus('Guardados: ' + data.files.join(', '), 'ok');
    else setStatus('Guardado: ' + data.file, 'ok');
  } catch (e) {
    setStatus(String(e.message || e), 'err');
  }
}

document.getElementById('nbins').addEventListener('change', () => { if (view === 'hist') drawHist(); });
document.getElementById('qfilter').addEventListener('change', () => { if (view === 'hist') drawHist(); });
document.getElementById('thresh').addEventListener('change', () => { if (view === 'hist') drawHist(); });
document.getElementById('ch').addEventListener('change', () => { if (view === 'overlay') refresh(); });
document.getElementById('n_max').addEventListener('change', () => { if (view === 'overlay') refresh(); });
document.getElementById('epoch').addEventListener('change', () => { if (view === 'stacked') refresh(); });
document.getElementById('chgrid').addEventListener('change', () => { if (view === 'stacked') drawStacked(); });
window.addEventListener('resize', () => {
  if (view === 'hist') drawHist();
  else if (view === 'overlay') drawOverlay();
  else drawStacked();
});

setView('hist');
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

"""Intenta PORT; si está ocupado, prueba PORT+1 … PORT+9."""
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
    error("No hay puerto libre en $(port)-$(port + n_try - 1). Cierra la instancia previa (p. ej. pkill -f plot_epochs.jl). Último error: $last_err")
end

function main()
    println("Cargando épocas desde $ROOT …")
    store = load_store()
    n_valid = count(String.(store.segments.status) .== "valid")
    n_rej = count(String.(store.segments.status) .== "rejected")
    println("  $(nrow(store.segments)) épocas · $n_valid válidas · $n_rej rechazada(s)")
    println("  Señal limpia: $(length(store.channels)) ch × $(size(store.clean, 2)) samp ($(round(store.duration_s; digits=1)) s)")
    println("  ICs rechazados: $(isempty(store.rejected_ics) ? "ninguno" : join(store.rejected_ics, ", "))")

    # Generar las 3 figuras canónicas al arrancar
    println("  Generando PNG canónicos…")
    for out in save_all_pngs(store)
        println("    → $(basename(out))")
    end

    server, port = _listen_available(HOST, PORT)
    url = "http://$HOST:$port/"
    println()
    println("UI épocas → $url")
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
