#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Histograma de amplitud interactivo (GUI navegador)
# ═══════════════════════════════════════════════════════════════
#
#  Rejilla de histogramas por canal (estética tipo informe QC):
#  barras coloreadas, media (línea roja) y μ ± σ (naranja discontinua).
#  El layout se escala con el nº de canales seleccionados.
#  Modo «combinado»: densidad de todos los canales + μ / ±2σ / ±3σ.
#
#  Entrada:  results/subjects/sub-M05/ses-T2/eyesclosed/tables/raw_signal.csv
#  Salida:   hist_<canales>[_t0-t1].png  |  hist_combined[_t0-t1].png (en este directorio)
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_raw_histogram.jl
#    # → http://127.0.0.1:8768/  (Ctrl+C para salir)
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_raw_histogram.jl
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
const PORT         = 8768   # distinto de plot_raw (8765), butterfly (8766), PSD (8767)

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK_ID = "EC"
const FS_HZ   = 500.0

# Paleta cualitativa (misma familia que butterfly / informe)
const BAR_COLORS = [
    RGBf(0.12, 0.47, 0.71), RGBf(0.20, 0.63, 0.17), RGBf(1.00, 0.50, 0.05),
    RGBf(0.89, 0.10, 0.11), RGBf(0.58, 0.40, 0.74), RGBf(0.55, 0.34, 0.29),
    RGBf(0.89, 0.47, 0.76), RGBf(0.50, 0.50, 0.50), RGBf(0.74, 0.74, 0.13),
    RGBf(0.09, 0.75, 0.81), RGBf(0.80, 0.20, 0.40), RGBf(0.15, 0.55, 0.35),
    RGBf(0.45, 0.25, 0.70), RGBf(0.90, 0.60, 0.10), RGBf(0.20, 0.40, 0.70),
    RGBf(0.55, 0.15, 0.15), RGBf(0.10, 0.60, 0.50), RGBf(0.70, 0.40, 0.10),
    RGBf(0.35, 0.35, 0.70), RGBf(0.60, 0.60, 0.20), RGBf(0.80, 0.30, 0.55),
    RGBf(0.25, 0.55, 0.75), RGBf(0.75, 0.25, 0.25), RGBf(0.40, 0.65, 0.30),
    RGBf(0.55, 0.30, 0.55), RGBf(0.85, 0.45, 0.15), RGBf(0.20, 0.50, 0.55),
    RGBf(0.65, 0.20, 0.45), RGBf(0.30, 0.45, 0.20), RGBf(0.50, 0.50, 0.75),
    RGBf(0.90, 0.35, 0.35),
]

const HEX_COLORS = [
    "#1f77b4","#2ca02c","#ff7f0e","#d62728","#9467bd","#8c564b",
    "#e377c2","#7f7f7f","#bcbd22","#17becf","#cc3366","#279158",
    "#7340b3","#e6991a","#3366b3","#8c2626","#1a9980","#b3661a",
    "#5959b3","#999933","#cc4d8c","#408cbc","#bf4040","#66a64d",
    "#8c4d8c","#d97326","#33808c","#a63373","#4d7333","#8080bf","#e65959",
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

# ── Histograma / stats ─────────────────────────────────────────

"""Bins equiespaciados + conteos (histograma clásico)."""
function _histogram(y::AbstractVector{<:Real}, nbins::Int)
    nbins = clamp(nbins, 5, 200)
    ymin, ymax = extrema(y)
    if !isfinite(ymin) || !isfinite(ymax) || ymax ≈ ymin
        c = isfinite(ymin) ? Float64(ymin) : 0.0
        edges = collect(range(c - 1.0, c + 1.0; length = nbins + 1))
        return edges, zeros(Int, nbins)
    end
    # margen 2 % para no cortar extremos en el borde del bin
    pad = 0.02 * (ymax - ymin)
    lo, hi = ymin - pad, ymax + pad
    edges = collect(range(lo, hi; length = nbins + 1))
    counts = zeros(Int, nbins)
    inv_w = nbins / (hi - lo)
    @inbounds for v in y
        isfinite(v) || continue
        k = clamp(floor(Int, (v - lo) * inv_w) + 1, 1, nbins)
        counts[k] += 1
    end
    return edges, counts
end

"""KDE gaussiana (Silverman); submuestrea si n es grande."""
function _kde(y::AbstractVector{<:Real}, xs::AbstractVector{<:Real};
              max_n::Int = 8000)
    vals = Float64[v for v in y if isfinite(v)]
    isempty(vals) && return zeros(length(xs))
    if length(vals) > max_n
        step = ceil(Int, length(vals) / max_n)
        vals = vals[1:step:end]
    end
    n = length(vals)
    σ = std(vals)
    σ <= 0 && (σ = 1.0)
    h = 1.06 * σ * n^(-0.2)
    dens = zeros(length(xs))
    inv = 1.0 / (h * sqrt(2π) * n)
    @inbounds for (i, xi) in enumerate(xs)
        s = 0.0
        for v in vals
            z = (v - xi) / h
            s += exp(-0.5 * z * z)
        end
        dens[i] = s * inv
    end
    return dens
end

function _channel_stats(y::AbstractVector{<:Real})
    vals = Float64[v for v in y if isfinite(v)]
    n = length(vals)
    n == 0 && return (n = 0, mean = NaN, std = NaN, var = NaN, min = NaN, max = NaN)
    μ = mean(vals)
    σ = std(vals; corrected = true)
    return (n = n, mean = μ, std = σ, var = σ^2, min = minimum(vals), max = maximum(vals))
end

"""Dimensiones de rejilla según nº de canales (máx. 4 columnas)."""
function _grid_dims(n::Int)::Tuple{Int,Int}
    n <= 0 && return (1, 1)
    n == 1 && return (1, 1)
    n == 2 && return (1, 2)
    n <= 4 && return (2, 2)
    n <= 6 && return (2, 3)
    n <= 9 && return (3, 3)
    n <= 12 && return (3, 4)
    n <= 16 && return (4, 4)
    ncols = min(4, ceil(Int, sqrt(n)))
    nrows = ceil(Int, n / ncols)
    return (nrows, ncols)
end

function _slice(store::RawStore, ch::AbstractString, t0::Float64, t1::Float64)
    mask = (store.t .>= t0) .& (store.t .<= t1)
    return store.data[String(ch)][mask]
end

# ── PNG (CairoMakie) ───────────────────────────────────────────

function save_hist_grid_png(
    store::RawStore,
    selected::Vector{String},
    t0::Float64,
    t1::Float64,
    nbins::Int,
    x_shared::Bool,
)::String
    isempty(selected) && error("Selecciona al menos un canal")
    t0 >= t1 && error("t0 debe ser < t1")
    for ch in selected
        haskey(store.data, ch) || error("Canal desconocido: $ch")
    end

    n = length(selected)
    nrows, ncols = _grid_dims(n)
    cell_w, cell_h = 280, 220
    fig_w = max(700, ncols * cell_w + 40)
    fig_h = max(420, nrows * cell_h + 60)

    # rango X compartido (opcional)
    shared_lo, shared_hi = Inf, -Inf
    if x_shared
        for ch in selected
            y = _slice(store, ch, t0, t1)
            isempty(y) && continue
            a, b = extrema(y)
            shared_lo = min(shared_lo, a)
            shared_hi = max(shared_hi, b)
        end
        if !isfinite(shared_lo)
            shared_lo, shared_hi = -1.0, 1.0
        else
            pad = 0.05 * max(shared_hi - shared_lo, 1.0)
            shared_lo -= pad
            shared_hi += pad
        end
    end

    fig = Figure(size = (fig_w, fig_h), fontsize = 12)
    Label(fig[0, 1:ncols],
          "Histograma de amplitud — $SUBJECT / $SESSION / task-$TASK_ID · $(round(t0; digits=1))–$(round(t1; digits=1)) s · $n canales";
          fontsize = 14, font = :bold, tellwidth = false)

    for (i, ch) in enumerate(selected)
        r = div(i - 1, ncols) + 1
        c = mod1(i, ncols)
        y = _slice(store, ch, t0, t1)
        st = _channel_stats(y)
        edges, counts = _histogram(y, nbins)
        centers = @. (edges[1:end-1] + edges[2:end]) / 2
        widths = diff(edges)
        col = BAR_COLORS[mod1(i, length(BAR_COLORS))]

        ax = Axis(fig[r, c];
            title  = "Canal $i: $ch",
            xlabel = (r == nrows || i + ncols > n) ? "Amplitud (µV)" : "",
            ylabel = (c == 1) ? "Frecuencia" : "",
            titlesize = 12,
        )
        ax.xgridvisible = true
        ax.ygridvisible = true

        if !isempty(counts) && maximum(counts) > 0
            barplot!(ax, centers, counts; width = widths, color = col, strokewidth = 0)
        end

        if isfinite(st.mean)
            vlines!(ax, [st.mean]; color = :red, linewidth = 1.8, label = "μ")
            if isfinite(st.std) && st.std > 0
                vlines!(ax, [st.mean - st.std, st.mean + st.std];
                        color = :orange, linestyle = :dash, linewidth = 1.4, label = "μ±σ")
            end
            # anotación en esquina superior derecha (legible)
            xmax = x_shared ? shared_hi : (isempty(edges) ? st.mean : edges[end])
            xmin = x_shared ? shared_lo : (isempty(edges) ? st.mean : edges[1])
            ymax = isempty(counts) ? 1.0 : float(maximum(counts))
            x_ann = xmin + 0.98 * (xmax - xmin)
            text!(ax, x_ann, ymax * 0.98;
                  text = "μ=$(round(st.mean; digits=1))\nσ=$(round(st.std; digits=1))\nσ²=$(round(st.var; digits=0))",
                  align = (:right, :top), fontsize = 10, color = :black,
                  font = :bold)
        end

        if x_shared
            xlims!(ax, shared_lo, shared_hi)
        end
    end

    tag_t = (isapprox(t0, store.t_min; atol=1e-3) && isapprox(t1, store.t_max; atol=1e-3)) ?
            "" : "_$(round(Int, t0))-$(round(Int, t1))s"
    ch_tag = n <= 4 ? join(selected, "_") : "$(n)ch"
    out_name = "hist_" * ch_tag * tag_t * ".png"
    out_path = joinpath(HERE, out_name)
    save(out_path, fig; px_per_unit = 2)
    return out_path
end

function save_hist_combined_png(
    store::RawStore,
    selected::Vector{String},
    t0::Float64,
    t1::Float64,
    nbins::Int,
)::String
    isempty(selected) && error("Selecciona al menos un canal")
    pooled = Float64[]
    for ch in selected
        append!(pooled, _slice(store, ch, t0, t1))
    end
    st = _channel_stats(pooled)
    edges, counts = _histogram(pooled, nbins)
    centers = @. (edges[1:end-1] + edges[2:end]) / 2
    widths = diff(edges)
    # densidad (área ≈ 1)
    bin_w = mean(widths)
    dens_hist = counts ./ (max(st.n, 1) * bin_w)
    xs = collect(range(edges[1], edges[end]; length = 200))
    dens_kde = _kde(pooled, xs)

    fig = Figure(size = (900, 520), fontsize = 13)
    ax = Axis(fig[1, 1];
        title  = "Distribución de amplitud ($(length(selected)) canales) — $SUBJECT / $SESSION",
        xlabel = "Amplitud (µV)",
        ylabel = "Densidad",
        yscale = log10,
    )
    ax.xgridvisible = true
    ax.ygridvisible = true

    # evitar ceros en log
    dens_plot = [d > 0 ? d : NaN for d in dens_hist]
    barplot!(ax, centers, dens_plot; width = widths,
             color = (RGBf(0.45, 0.65, 0.90), 0.85), strokewidth = 0)
    lines!(ax, xs, dens_kde; color = :black, linewidth = 1.8, label = "KDE")

    if isfinite(st.mean)
        vlines!(ax, [st.mean]; color = :blue, linestyle = :dash, linewidth = 1.6, label = "Media")
        if isfinite(st.std) && st.std > 0
            vlines!(ax, [st.mean - 2st.std, st.mean + 2st.std];
                    color = :green, linestyle = :dash, linewidth = 1.4, label = "±2 SD")
            vlines!(ax, [st.mean - 3st.std, st.mean + 3st.std];
                    color = :red, linestyle = :dash, linewidth = 1.4, label = "±3 SD")
        end
    end
    axislegend(ax; position = :rt, framevisible = true, labelsize = 11)

    # pie con stats
    Label(fig[2, 1],
          "μ = $(round(st.mean; digits=2)) µV · σ = $(round(st.std; digits=2)) · σ² = $(round(st.var; digits=1)) · n = $(st.n)";
          fontsize = 11, color = :gray, tellwidth = false)

    tag_t = (isapprox(t0, store.t_min; atol=1e-3) && isapprox(t1, store.t_max; atol=1e-3)) ?
            "" : "_$(round(Int, t0))-$(round(Int, t1))s"
    out_name = "hist_combined_$(length(selected))ch" * tag_t * ".png"
    out_path = joinpath(HERE, out_name)
    save(out_path, fig; px_per_unit = 2)
    return out_path
end

# ── HTTP mínimo ────────────────────────────────────────────────

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
        key = kv[1]
        val = replace(replace(kv[2], "%2C" => ","), "%20" => " ")
        qs[key] = val
    end
    return qs
end

"""Parseo mínimo JSON para /api/save."""
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
    _bool = function (key, default)
        m2 = match(Regex("\"" * key * "\"\\s*:\\s*(true|false)"), s)
        m2 === nothing && return Bool(default)
        return m2.captures[1] == "true"
    end
    return (
        channels = chs,
        t0 = _num("t0", 0.0),
        t1 = _num("t1", 10.0),
        nbins = Int(round(_num("nbins", 40.0))),
        mode = _str("mode", "grid"),
        x_shared = _bool("x_shared", false),
    )
end

function _hist_json_channel(y, nbins::Int, ch::String, idx::Int)::String
    st = _channel_stats(y)
    edges, counts = _histogram(y, nbins)
    centers = @. (edges[1:end-1] + edges[2:end]) / 2
    parts = [
        "\"channel\":\"$ch\"",
        "\"index\":$idx",
        "\"n\":$(st.n)",
        "\"mean\":$(isfinite(st.mean) ? round(st.mean; digits=4) : "null")",
        "\"std\":$(isfinite(st.std) ? round(st.std; digits=4) : "null")",
        "\"var\":$(isfinite(st.var) ? round(st.var; digits=2) : "null")",
        "\"min\":$(isfinite(st.min) ? round(st.min; digits=3) : "null")",
        "\"max\":$(isfinite(st.max) ? round(st.max; digits=3) : "null")",
        "\"centers\":[" * join(round.(centers; digits=3), ",") * "]",
        "\"counts\":[" * join(counts, ",") * "]",
        "\"edges\":[" * join(round.(edges; digits=3), ",") * "]",
    ]
    return "{" * join(parts, ",") * "}"
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
                   "\"n_samples\":$(length(store.t)),\"fs\":$(store.fs)," *
                   "\"channels\":[" * join(["\"$c\"" for c in store.channels], ",") * "]}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/hist"
            qs = _parse_qs(path)
            chs = String[String(c) for c in split(get(qs, "chs", ""), ',') if !isempty(c)]
            t0 = something(tryparse(Float64, get(qs, "t0", "0")), store.t_min)
            t1 = something(tryparse(Float64, get(qs, "t1", string(store.t_max))), store.t_max)
            nbins = something(tryparse(Int, get(qs, "nbins", "40")), 40)
            mode = String(get(qs, "mode", "grid"))
            x_shared = get(qs, "x_shared", "0") in ("1", "true", "True")

            isempty(chs) && error("chs vacío")
            nbins = clamp(nbins, 5, 200)

            if mode == "combined"
                pooled = Float64[]
                for ch in chs
                    haskey(store.data, ch) || continue
                    append!(pooled, _slice(store, ch, t0, t1))
                end
                st = _channel_stats(pooled)
                edges, counts = _histogram(pooled, nbins)
                centers = @. (edges[1:end-1] + edges[2:end]) / 2
                bin_w = mean(diff(edges))
                dens_hist = [c / (max(st.n, 1) * bin_w) for c in counts]
                xs = collect(range(edges[1], edges[end]; length = 160))
                dens_kde = _kde(pooled, xs)
                body = "{" *
                    "\"mode\":\"combined\"," *
                    "\"n_channels\":$(length(chs))," *
                    "\"n\":$(st.n)," *
                    "\"mean\":$(round(st.mean; digits=4))," *
                    "\"std\":$(round(st.std; digits=4))," *
                    "\"var\":$(round(st.var; digits=2))," *
                    "\"centers\":[" * join(round.(centers; digits=3), ",") * "]," *
                    "\"counts\":[" * join(counts, ",") * "]," *
                    "\"density\":[" * join(round.(dens_hist; digits=8), ",") * "]," *
                    "\"kde_x\":[" * join(round.(xs; digits=3), ",") * "]," *
                    "\"kde_y\":[" * join(round.(dens_kde; digits=8), ",") * "]" *
                    "}"
                _send(sock, 200, body; content_type = "application/json")
            else
                items = String[]
                glo, ghi = Inf, -Inf
                for (i, ch) in enumerate(chs)
                    haskey(store.data, ch) || continue
                    y = _slice(store, ch, t0, t1)
                    push!(items, _hist_json_channel(y, nbins, ch, i))
                    if !isempty(y)
                        a, b = extrema(y)
                        glo = min(glo, a); ghi = max(ghi, b)
                    end
                end
                if !isfinite(glo)
                    glo, ghi = -1.0, 1.0
                end
                body = "{" *
                    "\"mode\":\"grid\"," *
                    "\"x_shared\":$(x_shared ? "true" : "false")," *
                    "\"shared_min\":$(round(glo; digits=3))," *
                    "\"shared_max\":$(round(ghi; digits=3))," *
                    "\"items\":[" * join(items, ",") * "]" *
                    "}"
                _send(sock, 200, body; content_type = "application/json")
            end

        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_save_json(raw)
            isempty(req.channels) && error("channels vacío")
            out = if req.mode == "combined"
                save_hist_combined_png(store, req.channels, req.t0, req.t1, req.nbins)
            else
                save_hist_grid_png(store, req.channels, req.t0, req.t1, req.nbins, req.x_shared)
            end
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
    ch_opts = join([
        "<label class=\"ch\"><input type=\"checkbox\" value=\"$c\">$c</label>"
        for c in store.channels], "\n")
    hex_js = join(["'$c'" for c in HEX_COLORS], ",")
    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — Histograma de amplitud</title>
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
  .wrap { max-width:1280px; margin:16px auto; padding:0 16px; }
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
  .checkrow { display:flex; align-items:center; gap:8px; height:34px; font-size:13px; color:var(--muted); }
  .checkrow input { accent-color:var(--accent); }
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
  canvas { width:100%; display:block; background:#fff; border:1px solid var(--border); border-radius:8px; }
  .status { font-size:12px; color:var(--muted); margin-top:8px; min-height:1.2em; }
  .status.ok { color:#047857; }
  .status.err { color:#b91c1c; }
  .legend { font-size:12px; color:var(--muted); margin-top:6px; }
  .legend b { color:#dc2626; font-weight:600; }
  .legend i { color:#ea580c; font-style:normal; font-weight:600; }
</style>
</head>
<body>
<header>
  <h1>Histograma de amplitud</h1>
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
        <label>Bins</label>
        <select id="nbins">
          <option value="20">20</option>
          <option value="30">30</option>
          <option value="40" selected>40</option>
          <option value="50">50</option>
          <option value="80">80</option>
        </select>
      </div>
      <div class="field">
        <label>Vista</label>
        <select id="mode">
          <option value="grid" selected>Rejilla por canal</option>
          <option value="combined">Combinado (densidad)</option>
        </select>
      </div>
      <label class="checkrow">
        <input type="checkbox" id="x_shared"/> Eje X compartido
      </label>
      <button class="btn btn-primary" id="btn-plot" onclick="refreshPlot()">Actualizar</button>
      <button class="btn btn-save" id="btn-save" onclick="savePng()">Guardar PNG</button>
    </div>
  </div>

  <div class="card">
    <div style="display:flex;justify-content:space-between;align-items:baseline;margin-bottom:6px">
      <strong style="font-size:13px">Canales</strong>
      <span style="font-size:12px;color:var(--muted)" id="ch-count">0 seleccionados</span>
    </div>
    <div class="channels" id="channels">$ch_opts</div>
    <div class="quick">
      <button class="chip" type="button" onclick="selectAll(true)">Todos</button>
      <button class="chip" type="button" onclick="selectAll(false)">Ninguno</button>
      <button class="chip" type="button" onclick="preset(['Fz','F3','F4','Cz','C3','C4','Pz','O1','Oz','O2'])">10 centrales</button>
      <button class="chip" type="button" onclick="presetPage(0)">Lámina 1–16</button>
      <button class="chip" type="button" onclick="presetPage(1)">Lámina 17–31</button>
      <button class="chip" type="button" onclick="preset(['O1','Oz','O2'])">Occipitales</button>
    </div>
  </div>

  <div class="card">
    <canvas id="cv"></canvas>
    <div class="legend">
      Barras por canal · <b>línea roja: media (μ)</b> · <i>naranja discontinua: μ ± σ</i>
      · anotación: μ, σ y varianza (σ²)
    </div>
    <div class="status" id="status">Elige canales y pulsa Actualizar.</div>
  </div>
</div>

<script>
const T_MIN = $(store.t_min);
const T_MAX = $(store.t_max);
const ALL_CHS = $(replace(repr(store.channels), "\"" => "'"));
const COLORS = [$hex_js];

let _lastPlot = null;

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
function presetPage(page) {
  const start = page * 16;
  const slice = ALL_CHS.slice(start, start + 16);
  preset(slice);
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

function gridDims(n) {
  if (n <= 1) return [1, 1];
  if (n === 2) return [1, 2];
  if (n <= 4) return [2, 2];
  if (n <= 6) return [2, 3];
  if (n <= 9) return [3, 3];
  if (n <= 12) return [3, 4];
  if (n <= 16) return [4, 4];
  const cols = Math.min(4, Math.ceil(Math.sqrt(n)));
  const rows = Math.ceil(n / cols);
  return [rows, cols];
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

async function refreshPlot() {
  const chs = selectedChannels();
  if (!chs.length) { setStatus('Selecciona al menos un canal.', 'err'); return; }
  const [t0, t1] = timeRange();
  const nbins = document.getElementById('nbins').value;
  const mode = document.getElementById('mode').value;
  const xShared = document.getElementById('x_shared').checked ? '1' : '0';
  setStatus('Calculando histogramas…');
  const url = '/api/hist?chs=' + encodeURIComponent(chs.join(',')) +
              '&t0=' + t0 + '&t1=' + t1 +
              '&nbins=' + nbins + '&mode=' + mode + '&x_shared=' + xShared;
  const res = await fetch(url);
  const data = await res.json();
  if (!res.ok) { setStatus(data.error || 'Error', 'err'); return; }
  _lastPlot = { data, chs, t0, t1, mode };
  if (mode === 'combined') drawCombined(data, chs, t0, t1);
  else drawGrid(data, chs, t0, t1);
  const n = chs.length;
  setStatus(n + ' canal' + (n === 1 ? '' : 'es') + ' · ' +
            t0.toFixed(1) + '–' + t1.toFixed(1) + ' s · modo ' + mode);
}

function drawVLine(ctx, x, y0, y1, color, dash) {
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.5;
  ctx.setLineDash(dash || []);
  ctx.beginPath();
  ctx.moveTo(x, y0); ctx.lineTo(x, y1);
  ctx.stroke();
  ctx.restore();
}

function drawGrid(data, chs, t0, t1) {
  const items = data.items || [];
  const n = items.length;
  if (!n) return;
  const [nrows, ncols] = gridDims(n);
  const cellH = n <= 4 ? 200 : (n <= 9 ? 170 : (n <= 16 ? 150 : 130));
  const cssH = Math.max(360, nrows * cellH + 50);
  const { ctx, W, H } = setupHiDPICanvas(document.getElementById('cv'), cssH);

  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);

  const titleFs = 14;
  ctx.fillStyle = '#0f172a';
  ctx.font = '600 ' + titleFs + 'px "IBM Plex Sans", sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('Histograma de amplitud — ' + n + ' canales · ' +
               t0.toFixed(1) + '–' + t1.toFixed(1) + ' s', W / 2, 18);

  const gapX = 12, gapY = 10;
  const outerL = 8, outerR = 8, outerT = 28, outerB = 8;
  const cellW = (W - outerL - outerR - gapX * (ncols - 1)) / ncols;
  const usableH = H - outerT - outerB - gapY * (nrows - 1);
  const cH = usableH / nrows;
  const xShared = data.x_shared;
  const sMin = data.shared_min, sMax = data.shared_max;

  items.forEach((it, i) => {
    const r = Math.floor(i / ncols);
    const c = i % ncols;
    const x0 = outerL + c * (cellW + gapX);
    const y0 = outerT + r * (cH + gapY);
    drawOneHist(ctx, it, x0, y0, cellW, cH, COLORS[i % COLORS.length],
                xShared, sMin, sMax, i + 1);
  });
}

function drawOneHist(ctx, it, x0, y0, w, h, color, xShared, sMin, sMax, idx) {
  const pad = { l: 42, r: 8, t: 22, b: 28 };
  const pw = w - pad.l - pad.r;
  const ph = h - pad.t - pad.b;
  const left = x0 + pad.l, top = y0 + pad.t;

  // título
  ctx.fillStyle = '#1e293b';
  ctx.font = '600 11px "IBM Plex Sans", sans-serif';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'alphabetic';
  ctx.fillText('Canal ' + idx + ': ' + it.channel, x0 + 6, y0 + 14);

  const centers = it.centers || [];
  const counts = it.counts || [];
  if (!centers.length) return;

  let xmin, xmax;
  if (xShared && isFinite(sMin) && isFinite(sMax) && sMax > sMin) {
    xmin = sMin; xmax = sMax;
  } else {
    xmin = Math.min(...centers);
    xmax = Math.max(...centers);
    const edges = it.edges || [];
    if (edges.length >= 2) { xmin = edges[0]; xmax = edges[edges.length - 1]; }
  }
  const ymax = Math.max(1, Math.max(...counts));
  const xAt = v => left + ((v - xmin) / (xmax - xmin || 1)) * pw;
  const yAt = v => top + (1 - v / ymax) * ph;

  // grid
  ctx.strokeStyle = '#eef2f7';
  ctx.lineWidth = 1;
  for (let i = 0; i <= 3; i++) {
    const y = top + (i / 3) * ph;
    ctx.beginPath(); ctx.moveTo(left, y); ctx.lineTo(left + pw, y); ctx.stroke();
  }

  // barras
  const bw = pw / Math.max(centers.length, 1) * 0.92;
  ctx.fillStyle = color;
  for (let k = 0; k < centers.length; k++) {
    const x = xAt(centers[k]) - bw / 2;
    const y = yAt(counts[k]);
    const bh = top + ph - y;
    if (bh > 0) ctx.fillRect(x, y, bw, bh);
  }

  // frame
  ctx.strokeStyle = '#94a3b8';
  ctx.lineWidth = 1;
  ctx.strokeRect(left, top, pw, ph);

  // media / σ
  if (it.mean != null && isFinite(it.mean)) {
    drawVLine(ctx, xAt(it.mean), top, top + ph, '#dc2626', []);
    if (it.std != null && it.std > 0) {
      drawVLine(ctx, xAt(it.mean - it.std), top, top + ph, '#ea580c', [4, 3]);
      drawVLine(ctx, xAt(it.mean + it.std), top, top + ph, '#ea580c', [4, 3]);
    }
    // caja blanca + anotación (esquina superior derecha, legible sobre barras)
    const lines = [
      'μ = ' + it.mean.toFixed(1),
      'σ = ' + (it.std != null ? it.std.toFixed(1) : '—'),
      'σ² = ' + (it.var != null ? Math.round(it.var) : '—'),
    ];
    ctx.font = '600 11px "IBM Plex Sans", sans-serif';
    const lineH = 14;
    const padBox = 5;
    let tw = 0;
    lines.forEach(s => { tw = Math.max(tw, ctx.measureText(s).width); });
    const boxW = tw + padBox * 2;
    const boxH = lines.length * lineH + padBox * 2 - 2;
    const bx = left + pw - boxW - 4;
    const by = top + 4;
    ctx.fillStyle = 'rgba(255,255,255,0.92)';
    ctx.strokeStyle = '#cbd5e1';
    ctx.lineWidth = 1;
    ctx.fillRect(bx, by, boxW, boxH);
    ctx.strokeRect(bx, by, boxW, boxH);
    ctx.fillStyle = '#0f172a';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'top';
    lines.forEach((s, i) => {
      ctx.fillText(s, bx + padBox, by + padBox + i * lineH);
    });
  }

  // ticks X
  ctx.fillStyle = '#64748b';
  ctx.font = '10px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText(xmin.toFixed(0), left, top + ph + 14);
  ctx.fillText(xmax.toFixed(0), left + pw, top + ph + 14);
  ctx.fillText('µV', left + pw / 2, top + ph + 24);

  // ticks Y
  ctx.textAlign = 'right';
  ctx.fillText('0', left - 4, top + ph + 3);
  ctx.fillText(String(ymax), left - 4, top + 4);
}

function drawCombined(data, chs, t0, t1) {
  const cssH = 480;
  const { ctx, W, H } = setupHiDPICanvas(document.getElementById('cv'), cssH);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);

  const pad = { l: 64, r: 24, t: 36, b: 48 };
  const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;

  ctx.fillStyle = '#2563eb';
  ctx.font = '600 14px "IBM Plex Sans", sans-serif';
  ctx.textAlign = 'left';
  ctx.fillText('Distribución de amplitud (' + chs.length + ' canales)', pad.l, 22);

  const centers = data.centers || [];
  const density = data.density || [];
  const kdeX = data.kde_x || [];
  const kdeY = data.kde_y || [];
  if (!centers.length) return;

  const edges = [];
  // approx edges from centers spacing
  const dx = centers.length > 1 ? (centers[1] - centers[0]) : 1;
  let xmin = centers[0] - dx / 2, xmax = centers[centers.length - 1] + dx / 2;
  if (kdeX.length) { xmin = Math.min(xmin, kdeX[0]); xmax = Math.max(xmax, kdeX[kdeX.length - 1]); }

  // log density range
  const pos = density.filter(d => d > 0).concat(kdeY.filter(d => d > 0));
  let dmin = Math.min(...pos), dmax = Math.max(...pos);
  if (!(dmax > dmin)) { dmin = 1e-6; dmax = 1; }
  dmin = Math.max(dmin * 0.5, dmax * 1e-5);

  const xAt = v => pad.l + ((v - xmin) / (xmax - xmin || 1)) * pw;
  const yAt = v => {
    const lv = Math.log10(Math.max(v, dmin));
    const l0 = Math.log10(dmin), l1 = Math.log10(dmax);
    return pad.t + (1 - (lv - l0) / (l1 - l0 || 1)) * ph;
  };

  // grid
  ctx.strokeStyle = '#e2e8f0'; ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i++) {
    const y = pad.t + (i / 4) * ph;
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + pw, y); ctx.stroke();
  }
  for (let i = 0; i <= 4; i++) {
    const x = pad.l + (i / 4) * pw;
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
    const v = xmin + (i / 4) * (xmax - xmin);
    ctx.fillStyle = '#64748b'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
    ctx.fillText(v.toFixed(0), x, H - 28);
  }

  // hist bars
  const bw = pw / Math.max(centers.length, 1) * 0.9;
  ctx.fillStyle = 'rgba(100, 149, 237, 0.75)';
  for (let k = 0; k < centers.length; k++) {
    const d = density[k];
    if (!(d > 0)) continue;
    const x = xAt(centers[k]) - bw / 2;
    const y = yAt(d);
    const bh = pad.t + ph - y;
    if (bh > 0) ctx.fillRect(x, y, bw, bh);
  }

  // KDE
  if (kdeX.length) {
    ctx.strokeStyle = '#0f172a';
    ctx.lineWidth = 1.8;
    ctx.beginPath();
    let started = false;
    for (let k = 0; k < kdeX.length; k++) {
      if (!(kdeY[k] > 0)) continue;
      const x = xAt(kdeX[k]), y = yAt(kdeY[k]);
      if (!started) { ctx.moveTo(x, y); started = true; }
      else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }

  // stats lines
  const μ = data.mean, σ = data.std;
  if (μ != null && isFinite(μ)) {
    drawVLine(ctx, xAt(μ), pad.t, pad.t + ph, '#2563eb', [5, 4]);
    if (σ > 0) {
      drawVLine(ctx, xAt(μ - 2 * σ), pad.t, pad.t + ph, '#16a34a', [5, 4]);
      drawVLine(ctx, xAt(μ + 2 * σ), pad.t, pad.t + ph, '#16a34a', [5, 4]);
      drawVLine(ctx, xAt(μ - 3 * σ), pad.t, pad.t + ph, '#dc2626', [5, 4]);
      drawVLine(ctx, xAt(μ + 3 * σ), pad.t, pad.t + ph, '#dc2626', [5, 4]);
    }
  }

  ctx.strokeStyle = '#94a3b8';
  ctx.strokeRect(pad.l, pad.t, pw, ph);

  // legend
  const leg = [
    ['#2563eb', 'Media'],
    ['#16a34a', '±2 SD'],
    ['#dc2626', '±3 SD'],
    ['#0f172a', 'KDE'],
  ];
  let lx = pad.l + pw - 110, ly = pad.t + 14;
  ctx.font = '11px sans-serif';
  leg.forEach(([col, lab], i) => {
    ctx.strokeStyle = col; ctx.setLineDash([5, 4]); ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(lx, ly + i * 16); ctx.lineTo(lx + 18, ly + i * 16); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = '#334155'; ctx.textAlign = 'left';
    ctx.fillText(lab, lx + 24, ly + i * 16 + 3);
  });

  // axis labels + stats
  ctx.fillStyle = '#475569';
  ctx.font = '12px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('Amplitud (µV)', pad.l + pw / 2, H - 8);
  ctx.save();
  ctx.translate(16, pad.t + ph / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('Densidad (log)', 0, 0);
  ctx.restore();

  ctx.fillStyle = '#64748b';
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'left';
  ctx.fillText('μ=' + (μ != null ? μ.toFixed(2) : '—') +
               '  σ=' + (σ != null ? σ.toFixed(2) : '—') +
               '  σ²=' + (data.var != null ? Math.round(data.var) : '—') +
               '  n=' + data.n, pad.l, H - 30);
}

async function savePng() {
  const chs = selectedChannels();
  if (!chs.length) { setStatus('Selecciona al menos un canal.', 'err'); return; }
  const [t0, t1] = timeRange();
  const nbins = parseInt(document.getElementById('nbins').value, 10);
  const mode = document.getElementById('mode').value;
  const x_shared = document.getElementById('x_shared').checked;
  const btn = document.getElementById('btn-save');
  btn.disabled = true;
  setStatus('Generando PNG…');
  try {
    const res = await fetch('/api/save', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({channels: chs, t0, t1, nbins, mode, x_shared})
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
    if (!_lastPlot) return;
    if (_lastPlot.mode === 'combined')
      drawCombined(_lastPlot.data, _lastPlot.chs, _lastPlot.t0, _lastPlot.t1);
    else
      drawGrid(_lastPlot.data, _lastPlot.chs, _lastPlot.t0, _lastPlot.t1);
  }, 120);
});

document.getElementById('mode').addEventListener('change', () => {
  if (selectedChannels().length) refreshPlot();
});

// default: lámina 1 (16 primeros) o Cz si hay pocos
updateCount();
if (ALL_CHS.length >= 16) presetPage(0);
else preset(ALL_CHS.slice(0, Math.min(4, ALL_CHS.length)));
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
            "fs=$(round(store.fs; digits=1)) Hz · ",
            "$(round(store.t_min; digits=1))–$(round(store.t_max; digits=1)) s")

    server = listen(IPv4(HOST), PORT)
    url = "http://$HOST:$PORT/"
    println()
    println("UI histograma → $url")
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
