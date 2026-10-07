#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor interactivo de conectividad wPLI (v2)
# ═══════════════════════════════════════════════════════════════
#
#  Threshold XOR Top-N · heatmap (escala auto/0-1) · grafo 10-20
#  · hubs (strength/degree/clustering/betweenness) · interpretación
#
#  Fuente (paso [7/8]):
#    tables/connectivity/wpli_{BAND}.csv
#    tables/connectivity/connectivity_edges.csv
#    tables/connectivity/network_metrics.csv
#    json/connectivity_summary.json
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_connectivity.jl
#    # → http://127.0.0.1:8774/
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_connectivity.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      24-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics

const HERE         = @__DIR__
const RESULTS_UNIT = normpath(joinpath(HERE, "..", "..", "..", "..",
                     "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))
const CONN    = joinpath(RESULTS_UNIT, "tables", "connectivity")
const JSONDIR = joinpath(RESULTS_UNIT, "json")
const HOST    = "127.0.0.1"
const PORT    = 8774

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK_ID = "EC"

const BAND_ORDER = ["DELTA", "THETA", "ALPHA", "BETA_LOW", "BETA_MID", "BETA_HIGH", "GAMMA"]

const CH_POS = Dict{String,Tuple{Float64,Float64}}(
    "FZ"=>(0.00,0.72),"F3"=>(-0.35,0.55),"F4"=>(0.35,0.55),
    "F7"=>(-0.68,0.42),"F8"=>(0.68,0.42),"FT9"=>(-0.85,0.18),"FT10"=>(0.85,0.18),
    "FC5"=>(-0.50,0.28),"FC1"=>(-0.18,0.28),"FC2"=>(0.18,0.28),"FC6"=>(0.50,0.28),
    "C3"=>(-0.40,0.00),"CZ"=>(0.00,0.00),"C4"=>(0.40,0.00),
    "T7"=>(-0.80,0.00),"T8"=>(0.80,0.00),"TP9"=>(-0.85,-0.22),"TP10"=>(0.85,-0.22),
    "CP5"=>(-0.50,-0.28),"CP1"=>(-0.18,-0.28),"CP2"=>(0.18,-0.28),"CP6"=>(0.50,-0.28),
    "P3"=>(-0.35,-0.55),"PZ"=>(0.00,-0.55),"P4"=>(0.35,-0.55),
    "P7"=>(-0.68,-0.48),"P8"=>(0.68,-0.48),
    "O1"=>(-0.28,-0.82),"OZ"=>(0.00,-0.88),"O2"=>(0.28,-0.82),
)

const CH_REGION = Dict{String,String}(
    "FZ"=>"frontal","F3"=>"frontal","F4"=>"frontal","F7"=>"frontal","F8"=>"frontal",
    "FT9"=>"temporal","FT10"=>"temporal","FC5"=>"frontal","FC1"=>"frontal",
    "FC2"=>"frontal","FC6"=>"frontal","C3"=>"central","CZ"=>"central","C4"=>"central",
    "T7"=>"temporal","T8"=>"temporal","TP9"=>"temporal","TP10"=>"temporal",
    "CP5"=>"parietal","CP1"=>"parietal","CP2"=>"parietal","CP6"=>"parietal",
    "P3"=>"parietal","PZ"=>"parietal","P4"=>"parietal","P7"=>"parietal","P8"=>"parietal",
    "O1"=>"occipital","OZ"=>"occipital","O2"=>"occipital",
)

const REGION_LABEL = Dict(
    "frontal"=>"Frontal","central"=>"Central","parietal"=>"Parietal",
    "occipital"=>"Occipital","temporal"=>"Temporal",
)

mutable struct ConnStore
    channels::Vector{String}
    bands::Vector{String}
    matrices::Dict{String,Matrix{Float64}}
    edges::DataFrame
    metrics::DataFrame
    summary::Dict{String,Any}
    thr_default::Float64
    n_epochs::Int
end

# ── IO ─────────────────────────────────────────────────────────

function _parse_summary_json(path::String)::Dict{String,Any}
    out = Dict{String,Any}()
    isfile(path) || return out
    txt = read(path, String)
    for m in eachmatch(r"\"([^\"]+)\"\s*:\s*\"([^\"]*)\"", txt)
        out[String(m.captures[1])] = String(m.captures[2])
    end
    for m in eachmatch(r"\"([^\"]+)\"\s*:\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", txt)
        k = String(m.captures[1])
        haskey(out, k) && continue
        v = tryparse(Float64, m.captures[2])
        v === nothing || (out[k] = v)
    end
    return out
end

function _read_wpli_matrix(path::String)::Tuple{Vector{String},Matrix{Float64}}
    df = CSV.read(path, DataFrame)
    ch = String.(df[!, 1])
    n = length(ch)
    M = zeros(Float64, n, n)
    for (j, name) in enumerate(ch)
        col = Symbol(name)
        hasproperty(df, col) || error("Columna faltante $name en $path")
        M[:, j] = Float64.(df[!, col])
    end
    return ch, M
end

function load_store()::ConnStore
    isdir(CONN) || error("No encontrado: $CONN")
    edges_path = joinpath(CONN, "connectivity_edges.csv")
    metrics_path = joinpath(CONN, "network_metrics.csv")
    isfile(edges_path) || error("No encontrado: $edges_path")

    matrices = Dict{String,Matrix{Float64}}()
    channels = String[]
    bands = String[]
    for b in BAND_ORDER
        p = joinpath(CONN, "wpli_$(b).csv")
        isfile(p) || continue
        ch, M = _read_wpli_matrix(p)
        if isempty(channels)
            channels = ch
        else
            ch == channels || error("Canales distintos en $p")
        end
        matrices[b] = M
        push!(bands, b)
    end
    isempty(bands) && error("Ninguna matriz wpli_*.csv en $CONN")

    edges = CSV.read(edges_path, DataFrame)
    edges.ch_a = String.(edges.ch_a)
    edges.ch_b = String.(edges.ch_b)
    edges.band = String.(edges.band)

    metrics = isfile(metrics_path) ? CSV.read(metrics_path, DataFrame) : DataFrame()
    if nrow(metrics) > 0 && hasproperty(metrics, :channel)
        metrics.channel = String.(metrics.channel)
    end

    summary = _parse_summary_json(joinpath(JSONDIR, "connectivity_summary.json"))
    thr = Float64(get(summary, "threshold", 0.1))
    n_ep = Int(round(Float64(get(summary, "n_epochs_used", 0))))
    return ConnStore(channels, bands, matrices, edges, metrics, summary, thr, n_ep)
end

# ── Ciencia / métricas ─────────────────────────────────────────

function _band_strength(M::Matrix{Float64})::Vector{Float64}
    n = size(M, 1)
    s = zeros(Float64, n)
    n <= 1 && return s
    for i in 1:n
        s[i] = sum(M[i, j] for j in 1:n if j != i) / (n - 1)
    end
    return s
end

function _fisher_z(w::Float64)::Float64
    atanh(clamp(w, -0.999999, 0.999999))
end

function _upper_values(M::Matrix{Float64})::Vector{Float64}
    n = size(M, 1)
    out = Float64[]
    sizehint!(out, n * (n - 1) ÷ 2)
    for i in 1:n-1, j in i+1:n
        push!(out, M[i, j])
    end
    return out
end

function _percentile(v::Vector{Float64}, p::Float64)::Float64
    isempty(v) && return NaN
    s = sort(v)
    n = length(s)
    x = clamp(p, 0.0, 1.0) * (n - 1) + 1
    i = floor(Int, x)
    f = x - i
    i >= n && return s[end]
    i < 1 && return s[1]
    return s[i] * (1 - f) + s[min(i + 1, n)] * f
end

function _region(ch::String)::String
    get(CH_REGION, uppercase(ch), "other")
end

function _pos_xy(ch::String)::Tuple{Float64,Float64}
    get(CH_POS, uppercase(ch), (NaN, NaN))
end

function _edges_for_band(store::ConnStore, band::String)
    store.edges[store.edges.band .== band, :]
end

"""Select graph edges according to mode=threshold|topn."""
function _select_graph_edges(edf::DataFrame, mode::String, thr::Float64, topn::Int)::DataFrame
    edf = sort(edf, :rank)
    if mode == "topn"
        return first(edf, min(topn, nrow(edf)))
    else
        keep = edf[edf.wpli .>= thr, :]
        nrow(keep) > 400 && (keep = first(sort(keep, :wpli; rev=true), 400))
        return keep
    end
end

"""Binary adjacency from edge list (1-based indices)."""
function _binary_adj(channels::Vector{String}, graph_df::DataFrame)::BitMatrix
    n = length(channels)
    idx = Dict(c => i for (i, c) in enumerate(channels))
    A = falses(n, n)
    for r in eachrow(graph_df)
        i = get(idx, String(r.ch_a), 0)
        j = get(idx, String(r.ch_b), 0)
        (i == 0 || j == 0) && continue
        A[i, j] = true
        A[j, i] = true
    end
    return A
end

function _degree_from_adj(A::BitMatrix)::Vector{Int}
    [count(A[i, :]) for i in 1:size(A, 1)]
end

"""Local clustering coefficient (unweighted)."""
function _clustering(A::BitMatrix)::Vector{Float64}
    n = size(A, 1)
    c = fill(0.0, n)
    for i in 1:n
        nbrs = findall(A[i, :])
        k = length(nbrs)
        k < 2 && continue
        links = 0
        for a in 1:k-1, b in a+1:k
            A[nbrs[a], nbrs[b]] && (links += 1)
        end
        c[i] = 2 * links / (k * (k - 1))
    end
    return c
end

"""Brandes betweenness, normalized by (n-1)(n-2)."""
function _betweenness(A::BitMatrix)::Vector{Float64}
    n = size(A, 1)
    n <= 2 && return zeros(n)
    CB = zeros(Float64, n)
    for s in 1:n
        S = Int[]
        P = [Int[] for _ in 1:n]
        sigma = zeros(Float64, n); sigma[s] = 1.0
        dist = fill(-1, n); dist[s] = 0
        Q = Int[s]
        while !isempty(Q)
            v = popfirst!(Q)
            push!(S, v)
            for w in findall(A[v, :])
                if dist[w] < 0
                    push!(Q, w)
                    dist[w] = dist[v] + 1
                end
                if dist[w] == dist[v] + 1
                    sigma[w] += sigma[v]
                    push!(P[w], v)
                end
            end
        end
        delta = zeros(Float64, n)
        while !isempty(S)
            w = pop!(S)
            for v in P[w]
                delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            end
            w != s && (CB[w] += delta[w])
        end
    end
    # undirected: Brandes accumulates both directions; divide by 2
    CB .*= 0.5
    denom = (n - 1) * (n - 2)
    denom > 0 && (CB ./= denom)
    return CB
end

function _wpli_color(t::Float64)::RGBf
    # light blue → blue → violet
    t = clamp(t, 0.0, 1.0)
    if t < 0.5
        u = t / 0.5
        return RGBf(0.55 + (0.15 - 0.55) * u, 0.75 + (0.35 - 0.75) * u, 0.95 + (0.85 - 0.95) * u)
    else
        u = (t - 0.5) / 0.5
        return RGBf(0.15 + (0.50 - 0.15) * u, 0.35 + (0.15 - 0.35) * u, 0.85 + (0.70 - 0.85) * u)
    end
end

function _density_label(d::Float64)::String
    d < 0.25 && return "baja"
    d < 0.55 && return "moderada"
    return "alta"
end

function _interpret(
    band::String,
    dens::Float64,
    dens_all::Dict{String,Float64},
    hubs::AbstractVector,
    top_edges::DataFrame,
)::String
    lines = String[]
    # dens must be thr_default density (same basis as dens_all / comparativa)
    dens_pct = round(100 * dens; digits=1)
    push!(lines, "La banda $band presenta una densidad $(_density_label(dens)) ($(dens_pct)% de pares por encima del umbral de referencia).")

    # relative to other bands (all at thr_default — comparable across bands)
    if !isempty(dens_all)
        others = [dens_all[b] for b in keys(dens_all) if b != band]
        if !isempty(others)
            mu = mean(others)
            if dens > mu * 1.25
                push!(lines, "Su densidad es claramente superior a la media del resto de bandas.")
            elseif dens < mu * 0.75
                push!(lines, "Su densidad queda por debajo de la media del resto de bandas.")
            else
                push!(lines, "Su densidad es comparable a la del resto de bandas.")
            end
        end
    end

    if length(hubs) >= 2
        s1, s2 = hubs[1].strength, hubs[2].strength
        if s2 > 0 && s1 / s2 > 1.15
            push!(lines, "El nodo con mayor strength es $(hubs[1].channel) ($(round(s1; digits=3))); destaca frente a $(hubs[2].channel).")
        else
            push!(lines, "No se observa un hub claramente dominante (top: $(hubs[1].channel), $(hubs[2].channel)).")
        end
    elseif length(hubs) == 1
        push!(lines, "El nodo con mayor strength es $(hubs[1].channel).")
    end

    if nrow(top_edges) > 0
        nshow = min(5, nrow(top_edges))
        counts = Dict{String,Int}()
        for r in eachrow(first(top_edges, nshow))
            for ch in (String(r.ch_a), String(r.ch_b))
                reg = _region(ch)
                counts[reg] = get(counts, reg, 0) + 1
            end
        end
        if !isempty(counts)
            best = argmax(counts)
            push!(lines, "Las conexiones más fuertes se concentran en región $(get(REGION_LABEL, best, best)).")
        end
    end

    band == "DELTA" && push!(lines, "Aviso: DELTA puede ser menos fiable con epochs cortos (pocos ciclos por segmento).")
    return join(lines, " ")
end

# ── PNG ────────────────────────────────────────────────────────

function save_heatmap_png(store::ConnStore, band::String, scale::String)::String
    haskey(store.matrices, band) || error("Banda desconocida: $band")
    M0 = store.matrices[band]
    ch = store.channels
    n = length(ch)
    M = copy(M0)
    for i in 1:n
        M[i, i] = 0.0  # placeholder; grey overlay drawn separately
    end
    off = Float64[M0[i, j] for i in 1:n for j in 1:n if i != j]
    vmax = scale == "fixed" ? 1.0 : (isempty(off) ? 1.0 : maximum(off))
    vmax <= 0 && (vmax = 1.0)

    fig = Figure(size = (820, 720), fontsize = 12)
    ax = Axis(fig[1, 1];
        title = "wPLI $band — $SUBJECT/$SESSION/$TASK_ID · escala=$(scale == "fixed" ? "0-1" : "auto")",
        xlabel = "Canal", ylabel = "Canal",
        xticks = (1:n, ch), yticks = (1:n, ch),
        xticklabelrotation = π/3, xticklabelsize = 8, yticklabelsize = 8,
    )
    hm = heatmap!(ax, M; colormap = :viridis, colorrange = (0.0, vmax))
    for i in 1:n
        poly!(ax, Point2f[(i - 0.5, i - 0.5), (i + 0.5, i - 0.5), (i + 0.5, i + 0.5), (i - 0.5, i + 0.5)];
              color = RGBf(0.88, 0.90, 0.93), strokewidth = 0)
    end
    Colorbar(fig[1, 2], hm; label = "wPLI")
    out = joinpath(HERE, "connectivity_heatmap_$(band).png")
    save(out, fig; px_per_unit = 3)
    return out
end

function save_topo_png(store::ConnStore, band::String, mode::String, thr::Float64, topn::Int)::String
    haskey(store.matrices, band) || error("Banda desconocida: $band")
    M = store.matrices[band]
    ch = store.channels
    strength = _band_strength(M)
    smax = max(maximum(strength), 1e-9)
    edf = _edges_for_band(store, band)
    keep = _select_graph_edges(edf, mode, thr, topn)

    fig = Figure(size = (780, 780), fontsize = 12)
    mode_lbl = mode == "topn" ? "top-$topn" : "thr≥$(round(thr; digits=2))"
    ax = Axis(fig[1, 1];
        title = "Grafo wPLI $band · $mode_lbl — $SUBJECT/$SESSION/$TASK_ID",
        aspect = DataAspect(),
    )
    hidedecorations!(ax); hidespines!(ax)
    xlims!(ax, -1.25, 1.25); ylims!(ax, -1.25, 1.35)

    θ = range(0, 2π; length=256)
    lines!(ax, sin.(θ), cos.(θ); color = :black, linewidth = 1.8)
    θn = range(-π/6, π/6; length=30)
    lines!(ax, 0.12 .* sin.(θn), 1.0 .+ 0.10 .* cos.(θn); color = :black, linewidth = 1.6)

    # region labels
    for (lab, (x, y)) in (("Frontal", (0.0, 1.12)), ("Central", (0.0, 0.08)),
                          ("Parietal", (0.0, -0.62)), ("Occipital", (0.0, -1.05)),
                          ("Temporal L", (-1.12, 0.0)), ("Temporal R", (1.12, 0.0)))
        text!(ax, x, y; text = lab, align = (:center, :center),
              fontsize = 11, color = (RGBf(0.55, 0.60, 0.68), 0.55))
    end

    wmax = nrow(keep) > 0 ? maximum(Float64.(keep.wpli)) : 1.0
    for r in eachrow(keep)
        xa, ya = _pos_xy(String(r.ch_a))
        xb, yb = _pos_xy(String(r.ch_b))
        (isnan(xa) || isnan(xb)) && continue
        t = Float64(r.wpli) / wmax
        col = _wpli_color(t)
        lines!(ax, [xa, xb], [ya, yb]; color = (col, 0.20 + 0.70 * t), linewidth = 1.4)
    end

    xs = Float64[]; ys = Float64[]; ms = Float64[]; labs = String[]
    for (i, c) in enumerate(ch)
        x, y = _pos_xy(c)
        isnan(x) && continue
        push!(xs, x); push!(ys, y)
        push!(ms, 8 + 18 * (strength[i] / smax))
        push!(labs, c)
    end
    scatter!(ax, xs, ys; color = :white, markersize = ms,
             strokecolor = RGBf(0.05, 0.45, 0.45), strokewidth = 2)
    for (x, y, lab) in zip(xs, ys, labs)
        text!(ax, x, y + 0.06; text = lab, align = (:center, :bottom),
              fontsize = 9, color = (:black, 0.85))
    end

    tag = mode == "topn" ? "top$(topn)" : "thr$(round(Int, thr * 100))"
    out = joinpath(HERE, "connectivity_topo_$(band)_$(tag).png")
    save(out, fig; px_per_unit = 3)
    return out
end

function save_hubs_png(store::ConnStore, band::String)::String
    haskey(store.matrices, band) || error("Banda desconocida: $band")
    M = store.matrices[band]
    ch = store.channels
    strength = _band_strength(M)
    order = sortperm(strength; rev=true)
    names = ch[order]
    vals = strength[order]
    fig = Figure(size = (900, 520), fontsize = 12)
    ax = Axis(fig[1, 1];
        title = "Hubs (strength = media wPLI del nodo) — $band · $SUBJECT/$SESSION/$TASK_ID",
        xlabel = "Canal", ylabel = "Strength",
        xticks = (1:length(names), names),
        xticklabelrotation = π/3, xticklabelsize = 10,
    )
    barplot!(ax, 1:length(vals), vals; color = RGBf(0.15, 0.55, 0.55))
    out = joinpath(HERE, "connectivity_hubs_$(band).png")
    save(out, fig; px_per_unit = 3)
    return out
end

function save_connectivity_png(
    store::ConnStore, view::String, band::String, mode::String,
    thr::Float64, topn::Int, scale::String,
)::String
    view == "heatmap" && return save_heatmap_png(store, band, scale)
    view == "topo" && return save_topo_png(store, band, mode, thr, topn)
    view == "hubs" && return save_hubs_png(store, band)
    error("Vista desconocida: $view")
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

function _read_body(sock, headers::Dict{String,String})::String
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
    write(sock,
        "HTTP/1.1 $status\r\nContent-Type: $content_type\r\n" *
        "Content-Length: $(sizeof(body))\r\nConnection: close\r\n" *
        "Access-Control-Allow-Origin: *\r\n\r\n" * body)
end

function _json_escape(s::AbstractString)
    t = replace(String(s), "\\" => "\\\\")
    t = replace(t, "\"" => "\\\"")
    t = replace(t, "\n" => "\\n")
    t = replace(t, "\r" => "\\r")
    return t
end

function _parse_qs(path::AbstractString)::Dict{String,String}
    qs = Dict{String,String}()
    occursin('?', path) || return qs
    for pair in split(split(path, '?', limit=2)[2], '&')
        kv = split(pair, '=', limit=2)
        length(kv) == 2 && (qs[kv[1]] = kv[2])
    end
    return qs
end

function _parse_json(s::String)
    _str = (key, default) -> begin
        m = match(Regex("\"" * key * "\"\\s*:\\s*\"([^\"]*)\""), s)
        m === nothing ? String(default) : String(m.captures[1])
    end
    _num = (key, default) -> begin
        m = match(Regex("\"" * key * "\"\\s*:\\s*([-\\d.eE]+)"), s)
        m === nothing && return Float64(default)
        v = tryparse(Float64, m.captures[1])
        v === nothing ? Float64(default) : v
    end
    return (
        view  = _str("view", "heatmap"),
        band  = _str("band", "ALPHA"),
        mode  = _str("mode", "threshold"),
        thr   = _num("thr", 0.3),
        topn  = Int(round(_num("topn", 20))),
        scale = _str("scale", "auto"),
    )
end

# ── API payloads ───────────────────────────────────────────────

function _matrix_json(store::ConnStore, band::String)::String
    haskey(store.matrices, band) || error("Banda desconocida: $band")
    M = store.matrices[band]
    ch = store.channels
    n = length(ch)
    strength = _band_strength(M)
    flat = Float64[]
    sizehint!(flat, n * n)
    for j in 1:n, i in 1:n
        push!(flat, M[i, j])
    end
    ch_j = join(["\"$c\"" for c in ch], ",")
    return "{\"band\":\"$band\",\"n\":$n,\"channels\":[$ch_j]," *
           "\"matrix\":[" * join(round.(flat; digits=6), ",") * "]," *
           "\"strength\":[" * join(round.(strength; digits=6), ",") * "]," *
           "\"vmax\":$(round(maximum(M); digits=6))," *
           "\"n_edges\":$(n * (n - 1) ÷ 2)}"
end

function _edges_json(store::ConnStore, band::String, mode::String, thr::Float64, topn::Int)::String
    haskey(store.matrices, band) || error("Banda desconocida: $band")
    mode in ("threshold", "topn") || error("Modo desconocido: $mode")
    M = store.matrices[band]
    ch = store.channels
    n = length(ch)
    edf = sort(_edges_for_band(store, band), :rank)
    vals = _upper_values(M)
    strength = _band_strength(M)

    graph_df = _select_graph_edges(edf, mode, thr, topn)
    A = _binary_adj(ch, graph_df)
    degree = _degree_from_adj(A)
    clust = _clustering(A)
    betw = _betweenness(A)

    # KPIs on full distribution. Graph size (n_edges) depends on mode;
    # density for interpretation/KPI uses thr_default so bands are comparable
    # (topn density = topn/n_pairs is nearly constant and not informative).
    mean_v = isempty(vals) ? 0.0 : mean(vals)
    med_v  = isempty(vals) ? 0.0 : median(vals)
    std_v  = isempty(vals) ? 0.0 : std(vals)
    max_v  = isempty(vals) ? 0.0 : maximum(vals)
    p90_v  = _percentile(vals, 0.90)
    n_graph = nrow(graph_df)
    n_total = length(vals)
    n_above = mode == "topn" ? n_graph : count(>=(thr), vals)
    mean_str = mean(strength)

    hubs = [(channel=ch[i], strength=strength[i], degree=degree[i],
             clustering=clust[i], betweenness=betw[i]) for i in sortperm(strength; rev=true)]

    hubs_j = join([
        "{\"channel\":\"$(h.channel)\",\"strength\":$(round(h.strength; digits=6))," *
        "\"degree\":$(h.degree),\"clustering\":$(round(h.clustering; digits=4))," *
        "\"betweenness\":$(round(h.betweenness; digits=4))}"
        for h in hubs
    ], ",")

    # table: top connections always by rank (for interpretation), capped
    top_tbl = first(edf, min(max(topn, 20), nrow(edf)))
    edges_j = join([
        "{\"ch_a\":\"$(r.ch_a)\",\"ch_b\":\"$(r.ch_b)\",\"wpli\":$(round(Float64(r.wpli); digits=6))," *
        "\"rank\":$(Int(r.rank)),\"fisher_z\":$(round(_fisher_z(Float64(r.wpli)); digits=4))}"
        for r in eachrow(top_tbl)
    ], ",")

    graph_j = join([
        "{\"ch_a\":\"$(r.ch_a)\",\"ch_b\":\"$(r.ch_b)\",\"wpli\":$(round(Float64(r.wpli); digits=6))," *
        "\"rank\":$(Int(r.rank))}"
        for r in eachrow(graph_df)
    ], ",")

    # hist
    nb = 24
    lo, hi = 0.0, max(max_v, 1e-6)
    bw = (hi - lo) / nb
    counts = zeros(Int, nb)
    for v in vals
        k = min(nb, max(1, Int(floor((v - lo) / bw)) + 1))
        v >= hi && (k = nb)
        counts[k] += 1
    end

    # densities at thr_default for all bands (comparativa + interpretación)
    dens_all = Dict{String,Float64}()
    for b in store.bands
        vv = _upper_values(store.matrices[b])
        dens_all[b] = isempty(vv) ? 0.0 : count(>=(store.thr_default), vv) / length(vv)
    end
    dens = get(dens_all, band, 0.0)

    interp = _interpret(band, dens, dens_all, hubs, top_tbl)

    return "{\"band\":\"$band\",\"mode\":\"$mode\",\"thr\":$thr,\"topn\":$topn," *
           "\"n_nodes\":$n,\"n_edges\":$n_graph,\"n_total\":$n_total,\"n_above\":$n_above," *
           "\"mean\":$(round(mean_v; digits=6)),\"median\":$(round(med_v; digits=6))," *
           "\"std\":$(round(std_v; digits=6)),\"max\":$(round(max_v; digits=6))," *
           "\"p90\":$(round(p90_v; digits=6)),\"density\":$(round(dens; digits=6))," *
           "\"mean_strength\":$(round(mean_str; digits=6))," *
           "\"edges\":[$edges_j],\"graph_edges\":[$graph_j],\"hubs\":[$hubs_j]," *
           "\"hist\":{\"lo\":$lo,\"hi\":$hi,\"nbins\":$nb,\"counts\":[$(join(string.(counts), ","))]}," *
           "\"interpretation\":\"$(_json_escape(interp))\"}"
end

function _summary_json(store::ConnStore)::String
    parts = String[]
    thr = store.thr_default
    for b in store.bands
        M = store.matrices[b]
        vals = _upper_values(M)
        strength = _band_strength(M)
        dens = isempty(vals) ? 0.0 : count(>=(thr), vals) / length(vals)
        push!(parts,
            "\"$b\":{\"mean\":$(round(mean(vals); digits=6))," *
            "\"density\":$(round(dens; digits=6))," *
            "\"mean_strength\":$(round(mean(strength); digits=6))," *
            "\"max\":$(round(maximum(vals); digits=6))," *
            "\"n_edges\":$(length(vals))}")
    end
    method = string(get(store.summary, "method", "wpli"))
    return "{\"method\":\"$(_json_escape(method))\"," *
           "\"n_channels\":$(length(store.channels))," *
           "\"n_epochs\":$(store.n_epochs)," *
           "\"threshold\":$(store.thr_default)," *
           "\"bands\":{" * join(parts, ",") * "}}"
end

function handle_request(sock, store::ConnStore)
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
            ch_j = join(["\"$c\"" for c in store.channels], ",")
            b_j = join(["\"$b\"" for b in store.bands], ",")
            pos_parts = ["\"$c\":[$(round(_pos_xy(c)[1]; digits=4)),$(round(_pos_xy(c)[2]; digits=4))]" for c in store.channels]
            reg_parts = ["\"$c\":\"$(_region(c))\"" for c in store.channels]
            body = "{\"channels\":[$ch_j],\"bands\":[$b_j]," *
                   "\"thr_default\":$(store.thr_default)," *
                   "\"n_epochs\":$(store.n_epochs)," *
                   "\"n_channels\":$(length(store.channels))," *
                   "\"subject\":\"$SUBJECT\",\"session\":\"$SESSION\",\"task\":\"$TASK_ID\"," *
                   "\"positions\":{" * join(pos_parts, ",") * "}," *
                   "\"regions\":{" * join(reg_parts, ",") * "}}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/matrix"
            qs = _parse_qs(path)
            _send(sock, 200, _matrix_json(store, String(get(qs, "band", "ALPHA"))); content_type = "application/json")

        elseif method == "GET" && path_only == "/api/edges"
            qs = _parse_qs(path)
            band = String(get(qs, "band", "ALPHA"))
            mode = String(get(qs, "mode", "threshold"))
            thr = something(tryparse(Float64, get(qs, "thr", "0.3")), 0.3)
            topn = something(tryparse(Int, get(qs, "topn", "20")), 20)
            _send(sock, 200, _edges_json(store, band, mode, thr, topn); content_type = "application/json")

        elseif method == "GET" && path_only == "/api/summary"
            _send(sock, 200, _summary_json(store); content_type = "application/json")

        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_json(raw)
            out = save_connectivity_png(store, req.view, req.band, req.mode, req.thr, req.topn, req.scale)
            body = "{\"ok\":true,\"path\":\"$(_json_escape(out))\",\"file\":\"$(_json_escape(basename(out)))\"}"
            _send(sock, 200, body; content_type = "application/json")
            println("[$(Dates.format(now(), "HH:MM:SS"))] PNG → $out")
        else
            _send(sock, 404, "{\"error\":\"not found\"}"; content_type = "application/json")
        end
    catch e
        msg = sprint(showerror, e)
        @warn "Request error" exception = e
        _send(sock, 400, "{\"ok\":false,\"error\":\"$(_json_escape(msg))\"}"; content_type = "application/json")
    end
end

# ── HTML ───────────────────────────────────────────────────────

function html_page(store::ConnStore)::String
    band_opts = join([
        "<option value=\"$b\"" * (b == "ALPHA" ? " selected" : "") * ">$b</option>"
        for b in store.bands
    ], "\n")
    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — Conectividad wPLI v2</title>
<style>
  :root {
    --bg:#f4f6f8; --card:#fff; --border:#d8dee6; --text:#1e293b;
    --muted:#64748b; --accent:#2563eb; --accent2:#0f766e; --warn:#b45309;
  }
  * { box-sizing: border-box; }
  body { margin:0; font-family:"IBM Plex Sans","Segoe UI",sans-serif; background:var(--bg); color:var(--text); }
  header { padding:14px 20px; background:#0f172a; color:#e2e8f0; display:flex; align-items:baseline; gap:14px; flex-wrap:wrap; }
  header h1 { margin:0; font-size:16px; font-weight:600; }
  header span { font-size:12px; color:#94a3b8; }
  .wrap { max-width:1440px; margin:16px auto; padding:0 16px; }
  .card { background:var(--card); border:1px solid var(--border); border-radius:10px; padding:14px 16px; margin-bottom:14px; }
  .toolbar { display:flex; flex-wrap:wrap; gap:12px 18px; align-items:end; }
  .field { display:flex; flex-direction:column; gap:4px; }
  .field label { font-size:11px; color:var(--muted); text-transform:uppercase; letter-spacing:.04em; }
  select, input[type=range] {
    height:34px; padding:0 10px; border:1px solid var(--border);
    border-radius:6px; background:#fff; font-size:13px; min-width:110px;
  }
  input[type=range] { width:150px; padding:0; }
  select:disabled, input:disabled { opacity:.45; background:#f1f5f9; }
  .modebox { display:flex; gap:14px; align-items:center; height:34px; font-size:13px; }
  .modebox label { display:flex; align-items:center; gap:6px; cursor:pointer; }
  .btn { height:34px; padding:0 16px; border:none; border-radius:6px; font-size:13px; font-weight:600; cursor:pointer; }
  .btn-primary { background:var(--accent); color:#fff; }
  .btn-save { background:var(--accent2); color:#fff; }
  .btn:disabled { opacity:.5; cursor:wait; }
  .src { font-size:12.5px; line-height:1.45; color:#334155; background:#f0f9ff; border:1px solid #bae6fd; border-radius:8px; padding:10px 12px; margin-bottom:12px; }
  .src code { font-size:12px; background:#e0f2fe; padding:1px 5px; border-radius:4px; }
  .warn { color:var(--warn); font-size:12px; margin-top:6px; }
  .grid2 { display:grid; grid-template-columns:1fr 1fr; gap:14px; }
  @media (max-width:1000px) { .grid2 { grid-template-columns:1fr; } }
  .panel-title { font-size:13px; font-weight:600; margin:0 0 8px; }
  canvas { width:100%; display:block; background:#fff; border:1px solid var(--border); border-radius:8px; cursor:crosshair; }
  .status { font-size:12px; color:var(--muted); margin-top:8px; min-height:1.2em; }
  .status.ok { color:#047857; } .status.err { color:#b91c1c; }
  .kpis { display:flex; flex-wrap:wrap; gap:8px; margin-bottom:12px; }
  .kpi { background:#f8fafc; border:1px solid var(--border); border-radius:8px; padding:8px 10px; min-width:88px; }
  .kpi .v { font-size:16px; font-weight:700; font-variant-numeric:tabular-nums; }
  .kpi .l { font-size:10px; color:var(--muted); text-transform:uppercase; }
  .interp { background:#f8fafc; border-left:3px solid var(--accent2); padding:10px 12px; margin:0 0 14px; font-size:13px; line-height:1.5; color:#334155; }
  .interp strong { display:block; margin-bottom:4px; font-size:12px; text-transform:uppercase; color:var(--muted); letter-spacing:.04em; }
  .node-card { background:#fffbeb; border:1px solid #fcd34d; border-radius:8px; padding:10px 12px; margin:0 0 12px; font-size:12.5px; display:none; }
  .node-card h4 { margin:0 0 6px; font-size:14px; }
  .node-card .row { display:flex; justify-content:space-between; gap:12px; padding:2px 0; }
  .node-card .k { color:var(--muted); }
  table.tbl { width:100%; border-collapse:collapse; font-size:12.5px; }
  table.tbl th, table.tbl td { padding:7px 8px; border-bottom:1px solid var(--border); text-align:right; }
  table.tbl th:first-child, table.tbl td:first-child { text-align:left; }
  table.tbl th { background:#f8fafc; color:#475569; font-size:11px; text-transform:uppercase; cursor:help; }
  table.tbl tr.sel td { background:#bfdbfe; outline:2px solid #2563eb; }
  table.tbl tr:hover td { background:#f1f5f9; cursor:pointer; }
  .cols { display:grid; grid-template-columns:1.1fr 1fr 1.1fr; gap:14px; }
  @media (max-width:1100px) { .cols { grid-template-columns:1fr; } }
  #tooltip { position:fixed; pointer-events:none; z-index:20; background:#0f172a; color:#e2e8f0; font-size:12px; padding:6px 10px; border-radius:6px; display:none; max-width:280px; }
  table.bandtbl th { font-size:10px; }
  table.bandtbl tr.active td { background:#dbeafe; font-weight:600; }
</style>
</head>
<body>
<header>
  <h1>Conectividad funcional — wPLI v2</h1>
  <span>$SUBJECT / $SESSION / task-$TASK_ID · $(length(store.channels)) canales · $(store.n_epochs) epochs</span>
</header>
<div class="wrap">
  <div class="card">
    <div class="src">
      <strong>Fuente (paso [7/8]).</strong>
      <code>wpli_*.csv</code> · <code>connectivity_edges.csv</code> · Hilbert, sensor space, sin CSD.
      <em>Strength</em> = media wPLI del nodo en la banda.
      <em>Degree</em> = nº de vecinos en el grafo del modo activo (threshold o top-N).
      Clustering / betweenness sobre el grafo binario del modo.
    </div>
    <div class="toolbar">
      <div class="field">
        <label>Banda</label>
        <select id="band">$band_opts</select>
      </div>
      <div class="field">
        <label>Modo grafo</label>
        <div class="modebox">
          <label><input type="radio" name="mode" value="threshold" checked> Threshold</label>
          <label><input type="radio" name="mode" value="topn"> Top-N</label>
        </div>
      </div>
      <div class="field" id="field-thr">
        <label>Umbral wPLI <span id="thr-val">0.30</span></label>
        <input type="range" id="thr" min="0" max="0.9" step="0.01" value="0.30"/>
      </div>
      <div class="field" id="field-topn">
        <label>Top-N</label>
        <select id="topn" disabled>
          <option value="10">10</option>
          <option value="20" selected>20</option>
          <option value="50">50</option>
        </select>
      </div>
      <div class="field">
        <label>Escala heatmap</label>
        <select id="scale">
          <option value="auto" selected>Auto (0–max)</option>
          <option value="fixed">Fija (0–1)</option>
        </select>
      </div>
      <div class="field">
        <label>Guardar PNG</label>
        <select id="save_view">
          <option value="heatmap">Heatmap</option>
          <option value="topo" selected>Grafo topo</option>
          <option value="hubs">Hubs</option>
        </select>
      </div>
      <button class="btn btn-primary" onclick="refreshAll()">Actualizar</button>
      <button class="btn btn-save" id="btn-save" onclick="savePng()">Guardar PNG</button>
    </div>
    <div class="warn" id="delta-warn" style="display:none">Aviso: DELTA puede ser menos fiable con epochs cortos.</div>
  </div>

  <div class="grid2">
    <div class="card">
      <div class="panel-title">A · Heatmap wPLI</div>
      <canvas id="cv-hm"></canvas>
      <div class="status" id="st-hm">Cargando…</div>
    </div>
    <div class="card">
      <div class="panel-title">B · Grafo topográfico 10-20</div>
      <canvas id="cv-topo"></canvas>
      <div class="status" id="st-topo">Cargando…</div>
    </div>
  </div>

  <div class="card">
    <div class="panel-title">C · Interpretación</div>
    <div class="interp" id="interp"><strong>Interpretación automática</strong><span id="interp-txt">—</span></div>
    <div class="node-card" id="node-card">
      <h4 id="nc-name">—</h4>
      <div class="row"><span class="k" title="Media de wPLI del nodo en la banda">strength</span><span id="nc-str">—</span></div>
      <div class="row"><span class="k" title="Nº de vecinos en el grafo del modo activo">degree</span><span id="nc-deg">—</span></div>
      <div class="row"><span class="k" title="Coeficiente de clustering local (grafo binario)">clustering</span><span id="nc-cl">—</span></div>
      <div class="row"><span class="k" title="Betweenness normalizado (Brandes)">betweenness</span><span id="nc-bt">—</span></div>
    </div>
    <div class="kpis" id="kpis"></div>
    <div class="cols">
      <div>
        <strong style="font-size:12.5px">Top conexiones</strong>
        <div style="overflow:auto;max-height:300px;margin-top:8px">
          <table class="tbl" id="edge-table"><thead></thead><tbody></tbody></table>
        </div>
      </div>
      <div>
        <strong style="font-size:12.5px">Hubs</strong>
        <div style="overflow:auto;max-height:300px;margin-top:8px">
          <table class="tbl" id="hub-table"><thead></thead><tbody></tbody></table>
        </div>
      </div>
      <div>
        <strong style="font-size:12.5px">Distribución wPLI</strong>
        <canvas id="cv-hist" style="margin-top:8px;height:170px"></canvas>
        <strong style="font-size:12.5px;display:block;margin-top:16px">Comparativa entre bandas</strong>
        <p style="font-size:11px;color:var(--muted);margin:4px 0 8px">Density con thr_default del summary (comparación justa).</p>
        <table class="tbl bandtbl" id="band-table"><thead></thead><tbody></tbody></table>
      </div>
    </div>
    <div class="status" id="st-sel" style="margin-top:12px"></div>
  </div>
</div>
<div id="tooltip"></div>

<script>
const SUBJECT = '$SUBJECT / $SESSION / $TASK_ID';
let POS = {}, REGIONS = {};
let matrixData = null, edgesData = null, summaryData = null;
let selEdge = null, selNode = null, flashUntil = 0;

function \$(id) { return document.getElementById(id); }
function setStatus(id, msg, kind) {
  const el = \$(id); el.textContent = msg;
  el.className = 'status' + (kind ? ' ' + kind : '');
}
function modeValue() {
  const r = document.querySelector('input[name=mode]:checked');
  return r ? r.value : 'threshold';
}
function thrValue() { return parseFloat(\$('thr').value); }
function bandValue() { return \$('band').value; }
function topnValue() { return parseInt(\$('topn').value, 10); }
function scaleValue() { return \$('scale').value; }

function updateModeUI() {
  const m = modeValue();
  const thrOn = m === 'threshold';
  \$('thr').disabled = !thrOn;
  \$('topn').disabled = thrOn;
  \$('field-thr').style.opacity = thrOn ? '1' : '0.45';
  \$('field-topn').style.opacity = thrOn ? '0.45' : '1';
  \$('thr-val').textContent = thrValue().toFixed(2);
  \$('delta-warn').style.display = bandValue() === 'DELTA' ? 'block' : 'none';
}

function setupHiDPI(canvas, cssH) {
  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || 560;
  canvas.style.height = cssH + 'px';
  canvas.width = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, W: cssW, H: cssH };
}

function viridis(t) {
  t = Math.max(0, Math.min(1, t));
  const stops = [[68,1,84],[72,35,116],[64,67,135],[52,94,141],[41,120,142],
    [32,144,140],[34,167,132],[68,190,112],[121,209,81],[189,222,38],[253,231,37]];
  const x = t * (stops.length - 1), i = Math.floor(x), f = x - i;
  const a = stops[i], b = stops[Math.min(i+1, stops.length-1)];
  return 'rgb(' + [0,1,2].map(k => Math.round(a[k]+(b[k]-a[k])*f)).join(',') + ')';
}
function wpliStroke(t) {
  t = Math.max(0, Math.min(1, t));
  let r,g,b;
  if (t < 0.5) {
    const u = t/0.5;
    r = 0.55+(0.15-0.55)*u; g = 0.75+(0.35-0.75)*u; b = 0.95+(0.85-0.95)*u;
  } else {
    const u = (t-0.5)/0.5;
    r = 0.15+(0.50-0.15)*u; g = 0.35+(0.15-0.35)*u; b = 0.85+(0.70-0.85)*u;
  }
  return [Math.round(r*255), Math.round(g*255), Math.round(b*255)];
}
function showTip(x,y,text) {
  const el = \$('tooltip');
  el.style.display = 'block'; el.style.left = (x+12)+'px'; el.style.top = (y+12)+'px';
  el.textContent = text;
}
function hideTip() { \$('tooltip').style.display = 'none'; }
function edgeKey(a,b) { return a < b ? a+'|'+b : b+'|'+a; }
function isSelEdge(a,b) { return selEdge && edgeKey(a,b) === edgeKey(selEdge.a, selEdge.b); }

async function refreshAll() {
  updateModeUI();
  const band = bandValue(), thr = thrValue(), topn = topnValue(), mode = modeValue();
  setStatus('st-hm', 'Cargando…'); setStatus('st-topo', 'Cargando…');
  const [mRes, eRes] = await Promise.all([
    fetch('/api/matrix?band=' + encodeURIComponent(band)),
    fetch('/api/edges?band=' + encodeURIComponent(band) + '&mode=' + mode + '&thr=' + thr + '&topn=' + topn)
  ]);
  matrixData = await mRes.json();
  edgesData = await eRes.json();
  if (!mRes.ok) { setStatus('st-hm', matrixData.error || 'Error', 'err'); return; }
  if (!eRes.ok) { setStatus('st-topo', edgesData.error || 'Error', 'err'); return; }
  drawHeatmap(); drawTopo(); fillInterpretation();
  setStatus('st-hm', band + ' · escala ' + scaleValue() + ' · click celda / diagonal = autoconectividad');
  setStatus('st-topo', 'modo ' + mode + ' · ' + (edgesData.graph_edges||[]).length + ' aristas · nodos ∝ strength');
  updateSelectionStatus();
}

function fillInterpretation() {
  const d = edgesData;
  \$('interp-txt').textContent = d.interpretation || '—';
  \$('kpis').innerHTML = [
    ['banda', d.band], ['modo', d.mode], ['nodos', d.n_nodes],
    ['edges', d.n_edges], ['mean', d.mean.toFixed(3)], ['median', d.median.toFixed(3)],
    ['std', d.std.toFixed(3)], ['max', d.max.toFixed(3)],
    ['density', (100*d.density).toFixed(1)+'%'], ['mean str', d.mean_strength.toFixed(3)]
  ].map(([l,v]) => '<div class="kpi"><div class="v">'+v+'</div><div class="l">'+l+'</div></div>').join('');

  const eth = \$('edge-table').querySelector('thead');
  const etb = \$('edge-table').querySelector('tbody');
  eth.innerHTML = '<tr><th>#</th><th>Par</th><th>wPLI</th><th>Z</th></tr>';
  etb.innerHTML = (d.edges||[]).map(e => {
    const sel = isSelEdge(e.ch_a,e.ch_b) ? ' class="sel"' : '';
    return '<tr data-a="'+e.ch_a+'" data-b="'+e.ch_b+'"'+sel+'>'+
      '<td>'+e.rank+'</td><td>'+e.ch_a+'–'+e.ch_b+'</td>'+
      '<td>'+e.wpli.toFixed(4)+'</td><td>'+e.fisher_z.toFixed(3)+'</td></tr>';
  }).join('');
  etb.querySelectorAll('tr').forEach(tr => tr.onclick = () => selectEdge(tr.dataset.a, tr.dataset.b));

  const hth = \$('hub-table').querySelector('thead');
  const htb = \$('hub-table').querySelector('tbody');
  hth.innerHTML = '<tr>'+
    '<th>#</th><th>Canal</th>'+
    '<th title="Media de wPLI del nodo en la banda">Str</th>'+
    '<th title="Nº de vecinos en el grafo del modo activo">Deg</th>'+
    '<th title="Clustering local (grafo binario)">Clust</th>'+
    '<th title="Betweenness normalizado">Betw</th></tr>';
  htb.innerHTML = (d.hubs||[]).slice(0,15).map((h,i) => {
    const sel = selNode === h.channel ? ' class="sel"' : '';
    return '<tr data-ch="'+h.channel+'"'+sel+'>'+
      '<td>'+(i+1)+'</td><td>'+h.channel+'</td>'+
      '<td>'+h.strength.toFixed(3)+'</td><td>'+h.degree+'</td>'+
      '<td>'+h.clustering.toFixed(2)+'</td><td>'+h.betweenness.toFixed(3)+'</td></tr>';
  }).join('');
  htb.querySelectorAll('tr').forEach(tr => tr.onclick = () => selectNode(tr.dataset.ch));

  drawHist(d);
  renderBandTable();
  updateNodeCard();
}

function updateNodeCard() {
  const card = \$('node-card');
  if (!selNode || !edgesData) { card.style.display = 'none'; return; }
  const h = (edgesData.hubs||[]).find(x => x.channel === selNode);
  if (!h) { card.style.display = 'none'; return; }
  card.style.display = 'block';
  \$('nc-name').textContent = h.channel;
  \$('nc-str').textContent = h.strength.toFixed(4);
  \$('nc-deg').textContent = String(h.degree);
  \$('nc-cl').textContent = h.clustering.toFixed(3);
  \$('nc-bt').textContent = h.betweenness.toFixed(4);
}

function drawHist(d) {
  const hist = d.hist; if (!hist) return;
  const canvas = \$('cv-hist');
  const { ctx, W, H } = setupHiDPI(canvas, 170);
  ctx.clearRect(0,0,W,H); ctx.fillStyle='#fff'; ctx.fillRect(0,0,W,H);
  const pad = {l:30,r:10,t:12,b:28};
  const pw = W-pad.l-pad.r, ph = H-pad.t-pad.b;
  const counts = hist.counts||[];
  const maxC = Math.max(1, ...counts);
  const bw = pw / counts.length;
  counts.forEach((c,i) => {
    const h = (c/maxC)*ph;
    ctx.fillStyle = '#93c5fd';
    ctx.fillRect(pad.l+i*bw+1, pad.t+ph-h, Math.max(1,bw-2), h);
  });
  const xAt = v => pad.l + ((v-hist.lo)/(hist.hi-hist.lo||1))*pw;
  function vline(v, color, label) {
    if (!Number.isFinite(v)) return;
    const x = xAt(v);
    ctx.strokeStyle = color; ctx.lineWidth = 1.4; ctx.setLineDash(label==='thr'?[4,3]:[]);
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t+ph); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = color; ctx.font = '10px sans-serif'; ctx.textAlign = 'center';
    ctx.fillText(label, x, pad.t+10);
  }
  vline(d.mean, '#0f766e', 'μ');
  vline(d.median, '#1d4ed8', 'med');
  vline(d.p90, '#7c3aed', 'P90');
  if (d.mode === 'threshold') vline(d.thr, '#dc2626', 'thr');
  ctx.fillStyle='#64748b'; ctx.font='11px sans-serif'; ctx.textAlign='center';
  ctx.fillText('0', pad.l, H-6);
  ctx.fillText(hist.hi.toFixed(2), pad.l+pw, H-6);
}

async function renderBandTable() {
  if (!summaryData) summaryData = await (await fetch('/api/summary')).json();
  const bands = summaryData.bands || {};
  const keys = Object.keys(bands);
  const thead = \$('band-table').querySelector('thead');
  const tbody = \$('band-table').querySelector('tbody');
  thead.innerHTML = '<tr><th>Banda</th><th>Mean</th><th>Density</th><th>Mean str</th></tr>';
  const active = bandValue();
  tbody.innerHTML = keys.map(k => {
    const b = bands[k];
    const cls = k === active ? ' class="active"' : '';
    return '<tr data-band="'+k+'"'+cls+'>'+
      '<td>'+k+'</td><td>'+b.mean.toFixed(3)+'</td>'+
      '<td>'+(100*b.density).toFixed(1)+'%</td>'+
      '<td>'+b.mean_strength.toFixed(3)+'</td></tr>';
  }).join('');
  tbody.querySelectorAll('tr').forEach(tr => {
    tr.onclick = () => { \$('band').value = tr.dataset.band; selEdge=null; selNode=null; refreshAll(); };
  });
}

function drawHeatmap() {
  if (!matrixData) return;
  const canvas = \$('cv-hm');
  const { ctx, W, H } = setupHiDPI(canvas, 480);
  const ch = matrixData.channels, n = matrixData.n, M = matrixData.matrix;
  const scale = scaleValue();
  const vmax = scale === 'fixed' ? 1.0 : (matrixData.vmax || 1);
  const pad = {l:52,r:56,t:28,b:52};
  const pw = W-pad.l-pad.r, ph = H-pad.t-pad.b;
  const cell = Math.min(pw, ph) / n;
  const ox = pad.l + (pw - cell*n)/2;
  const oy = pad.t + (ph - cell*n)/2;
  ctx.clearRect(0,0,W,H); ctx.fillStyle='#fff'; ctx.fillRect(0,0,W,H);
  ctx.fillStyle='#0f172a'; ctx.font='600 13px "IBM Plex Sans",sans-serif'; ctx.textAlign='center';
  ctx.fillText('wPLI '+matrixData.band+' — '+SUBJECT+' · '+scale, W/2, 16);

  for (let j=0;j<n;j++) for (let i=0;i<n;i++) {
    if (i===j) { ctx.fillStyle='#e2e8f0'; }
    else { ctx.fillStyle = viridis(Math.min(1, M[j*n+i]/vmax)); }
    ctx.fillRect(ox+j*cell, oy+i*cell, cell+0.5, cell+0.5);
  }
  const flash = Date.now() < flashUntil;
  if (selEdge) {
    const ia = ch.indexOf(selEdge.a), ib = ch.indexOf(selEdge.b);
    if (ia>=0 && ib>=0) {
      ctx.strokeStyle = flash ? '#fbbf24' : '#e11d48';
      ctx.lineWidth = flash ? 4 : 3;
      ctx.strokeRect(ox+ib*cell, oy+ia*cell, cell, cell);
      ctx.strokeRect(ox+ia*cell, oy+ib*cell, cell, cell);
    }
  }
  if (selNode) {
    const i = ch.indexOf(selNode);
    if (i>=0) {
      ctx.strokeStyle='#0f766e'; ctx.lineWidth=2.5;
      ctx.strokeRect(ox, oy+i*cell, cell*n, cell);
      ctx.strokeRect(ox+i*cell, oy, cell, cell*n);
    }
  }
  const cbX = ox+n*cell+12, cbY = oy, cbH = n*cell, cbW=12;
  for (let y=0;y<cbH;y++) { ctx.fillStyle=viridis(1-y/cbH); ctx.fillRect(cbX, cbY+y, cbW, 1); }
  ctx.strokeStyle='#64748b'; ctx.strokeRect(cbX,cbY,cbW,cbH);
  ctx.fillStyle='#475569'; ctx.font='11px sans-serif'; ctx.textAlign='left';
  ctx.fillText(vmax.toFixed(2), cbX+16, cbY+8);
  ctx.fillText('0', cbX+16, cbY+cbH);
  canvas._hm = {ox,oy,cell,n,ch,M,vmax};
}

function drawTopo() {
  if (!matrixData || !edgesData) return;
  const canvas = \$('cv-topo');
  const { ctx, W, H } = setupHiDPI(canvas, 480);
  const pad = 28;
  const size = Math.min(W,H)-2*pad;
  const cx = W/2, cy = H/2+8, R = size/2;
  ctx.clearRect(0,0,W,H); ctx.fillStyle='#fff'; ctx.fillRect(0,0,W,H);
  ctx.fillStyle='#0f172a'; ctx.font='600 13px "IBM Plex Sans",sans-serif'; ctx.textAlign='center';
  ctx.fillText('Grafo 10-20 · '+matrixData.band+' · '+edgesData.mode, W/2, 16);

  ctx.strokeStyle='#0f172a'; ctx.lineWidth=2;
  ctx.beginPath(); ctx.arc(cx,cy,R,0,Math.PI*2); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(cx-10,cy-R); ctx.quadraticCurveTo(cx,cy-R-18,cx+10,cy-R); ctx.stroke();

  // region labels
  ctx.fillStyle='rgba(100,116,139,0.55)'; ctx.font='11px "IBM Plex Sans",sans-serif';
  ctx.fillText('Frontal', cx, cy-R*0.92);
  ctx.fillText('Central', cx, cy+4);
  ctx.fillText('Parietal', cx, cy+R*0.55);
  ctx.fillText('Occipital', cx, cy+R*0.92);
  ctx.fillText('Temporal', cx-R*0.95, cy);
  ctx.fillText('Temporal', cx+R*0.95, cy);

  const ch = matrixData.channels;
  const strength = matrixData.strength||[];
  const smax = Math.max(1e-9, ...strength);
  const xy = {};
  ch.forEach((c,i) => {
    const p = POS[c]||[0,0];
    xy[c] = {x: cx+p[0]*R*0.92, y: cy-p[1]*R*0.92, i};
  });

  const gEdges = edgesData.graph_edges||[];
  const wmax = Math.max(1e-9, ...gEdges.map(e=>e.wpli));
  gEdges.forEach(e => {
    const A = xy[e.ch_a], B = xy[e.ch_b]; if (!A||!B) return;
    const t = e.wpli/wmax;
    const sel = isSelEdge(e.ch_a,e.ch_b);
    const nodeHit = selNode && (e.ch_a===selNode||e.ch_b===selNode);
    const [r,g,b] = wpliStroke(t);
    const alpha = sel ? 1 : (0.18 + 0.72*t);
    ctx.strokeStyle = sel ? (Date.now()<flashUntil ? '#fbbf24' : '#e11d48')
      : (nodeHit ? '#0f766e' : 'rgba('+r+','+g+','+b+','+alpha+')');
    ctx.lineWidth = sel ? 3.5 : 1.35;
    ctx.beginPath(); ctx.moveTo(A.x,A.y); ctx.lineTo(B.x,B.y); ctx.stroke();
  });

  ch.forEach((c,i) => {
    const p = xy[c]; if (!p) return;
    const rad = 5 + 10*(strength[i]/smax);
    const active = selNode===c || (selEdge && (selEdge.a===c||selEdge.b===c));
    ctx.beginPath(); ctx.arc(p.x,p.y,rad,0,Math.PI*2);
    ctx.fillStyle = active ? '#fef3c7' : '#fff'; ctx.fill();
    ctx.strokeStyle = active ? '#b45309' : '#0f766e';
    ctx.lineWidth = active ? 2.6 : 1.6; ctx.stroke();
    ctx.fillStyle='#0f172a'; ctx.font='10px "IBM Plex Sans",sans-serif'; ctx.textAlign='center';
    ctx.fillText(c, p.x, p.y-rad-3);
  });
  canvas._topo = {xy, gEdges};
}

function selectEdge(a,b) {
  selEdge = {a,b}; selNode = null; flashUntil = Date.now()+450;
  drawHeatmap(); drawTopo(); fillInterpretation(); updateSelectionStatus();
  setTimeout(() => { drawHeatmap(); drawTopo(); }, 500);
}
function selectNode(ch) {
  selNode = ch; selEdge = null;
  drawHeatmap(); drawTopo(); fillInterpretation(); updateSelectionStatus();
}
function updateSelectionStatus() {
  if (selEdge && matrixData) {
    const ia = matrixData.channels.indexOf(selEdge.a);
    const ib = matrixData.channels.indexOf(selEdge.b);
    const w = (ia>=0&&ib>=0) ? matrixData.matrix[ib*matrixData.n+ia] : NaN;
    const ed = (edgesData.edges||[]).find(e => edgeKey(e.ch_a,e.ch_b)===edgeKey(selEdge.a,selEdge.b));
    setStatus('st-sel', 'Arista '+selEdge.a+'–'+selEdge.b+' · wPLI='+(Number.isFinite(w)?w.toFixed(4):'—')+
      (ed?(' · rank='+ed.rank+' · Z='+ed.fisher_z.toFixed(3)):''));
  } else if (selNode) {
    setStatus('st-sel', 'Nodo '+selNode+' — ver tarjeta de métricas arriba.');
  } else {
    setStatus('st-sel', 'Sin selección — click en heatmap, grafo o tablas.');
  }
}

function hitHeatmap(ev) {
  const info = \$('cv-hm')._hm; if (!info) return null;
  const rect = \$('cv-hm').getBoundingClientRect();
  const x = ev.clientX-rect.left, y = ev.clientY-rect.top;
  const j = Math.floor((x-info.ox)/info.cell);
  const i = Math.floor((y-info.oy)/info.cell);
  if (i<0||j<0||i>=info.n||j>=info.n) return null;
  return {a:info.ch[i], b:info.ch[j], v:info.M[j*info.n+i], i, j};
}
\$('cv-hm').addEventListener('mousemove', ev => {
  const h = hitHeatmap(ev); if (!h) { hideTip(); return; }
  if (h.i===h.j) showTip(ev.clientX,ev.clientY, h.a+' · autoconectividad (no interpretada)');
  else showTip(ev.clientX,ev.clientY, h.a+'–'+h.b+' = '+h.v.toFixed(4));
});
\$('cv-hm').addEventListener('mouseleave', hideTip);
\$('cv-hm').addEventListener('click', ev => {
  const h = hitHeatmap(ev); if (!h || h.i===h.j) return; selectEdge(h.a,h.b);
});

\$('cv-topo').addEventListener('mousemove', ev => {
  const info = \$('cv-topo')._topo; if (!info) return;
  const rect = \$('cv-topo').getBoundingClientRect();
  const x = ev.clientX-rect.left, y = ev.clientY-rect.top;
  let best=null, bestD=14;
  for (const [c,p] of Object.entries(info.xy)) {
    const d = Math.hypot(p.x-x,p.y-y); if (d<bestD) {bestD=d; best=c;}
  }
  if (best) showTip(ev.clientX,ev.clientY, best+(REGIONS[best]?' ('+REGIONS[best]+')':''));
  else hideTip();
});
\$('cv-topo').addEventListener('mouseleave', hideTip);
\$('cv-topo').addEventListener('click', ev => {
  const info = \$('cv-topo')._topo; if (!info) return;
  const rect = \$('cv-topo').getBoundingClientRect();
  const x = ev.clientX-rect.left, y = ev.clientY-rect.top;
  let best=null, bestD=14;
  for (const [c,p] of Object.entries(info.xy)) {
    const d = Math.hypot(p.x-x,p.y-y); if (d<bestD) {bestD=d; best=c;}
  }
  if (best) selectNode(best);
});

async function savePng() {
  const btn = \$('btn-save'); btn.disabled = true;
  setStatus('st-sel', 'Generando PNG…');
  try {
    const res = await fetch('/api/save', {
      method:'POST', headers:{'Content-Type':'application/json'},
      body: JSON.stringify({
        view: \$('save_view').value, band: bandValue(), mode: modeValue(),
        thr: thrValue(), topn: topnValue(), scale: scaleValue()
      })
    });
    const data = await res.json();
    if (!res.ok || !data.ok) throw new Error(data.error || 'Error');
    setStatus('st-sel', 'Guardado: '+data.file, 'ok');
  } catch (e) { setStatus('st-sel', String(e.message||e), 'err'); }
  finally { btn.disabled = false; }
}

let _rt = null;
window.addEventListener('resize', () => {
  clearTimeout(_rt);
  _rt = setTimeout(() => { if (matrixData) { drawHeatmap(); drawTopo(); if (edgesData) drawHist(edgesData); } }, 120);
});

document.querySelectorAll('input[name=mode]').forEach(r => r.addEventListener('change', () => { updateModeUI(); refreshAll(); }));
\$('band').addEventListener('change', () => { selEdge=null; selNode=null; refreshAll(); });
\$('thr').addEventListener('input', () => { \$('thr-val').textContent = thrValue().toFixed(2); });
\$('thr').addEventListener('change', refreshAll);
\$('topn').addEventListener('change', refreshAll);
\$('scale').addEventListener('change', () => { drawHeatmap(); });

(async function init() {
  const meta = await (await fetch('/api/meta')).json();
  POS = meta.positions || {};
  REGIONS = meta.regions || {};
  updateModeUI();
  await refreshAll();
})();
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

function _listen_available(host::String, preferred::Int)
    last_err = nothing
    for p in preferred:(preferred + 20)
        try
            return listen(IPv4(host), p), p
        catch e
            last_err = e
            isa(e, Base.IOError) && occursin("EADDRINUSE", sprint(showerror, e)) && continue
            rethrow(e)
        end
    end
    error("No hay puerto libre en $(preferred)-$(preferred + 20): $last_err")
end

function main()
    println("Cargando conectividad desde $CONN …")
    store = load_store()
    println("  $(length(store.channels)) canales · bandas: ", join(store.bands, ", "))
    println("  edges=$(nrow(store.edges)) · thr_default=$(store.thr_default) · epochs=$(store.n_epochs)")

    server, port = _listen_available(HOST, PORT)
    url = "http://$HOST:$port/"
    println()
    println("UI conectividad wPLI v2 → $url")
    port != PORT && println("(puerto $PORT ocupado; usando $port)")
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
