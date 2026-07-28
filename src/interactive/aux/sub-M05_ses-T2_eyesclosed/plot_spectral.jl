#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor interactivo análisis espectral (4 modos)
# ═══════════════════════════════════════════════════════════════
#
#  Modos de salida:
#    1. bands_abs  — Potencia por bandas (µV²)
#    2. psd        — Densidad espectral de potencia absoluta (µV²/Hz)
#    3. bands_rel  — Potencia relativa (% sobre Σ bandas NeuroMIND)
#    4. asd        — Densidad espectral de amplitud √PSD (µV/√Hz)
#
#  ⚠  No es señal cruda: el PSD sale del paso [6/8] del pipeline
#     (filtrado → ICA → segmentación → AR → FFT Hamming), guardado en:
#       tables/psd_by_channel.csv      ← curvas PSD / ASD
#       tables/band_power_summary.csv  ← potencia por banda
#       tables/spectral_indices.csv    ← cocientes α/θ, pico α, …
#       tables/regional_psd.csv        ← media regional (referencia)
#     Los valles a 50/100 Hz son el notch de red del preprocesado.
#
#  Salida PNG:
#    spectral_bands_abs_<canal>.png
#    spectral_psd_<canal>_fmax<N>.png
#    spectral_bands_rel_<canal>.png
#    spectral_asd_<canal>_fmax<N>.png
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_spectral.jl
#    # → http://127.0.0.1:8767/  (Ctrl+C para salir)
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_spectral.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      24-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics

const HERE         = @__DIR__
const RESULTS_UNIT = normpath(joinpath(HERE, "..", "..", "..", "..",
                     "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))
const TABLES = joinpath(RESULTS_UNIT, "tables")
const HOST   = "127.0.0.1"
const PORT   = 8767

const SUBJECT = "sub-M05"
const SESSION = "ses-T2"
const TASK_ID = "EC"
const FS_HZ   = 500.0

const MODES = ["bands_abs", "psd", "bands_rel", "asd"]
const CURVE_MODES = Set(["psd", "asd"])
const BAR_MODES   = Set(["bands_abs", "bands_rel"])

# Bandas NeuroMIND (config/pipeline.toml [bands])
const BAND_ORDER = ["DELTA", "THETA", "ALPHA", "BETA_LOW", "BETA_MID", "BETA_HIGH", "GAMMA"]
const BAND_HZ = Dict(
    "DELTA"     => (0.5,  4.0),
    "THETA"     => (4.0,  8.0),
    "ALPHA"     => (7.8, 11.7),
    "BETA_LOW"  => (12.0, 15.0),
    "BETA_MID"  => (15.0, 18.0),
    "BETA_HIGH" => (18.0, 30.0),
    "GAMMA"     => (30.0, 50.0),
)
const BAND_RGB = Dict(
    "DELTA"     => (0.996, 0.843, 0.667),
    "THETA"     => (0.996, 0.941, 0.541),
    "ALPHA"     => (0.733, 0.969, 0.816),
    "BETA_LOW"  => (0.749, 0.859, 0.996),
    "BETA_MID"  => (0.780, 0.824, 0.996),
    "BETA_HIGH" => (0.867, 0.839, 0.996),
    "GAMMA"     => (0.898, 0.906, 0.922),
)

mutable struct SpectralStore
    channels::Vector{String}
    freqs::Vector{Float64}
    psd::Dict{String,Vector{Float64}}
    mean_psd::Vector{Float64}
    asd::Dict{String,Vector{Float64}}
    mean_asd::Vector{Float64}
    band_power::DataFrame
    indices::DataFrame
    regional::DataFrame
    fs::Float64
end

function load_store()::SpectralStore
    psd_path = joinpath(TABLES, "psd_by_channel.csv")
    bp_path  = joinpath(TABLES, "band_power_summary.csv")
    idx_path = joinpath(TABLES, "spectral_indices.csv")
    reg_path = joinpath(TABLES, "regional_psd.csv")
    isfile(psd_path) || error("No encontrado: $psd_path")
    isfile(bp_path)  || error("No encontrado: $bp_path")

    df = CSV.read(psd_path, DataFrame)
    df.channel = String.(df.channel)
    channels = unique(String.(df.channel))
    ch0 = channels[1]
    sub0 = sort(df[df.channel .== ch0, :], :freq_hz)
    freqs = Float64.(sub0.freq_hz)
    n_f = length(freqs)

    psd = Dict{String,Vector{Float64}}()
    asd = Dict{String,Vector{Float64}}()
    for ch in channels
        sub = sort(df[df.channel .== ch, :], :freq_hz)
        length(sub.freq_hz) == n_f || error("Grid de frecuencias distinto en $ch")
        y = Float64.(sub.power_uv2)
        psd[ch] = y
        asd[ch] = sqrt.(max.(y, 0.0))
    end
    mean_psd = vec(mean(hcat([psd[ch] for ch in channels]...); dims=2))
    mean_asd = sqrt.(max.(mean_psd, 0.0))

    bp = CSV.read(bp_path, DataFrame)
    bp.channel = String.(bp.channel)
    idx = isfile(idx_path) ? CSV.read(idx_path, DataFrame) : DataFrame()
    if nrow(idx) > 0 && hasproperty(idx, :channel)
        idx.channel = String.(idx.channel)
    end
    reg = isfile(reg_path) ? CSV.read(reg_path, DataFrame) : DataFrame()

    return SpectralStore(channels, freqs, psd, mean_psd, asd, mean_asd, bp, idx, reg, FS_HZ)
end

function _band_row(store::SpectralStore, channel::String)::Dict{String,Float64}
    rows = store.band_power[store.band_power.channel .== channel, :]
    isempty(rows) && return Dict{String,Float64}()
    r = rows[1, :]
    out = Dict{String,Float64}()
    for b in BAND_ORDER
        sym = Symbol(b)
        hasproperty(r, sym) || continue
        out[b] = Float64(r[sym])
    end
    return out
end

function _mean_band_power(store::SpectralStore)::Dict{String,Float64}
    out = Dict{String,Float64}()
    for b in BAND_ORDER
        col = Symbol(b)
        hasproperty(store.band_power, col) || continue
        out[b] = mean(Float64.(store.band_power[!, col]))
    end
    return out
end

function _rel_pct(powers::Dict{String,Float64})::Dict{String,Float64}
    tot = sum(values(powers); init=0.0)
    tot <= 0 && return Dict(b => 0.0 for b in keys(powers))
    return Dict(b => 100.0 * v / tot for (b, v) in powers)
end

function _index_row(store::SpectralStore, channel::String)::Dict{String,Float64}
    nrow(store.indices) == 0 && return Dict{String,Float64}()
    rows = store.indices[store.indices.channel .== channel, :]
    isempty(rows) && return Dict{String,Float64}()
    r = rows[1, :]
    out = Dict{String,Float64}()
    for n in names(r)
        n == "channel" && continue
        v = r[n]
        v isa Number && (out[String(n)] = Float64(v))
    end
    return out
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

function _mode_meta(mode::String)
    if mode == "bands_abs"
        return (kind="bars", title="Potencia por bandas", ylabel="Potencia (µV²)", unit="µV²")
    elseif mode == "bands_rel"
        return (kind="bars", title="Potencia relativa", ylabel="Potencia relativa (%)", unit="%")
    elseif mode == "psd"
        return (kind="curve", title="PSD absoluta", ylabel="PSD (µV²/Hz)", unit="µV²/Hz")
    elseif mode == "asd"
        return (kind="curve", title="ASD", ylabel="ASD (µV/√Hz)", unit="µV/√Hz")
    else
        error("Modo desconocido: $mode")
    end
end

function _curve_series(store::SpectralStore, mode::String, channel::String)
    if mode == "psd"
        return store.psd[channel], store.mean_psd
    elseif mode == "asd"
        return store.asd[channel], store.mean_asd
    else
        error("Modo de curva inválido: $mode")
    end
end

function _band_payload(store::SpectralStore, channel::String, mode::String)
    bp = _band_row(store, channel)
    bp_mu = _mean_band_power(store)
    pct = _rel_pct(bp)
    pct_mu = _rel_pct(bp_mu)
    use_rel = mode == "bands_rel"
    bands_j = String[]
    for b in BAND_ORDER
        f1, f2 = BAND_HZ[b]
        r, g, bl = BAND_RGB[b]
        val = use_rel ? get(pct, b, NaN) : get(bp, b, NaN)
        mval = use_rel ? get(pct_mu, b, NaN) : get(bp_mu, b, NaN)
        push!(bands_j,
            "{\"name\":\"$b\",\"f1\":$f1,\"f2\":$f2," *
            "\"value\":$(round(val; digits=6)),\"mean_value\":$(round(mval; digits=6))," *
            "\"color\":\"rgb($(round(Int,255*r)),$(round(Int,255*g)),$(round(Int,255*bl)))\"}")
    end
    return join(bands_j, ","), bp, bp_mu, pct, pct_mu
end

# ── PNG ────────────────────────────────────────────────────────

function _save_curve_png(
    store::SpectralStore,
    mode::String,
    channel::String,
    fmax::Float64,
    yscale::String,
    show_mean::Bool,
)::String
    haskey(store.psd, channel) || error("Canal desconocido: $channel")
    meta = _mode_meta(mode)
    y_ch0, y_mu0 = _curve_series(store, mode, channel)
    mask = (store.freqs .>= 0.5) .& (store.freqs .<= fmax)
    f = store.freqs[mask]
    y_ch = copy(y_ch0[mask])
    y_mu = copy(y_mu0[mask])
    y_ch .= max.(y_ch, 1e-12)
    y_mu .= max.(y_mu, 1e-12)

    fig = Figure(size = (1100, 540), fontsize = 13)
    ylab = if yscale == "log10"
        mode == "asd" ? "log₁₀ ASD (µV/√Hz)" : "log₁₀ PSD (µV²/Hz)"
    else
        meta.ylabel
    end
    ax = Axis(fig[1, 1];
        title = "$(meta.title) — $SUBJECT/$SESSION/$TASK_ID · $channel" *
                (show_mean ? " + promedio ($(length(store.channels)) ch)" : "") *
                " · fs=$(round(store.fs; digits=1)) Hz",
        xlabel = "Frecuencia (Hz)",
        ylabel = ylab,
        titlesize = 14,
        xlabelsize = 13,
        ylabelsize = 13,
        xticklabelsize = 12,
        yticklabelsize = 12,
    )
    ax.xgridvisible = true
    ax.ygridvisible = true

    for b in BAND_ORDER
        f1, f2 = BAND_HZ[b]
        f1 = max(f1, 0.5)
        f2 = min(f2, fmax)
        f1 >= f2 && continue
        r, g, bl = BAND_RGB[b]
        vspan!(ax, f1, f2; color = (RGBf(r, g, bl), 0.40))
    end
    if fmax >= 50
        vlines!(ax, [50.0]; color = (:orange, 0.85), linestyle = :dash, linewidth = 1.4)
    end
    if fmax >= 100
        vlines!(ax, [100.0]; color = (:red, 0.75), linestyle = :dash, linewidth = 1.4)
    end

    y_plot_ch = yscale == "log10" ? log10.(y_ch) : y_ch
    y_plot_mu = yscale == "log10" ? log10.(y_mu) : y_mu
    if show_mean
        lines!(ax, f, y_plot_mu; color = RGBf(0.05, 0.55, 0.55), linewidth = 2.0, label = "Promedio")
    end
    lines!(ax, f, y_plot_ch; color = RGBf(0.15, 0.35, 0.85), linewidth = 1.8, label = channel)
    axislegend(ax; position = :rt, framevisible = true, labelsize = 12)
    xlims!(ax, 0.5, fmax)

    bp = _band_row(store, channel)
    pct = _rel_pct(bp)
    vals = join(["$(b)=$(round(bp[b]; digits=3)) ($(round(pct[b]; digits=1))%)"
                 for b in BAND_ORDER if haskey(bp, b)], "   ")
    Label(fig[2, 1], "Potencia $channel (µV²):  " * vals;
          fontsize = 10, tellwidth = false, halign = :left)

    tag = mode == "asd" ? "asd" : "psd"
    out = joinpath(HERE, "spectral_$(tag)_$(channel)_fmax$(round(Int, fmax)).png")
    save(out, fig; px_per_unit = 3)
    return out
end

function _save_bars_png(
    store::SpectralStore,
    mode::String,
    channel::String,
    show_mean::Bool,
)::String
    haskey(store.psd, channel) || error("Canal desconocido: $channel")
    meta = _mode_meta(mode)
    bp = _band_row(store, channel)
    bp_mu = _mean_band_power(store)
    pct = _rel_pct(bp)
    pct_mu = _rel_pct(bp_mu)
    use_rel = mode == "bands_rel"

    names = [b for b in BAND_ORDER if haskey(bp, b)]
    isempty(names) && error("Sin potencias de banda para $channel")
    vals = [use_rel ? pct[b] : bp[b] for b in names]
    mvals = [use_rel ? get(pct_mu, b, 0.0) : get(bp_mu, b, 0.0) for b in names]
    colors = [RGBf(BAND_RGB[b]...) for b in names]
    xs = 1:length(names)

    fig = Figure(size = (1100, 540), fontsize = 13)
    ax = Axis(fig[1, 1];
        title = "$(meta.title) — $SUBJECT/$SESSION/$TASK_ID · $channel" *
                (show_mean ? " + promedio ($(length(store.channels)) ch)" : "") *
                " · fs=$(round(store.fs; digits=1)) Hz",
        xlabel = "Banda",
        ylabel = meta.ylabel,
        titlesize = 14,
        xticks = (collect(xs), names),
        xticklabelrotation = π / 6,
    )
    ax.xgridvisible = false
    ax.ygridvisible = true

    barplot!(ax, xs, vals; color = colors, width = show_mean ? 0.38 : 0.62,
             dodge_gap = 0.0, label = channel)
    if show_mean
        barplot!(ax, xs .+ 0.42, mvals; color = RGBf(0.05, 0.55, 0.55), width = 0.38,
                 label = "Promedio")
        axislegend(ax; position = :rt, framevisible = true, labelsize = 12)
    end
    if use_rel
        ylims!(ax, 0, max(100.0, maximum(vcat(vals, show_mean ? mvals : Float64[])) * 1.08))
    else
        ymax = maximum(vcat(vals, show_mean ? mvals : Float64[]))
        ylims!(ax, 0, ymax * 1.12)
    end

    foot = if use_rel
        join(["$(b)=$(round(pct[b]; digits=1))%" for b in names], "   ")
    else
        join(["$(b)=$(round(bp[b]; digits=3))" for b in names], "   ")
    end
    Label(fig[2, 1], "$(meta.title) $channel ($(meta.unit)):  " * foot;
          fontsize = 10, tellwidth = false, halign = :left)

    tag = use_rel ? "bands_rel" : "bands_abs"
    out = joinpath(HERE, "spectral_$(tag)_$(channel).png")
    save(out, fig; px_per_unit = 3)
    return out
end

function save_spectral_png(
    store::SpectralStore,
    mode::String,
    channel::String,
    fmax::Float64,
    yscale::String,
    show_mean::Bool,
)::String
    mode in MODES || error("Modo desconocido: $mode")
    if mode in CURVE_MODES
        return _save_curve_png(store, mode, channel, fmax, yscale, show_mean)
    else
        return _save_bars_png(store, mode, channel, show_mean)
    end
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

function _parse_json(s::String)
    _str = function (key, default)
        m = match(Regex("\"" * key * "\"\\s*:\\s*\"([^\"]*)\""), s)
        return m === nothing ? String(default) : String(m.captures[1])
    end
    _num = function (key, default)
        m = match(Regex("\"" * key * "\"\\s*:\\s*([-\\d.eE]+)"), s)
        m === nothing && return Float64(default)
        v = tryparse(Float64, m.captures[1])
        return v === nothing ? Float64(default) : v
    end
    _bool = function (key, default)
        m = match(Regex("\"" * key * "\"\\s*:\\s*(true|false)"), s)
        m === nothing && return default
        return m.captures[1] == "true"
    end
    return (
        mode      = _str("mode", "psd"),
        channel   = _str("channel", "Cz"),
        fmax      = _num("fmax", 150.0),
        yscale    = _str("yscale", "log10"),
        show_mean = _bool("show_mean", true),
    )
end

function _json_escape(s::AbstractString)
    replace(replace(s, "\\" => "\\\\"), "\"" => "\\\"")
end

function _spectrum_json(store::SpectralStore, mode::String, channel::String, fmax::Float64)::String
    mode in MODES || error("Modo desconocido: $mode")
    haskey(store.psd, channel) || error("Canal desconocido: $channel")
    meta = _mode_meta(mode)
    bands_arr, bp, bp_mu, pct, pct_mu = _band_payload(store, channel, mode)
    idx = _index_row(store, channel)

    bands_abs = join(["\"$b\":$(round(get(bp, b, NaN); digits=6))" for b in BAND_ORDER], ",")
    bands_pct = join(["\"$b\":$(round(get(pct, b, NaN); digits=2))" for b in BAND_ORDER], ",")
    bands_mu  = join(["\"$b\":$(round(get(bp_mu, b, NaN); digits=6))" for b in BAND_ORDER], ",")
    bands_mup = join(["\"$b\":$(round(get(pct_mu, b, NaN); digits=2))" for b in BAND_ORDER], ",")
    idx_j = join(["\"$k\":$(round(v; digits=4))" for (k, v) in idx], ",")

    common = "\"mode\":\"$mode\",\"kind\":\"$(meta.kind)\"," *
             "\"title\":\"$(_json_escape(meta.title))\"," *
             "\"ylabel\":\"$(_json_escape(meta.ylabel))\"," *
             "\"unit\":\"$(_json_escape(meta.unit))\"," *
             "\"channel\":\"$channel\"," *
             "\"bands\":[$bands_arr]," *
             "\"band_power\":{$bands_abs}," *
             "\"band_pct\":{$bands_pct}," *
             "\"band_power_mean\":{$bands_mu}," *
             "\"band_pct_mean\":{$bands_mup}," *
             "\"indices\":{$idx_j}"

    if mode in CURVE_MODES
        mask = findall(i -> store.freqs[i] >= 0.5 && store.freqs[i] <= fmax, eachindex(store.freqs))
        if length(mask) > 800
            step = ceil(Int, length(mask) / 800)
            mask = mask[1:step:end]
        end
        f = store.freqs[mask]
        ych, ymu = _curve_series(store, mode, channel)
        ych = ych[mask]
        ymu = ymu[mask]
        return "{" * common * "," *
               "\"f\":[" * join(round.(f; digits=4), ",") * "]," *
               "\"y\":[" * join(round.(ych; digits=8), ",") * "]," *
               "\"y_mean\":[" * join(round.(ymu; digits=8), ",") * "]}"
    else
        return "{" * common * "}"
    end
end

function handle_request(sock, store::SpectralStore)
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
            bands_j = join(["{\"name\":\"$b\",\"f1\":$(BAND_HZ[b][1]),\"f2\":$(BAND_HZ[b][2])," *
                            "\"color\":\"rgb($(round(Int,255*BAND_RGB[b][1])),$(round(Int,255*BAND_RGB[b][2])),$(round(Int,255*BAND_RGB[b][3])))\"}"
                            for b in BAND_ORDER], ",")
            modes_j = join(["\"$m\"" for m in MODES], ",")
            body = "{\"channels\":[" * join(["\"$c\"" for c in store.channels], ",") * "]," *
                   "\"fs\":$(store.fs),\"f_max\":$(maximum(store.freqs))," *
                   "\"n_channels\":$(length(store.channels))," *
                   "\"modes\":[$modes_j],\"bands\":[$bands_j]}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/spectrum"
            qs = _parse_qs(path)
            mode = String(get(qs, "mode", "psd"))
            ch = String(get(qs, "channel", "Cz"))
            fmax = something(tryparse(Float64, get(qs, "fmax", "150")), 150.0)
            body = _spectrum_json(store, mode, ch, fmax)
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            req = _parse_json(raw)
            out = save_spectral_png(store, req.mode, req.channel, req.fmax, req.yscale, req.show_mean)
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

function html_page(store::SpectralStore)::String
    ch_opts = join(["<option value=\"$c\"" * (c == "Cz" ? " selected" : "") * ">$c</option>"
                    for c in store.channels], "\n")
    n_ch = length(store.channels)
    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND — Análisis espectral (4 modos)</title>
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
    border-radius:6px; background:#fff; font-size:13px; min-width:110px;
  }
  select:disabled { opacity:.45; background:#f1f5f9; }
  .chk { display:flex; align-items:center; gap:8px; height:34px; font-size:13px; }
  .btn { height:34px; padding:0 16px; border:none; border-radius:6px; font-size:13px; font-weight:600; cursor:pointer; }
  .btn-primary { background:var(--accent); color:#fff; }
  .btn-primary:hover { background:#1d4ed8; }
  .btn-save { background:var(--accent2); color:#fff; }
  .btn-save:hover { background:#0d9488; }
  .btn:disabled { opacity:.5; cursor:wait; }
  canvas { width:100%; display:block; background:#fff; border:1px solid var(--border); border-radius:8px; }
  .status { font-size:12px; color:var(--muted); margin-top:8px; min-height:1.2em; }
  .status.ok { color:#047857; } .status.err { color:#b91c1c; }
  table.bands { width:100%; border-collapse:collapse; font-size:12.5px; }
  table.bands th, table.bands td { padding:8px 10px; border-bottom:1px solid var(--border); text-align:right; }
  table.bands th:first-child, table.bands td:first-child { text-align:left; }
  table.bands th { background:#f8fafc; color:#475569; font-weight:600; font-size:11px; text-transform:uppercase; }
  table.bands tr:last-child td { border-bottom:none; }
  table.bands.emphasize-pct td span.pct { color:#0f766e; font-weight:600; }
  .src {
    font-size:12.5px; line-height:1.45; color:#334155;
    background:#f0f9ff; border:1px solid #bae6fd; border-radius:8px;
    padding:10px 12px; margin-bottom:12px;
  }
  .src code { font-size:12px; background:#e0f2fe; padding:1px 5px; border-radius:4px; }
  .band-legend { display:flex; flex-wrap:wrap; gap:8px; margin:0 0 12px; }
  .band-chip {
    display:inline-flex; align-items:center; gap:6px;
    font-size:12px; font-weight:600; color:#1e293b;
    padding:4px 10px; border-radius:999px; border:1px solid var(--border);
    background:#fff;
  }
  .band-chip i {
    width:12px; height:12px; border-radius:3px; display:inline-block;
    border:1px solid rgba(0,0,0,.12);
  }
  .idx-wrap { margin-top:14px; }
  .idx-wrap h3 { margin:0 0 6px; font-size:13px; }
  .idx-wrap p.hint { margin:0 0 10px; font-size:12px; color:var(--muted); }
  table.idx { width:100%; border-collapse:collapse; font-size:12.5px; }
  table.idx th, table.idx td {
    padding:8px 10px; border-bottom:1px solid var(--border); text-align:left; vertical-align:top;
  }
  table.idx th { background:#f8fafc; color:#475569; font-size:11px; text-transform:uppercase; }
  table.idx td.val { text-align:right; font-variant-numeric:tabular-nums; font-weight:600; white-space:nowrap; }
  table.idx td.key { color:#64748b; font-family:ui-monospace,monospace; font-size:11.5px; }
</style>
</head>
<body>
<header>
  <h1>Análisis espectral — 4 modos</h1>
  <span>$SUBJECT / $SESSION / task-$TASK_ID · fs=$(round(Int, store.fs)) Hz · $n_ch canales</span>
</header>
<div class="wrap">
  <div class="card">
    <div class="src">
      <strong>Fuente de datos (no es señal cruda).</strong>
      Curvas: <code>tables/psd_by_channel.csv</code> (PSD) → ASD = √PSD.
      Potencias: <code>band_power_summary.csv</code> · relativa = banda / Σ bandas NeuroMIND × 100.
      Índices: <code>spectral_indices.csv</code>.
      Paso [6/8]: filtrado → ICA → segmentación → AR → FFT Hamming.
    </div>
    <div class="toolbar">
      <div class="field">
        <label>Modo</label>
        <select id="mode">
          <option value="bands_abs">1. Potencia por bandas</option>
          <option value="psd" selected>2. PSD absoluta</option>
          <option value="bands_rel">3. Potencia relativa (%)</option>
          <option value="asd">4. ASD (√PSD)</option>
        </select>
      </div>
      <div class="field">
        <label>Canal</label>
        <select id="channel">$ch_opts</select>
      </div>
      <div class="field" id="field-fmax">
        <label>Frecuencia máx.</label>
        <select id="fmax">
          <option value="50">50 Hz</option>
          <option value="100">100 Hz</option>
          <option value="150" selected>150 Hz</option>
          <option value="250">250 Hz (Nyquist)</option>
        </select>
      </div>
      <div class="field" id="field-yscale">
        <label>Eje Y</label>
        <select id="yscale">
          <option value="log10" selected>log₁₀</option>
          <option value="linear">lineal</option>
        </select>
      </div>
      <label class="chk"><input type="checkbox" id="show_mean" checked> Mostrar promedio</label>
      <button class="btn btn-primary" onclick="refreshPlot()">Actualizar</button>
      <button class="btn btn-save" id="btn-save" onclick="savePng()">Guardar PNG</button>
    </div>
  </div>

  <div class="card">
    <div class="band-legend" id="band-legend"></div>
    <canvas id="cv"></canvas>
    <div class="status" id="status">Cargando…</div>
  </div>

  <div class="card">
    <strong style="font-size:13px" id="band-card-title">Potencia por bandas (µV²)</strong>
    <p style="margin:4px 0 0;font-size:12px;color:var(--muted)" id="band-card-hint">
      Valores de <code>band_power_summary.csv</code> · % relativo a la suma de bandas NeuroMIND
    </p>
    <div style="overflow-x:auto;margin-top:10px">
      <table class="bands" id="band-table">
        <thead></thead><tbody></tbody>
      </table>
    </div>
    <div class="idx-wrap">
      <h3>Índices espectrales del canal</h3>
      <p class="hint">Fuente: <code>spectral_indices.csv</code> — cocientes de potencia entre bandas y pico α del canal seleccionado.</p>
      <table class="idx" id="idx-table">
        <thead>
          <tr><th>Índice</th><th>Significado</th><th>Valor</th></tr>
        </thead>
        <tbody></tbody>
      </table>
    </div>
  </div>
</div>

<script>
const BANDS = $(
    "[" * join(["{name:\"$b\",f1:$(BAND_HZ[b][1]),f2:$(BAND_HZ[b][2])," *
                "rgb:[$(BAND_RGB[b][1]),$(BAND_RGB[b][2]),$(BAND_RGB[b][3])]}"
                for b in BAND_ORDER], ",") * "]"
);
const FS = $(store.fs);
const N_CH = $(n_ch);
const INDEX_INFO = {
  alpha_theta:   { label: 'α / θ',           meaning: 'Potencia Alpha dividida por Theta (alerta vs somnolencia)' },
  beta_alpha:    { label: 'β / α',           meaning: 'Potencia Beta (media low/mid/high) dividida por Alpha' },
  theta_beta:    { label: 'θ / β',           meaning: 'Potencia Theta dividida por Beta (índice atencional clásico)' },
  gamma_alpha:   { label: 'γ / α',           meaning: 'Potencia Gamma dividida por Alpha' },
  peak_alpha_hz: { label: 'Pico α (Hz)',     meaning: 'Frecuencia del máximo de PSD dentro de la banda Alpha' },
  peak_alpha_uv2:{ label: 'Pico α (µV²/Hz)', meaning: 'Amplitud PSD en ese pico Alpha' },
};
let _last = null;

function isCurveMode(mode) { return mode === 'psd' || mode === 'asd'; }

function bandColorCSS(b) {
  const [r,g,bl] = b.rgb;
  return 'rgb(' + Math.round(r*255) + ',' + Math.round(g*255) + ',' + Math.round(bl*255) + ')';
}

function renderBandLegend() {
  document.getElementById('band-legend').innerHTML = BANDS.map(b =>
    '<span class="band-chip"><i style="background:' + bandColorCSS(b) + '"></i>' +
    b.name + ' <span style="font-weight:400;color:#64748b">(' + b.f1 + '–' + b.f2 + ' Hz)</span></span>'
  ).join('');
}

function setStatus(msg, kind) {
  const el = document.getElementById('status');
  el.textContent = msg;
  el.className = 'status' + (kind ? ' ' + kind : '');
}

function updateControls() {
  const mode = document.getElementById('mode').value;
  const curve = isCurveMode(mode);
  document.getElementById('fmax').disabled = !curve;
  document.getElementById('yscale').disabled = !curve;
  document.getElementById('field-fmax').style.opacity = curve ? '1' : '0.45';
  document.getElementById('field-yscale').style.opacity = curve ? '1' : '0.45';
  const title = document.getElementById('band-card-title');
  const hint = document.getElementById('band-card-hint');
  const tbl = document.getElementById('band-table');
  if (mode === 'bands_rel') {
    title.textContent = 'Potencia relativa (%)';
    hint.innerHTML = '% = banda / Σ bandas NeuroMIND × 100 · fuente <code>band_power_summary.csv</code>';
    tbl.classList.add('emphasize-pct');
  } else {
    title.textContent = 'Potencia por bandas (µV²)';
    hint.innerHTML = 'Valores de <code>band_power_summary.csv</code> · % relativo a la suma de bandas NeuroMIND';
    tbl.classList.toggle('emphasize-pct', false);
  }
}

function setupHiDPICanvas(canvas, cssH) {
  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.clientWidth || 1100;
  canvas.style.height = cssH + 'px';
  canvas.width  = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, W: cssW, H: cssH };
}

async function refreshPlot() {
  updateControls();
  const mode = document.getElementById('mode').value;
  const channel = document.getElementById('channel').value;
  const fmax = parseFloat(document.getElementById('fmax').value);
  const yscale = document.getElementById('yscale').value;
  const showMean = document.getElementById('show_mean').checked;
  setStatus('Cargando…');
  const url = '/api/spectrum?mode=' + encodeURIComponent(mode) +
              '&channel=' + encodeURIComponent(channel) + '&fmax=' + fmax;
  const res = await fetch(url);
  const data = await res.json();
  if (!res.ok) { setStatus(data.error || 'Error', 'err'); return; }
  _last = { data, mode, fmax, yscale, showMean };
  drawSpectrum(data, mode, fmax, yscale, showMean);
  fillTable(data, mode);
  const src = isCurveMode(mode)
    ? (mode === 'asd' ? 'ASD = √PSD desde psd_by_channel.csv' : 'psd_by_channel.csv (post-ICA/AR)')
    : 'band_power_summary.csv';
  setStatus(data.title + ' · ' + channel + (isCurveMode(mode) ? (' · 0.5–' + fmax + ' Hz') : '') + ' · ' + src);
}

function fillTable(data, mode) {
  const thead = document.querySelector('#band-table thead');
  const tbody = document.querySelector('#band-table tbody');
  const names = BANDS.map(b => b.name);
  thead.innerHTML = '<tr><th></th>' + names.map(n => {
    const b = BANDS.find(x => x.name === n);
    return '<th style="background:' + bandColorCSS(b) + ';color:#1e293b">' + n +
           '<br><span style="font-weight:400;text-transform:none">(' +
           b.f1 + '–' + b.f2 + ')</span></th>';
  }).join('') + '</tr>';

  const emphasizeRel = mode === 'bands_rel';
  const row = (label, abs, pct) => '<tr><td><strong>' + label + '</strong></td>' +
    names.map(n => {
      const a = abs[n], p = pct[n];
      if (a == null || Number.isNaN(a)) return '<td>—</td>';
      if (emphasizeRel) {
        return '<td><span class="pct">' + p.toFixed(1) + '%</span><br>' +
               '<span style="color:#64748b">' + a.toFixed(3) + ' µV²</span></td>';
      }
      return '<td>' + a.toFixed(3) + '<br><span class="pct" style="color:#64748b">' +
             p.toFixed(1) + '%</span></td>';
    }).join('') + '</tr>';

  tbody.innerHTML =
    row(data.channel + (emphasizeRel ? ' (%)' : ' (µV²)'), data.band_power, data.band_pct) +
    row('Promedio' + (emphasizeRel ? ' (%)' : ' (µV²)'), data.band_power_mean, data.band_pct_mean);

  const idx = data.indices || {};
  const order = ['alpha_theta','beta_alpha','theta_beta','gamma_alpha','peak_alpha_hz','peak_alpha_uv2'];
  const keys = order.filter(k => idx[k] != null).concat(Object.keys(idx).filter(k => !order.includes(k)));
  const tb = document.querySelector('#idx-table tbody');
  tb.innerHTML = keys.map(k => {
    const info = INDEX_INFO[k] || { label: k, meaning: 'Índice de spectral_indices.csv' };
    const v = idx[k];
    const vs = (typeof v === 'number') ? (Math.abs(v) >= 100 ? v.toFixed(2) : v.toFixed(3)) : String(v);
    return '<tr><td class="key">' + k + '</td><td><strong>' + info.label + '</strong> — ' +
           info.meaning + '</td><td class="val">' + vs + '</td></tr>';
  }).join('');
}

function drawSpectrum(data, mode, fmax, yscale, showMean) {
  if (isCurveMode(mode)) drawCurve(data, mode, fmax, yscale, showMean);
  else drawBars(data, mode, showMean);
}

function drawBars(data, mode, showMean) {
  const canvas = document.getElementById('cv');
  const { ctx, W, H } = setupHiDPICanvas(canvas, 460);
  const pad = { l: 68, r: 20, t: 44, b: 64 };
  const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);

  const bands = data.bands || [];
  const n = bands.length;
  if (!n) {
    setStatus('Sin datos de bandas', 'err');
    return;
  }
  const vals = bands.map(b => b.value);
  const mvals = bands.map(b => b.mean_value);
  let ymax = Math.max(...vals, ...(showMean ? mvals : [0]), 1e-9);
  if (mode === 'bands_rel') ymax = Math.max(100, ymax);
  ymax *= 1.12;
  const ymin = 0;
  const yAt = v => pad.t + (1 - (v - ymin) / (ymax - ymin)) * ph;
  const groupW = pw / n;
  const barW = showMean ? groupW * 0.36 : groupW * 0.55;

  // grid Y
  const nY = 5;
  for (let i = 0; i <= nY; i++) {
    const v = ymin + (i / nY) * (ymax - ymin);
    const y = yAt(v);
    ctx.strokeStyle = '#e2e8f0';
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + pw, y); ctx.stroke();
    ctx.fillStyle = '#475569'; ctx.font = '12px "IBM Plex Sans",sans-serif';
    ctx.textAlign = 'right';
    ctx.fillText(mode === 'bands_rel' ? v.toFixed(0) : v.toExponential(1), pad.l - 8, y + 4);
  }
  ctx.strokeStyle = '#64748b'; ctx.lineWidth = 1.2;
  ctx.strokeRect(pad.l, pad.t, pw, ph);

  for (let i = 0; i < n; i++) {
    const b = bands[i];
    const cx = pad.l + (i + 0.5) * groupW;
    const h = Math.max(1, (b.value / ymax) * ph);
    const x0 = showMean ? cx - barW - 2 : cx - barW / 2;
    ctx.fillStyle = b.color || '#93c5fd';
    ctx.fillRect(x0, pad.t + ph - h, barW, h);
    if (showMean) {
      const hm = Math.max(1, (b.mean_value / ymax) * ph);
      ctx.fillStyle = '#0d9488';
      ctx.fillRect(cx + 2, pad.t + ph - hm, barW, hm);
    }
    ctx.fillStyle = '#334155';
    ctx.font = '11px "IBM Plex Sans",sans-serif';
    ctx.textAlign = 'center';
    const label = mode === 'bands_rel' ? b.value.toFixed(1) + '%' : b.value.toExponential(1);
    ctx.fillText(label, showMean ? cx - barW / 2 : cx, pad.t + ph - h - 6);
    ctx.fillStyle = '#475569';
    ctx.font = '11px "IBM Plex Sans",sans-serif';
    ctx.fillText(b.name, cx, H - 28);
    ctx.fillStyle = '#94a3b8';
    ctx.font = '10px "IBM Plex Sans",sans-serif';
    ctx.fillText(b.f1 + '–' + b.f2, cx, H - 12);
  }

  // leyenda
  let lx = pad.l + 10, ly = pad.t + 16;
  ctx.font = '12px "IBM Plex Sans",sans-serif';
  ctx.fillStyle = '#2563eb';
  ctx.fillRect(lx, ly - 6, 14, 10);
  ctx.fillStyle = '#0f172a'; ctx.textAlign = 'left';
  ctx.fillText(data.channel, lx + 20, ly + 2);
  if (showMean) {
    lx += 90;
    ctx.fillStyle = '#0d9488';
    ctx.fillRect(lx, ly - 6, 14, 10);
    ctx.fillStyle = '#0f172a';
    ctx.fillText('Promedio (' + N_CH + ' ch)', lx + 20, ly + 2);
  }

  ctx.fillStyle = '#0f172a';
  ctx.font = '600 14px "IBM Plex Sans",sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText(data.title + ' — ' + SUBJECT_LABEL() + ' · fs=' + FS.toFixed(1) + ' Hz',
               pad.l + pw / 2, 18);
  ctx.save();
  ctx.translate(16, pad.t + ph / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = '#334155';
  ctx.font = '13px "IBM Plex Sans",sans-serif';
  ctx.fillText(data.ylabel, 0, 0);
  ctx.restore();
}

function SUBJECT_LABEL() { return 'sub-M05/ses-T2/EC'; }

function drawCurve(data, mode, fmax, yscale, showMean) {
  const canvas = document.getElementById('cv');
  const { ctx, W, H } = setupHiDPICanvas(canvas, 460);
  const pad = { l: 68, r: 20, t: 40, b: 48 };
  const pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;

  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);

  const f = data.f || [];
  let yCh = (data.y || []).map(v => Math.max(v, 1e-12));
  let yMu = (data.y_mean || []).map(v => Math.max(v, 1e-12));
  if (yscale === 'log10') {
    yCh = yCh.map(Math.log10);
    yMu = yMu.map(Math.log10);
  }

  let ymin = Infinity, ymax = -Infinity;
  const consider = showMean ? yCh.concat(yMu) : yCh;
  for (const v of consider) { if (v < ymin) ymin = v; if (v > ymax) ymax = v; }
  const ypad = 0.06 * Math.max(ymax - ymin, 1e-6);
  ymin -= ypad; ymax += ypad;

  const x0 = 0.5, x1 = fmax;
  const xAt = u => pad.l + ((u - x0) / (x1 - x0)) * pw;
  const yAt = v => pad.t + (1 - (v - ymin) / (ymax - ymin)) * ph;

  BANDS.forEach(b => {
    const f1 = Math.max(b.f1, x0), f2 = Math.min(b.f2, x1);
    if (f1 >= f2) return;
    const [r,g,bl] = b.rgb;
    ctx.fillStyle = 'rgba(' + Math.round(r*255) + ',' + Math.round(g*255) + ',' + Math.round(bl*255) + ',0.40)';
    ctx.fillRect(xAt(f1), pad.t, xAt(f2) - xAt(f1), ph);
  });

  ctx.strokeStyle = '#e2e8f0'; ctx.lineWidth = 1;
  const xTicks = [1, 10, 50, 100, 150, 200, 250].filter(u => u >= x0 && u <= x1);
  if (!xTicks.includes(x0)) xTicks.unshift(x0);
  if (!xTicks.includes(x1)) xTicks.push(x1);
  xTicks.forEach(u => {
    const x = xAt(u);
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
    ctx.fillStyle = '#475569'; ctx.font = '12px "IBM Plex Sans",sans-serif';
    ctx.textAlign = 'center'; ctx.fillText(String(u), x, H - 18);
  });
  const nY = 5;
  for (let i = 0; i <= nY; i++) {
    const v = ymin + (i / nY) * (ymax - ymin);
    const y = yAt(v);
    ctx.strokeStyle = '#e2e8f0';
    ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + pw, y); ctx.stroke();
    ctx.fillStyle = '#475569'; ctx.font = '12px "IBM Plex Sans",sans-serif';
    ctx.textAlign = 'right';
    ctx.fillText(yscale === 'log10' ? v.toFixed(1) : v.toExponential(1), pad.l - 8, y + 4);
  }

  function markLine(freq, color, label) {
    if (freq < x0 || freq > x1) return;
    const x = xAt(freq);
    ctx.strokeStyle = color; ctx.setLineDash([5, 4]); ctx.lineWidth = 1.4;
    ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + ph); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = color; ctx.font = '11px "IBM Plex Sans",sans-serif';
    ctx.textAlign = 'left'; ctx.fillText(label, x + 4, pad.t + ph - 8);
  }
  markLine(50, '#ea580c', '50 Hz (notch)');
  markLine(100, '#dc2626', '100 Hz');

  ctx.strokeStyle = '#64748b'; ctx.lineWidth = 1.2;
  ctx.strokeRect(pad.l, pad.t, pw, ph);

  ctx.save();
  ctx.beginPath(); ctx.rect(pad.l, pad.t, pw, ph); ctx.clip();
  function strokeSeries(yy, color, lw) {
    ctx.strokeStyle = color; ctx.lineWidth = lw;
    ctx.beginPath();
    for (let i = 0; i < f.length; i++) {
      const x = xAt(f[i]), y = yAt(yy[i]);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();
  }
  if (showMean) strokeSeries(yMu, '#0d9488', 2.0);
  strokeSeries(yCh, '#2563eb', 1.8);
  ctx.restore();

  let lx = pad.l + 10, ly = pad.t + 18;
  ctx.font = '12px "IBM Plex Sans",sans-serif';
  if (showMean) {
    ctx.strokeStyle = '#0d9488'; ctx.lineWidth = 2; ctx.beginPath();
    ctx.moveTo(lx, ly); ctx.lineTo(lx + 18, ly); ctx.stroke();
    ctx.fillStyle = '#0f172a'; ctx.textAlign = 'left';
    ctx.fillText('Promedio (' + N_CH + ' ch)', lx + 24, ly + 4);
    lx += 150;
  }
  ctx.strokeStyle = '#2563eb'; ctx.lineWidth = 2; ctx.beginPath();
  ctx.moveTo(lx, ly); ctx.lineTo(lx + 18, ly); ctx.stroke();
  ctx.fillStyle = '#0f172a'; ctx.textAlign = 'left';
  ctx.fillText(data.channel, lx + 24, ly + 4);

  ctx.fillStyle = '#0f172a';
  ctx.font = '600 14px "IBM Plex Sans",sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText(data.title + ' — ' + SUBJECT_LABEL() + ' · fs=' + FS.toFixed(1) + ' Hz',
               pad.l + pw / 2, 18);
  ctx.fillStyle = '#334155';
  ctx.font = '13px "IBM Plex Sans",sans-serif';
  ctx.fillText('Frecuencia (Hz)', pad.l + pw / 2, H - 4);
  ctx.save();
  ctx.translate(16, pad.t + ph / 2); ctx.rotate(-Math.PI / 2);
  const ylab = yscale === 'log10'
    ? (mode === 'asd' ? 'log₁₀ ASD (µV/√Hz)' : 'log₁₀ PSD (µV²/Hz)')
    : data.ylabel;
  ctx.fillText(ylab, 0, 0);
  ctx.restore();
}

async function savePng() {
  const mode = document.getElementById('mode').value;
  const channel = document.getElementById('channel').value;
  const fmax = parseFloat(document.getElementById('fmax').value);
  const yscale = document.getElementById('yscale').value;
  const show_mean = document.getElementById('show_mean').checked;
  const btn = document.getElementById('btn-save');
  btn.disabled = true;
  setStatus('Generando PNG…');
  try {
    const res = await fetch('/api/save', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({mode, channel, fmax, yscale, show_mean})
    });
    const data = await res.json();
    if (!res.ok || !data.ok) throw new Error(data.error || 'Error');
    setStatus('Guardado: ' + data.file, 'ok');
  } catch (e) {
    setStatus(String(e.message || e), 'err');
  } finally {
    btn.disabled = false;
  }
}

let _rt = null;
window.addEventListener('resize', () => {
  clearTimeout(_rt);
  _rt = setTimeout(() => {
    if (_last) drawSpectrum(_last.data, _last.mode, _last.fmax, _last.yscale, _last.showMean);
  }, 120);
});

document.getElementById('mode').addEventListener('change', () => { updateControls(); refreshPlot(); });
document.getElementById('channel').addEventListener('change', refreshPlot);
document.getElementById('fmax').addEventListener('change', refreshPlot);
document.getElementById('yscale').addEventListener('change', refreshPlot);
document.getElementById('show_mean').addEventListener('change', refreshPlot);

renderBandLegend();
updateControls();
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
    println("Cargando tablas espectrales desde $TABLES …")
    store = load_store()
    println("  $(length(store.channels)) canales · $(length(store.freqs)) freqs · ",
            "0–$(round(maximum(store.freqs); digits=1)) Hz")
    println("  Modos: ", join(MODES, ", "))

    server = listen(IPv4(HOST), PORT)
    url = "http://$HOST:$PORT/"
    println()
    println("UI análisis espectral → $url")
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
