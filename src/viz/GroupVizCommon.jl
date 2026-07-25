# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — GroupVizCommon
#  Helpers compartidos para visores y figuras de cohorte
#  (longitudinal T1→T2 · transversal MS vs Control)
# ═══════════════════════════════════════════════════════════════
#
#  · IO CSV/JSON · posiciones 10-20 · orden fisiológico de bandas
#  · Heatmaps con escala emparejada + RdBu_r (positivo = rojo)
#  · Redes topo · paired/group means · topo Δ band-power
#
#  Uso (scripts / plotters standalone):
#    include(joinpath(@__DIR__, "..", "viz", "GroupVizCommon.jl"))
#    using .GroupVizCommon
#
# ───────────────────────────────────────────────────────────────

module GroupVizCommon

using CSV, DataFrames, CairoMakie, Statistics, TOML

export BAND_ORDER, CH_POS, DIVERGING_CMAP
export parse_summary_json, read_mat_csv, safe_csv, resolve_results_root
export json_escape, mat_flat, df_rows_json, pos_json
export sort_bands_physio, honest_best_band, count_fdr, count_uncorr
export save_heatmap, save_heatmap_triplet, save_topo_network
export save_paired_means, save_group_means, save_topo_delta

const BAND_ORDER = ["DELTA", "THETA", "ALPHA", "BETA_LOW", "BETA_MID", "BETA_HIGH", "GAMMA"]

# Positivo (MS−Ctrl / T2−T1) = rojo; negativo = azul. Coincide con el JS del visor.
# Makie no expone :RdBu_r; Reverse(:RdBu) es el equivalente (ColorBrewer RdBu invertido).
const DIVERGING_CMAP = Reverse(:RdBu)

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
    "FP1"=>(-0.25,0.88),"FP2"=>(0.25,0.88),
)

# ── IO ─────────────────────────────────────────────────────────

function parse_summary_json(path::String)::Dict{String,Any}
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
    for m in eachmatch(r"\"([^\"]+)\"\s*:\s*(true|false)", txt)
        k = String(m.captures[1])
        haskey(out, k) && continue
        out[k] = m.captures[2] == "true"
    end
    return out
end

function read_mat_csv(path::String)
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    isempty(df) && return nothing
    ch = String.(df[!, 1])
    n = length(ch)
    M = zeros(Float64, n, n)
    for (j, name) in enumerate(ch)
        col = Symbol(name)
        hasproperty(df, col) || return nothing
        M[:, j] = Float64.(df[!, col])
    end
    return (ch, M)
end

function safe_csv(path::String)::DataFrame
    isfile(path) || return DataFrame()
    try
        return CSV.read(path, DataFrame)
    catch
        return DataFrame()
    end
end

function resolve_results_root(proj::String, explicit::Union{Nothing,String}=nothing)::String
    explicit !== nothing && return abspath(explicit)
    cfg_p = joinpath(proj, "config", "pipeline.toml")
    if isfile(cfg_p)
        cfg = TOML.parsefile(cfg_p)
        r = get(get(cfg, "paths", Dict()), "results", "results")
        return isabspath(r) ? r : joinpath(proj, r)
    end
    return joinpath(proj, "results")
end

# ── JSON helpers ───────────────────────────────────────────────

json_escape(s::AbstractString) =
    replace(replace(replace(String(s), "\\" => "\\\\"), "\"" => "\\\""), "\n" => "\\n")

function mat_flat(M::Matrix{Float64})::String
    n = size(M, 1)
    parts = String[]
    sizehint!(parts, n * n)
    for j in 1:n, i in 1:n
        push!(parts, string(round(M[i, j]; digits=6)))
    end
    return "[" * join(parts, ",") * "]"
end

function df_rows_json(df::DataFrame, cols::Vector{Symbol}; limit::Int=0)::String
    nrow(df) == 0 && return "[]"
    n = limit > 0 ? min(limit, nrow(df)) : nrow(df)
    rows = String[]
    for i in 1:n
        fields = String[]
        for c in cols
            hasproperty(df, c) || continue
            v = df[i, c]
            if v isa AbstractString
                push!(fields, "\"$c\":\"$(json_escape(string(v)))\"")
            elseif v isa Bool
                push!(fields, "\"$c\":$(v ? "true" : "false")")
            elseif v isa Integer
                push!(fields, "\"$c\":$v")
            elseif v isa Real
                x = Float64(v)
                push!(fields, "\"$c\":$(isnan(x) || isinf(x) ? "null" : round(x; digits=6))")
            else
                push!(fields, "\"$c\":\"$(json_escape(string(v)))\"")
            end
        end
        push!(rows, "{" * join(fields, ",") * "}")
    end
    return "[" * join(rows, ",") * "]"
end

function pos_json()::String
    parts = ["\"$k\":[$(v[1]),$(v[2])]" for (k, v) in CH_POS]
    return "{" * join(parts, ",") * "}"
end

# ── Band / KPI helpers ─────────────────────────────────────────

function sort_bands_physio(bands)::Vector{String}
    bs = String.(collect(bands))
    order = Dict(b => i for (i, b) in enumerate(BAND_ORDER))
    return sort(bs; by = b -> get(order, b, 1000 + hash(b) % 100))
end

"""
    honest_best_band(band_stats_rows; n_sig_key=:n_sig, diff_key=:diff_mean)

Si hay FDR, elige la banda con más edges q<0.05.
Si n_total_sig=0, devuelve "" (el visor mostrará '— (sin FDR)').
"""
function honest_best_band(rows; n_sig_field=:n_sig, diff_field=:diff_mean)::String
    isempty(rows) && return ""
    n_sigs = [Int(getfield(r, n_sig_field)) for r in rows]
    total = sum(n_sigs)
    total == 0 && return ""
    # prefer more FDR hits; tie-break by |diff|
    scores = [n_sigs[i] * 1000.0 + abs(Float64(getfield(rows[i], diff_field))) for i in eachindex(rows)]
    return string(getfield(rows[argmax(scores)], :band))
end

function count_fdr(edf::DataFrame; qcol=:q_value)::Int
    nrow(edf) == 0 && return 0
    hasproperty(edf, qcol) || return 0
    return count(r -> Float64(r[qcol]) < 0.05, eachrow(edf))
end

function count_uncorr(edf::DataFrame; pcol=:p_value)::Int
    nrow(edf) == 0 && return 0
    hasproperty(edf, pcol) || return 0
    return count(r -> Float64(r[pcol]) < 0.05, eachrow(edf))
end

function _ch_xy(ch::AbstractString)
    k = uppercase(String(ch))
    return get(CH_POS, k, (0.0, 0.0))
end

# ── Static figures (CairoMakie) ────────────────────────────────

"""
    save_heatmap(path, W, ch, title; cmap=:viridis, diverging=false, colorrange=nothing)

Diverging usa `DIVERGING_CMAP` (= `:RdBu_r`) para que positivo = rojo.
"""
function save_heatmap(path, W, ch, title; cmap=:viridis, diverging=false,
                      colorrange=nothing, colorbar_label="")
    n = length(ch)
    fig = Figure(size=(720, 640), fontsize=11)
    ax = Axis(fig[1, 1]; title=title, xlabel="Canal", ylabel="Canal",
              xticks=(1:n, ch), yticks=(1:n, ch),
              xticklabelrotation=π/3, xticklabelsize=7, yticklabelsize=7)
    if diverging
        lim = colorrange === nothing ? maximum(abs, W) : maximum(abs, colorrange)
        lim < 1e-12 && (lim = 1.0)
        hm = heatmap!(ax, W; colormap=DIVERGING_CMAP, colorrange=(-lim, lim))
        Colorbar(fig[1, 2], hm; label=isempty(colorbar_label) ? "Δ" : colorbar_label)
    else
        cr = colorrange
        hm = cr === nothing ? heatmap!(ax, W; colormap=cmap) :
             heatmap!(ax, W; colormap=cmap, colorrange=cr)
        Colorbar(fig[1, 2], hm; label=isempty(colorbar_label) ? "wPLI" : colorbar_label)
    end
    mkpath(dirname(path))
    save(path, fig)
    return path
end

"""
    save_heatmap_triplet(path, Wa, Wb, Wd, ch, titles; shared_abs=true)

Tres paneles A | B | Δ con escala compartida en A/B y Δ divergente simétrica.
`titles` = (title_a, title_b, title_d).
"""
function save_heatmap_triplet(path, Wa, Wb, Wd, ch, titles;
                              label_a="wPLI", label_d="Δ wPLI")
    n = length(ch)
    lim_ab = max(maximum(Wa), maximum(Wb), 1e-12)
    lim_d = maximum(abs, Wd); lim_d < 1e-12 && (lim_d = 1.0)
    fig = Figure(size=(1400, 480), fontsize=10)
    for (col, (mat, ttl, diverging)) in enumerate([
        (Wa, titles[1], false),
        (Wb, titles[2], false),
        (Wd, titles[3], true),
    ])
        ax = Axis(fig[1, col]; title=ttl,
                  xticks=(1:n, ch), yticks=(1:n, ch),
                  xticklabelrotation=π/3, xticklabelsize=6, yticklabelsize=6)
        if diverging
            hm = heatmap!(ax, mat; colormap=DIVERGING_CMAP, colorrange=(-lim_d, lim_d))
            Colorbar(fig[2, col], hm; vertical=false, label=label_d, flipaxis=false)
        else
            hm = heatmap!(ax, mat; colormap=:viridis, colorrange=(0.0, lim_ab))
            Colorbar(fig[2, col], hm; vertical=false, label=label_a, flipaxis=false)
        end
    end
    mkpath(dirname(path))
    save(path, fig)
    return path
end

"""
    save_topo_network(path, ch, edge_rows, title; topn=20)

Grafo en layout 10–20. `edge_rows` iterable de NamedTuples/rows con ch_a, ch_b, diff.
Si vacío, no escribe (caller genera explore top-N).
"""
function save_topo_network(path, ch, edge_rows, title; linewidth_base=1.2)
    isempty(edge_rows) && return nothing
    fig = Figure(size=(720, 700), fontsize=11)
    ax = Axis(fig[1, 1]; title=title, aspect=DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    # head outline
    θs = range(0, 2π; length=120)
    lines!(ax, 1.05 .* cos.(θs), 1.05 .* sin.(θs); color=:gray60, linewidth=1)
    xy = Dict(c => _ch_xy(c) for c in ch)
    diffs = [Float64(r.diff) for r in edge_rows]
    dmax = max(maximum(abs, diffs), 1e-12)
    for r in edge_rows
        a = string(r.ch_a); b = string(r.ch_b)
        (haskey(xy, a) && haskey(xy, b)) || continue
        d = Float64(r.diff)
        col = d >= 0 ? (:firebrick) : (:steelblue)
        lw = linewidth_base + 2.5 * abs(d) / dmax
        lines!(ax, [xy[a][1], xy[b][1]], [xy[a][2], xy[b][2]];
               color=col, linewidth=lw)
    end
    xs = [xy[c][1] for c in ch]; ys = [xy[c][2] for c in ch]
    scatter!(ax, xs, ys; color=:white, strokecolor=:gray20, strokewidth=1.2, markersize=14)
    for c in ch
        text!(ax, xy[c][1], xy[c][2] + 0.06, text=c; fontsize=8, align=(:center, :bottom))
    end
    xlims!(ax, -1.25, 1.25); ylims!(ax, -1.25, 1.25)
    # legend
    Label(fig[2, 1], "Rojo: Δ>0 (↑)   Azul: Δ<0 (↓)   Grosor ∝ |Δ|";
          fontsize=10, color=:gray40, tellwidth=false)
    mkpath(dirname(path))
    save(path, fig)
    return path
end

function save_paired_means(path, subject_means_df; band_order=BAND_ORDER)
    isempty(subject_means_df) && return nothing
    raw = unique(string.(subject_means_df.band))
    bands = sort_bands_physio(intersect(band_order, raw))
    isempty(bands) && (bands = sort_bands_physio(raw))
    fig = Figure(size=(900, 420), fontsize=11)
    ax = Axis(fig[1, 1]; title="Mean wPLI pareado T1 vs T2",
              xlabel="Banda", ylabel="Mean wPLI (triángulo superior)",
              xticks=(1:length(bands), bands))
    for (bi, b) in enumerate(bands)
        sub = filter(r -> string(r.band) == b, eachrow(subject_means_df))
        by_s = Dict{String, Dict{String,Float64}}()
        for r in sub
            sid = string(r.subject_id)
            haskey(by_s, sid) || (by_s[sid] = Dict{String,Float64}())
            by_s[sid][string(r.timepoint)] = Float64(r.mean_wpli)
        end
        for (_, tp) in by_s
            haskey(tp, "T1") && haskey(tp, "T2") || continue
            lines!(ax, [bi - 0.15, bi + 0.15], [tp["T1"], tp["T2"]];
                   color=(:gray60, 0.45), linewidth=1)
            scatter!(ax, [bi - 0.15], [tp["T1"]]; color=:steelblue, markersize=8)
            scatter!(ax, [bi + 0.15], [tp["T2"]]; color=:darkorange, markersize=8)
        end
    end
    mkpath(dirname(path))
    save(path, fig)
    return path
end

function save_group_means(path, subject_means_df; band_order=BAND_ORDER)
    isempty(subject_means_df) && return nothing
    raw = unique(string.(subject_means_df.band))
    bands = sort_bands_physio(intersect(band_order, raw))
    isempty(bands) && (bands = sort_bands_physio(raw))
    fig = Figure(size=(900, 420), fontsize=11)
    ax = Axis(fig[1, 1]; title="Mean wPLI por grupo (MS vs Control)",
              xlabel="Banda", ylabel="Mean wPLI (triángulo superior)",
              xticks=(1:length(bands), bands))
    for (bi, b) in enumerate(bands)
        sub = filter(r -> string(r.band) == b, eachrow(subject_means_df))
        ms_v = [Float64(r.mean_wpli) for r in sub if lowercase(string(r.group)) in ("ms", "em", "patient")]
        ct_v = [Float64(r.mean_wpli) for r in sub if lowercase(string(r.group)) in ("control", "ctrl", "hc")]
        for v in ms_v
            scatter!(ax, [bi - 0.15], [v]; color=(:firebrick, 0.55), markersize=7)
        end
        for v in ct_v
            scatter!(ax, [bi + 0.15], [v]; color=(:steelblue, 0.55), markersize=7)
        end
        if !isempty(ms_v)
            m = mean(ms_v); s = length(ms_v) > 1 ? std(ms_v) / sqrt(length(ms_v)) : 0.0
            lines!(ax, [bi - 0.22, bi - 0.08], [m, m]; color=:firebrick, linewidth=2.5)
            lines!(ax, [bi - 0.15, bi - 0.15], [m - s, m + s]; color=:firebrick, linewidth=1.5)
        end
        if !isempty(ct_v)
            m = mean(ct_v); s = length(ct_v) > 1 ? std(ct_v) / sqrt(length(ct_v)) : 0.0
            lines!(ax, [bi + 0.08, bi + 0.22], [m, m]; color=:steelblue, linewidth=2.5)
            lines!(ax, [bi + 0.15, bi + 0.15], [m - s, m + s]; color=:steelblue, linewidth=1.5)
        end
    end
    mkpath(dirname(path))
    save(path, fig)
    return path
end

function save_topo_delta(path, ch, delta_vals, xy, title; colorbar_label="Δ")
    fig = Figure(size=(560, 520), fontsize=11)
    ax = Axis(fig[1, 1]; title=title, aspect=DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    θs = range(0, 2π; length=120)
    lines!(ax, 1.05 .* cos.(θs), 1.05 .* sin.(θs); color=:gray60, linewidth=1)
    xs = Float64[]; ys = Float64[]; vs = Float64[]; labs = String[]
    for (i, c) in enumerate(ch)
        if haskey(xy, c)
            # proyectar x,y BIDS a [-1,1] aprox. usando solo x,y
            push!(xs, xy[c][1]); push!(ys, xy[c][2])
        else
            p = _ch_xy(c)
            push!(xs, p[1]); push!(ys, p[2])
        end
        push!(vs, delta_vals[i]); push!(labs, c)
    end
    # normalizar si coords BIDS fuera de [-1.2,1.2]
    if !isempty(xs)
        mx = max(maximum(abs, xs), maximum(abs, ys), 1e-9)
        if mx > 1.2
            xs ./= mx; ys ./= mx
        end
    end
    lim = maximum(abs, vs); lim < 1e-12 && (lim = 1.0)
    sc = scatter!(ax, xs, ys; color=vs, colormap=DIVERGING_CMAP, colorrange=(-lim, lim),
                  markersize=28)
    for (i, lab) in enumerate(labs)
        text!(ax, xs[i], ys[i] + 0.05, text=lab; fontsize=7, align=(:center, :bottom))
    end
    Colorbar(fig[1, 2], sc; label=colorbar_label)
    mkpath(dirname(path))
    save(path, fig)
    return path
end

end # module
