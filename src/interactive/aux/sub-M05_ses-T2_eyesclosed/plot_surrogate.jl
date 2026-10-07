#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor interactivo de surrogates / inferencia wPLI
# ═══════════════════════════════════════════════════════════════
#
#  Overview · heatmaps p/q/z/sig · contraste obs vs nula
#  · grafo FDR-significativo · volcano + tabla
#
#  Fuente (paso [SUR]):
#    tables/surrogate/*
#    json/surrogate_summary.json
#    cache/surrogate/null_distribution_{band}.jls  (opcional, hist empírico)
#
#  Uso:
#    julia --project=. src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_surrogate.jl
#    # → http://127.0.0.1:8775/
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/aux/sub-M05_ses-T2_eyesclosed/plot_surrogate.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      25-07-2026
#  Modificado  25-07-2026 — versionado fuera de results/ (antes figures/aux/)
# ───────────────────────────────────────────────────────────────

using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics, Serialization

const HERE         = @__DIR__
const RESULTS_UNIT = normpath(joinpath(HERE, "..", "..", "..", "..",
                     "results", "subjects", "sub-M05", "ses-T2", "eyesclosed"))
const SURR_DIR  = joinpath(RESULTS_UNIT, "tables", "surrogate")
const JSONDIR   = joinpath(RESULTS_UNIT, "json")
const CACHE_DIR = joinpath(RESULTS_UNIT, "cache", "surrogate")
const FIG_DIR   = joinpath(RESULTS_UNIT, "figures", "surrogate")
const HOST      = "127.0.0.1"
const PORT      = 8775

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
    "occipital"=>"Occipital","temporal"=>"Temporal","other"=>"Otra",
)

mutable struct BandPack
    observed::Matrix{Float64}
    pvalues::Matrix{Float64}
    qvalues::Matrix{Float64}
    significant::Matrix{Float64}
    stats::DataFrame                 # null_stats rows
    null_dist::Union{Nothing,Array{Float64,3}}
    fdr_thr::Float64
    n_surrogates::Int
end

mutable struct SurrStore
    channels::Vector{String}
    bands::Vector{String}
    packs::Dict{String,BandPack}
    significant::DataFrame           # all FDR-sig edges
    quality::DataFrame
    summary::Dict{String,Any}
    alpha::Float64
    n_surrogates::Int
    method::String
    fdr_method::String
    seed::Int
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

function _read_matrix(path::String)::Tuple{Vector{String},Matrix{Float64}}
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

"""Lee clave de Dict cache admitiendo String o Symbol."""
function _cache_get(cache, key::String, default=nothing)
    haskey(cache, key) && return cache[key]
    sk = Symbol(key)
    haskey(cache, sk) && return cache[sk]
    return default
end

function _try_load_null_cache(band::String)
    p = joinpath(CACHE_DIR, "null_distribution_$(band).jls")
    isfile(p) || return nothing
    try
        return Serialization.deserialize(p)
    catch e
        @warn "No se pudo leer cache nulo $band: $e"
        return nothing
    end
end

function _extract_null_dist(cache, n_ch::Int)::Union{Nothing,Array{Float64,3}}
    cache === nothing && return nothing
    raw = _cache_get(cache, "null_distribution")
    raw === nothing && return nothing
    nd = raw isa Array{Float64,3} ? raw : Array{Float64,3}(raw)
    size(nd, 1) == n_ch && size(nd, 2) == n_ch || return nothing
    size(nd, 3) >= 1 || return nothing
    return nd
end

function load_store()::SurrStore
    isdir(SURR_DIR) || error("No encontrado: $SURR_DIR")
    summary = _parse_summary_json(joinpath(JSONDIR, "surrogate_summary.json"))
    quality_path = joinpath(SURR_DIR, "surrogate_quality.csv")
    sig_path = joinpath(SURR_DIR, "significant_connections.csv")
    isfile(quality_path) || error("No encontrado: $quality_path")

    quality = CSV.read(quality_path, DataFrame)
    quality.band = String.(quality.band)

    significant = isfile(sig_path) ? CSV.read(sig_path, DataFrame) : DataFrame(
        ch_a=String[], ch_b=String[], band=String[],
        wpli_obs=Float64[], p_value=Float64[], q_value=Float64[], z_score=Float64[])
    if nrow(significant) > 0
        significant.ch_a = String.(significant.ch_a)
        significant.ch_b = String.(significant.ch_b)
        significant.band = String.(significant.band)
    end

    channels = String[]
    bands = String[]
    packs = Dict{String,BandPack}()
    n_cache_ok = 0

    for b in BAND_ORDER
        obs_p = joinpath(SURR_DIR, "wpli_observed_$(b).csv")
        p_p   = joinpath(SURR_DIR, "wpli_pvalues_$(b).csv")
        q_p   = joinpath(SURR_DIR, "wpli_qvalues_$(b).csv")
        s_p   = joinpath(SURR_DIR, "wpli_significant_$(b).csv")
        st_p  = joinpath(SURR_DIR, "surrogate_null_stats_$(b).csv")
        all(isfile, (obs_p, p_p, q_p, s_p, st_p)) || continue

        ch, Mobs = _read_matrix(obs_p)
        _, Mp = _read_matrix(p_p)
        _, Mq = _read_matrix(q_p)
        _, Ms = _read_matrix(s_p)
        if isempty(channels)
            channels = ch
        else
            ch == channels || error("Canales distintos en $obs_p")
        end
        stats = CSV.read(st_p, DataFrame)
        stats.ch_a = String.(stats.ch_a)
        stats.ch_b = String.(stats.ch_b)

        fdr = 0.0
        nsur = Int(round(Float64(get(summary, "n_surrogates", 500))))
        qrow = quality[quality.band .== b, :]
        if nrow(qrow) > 0
            fdr = Float64(qrow.fdr_threshold[1])
            nsur = Int(qrow.n_surrogates[1])
        else
            fdr = Float64(get(summary, "$(b)_fdr_thr", 0.0))
        end

        cache = _try_load_null_cache(b)
        null_d = _extract_null_dist(cache, length(ch))
        if null_d !== nothing
            n_cache_ok += 1
            nsur = Int(_cache_get(cache, "n_surrogates", nsur))
            fdr = Float64(_cache_get(cache, "fdr_threshold", fdr))
            # Si el cache trae matrices observadas/p, preferirlas (misma corrida)
            obs_c = _cache_get(cache, "observed")
            p_c = _cache_get(cache, "p_values")
            if obs_c isa AbstractMatrix && size(obs_c) == size(Mobs)
                Mobs = Array{Float64}(obs_c)
            end
            if p_c isa AbstractMatrix && size(p_c) == size(Mp)
                Mp = Array{Float64}(p_c)
            end
        end

        packs[b] = BandPack(Mobs, Mp, Mq, Ms, stats, null_d, fdr, nsur)
        push!(bands, b)
    end
    isempty(bands) && error("Ninguna banda surrogate completa en $SURR_DIR")

    alpha = Float64(get(summary, "alpha", 0.05))
    nsur = Int(round(Float64(get(summary, "n_surrogates", packs[bands[1]].n_surrogates))))
    method = string(get(summary, "method", "circular_shift"))
    fdr_m = string(get(summary, "fdr_method", "bh"))
    seed = Int(round(Float64(get(summary, "seed", 42))))

    println("  tablas: $SURR_DIR")
    println("  cache:  $CACHE_DIR  ($n_cache_ok/$(length(bands)) bandas con null_distribution)")

    return SurrStore(channels, bands, packs, significant, quality, summary,
                     alpha, nsur, method, fdr_m, seed)
end

# ── Helpers ────────────────────────────────────────────────────

function _json_escape(s::AbstractString)
    t = replace(String(s), "\\" => "\\\\")
    t = replace(t, "\"" => "\\\"")
    t = replace(t, "\n" => "\\n")
    t = replace(t, "\r" => "\\r")
    return t
end

function _region(ch::String)::String
    get(CH_REGION, uppercase(ch), "other")
end

function _pos_xy(ch::String)::Tuple{Float64,Float64}
    get(CH_POS, uppercase(ch), (NaN, NaN))
end

function _resolve_channel(store::SurrStore, ch::String)::String
    ch in store.channels && return ch
    i = findfirst(c -> uppercase(c) == uppercase(ch), store.channels)
    i === nothing && error("Canal desconocido: $ch")
    return store.channels[i]
end

function _ch_index(store::SurrStore, ch::String)::Int
    c = _resolve_channel(store, ch)
    return something(findfirst(==(c), store.channels))
end

function _p_floor(n_sur::Int)::Float64
    1.0 / (n_sur + 1)
end

function _pct_p_floor(stats::DataFrame, n_sur::Int)::Float64
    nrow(stats) == 0 && return 0.0
    pf = _p_floor(n_sur)
    count(p -> abs(Float64(p) - pf) < 1e-9 || Float64(p) <= pf + 1e-12,
          stats.p_value) / nrow(stats) * 100
end

function _matrix_layer(pack::BandPack, layer::String)::Matrix{Float64}
    layer == "obs"  && return pack.observed
    layer == "p"    && return pack.pvalues
    layer == "q"    && return pack.qvalues
    layer == "sig"  && return pack.significant
    layer == "z"    && return _z_matrix(pack)
    layer == "neglogp" && return .-log10.(clamp.(pack.pvalues, 1e-12, 1.0))
    error("Capa desconocida: $layer")
end

function _z_matrix(pack::BandPack)::Matrix{Float64}
    n = size(pack.observed, 1)
    Z = zeros(Float64, n, n)
    ch = nothing
    # build from stats via channel names in pack.stats
    # We need channel order — use store channels passed separately; rebuild from stats pairs
    # Fallback: compute from obs/null if columns present
    for r in eachrow(pack.stats)
        # filled later with channel index map in API
        nothing
    end
    return Z
end

function _z_matrix(store::SurrStore, band::String)::Matrix{Float64}
    pack = store.packs[band]
    n = length(store.channels)
    Z = zeros(Float64, n, n)
    idx = Dict(c => i for (i, c) in enumerate(store.channels))
    for r in eachrow(pack.stats)
        i = get(idx, String(r.ch_a), 0)
        j = get(idx, String(r.ch_b), 0)
        (i == 0 || j == 0) && continue
        z = Float64(r.z_score)
        Z[i, j] = z
        Z[j, i] = z
    end
    return Z
end

function _layer_matrix(store::SurrStore, band::String, layer::String)::Matrix{Float64}
    pack = store.packs[band]
    layer == "z" && return _z_matrix(store, band)
    layer == "neglogp" && return .-log10.(clamp.(pack.pvalues, 1e-12, 1.0))
    layer == "obs" && return pack.observed
    layer == "p" && return pack.pvalues
    layer == "q" && return pack.qvalues
    layer == "sig" && return pack.significant
    error("Capa desconocida: $layer")
end

function _flat_upper(M::Matrix{Float64})::Vector{Float64}
    n = size(M, 1)
    Float64[M[i, j] for i in 1:n-1 for j in i+1:n]
end

function _edge_row(pack::BandPack, a::String, b::String)
    for r in eachrow(pack.stats)
        ca, cb = String(r.ch_a), String(r.ch_b)
        if (ca == a && cb == b) || (ca == b && cb == a)
            return r
        end
    end
    return nothing
end

function _null_samples(pack::BandPack, store::SurrStore, a::String, b::String)::Union{Nothing,Vector{Float64}}
    pack.null_dist === nothing && return nothing
    i = _ch_index(store, a)
    j = _ch_index(store, b)
    return vec(pack.null_dist[i, j, :])
end

function _normal_pdf_samples(μ::Float64, σ::Float64; n::Int = 400, lo=nothing, hi=nothing)
    σ <= 0 && (σ = 1e-6)
    xlo = lo === nothing ? μ - 4σ : Float64(lo)
    xhi = hi === nothing ? μ + 4σ : Float64(hi)
    xlo = min(xlo, 0.0)
    xs = range(xlo, xhi; length=n)
    ys = [exp(-0.5 * ((x - μ) / σ)^2) / (σ * sqrt(2π)) for x in xs]
    return collect(xs), ys
end

function _hist_counts(vals::Vector{Float64}; nbins::Int = 30)
    isempty(vals) && return (0.0, 1.0, zeros(Int, nbins))
    lo = min(0.0, minimum(vals))
    hi = max(maximum(vals), lo + 1e-6)
    bw = (hi - lo) / nbins
    counts = zeros(Int, nbins)
    for v in vals
        k = min(nbins, max(1, Int(floor((v - lo) / bw)) + 1))
        v >= hi && (k = nbins)
        counts[k] += 1
    end
    return lo, hi, counts
end

function _auto_interpret(store::SurrStore)::String
    lines = String[]
    push!(lines, "Test $(store.method) · N=$(store.n_surrogates) · FDR-$(store.fdr_method) α=$(store.alpha) · p_min≈$(round(_p_floor(store.n_surrogates); digits=5)).")
    # rank bands by pct_sig
    rows = [(String(r.band), Float64(r.pct_sig), Int(r.n_sig)) for r in eachrow(store.quality)]
    sort!(rows; by = x -> -x[2])
    sig = filter(r -> r[3] > 0, rows)
    zero = filter(r -> r[3] == 0, rows)
    if !isempty(sig)
        tops = join(["$(b) ($(round(p; digits=1))%)" for (b, p, _) in sig[1:min(3, length(sig))]], ", ")
        push!(lines, "Bandas con conexiones FDR-significativas: $tops.")
    else
        push!(lines, "Ninguna banda supera FDR a α=$(store.alpha).")
    end
    if !isempty(zero)
        push!(lines, "Sin pares FDR-sig: " * join(first.(zero), ", ") *
              " (pueden tener p<0.05 o z altos que no sobreviven a la corrección múltiple).")
    end
    # hubs
    if nrow(store.significant) > 0
        deg = Dict{String,Int}()
        for r in eachrow(store.significant)
            deg[String(r.ch_a)] = get(deg, String(r.ch_a), 0) + 1
            deg[String(r.ch_b)] = get(deg, String(r.ch_b), 0) + 1
        end
        hubs = sort(collect(deg); by = x -> -x[2])
        top = join(["$(h[1])($(h[2]))" for h in hubs[1:min(4, length(hubs))]], ", ")
        push!(lines, "Hubs en red significativa (degree): $top.")
    end
    return join(lines, " ")
end

# ── API JSON ───────────────────────────────────────────────────

function _overview_json(store::SurrStore)::String
    band_rows = String[]
    for b in store.bands
        pack = store.packs[b]
        qrow = store.quality[store.quality.band .== b, :]
        n_sig = nrow(qrow) > 0 ? Int(qrow.n_sig[1]) : 0
        pct = nrow(qrow) > 0 ? Float64(qrow.pct_sig[1]) : 0.0
        mean_p = nrow(qrow) > 0 ? Float64(qrow.mean_p[1]) : mean(pack.stats.p_value)
        obs_mean = nrow(qrow) > 0 ? Float64(qrow.obs_mean[1]) : mean(pack.stats.wpli_obs)
        obs_max = nrow(qrow) > 0 ? Float64(qrow.obs_max[1]) : maximum(pack.stats.wpli_obs)
        pf_pct = round(_pct_p_floor(pack.stats, pack.n_surrogates); digits=1)
        has_cache = pack.null_dist !== nothing
        # top edge
        top = pack.stats[argmax(pack.stats.wpli_obs), :]
        push!(band_rows,
            "{\"band\":\"$b\",\"n_sig\":$n_sig,\"pct_sig\":$(round(pct; digits=2))," *
            "\"fdr_thr\":$(round(pack.fdr_thr; digits=4)),\"mean_p\":$(round(mean_p; digits=4))," *
            "\"obs_mean\":$(round(obs_mean; digits=4)),\"obs_max\":$(round(obs_max; digits=4))," *
            "\"pct_p_floor\":$pf_pct,\"has_null_cache\":$(has_cache)," *
            "\"top_edge\":\"$(top.ch_a)–$(top.ch_b)\",\"top_wpli\":$(round(Float64(top.wpli_obs); digits=4))," *
            "\"top_z\":$(round(Float64(top.z_score); digits=3)),\"top_q\":$(round(Float64(top.q_value); digits=4))}")
    end
    interp = _auto_interpret(store)
    return "{\"subject\":\"$SUBJECT\",\"session\":\"$SESSION\",\"task\":\"$TASK_ID\"," *
           "\"method\":\"$(_json_escape(store.method))\"," *
           "\"n_surrogates\":$(store.n_surrogates)," *
           "\"alpha\":$(store.alpha)," *
           "\"fdr_method\":\"$(_json_escape(store.fdr_method))\"," *
           "\"seed\":$(store.seed)," *
           "\"p_min\":$(round(_p_floor(store.n_surrogates); digits=6))," *
           "\"n_channels\":$(length(store.channels))," *
           "\"n_pairs\":$(length(store.channels)*(length(store.channels)-1)÷2)," *
           "\"n_sig_total\":$(nrow(store.significant))," *
           "\"bands\":[$(join(band_rows, ","))]," *
           "\"interpretation\":\"$(_json_escape(interp))\"}"
end

function _band_json(store::SurrStore, band::String, layer::String)::String
    haskey(store.packs, band) || error("Banda desconocida: $band")
    pack = store.packs[band]
    M = _layer_matrix(store, band, layer)
    ch = store.channels
    n = length(ch)
    flat = Float64[]
    sizehint!(flat, n * n)
    for j in 1:n, i in 1:n
        push!(flat, M[i, j])
    end
    off = Float64[M[i, j] for i in 1:n for j in 1:n if i != j]
    vmax = isempty(off) ? 1.0 : maximum(abs.(off))
    vmin = layer in ("z",) ? -vmax : 0.0
    if layer == "sig"
        vmin, vmax = 0.0, 1.0
    elseif layer == "obs"
        vmin, vmax = 0.0, max(maximum(off), 1e-6)
    elseif layer in ("p", "q")
        vmin, vmax = 0.0, 1.0
    end

    # edges for volcano / table
    edges_j = join([
        "{\"ch_a\":\"$(r.ch_a)\",\"ch_b\":\"$(r.ch_b)\"," *
        "\"wpli\":$(round(Float64(r.wpli_obs); digits=6))," *
        "\"null_mean\":$(round(Float64(r.null_mean); digits=6))," *
        "\"null_std\":$(round(Float64(r.null_std); digits=6))," *
        "\"p\":$(round(Float64(r.p_value); digits=6))," *
        "\"q\":$(round(Float64(r.q_value); digits=6))," *
        "\"z\":$(round(Float64(r.z_score); digits=4))," *
        "\"sig\":$(Float64(r.q_value) < store.alpha)}"
        for r in eachrow(pack.stats)
    ], ",")

    ch_j = join(["\"$c\"" for c in ch], ",")
    return "{\"band\":\"$band\",\"layer\":\"$layer\",\"n\":$n,\"channels\":[$ch_j]," *
           "\"matrix\":[" * join(round.(flat; digits=6), ",") * "]," *
           "\"vmin\":$(round(vmin; digits=6)),\"vmax\":$(round(vmax; digits=6))," *
           "\"fdr_thr\":$(round(pack.fdr_thr; digits=6))," *
           "\"n_surrogates\":$(pack.n_surrogates)," *
           "\"has_null_cache\":$(pack.null_dist !== nothing)," *
           "\"alpha\":$(store.alpha)," *
           "\"edges\":[$edges_j]}"
end

function _contrast_json(store::SurrStore, band::String, a::String, b::String)::String
    haskey(store.packs, band) || error("Banda desconocida: $band")
    pack = store.packs[band]
    a = _resolve_channel(store, a)
    b = _resolve_channel(store, b)
    row = _edge_row(pack, a, b)
    row === nothing && error("Par no encontrado: $(a)–$(b)")
    obs = Float64(row.wpli_obs)
    μ = Float64(row.null_mean)
    σ = Float64(row.null_std)
    p = Float64(row.p_value)
    q = Float64(row.q_value)
    z = Float64(row.z_score)
    sig = q < store.alpha

    samples = _null_samples(pack, store, a, b)
    source = "empirical"
    hist_lo, hist_hi, counts = 0.0, 1.0, zeros(Int, 30)
    pdf_x = Float64[]; pdf_y = Float64[]
    if samples !== nothing && !isempty(samples)
        hist_lo, hist_hi, counts = _hist_counts(samples; nbins=32)
        hist_hi = max(hist_hi, obs * 1.05, μ + 3σ)
        # media empírica de las muestras (puede diferir levemente del CSV redondeado)
        μ = mean(samples)
    else
        source = "normal_approx"
        hist_hi = max(obs * 1.1, μ + 4σ, 0.2)
        pdf_x, pdf_y = _normal_pdf_samples(μ, σ; lo=0.0, hi=hist_hi)
        synth = clamp.(μ .+ σ .* randn(pack.n_surrogates), 0.0, 1.0)
        hist_lo, hist_hi, counts = _hist_counts(synth; nbins=32)
        hist_hi = max(hist_hi, obs * 1.05)
    end

    verdict = sig ? "conectividad significativa (q < α)" : "no significativa tras FDR"
    return "{\"band\":\"$band\",\"ch_a\":\"$a\",\"ch_b\":\"$b\"," *
           "\"wpli_obs\":$(round(obs; digits=6))," *
           "\"null_mean\":$(round(μ; digits=6)),\"null_std\":$(round(σ; digits=6))," *
           "\"p\":$(round(p; digits=6)),\"q\":$(round(q; digits=6)),\"z\":$(round(z; digits=4))," *
           "\"sig\":$sig,\"alpha\":$(store.alpha)," *
           "\"fdr_thr\":$(round(pack.fdr_thr; digits=6))," *
           "\"n_surrogates\":$(pack.n_surrogates)," *
           "\"p_min\":$(round(_p_floor(pack.n_surrogates); digits=6))," *
           "\"source\":\"$source\",\"verdict\":\"$(_json_escape(verdict))\"," *
           "\"hist\":{\"lo\":$hist_lo,\"hi\":$hist_hi,\"nbins\":$(length(counts))," *
           "\"counts\":[$(join(string.(counts), ","))]}," *
           "\"pdf\":{\"x\":[$(join(string.(round.(pdf_x; digits=5)), ","))]," *
           "\"y\":[$(join(string.(round.(pdf_y; digits=6)), ","))]}}"
end

function _topo_json(store::SurrStore, band::String)::String
    haskey(store.packs, band) || error("Banda desconocida: $band")
    edf = store.significant[store.significant.band .== band, :]
    deg = Dict{String,Int}(c => 0 for c in store.channels)
    str = Dict{String,Float64}(c => 0.0 for c in store.channels)
    for r in eachrow(edf)
        a, b = String(r.ch_a), String(r.ch_b)
        deg[a] = get(deg, a, 0) + 1
        deg[b] = get(deg, b, 0) + 1
        w = Float64(r.wpli_obs)
        str[a] = get(str, a, 0.0) + w
        str[b] = get(str, b, 0.0) + w
    end
    nodes_j = join([
        let xy = _pos_xy(c)
            "{\"channel\":\"$c\",\"x\":$(xy[1]),\"y\":$(xy[2])," *
            "\"degree\":$(deg[c]),\"strength\":$(round(str[c]; digits=4))," *
            "\"region\":\"$(_region(c))\"}"
        end for c in store.channels
    ], ",")
    edges_j = join([
        "{\"ch_a\":\"$(r.ch_a)\",\"ch_b\":\"$(r.ch_b)\"," *
        "\"wpli\":$(round(Float64(r.wpli_obs); digits=4))," *
        "\"q\":$(round(Float64(r.q_value); digits=4))," *
        "\"z\":$(round(Float64(r.z_score); digits=3))}"
        for r in eachrow(edf)
    ], ",")
    hubs = sort(collect(deg); by = x -> -x[2])
    hubs_j = join([
        "{\"channel\":\"$(h[1])\",\"degree\":$(h[2]),\"strength\":$(round(str[h[1]]; digits=4))}"
        for h in hubs if h[2] > 0
    ], ",")
    return "{\"band\":\"$band\",\"n_sig\":$(nrow(edf)),\"alpha\":$(store.alpha)," *
           "\"nodes\":[$nodes_j],\"edges\":[$edges_j],\"hubs\":[$hubs_j]}"
end

# ── PNG export ─────────────────────────────────────────────────

function save_overview_png(store::SurrStore)::String
    mkpath(FIG_DIR)
    bands = store.bands
    pct = Float64[let q = store.quality[store.quality.band .== b, :]
        nrow(q) > 0 ? Float64(q.pct_sig[1]) : 0.0
    end for b in bands]
    fig = Figure(size=(720, 420), fontsize=12)
    ax = Axis(fig[1, 1];
        title="$(SUBJECT)/$(SESSION)/$(TASK_ID) — % pares FDR-sig por banda",
        xlabel="Banda", ylabel="% significativos",
        xticks=(1:length(bands), bands), xticklabelrotation=π/4)
    barplot!(ax, 1:length(bands), pct; color=:steelblue)
    hlines!(ax, [store.alpha * 100]; color=:gray, linestyle=:dash, label="α×100 (ref)")
    out = joinpath(FIG_DIR, "surrogate_overview_pct_sig.png")
    save(out, fig; px_per_unit=2)
    return out
end

function save_contrast_png(store::SurrStore, band::String, a::String, b::String)::String
    mkpath(FIG_DIR)
    pack = store.packs[band]
    row = _edge_row(pack, a, b)
    row === nothing && error("Par no encontrado")
    obs = Float64(row.wpli_obs)
    μ = Float64(row.null_mean)
    samples = _null_samples(pack, store, a, b)
    fig = Figure(size=(700, 420), fontsize=12)
    ax = Axis(fig[1, 1];
        title="Contraste wPLI — $band · $a ↔ $b",
        xlabel="wPLI", ylabel="Frecuencia")
    if samples !== nothing
        hist!(ax, samples; bins=32, color=(:steelblue, 0.65), label="Valores subrogados")
    else
        σ = max(Float64(row.null_std), 1e-6)
        synth = clamp.(μ .+ σ .* randn(pack.n_surrogates), 0.0, 1.0)
        hist!(ax, synth; bins=32, color=(:steelblue, 0.45), label="Nula ≈ Normal(μ,σ)")
    end
    vlines!(ax, [μ]; color=:royalblue, linestyle=:dash, linewidth=2, label="Media subrogada")
    vlines!(ax, [obs]; color=:tomato, linewidth=2, label="wPLI observado")
    axislegend(ax; position=:rt, labelsize=9)
    Label(fig[2, 1],
        "p=$(round(Float64(row.p_value); digits=4)) · q=$(round(Float64(row.q_value); digits=4)) · " *
        "z=$(round(Float64(row.z_score); digits=2)) · " *
        (Float64(row.q_value) < store.alpha ? "FDR-sig" : "no sig");
        fontsize=11)
    out = joinpath(FIG_DIR, "surrogate_contrast_$(band)_$(a)_$(b).png")
    save(out, fig; px_per_unit=2)
    return out
end

function save_topo_png(store::SurrStore, band::String)::String
    mkpath(FIG_DIR)
    edf = store.significant[store.significant.band .== band, :]
    fig = Figure(size=(640, 640), fontsize=11)
    ax = Axis(fig[1, 1]; title="Red FDR-sig — $band (α=$(store.alpha))",
              aspect=DataAspect())
    hidespines!(ax); hidedecorations!(ax)
    θ = range(0, 2π; length=200)
    lines!(ax, cos.(θ), sin.(θ); color=:gray70, linewidth=1)
    # edges
    for r in eachrow(edf)
        xa, ya = _pos_xy(String(r.ch_a))
        xb, yb = _pos_xy(String(r.ch_b))
        (isnan(xa) || isnan(xb)) && continue
        w = Float64(r.wpli_obs)
        lines!(ax, [xa, xb], [ya, yb]; color=(:steelblue, 0.35 + 0.5*clamp(w, 0, 1)),
               linewidth=0.8 + 2.5*clamp(w, 0, 1))
    end
    deg = Dict{String,Int}(c => 0 for c in store.channels)
    for r in eachrow(edf)
        deg[String(r.ch_a)] = get(deg, String(r.ch_a), 0) + 1
        deg[String(r.ch_b)] = get(deg, String(r.ch_b), 0) + 1
    end
    dmax = max(maximum(values(deg)), 1)
    for c in store.channels
        x, y = _pos_xy(c)
        isnan(x) && continue
        r = 0.035 + 0.06 * deg[c] / dmax
        poly!(ax, Point2f[(x + r*cos(t), y + r*sin(t)) for t in range(0, 2π; length=24)];
              color=deg[c] > 0 ? :tomato : :gray75, strokecolor=:white, strokewidth=1)
        text!(ax, x, y + 0.08, text=c; align=(:center, :bottom), fontsize=8)
    end
    out = joinpath(FIG_DIR, "surrogate_topo_sig_$(band).png")
    save(out, fig; px_per_unit=2)
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

function _parse_qs(path::AbstractString)::Dict{String,String}
    qs = Dict{String,String}()
    occursin('?', path) || return qs
    for pair in split(split(path, '?', limit=2)[2], '&')
        kv = split(pair, '=', limit=2)
        length(kv) == 2 && (qs[String(kv[1])] = String(kv[2]))
    end
    return qs
end

function _urldecode(s::String)::String
    replace(s, "%20" => " ", "%2D" => "-", "%2d" => "-")
end

function handle_request(sock, store::SurrStore)
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
            cache_ok = Dict(b => store.packs[b].null_dist !== nothing for b in store.bands)
            cache_j = join(["\"$b\":$(cache_ok[b])" for b in store.bands], ",")
            body = "{\"channels\":[$ch_j],\"bands\":[$b_j]," *
                   "\"alpha\":$(store.alpha),\"n_surrogates\":$(store.n_surrogates)," *
                   "\"method\":\"$(_json_escape(store.method))\"," *
                   "\"subject\":\"$SUBJECT\",\"session\":\"$SESSION\",\"task\":\"$TASK_ID\"," *
                   "\"positions\":{" * join(pos_parts, ",") * "}," *
                   "\"null_cache\":{$cache_j}}"
            _send(sock, 200, body; content_type = "application/json")

        elseif method == "GET" && path_only == "/api/overview"
            _send(sock, 200, _overview_json(store); content_type = "application/json")

        elseif method == "GET" && path_only == "/api/band"
            qs = _parse_qs(path)
            band = _urldecode(String(get(qs, "band", "ALPHA")))
            layer = String(get(qs, "layer", "neglogp"))
            _send(sock, 200, _band_json(store, band, layer); content_type = "application/json")

        elseif method == "GET" && path_only == "/api/contrast"
            qs = _parse_qs(path)
            band = _urldecode(String(get(qs, "band", "ALPHA")))
            a = _urldecode(String(get(qs, "a", "")))
            b = _urldecode(String(get(qs, "b", "")))
            if isempty(a) || isempty(b)
                # default: top wpli of band
                pack = store.packs[band]
                top = pack.stats[argmax(pack.stats.wpli_obs), :]
                a, b = String(top.ch_a), String(top.ch_b)
            end
            _send(sock, 200, _contrast_json(store, band, a, b); content_type = "application/json")

        elseif method == "GET" && path_only == "/api/topo"
            qs = _parse_qs(path)
            band = _urldecode(String(get(qs, "band", "ALPHA")))
            _send(sock, 200, _topo_json(store, band); content_type = "application/json")

        elseif method == "POST" && path_only == "/api/save"
            raw = _read_body(sock, headers)
            view = let m = match(r"\"view\"\s*:\s*\"([^\"]*)\"", raw); m === nothing ? "overview" : m.captures[1] end
            band = let m = match(r"\"band\"\s*:\s*\"([^\"]*)\"", raw); m === nothing ? "ALPHA" : m.captures[1] end
            a = let m = match(r"\"a\"\s*:\s*\"([^\"]*)\"", raw); m === nothing ? "" : m.captures[1] end
            b = let m = match(r"\"b\"\s*:\s*\"([^\"]*)\"", raw); m === nothing ? "" : m.captures[1] end
            out = if view == "overview"
                save_overview_png(store)
            elseif view == "contrast"
                isempty(a) && (a = String(store.packs[band].stats[argmax(store.packs[band].stats.wpli_obs), :].ch_a))
                isempty(b) && (b = String(store.packs[band].stats[argmax(store.packs[band].stats.wpli_obs), :].ch_b))
                save_contrast_png(store, band, a, b)
            elseif view == "topo"
                save_topo_png(store, band)
            else
                error("Vista save desconocida: $view")
            end
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

function html_page(store::SurrStore)::String
    band_opts = join([
        "<option value=\"$b\"" * (b == "ALPHA" ? " selected" : "") * ">$b</option>"
        for b in store.bands
    ], "\n")
    """
<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<title>NeuroMIND · Surrogates</title>
<style>
  :root { --bg:#f4f6f8; --card:#fff; --ink:#1a2332; --muted:#5b6b7c; --line:#d5dde6;
          --accent:#1f6feb; --ok:#0a7a3e; --warn:#b45309; --bad:#b42318; }
  * { box-sizing: border-box; }
  body { margin:0; font-family: "IBM Plex Sans", "Segoe UI", sans-serif; background:var(--bg); color:var(--ink); }
  header { display:flex; flex-wrap:wrap; gap:12px; align-items:center; justify-content:space-between;
           padding:14px 18px; background:linear-gradient(120deg,#0f2744,#1f4e79); color:#fff; }
  header h1 { margin:0; font-size:18px; font-weight:600; letter-spacing:.02em; }
  header .meta { font-size:12px; opacity:.85; }
  .controls { display:flex; flex-wrap:wrap; gap:10px; align-items:center; padding:10px 18px; background:var(--card); border-bottom:1px solid var(--line); }
  label { font-size:12px; color:var(--muted); display:flex; flex-direction:column; gap:3px; }
  select, button { font:inherit; font-size:13px; padding:6px 10px; border:1px solid var(--line); border-radius:6px; background:#fff; }
  button { cursor:pointer; background:var(--accent); color:#fff; border-color:var(--accent); }
  button.secondary { background:#fff; color:var(--ink); }
  .grid { display:grid; grid-template-columns: 1.1fr 1fr; gap:14px; padding:14px 18px 28px; }
  .panel { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:12px 14px; }
  .panel-title { font-size:13px; font-weight:600; margin:0 0 8px; }
  .full { grid-column: 1 / -1; }
  .interp { background:#eef5ff; border:1px solid #c9dbf5; border-radius:8px; padding:10px 12px; font-size:13px; line-height:1.45; }
  .kpis { display:grid; grid-template-columns: repeat(auto-fit,minmax(90px,1fr)); gap:8px; margin:10px 0; }
  .kpi { background:#f7f9fc; border:1px solid var(--line); border-radius:8px; padding:8px; text-align:center; }
  .kpi .v { font-size:16px; font-weight:650; }
  .kpi .l { font-size:10px; color:var(--muted); text-transform:uppercase; letter-spacing:.04em; }
  canvas { width:100%; background:#fafbfc; border:1px solid var(--line); border-radius:8px; cursor:crosshair; }
  table { width:100%; border-collapse:collapse; font-size:12px; }
  th, td { padding:5px 6px; border-bottom:1px solid var(--line); text-align:left; }
  th { color:var(--muted); font-weight:600; }
  tr:hover { background:#f0f6ff; cursor:pointer; }
  tr.sel { background:#dbeafe; }
  .status { font-size:11px; color:var(--muted); margin-top:6px; }
  .sidebox { font-size:12px; line-height:1.45; background:#f8fafc; border:1px solid var(--line); border-radius:8px; padding:10px; margin-top:8px; }
  .sidebox strong { display:block; margin-bottom:4px; }
  .ok { color:var(--ok); } .bad { color:var(--bad); } .warn { color:var(--warn); }
  .twocol { display:grid; grid-template-columns: 1.2fr .8fr; gap:10px; }
  @media (max-width: 1100px) { .grid, .twocol { grid-template-columns: 1fr; } }
</style>
</head>
<body>
<header>
  <div>
    <h1>Surrogates · Inferencia wPLI</h1>
    <div class="meta">$SUBJECT / $SESSION / $TASK_ID · circular_shift · N=$(store.n_surrogates)</div>
  </div>
  <div class="meta">p_min ≈ $(round(_p_floor(store.n_surrogates); digits=5)) · FDR-$(store.fdr_method) α=$(store.alpha)</div>
</header>

<div class="controls">
  <label>Banda
    <select id="band">$band_opts</select>
  </label>
  <label>Capa heatmap
    <select id="layer">
      <option value="neglogp" selected>−log₁₀(p)</option>
      <option value="q">q-value</option>
      <option value="z">z-score</option>
      <option value="obs">wPLI obs</option>
      <option value="sig">máscara FDR</option>
      <option value="p">p-value</option>
    </select>
  </label>
  <button onclick="refreshAll()">Actualizar</button>
  <button class="secondary" onclick="saveView('overview')">PNG overview</button>
  <button class="secondary" onclick="saveView('contrast')">PNG contraste</button>
  <button class="secondary" onclick="saveView('topo')">PNG topo</button>
  <span class="status" id="st">—</span>
</div>

<div class="grid">
  <div class="panel full">
    <div class="panel-title">A · Overview del test</div>
    <div class="interp" id="interp">—</div>
    <div class="kpis" id="kpis"></div>
    <canvas id="cvBars" height="220"></canvas>
    <div class="status" id="stCache">Cache nulos: —</div>
  </div>

  <div class="panel">
    <div class="panel-title">B · Mapa de inferencia</div>
    <canvas id="cvHeat" height="420"></canvas>
    <div class="status" id="stHeat">Click en celda para seleccionar par</div>
  </div>

  <div class="panel">
    <div class="panel-title">C · Contraste obs vs nula</div>
    <div class="twocol">
      <canvas id="cvContrast" height="320"></canvas>
      <div class="sidebox" id="contrastBox">Selecciona un par…</div>
    </div>
  </div>

  <div class="panel">
    <div class="panel-title">D · Grafo FDR-significativo (10–20)</div>
    <canvas id="cvTopo" height="420"></canvas>
    <div class="status" id="stTopo">—</div>
  </div>

  <div class="panel">
    <div class="panel-title">E · Volcano / explorador</div>
    <canvas id="cvVolcano" height="280"></canvas>
    <div style="max-height:260px;overflow:auto;margin-top:8px">
      <table>
        <thead><tr><th>Par</th><th>wPLI</th><th>z</th><th>p</th><th>q</th><th>sig</th></tr></thead>
        <tbody id="edgeBody"></tbody>
      </table>
    </div>
  </div>
</div>

<script>
const \$ = id => document.getElementById(id);
let overview=null, bandData=null, contrast=null, topo=null, meta=null;
let sel = {a:null, b:null};

function bandValue(){ return \$('band').value; }
function layerValue(){ return \$('layer').value; }

function setStatus(msg, cls){ const el=\$('st'); el.textContent=msg; el.className='status '+(cls||''); }

async function refreshAll(){
  setStatus('Cargando…');
  const band = bandValue();
  const [oRes, bRes, tRes] = await Promise.all([
    fetch('/api/overview'),
    fetch('/api/band?band='+encodeURIComponent(band)+'&layer='+encodeURIComponent(layerValue())),
    fetch('/api/topo?band='+encodeURIComponent(band))
  ]);
  overview = await oRes.json();
  bandData = await bRes.json();
  topo = await tRes.json();
  fillOverview();
  drawBars();
  drawHeat();
  drawVolcano();
  fillTable();
  drawTopo();
  // contrast: keep selection or top edge
  if(!sel.a || !sel.b){
    const top = bandData.edges.slice().sort((x,y)=>y.wpli-x.wpli)[0];
    if(top){ sel.a=top.ch_a; sel.b=top.ch_b; }
  }
  await loadContrast();
  setStatus(band+' · capa '+layerValue()+' · '+
    (bandData.has_null_cache
      ? 'datos reales · cache nulo OK ('+bandData.n_surrogates+' surr)'
      : 'sin cache nulo (approx Normal)'),
    bandData.has_null_cache ? 'ok' : 'warn');
}

function fillOverview(){
  \$('interp').innerHTML = '<strong>Interpretación automática</strong><br>'+(overview.interpretation||'—');
  \$('kpis').innerHTML = [
    ['N surr', overview.n_surrogates],
    ['α FDR', overview.alpha],
    ['p_min', overview.p_min],
    ['pares', overview.n_pairs],
    ['sig total', overview.n_sig_total],
    ['método', overview.method],
    ['seed', overview.seed],
    ['canales', overview.n_channels]
  ].map(([l,v])=>'<div class="kpi"><div class="v">'+v+'</div><div class="l">'+l+'</div></div>').join('');
  const nOk = (overview.bands||[]).filter(b=>b.has_null_cache).length;
  const nB = (overview.bands||[]).length;
  const caches = (overview.bands||[]).map(b=>b.band+':'+(b.has_null_cache?'✓':'·')).join('  ');
  \$('stCache').textContent = 'Cache nulos empíricos: '+nOk+'/'+nB+' — '+caches+
    (nOk===nB ? ' · histograma C usa muestras reales (N surr)' : ' · bandas sin cache → Normal(μ,σ)');
}

function drawBars(){
  const cv=\$('cvBars'), ctx=cv.getContext('2d');
  const W=cv.width=cv.clientWidth*devicePixelRatio, H=cv.height=220*devicePixelRatio;
  ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);
  const w=cv.clientWidth, h=220;
  ctx.clearRect(0,0,w,h);
  const bands=overview.bands||[];
  if(!bands.length) return;
  const pad={l:40,r:12,t:16,b:48};
  const maxPct=Math.max(10, ...bands.map(b=>b.pct_sig));
  const bw=(w-pad.l-pad.r)/bands.length;
  bands.forEach((b,i)=>{
    const x=pad.l+i*bw+bw*0.15;
    const bh=((h-pad.t-pad.b)*b.pct_sig)/maxPct;
    const y=h-pad.b-bh;
    ctx.fillStyle = b.n_sig>0 ? '#1f6feb' : '#94a3b8';
    ctx.fillRect(x,y,bw*0.7,bh);
    ctx.fillStyle='#334155';
    ctx.font='11px sans-serif';
    ctx.textAlign='center';
    ctx.fillText(b.pct_sig.toFixed(1)+'%', x+bw*0.35, y-4);
    ctx.save();
    ctx.translate(x+bw*0.35, h-pad.b+14);
    ctx.rotate(-0.4);
    ctx.fillText(b.band, 0, 0);
    ctx.restore();
  });
  ctx.strokeStyle='#cbd5e1'; ctx.beginPath();
  ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, h-pad.b); ctx.lineTo(w-pad.r, h-pad.b); ctx.stroke();
}

function colorScale(t, layer){
  t=Math.max(0,Math.min(1,t));
  if(layer==='sig') return t>0.5 ? 'rgb(220,60,50)' : 'rgb(230,235,240)';
  if(layer==='z'){
    // diverging around 0 handled by vmin/vmax abs
    const r=Math.round(40+180*t), b=Math.round(200-160*t);
    return 'rgb('+r+',80,'+b+')';
  }
  // sequential blue-yellow-red
  const r=Math.round(30+220*t), g=Math.round(80+100*(1-Math.abs(t-0.5)*2)), bl=Math.round(200-180*t);
  return 'rgb('+r+','+g+','+bl+')';
}

function drawHeat(){
  const cv=\$('cvHeat'), ctx=cv.getContext('2d');
  const n=bandData.n, ch=bandData.channels, M=bandData.matrix;
  const W=cv.width=cv.clientWidth*devicePixelRatio, H=cv.height=420*devicePixelRatio;
  ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);
  const w=cv.clientWidth, h=420;
  ctx.clearRect(0,0,w,h);
  const pad=56, size=Math.min(w,h)-pad-20;
  const cell=size/n;
  const x0=(w-size)/2, y0=20;
  const vmin=bandData.vmin, vmax=bandData.vmax, layer=bandData.layer;
  for(let j=0;j<n;j++) for(let i=0;i<n;i++){
    let v=M[j*n+i];
    if(i===j){ ctx.fillStyle='#e8edf2'; }
    else {
      let t = vmax>vmin ? (v-vmin)/(vmax-vmin) : 0;
      if(layer==='z') t = (v+Math.abs(vmax))/(2*Math.abs(vmax)+1e-9);
      ctx.fillStyle=colorScale(t, layer);
    }
    ctx.fillRect(x0+j*cell, y0+i*cell, cell+0.5, cell+0.5);
  }
  // selection
  if(sel.a && sel.b){
    const ia=ch.indexOf(sel.a), ib=ch.indexOf(sel.b);
    if(ia>=0&&ib>=0){
      ctx.strokeStyle='#111'; ctx.lineWidth=2;
      ctx.strokeRect(x0+ib*cell, y0+ia*cell, cell, cell);
      ctx.strokeRect(x0+ia*cell, y0+ib*cell, cell, cell);
    }
  }
  ctx.fillStyle='#475569'; ctx.font='9px sans-serif';
  for(let i=0;i<n;i++){
    ctx.save();
    ctx.translate(x0+i*cell+cell/2, y0+size+10);
    ctx.rotate(-Math.PI/3);
    ctx.textAlign='right'; ctx.fillText(ch[i],0,0);
    ctx.restore();
    ctx.textAlign='right';
    ctx.fillText(ch[i], x0-4, y0+i*cell+cell*0.7);
  }
  cv.onclick = (ev)=>{
    const rect=cv.getBoundingClientRect();
    const x=(ev.clientX-rect.left), y=(ev.clientY-rect.top);
    const j=Math.floor((x-x0)/cell), i=Math.floor((y-y0)/cell);
    if(i<0||j<0||i>=n||j>=n||i===j) return;
    sel.a=ch[i]; sel.b=ch[j];
    loadContrast(); drawHeat(); fillTable();
  };
}

function drawVolcano(){
  const cv=\$('cvVolcano'), ctx=cv.getContext('2d');
  const W=cv.width=cv.clientWidth*devicePixelRatio, H=cv.height=280*devicePixelRatio;
  ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);
  const w=cv.clientWidth, h=280;
  ctx.clearRect(0,0,w,h);
  const edges=bandData.edges||[];
  const pad={l:44,r:12,t:16,b:36};
  let maxZ=1, maxNLP=1;
  edges.forEach(e=>{ maxZ=Math.max(maxZ, Math.abs(e.z)); maxNLP=Math.max(maxNLP, -Math.log10(Math.max(e.p,1e-12))); });
  function X(z){ return pad.l + (z+maxZ)/(2*maxZ)*(w-pad.l-pad.r); }
  function Y(nlp){ return h-pad.b - nlp/maxNLP*(h-pad.t-pad.b); }
  ctx.strokeStyle='#cbd5e1'; ctx.beginPath();
  ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,h-pad.b); ctx.lineTo(w-pad.r,h-pad.b); ctx.stroke();
  // alpha line approx for p=0.05
  const y05=Y(-Math.log10(0.05));
  ctx.strokeStyle='#94a3b8'; ctx.setLineDash([4,3]);
  ctx.beginPath(); ctx.moveTo(pad.l,y05); ctx.lineTo(w-pad.r,y05); ctx.stroke(); ctx.setLineDash([]);
  edges.forEach(e=>{
    const nlp=-Math.log10(Math.max(e.p,1e-12));
    ctx.fillStyle = e.sig ? '#dc2626' : '#93c5fd';
    ctx.beginPath(); ctx.arc(X(e.z), Y(nlp), e.sig?2.4:1.6, 0, Math.PI*2); ctx.fill();
  });
  ctx.fillStyle='#64748b'; ctx.font='11px sans-serif';
  ctx.fillText('z-score', w/2-20, h-10);
  ctx.save(); ctx.translate(14, h/2); ctx.rotate(-Math.PI/2); ctx.fillText('−log10(p)',0,0); ctx.restore();
}

function fillTable(){
  const edges=(bandData.edges||[]).slice().sort((a,b)=> (b.sig-a.sig) || (a.q-b.q) || (b.z-a.z));
  const top=edges.slice(0,80);
  \$('edgeBody').innerHTML = top.map(e=>{
    const selCls = (sel.a&&sel.b&&((e.ch_a===sel.a&&e.ch_b===sel.b)||(e.ch_a===sel.b&&e.ch_b===sel.a)))?' class="sel"':'';
    return '<tr data-a="'+e.ch_a+'" data-b="'+e.ch_b+'"'+selCls+'>'+
      '<td>'+e.ch_a+'–'+e.ch_b+'</td><td>'+e.wpli.toFixed(3)+'</td><td>'+e.z.toFixed(2)+'</td>'+
      '<td>'+e.p.toFixed(4)+'</td><td>'+e.q.toFixed(4)+'</td>'+
      '<td class="'+(e.sig?'ok':'bad')+'">'+(e.sig?'sí':'no')+'</td></tr>';
  }).join('');
  \$('edgeBody').querySelectorAll('tr').forEach(tr=>{
    tr.onclick=()=>{ sel.a=tr.dataset.a; sel.b=tr.dataset.b; loadContrast(); drawHeat(); fillTable(); };
  });
}

async function loadContrast(){
  if(!sel.a||!sel.b) return;
  const band=bandValue();
  const res=await fetch('/api/contrast?band='+encodeURIComponent(band)+'&a='+encodeURIComponent(sel.a)+'&b='+encodeURIComponent(sel.b));
  contrast=await res.json();
  drawContrast();
  const cls = contrast.sig ? 'ok' : 'bad';
  const srcLabel = contrast.source==='empirical'
    ? '<span class="ok">histograma empírico real (cache/surrogate · N='+contrast.n_surrogates+')</span>'
    : '<span class="warn">aproximación Normal(μ,σ) — falta null_distribution_*.jls</span>';
  \$('contrastBox').innerHTML =
    '<strong>'+contrast.band+' · '+contrast.ch_a+' ↔ '+contrast.ch_b+'</strong>'+
    'Fuente: '+srcLabel+'<br><br>'+
    'H₀: no hay acoplamiento de fase más allá del azar (circular-shift).<br>'+
    'p = proporción de surrogates con wPLI ≥ observado (+ corrección Monte Carlo).<br><br>'+
    'wPLI obs = <b>'+contrast.wpli_obs.toFixed(4)+'</b><br>'+
    'μ nulo = '+contrast.null_mean.toFixed(4)+' · σ = '+contrast.null_std.toFixed(4)+'<br>'+
    'p = <b>'+contrast.p.toFixed(4)+'</b> · q = <b>'+contrast.q.toFixed(4)+'</b> · z = <b>'+contrast.z.toFixed(2)+'</b><br>'+
    'N = '+contrast.n_surrogates+' · p_min = '+contrast.p_min+'<br><br>'+
    '<span class="'+cls+'"><b>Resultado:</b> '+contrast.verdict+'</span>';
}

function drawContrast(){
  const cv=\$('cvContrast'), ctx=cv.getContext('2d');
  if(!contrast) return;
  const W=cv.width=cv.clientWidth*devicePixelRatio, H=cv.height=320*devicePixelRatio;
  ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);
  const w=cv.clientWidth, h=320;
  ctx.clearRect(0,0,w,h);
  const pad={l:40,r:14,t:18,b:36};
  const hist=contrast.hist, counts=hist.counts, nb=hist.nbins;
  const lo=hist.lo, hi=hist.hi;
  const maxC=Math.max(1,...counts);
  const bw=(w-pad.l-pad.r)/nb;
  counts.forEach((c,i)=>{
    const bh=(h-pad.t-pad.b)*c/maxC;
    ctx.fillStyle='rgba(37,99,235,0.55)';
    ctx.fillRect(pad.l+i*bw, h-pad.b-bh, Math.max(1,bw-1), bh);
  });
  function X(v){ return pad.l + (v-lo)/Math.max(hi-lo,1e-9)*(w-pad.l-pad.r); }
  // mean null
  ctx.strokeStyle='#1d4ed8'; ctx.setLineDash([5,3]); ctx.lineWidth=2;
  ctx.beginPath(); ctx.moveTo(X(contrast.null_mean), pad.t); ctx.lineTo(X(contrast.null_mean), h-pad.b); ctx.stroke();
  // obs
  ctx.strokeStyle='#ea580c'; ctx.setLineDash([]); ctx.lineWidth=2.5;
  ctx.beginPath(); ctx.moveTo(X(contrast.wpli_obs), pad.t); ctx.lineTo(X(contrast.wpli_obs), h-pad.b); ctx.stroke();
  ctx.fillStyle='#334155'; ctx.font='11px sans-serif';
  ctx.fillText('media nula', X(contrast.null_mean)+4, pad.t+12);
  ctx.fillStyle='#c2410c';
  ctx.fillText('obs '+contrast.wpli_obs.toFixed(3), Math.min(w-90, X(contrast.wpli_obs)+4), pad.t+28);
  ctx.strokeStyle='#cbd5e1'; ctx.beginPath();
  ctx.moveTo(pad.l,pad.t); ctx.lineTo(pad.l,h-pad.b); ctx.lineTo(w-pad.r,h-pad.b); ctx.stroke();
  ctx.fillStyle='#64748b'; ctx.fillText('wPLI', w/2-10, h-10);
}

function drawTopo(){
  const cv=\$('cvTopo'), ctx=cv.getContext('2d');
  const W=cv.width=cv.clientWidth*devicePixelRatio, H=cv.height=420*devicePixelRatio;
  ctx.setTransform(devicePixelRatio,0,0,devicePixelRatio,0,0);
  const w=cv.clientWidth, h=420;
  ctx.clearRect(0,0,w,h);
  const cx=w/2, cy=h/2, R=Math.min(w,h)*0.38;
  ctx.strokeStyle='#94a3b8'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.arc(cx,cy,R,0,Math.PI*2); ctx.stroke();
  ctx.beginPath(); ctx.arc(cx,cy,R*0.08,0,Math.PI*2); ctx.stroke();
  function XY(n){ return [cx+n.x*R, cy-n.y*R]; }
  const nodes=topo.nodes||[], edges=topo.edges||[];
  const dmax=Math.max(1, ...nodes.map(n=>n.degree));
  edges.forEach(e=>{
    const na=nodes.find(n=>n.channel===e.ch_a), nb=nodes.find(n=>n.channel===e.ch_b);
    if(!na||!nb) return;
    const [x1,y1]=XY(na), [x2,y2]=XY(nb);
    ctx.strokeStyle='rgba(30,100,220,'+(0.2+0.6*Math.min(e.wpli,1))+')';
    ctx.lineWidth=0.7+2.2*Math.min(e.wpli,1);
    ctx.beginPath(); ctx.moveTo(x1,y1); ctx.lineTo(x2,y2); ctx.stroke();
  });
  nodes.forEach(n=>{
    const [x,y]=XY(n);
    const r=4+10*(n.degree/dmax);
    ctx.fillStyle = n.degree>0 ? '#ef4444' : '#cbd5e1';
    ctx.beginPath(); ctx.arc(x,y,r,0,Math.PI*2); ctx.fill();
    ctx.fillStyle='#1e293b'; ctx.font='10px sans-serif'; ctx.textAlign='center';
    ctx.fillText(n.channel, x, y-r-3);
  });
  \$('stTopo').textContent = topo.n_sig+' aristas FDR-sig · hubs: '+
    (topo.hubs||[]).slice(0,5).map(h=>h.channel+'('+h.degree+')').join(', ');
  cv.onclick=(ev)=>{
    const rect=cv.getBoundingClientRect();
    const x=ev.clientX-rect.left, y=ev.clientY-rect.top;
    let best=null, bd=1e9;
    nodes.forEach(n=>{
      const [nx,ny]=XY(n); const d=(x-nx)**2+(y-ny)**2;
      if(d<bd && d<14*14){ bd=d; best=n.channel; }
    });
    if(!best) return;
    // pick strongest sig edge involving node
    const cand=edges.filter(e=>e.ch_a===best||e.ch_b===best).sort((a,b)=>b.wpli-a.wpli)[0];
    if(!cand) return;
    sel.a=cand.ch_a; sel.b=cand.ch_b; loadContrast(); drawHeat(); fillTable();
  };
}

async function saveView(view){
  setStatus('Guardando PNG…');
  const body={view, band:bandValue(), a:sel.a||'', b:sel.b||''};
  const res=await fetch('/api/save',{method:'POST', headers:{'Content-Type':'application/json'}, body:JSON.stringify(body)});
  const j=await res.json();
  if(j.ok) setStatus('PNG → '+j.file, 'ok');
  else setStatus(j.error||'Error', 'bad');
}

\$('band').onchange = ()=>{ sel.a=null; sel.b=null; refreshAll(); };
\$('layer').onchange = refreshAll;

(async function init(){
  meta = await (await fetch('/api/meta')).json();
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
    println("Cargando surrogates desde $SURR_DIR …")
    store = load_store()
    n_cache = count(b -> store.packs[b].null_dist !== nothing, store.bands)
    println("  $(length(store.channels)) canales · bandas: ", join(store.bands, ", "))
    println("  N=$(store.n_surrogates) · α=$(store.alpha) · FDR=$(store.fdr_method) · sig_total=$(nrow(store.significant))")
    println("  cache nulos: $n_cache/$(length(store.bands)) bandas" *
            (n_cache == 0 ? "  (aún no generado — hist usará Normal(μ,σ))" : ""))

    server, port = _listen_available(HOST, PORT)
    url = "http://$HOST:$port/"
    println()
    println("UI surrogates → $url")
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

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
