# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Análisis longitudinal profesional (T1 → T2)
# ═══════════════════════════════════════════════════════════════
#
#  Diseño experimental (Fig. 3.1):
#    · Solo pacientes EM con par completo T1+T2 (N diseño = 30)
#    · Controles NO entran; pérdidas de seguimiento (sin T2) se excluyen
#      (no hay comparación posible — no se analizan)
#    · EC y EO en paralelo: mismo contraste T1→T2, sin pooling
#
#  Comparación intra-sujeto T1→T2 sobre resultados del pipeline:
#    · Inclusión por condición + QC (qc_decision_table)
#    · wPLI edge-wise: Wilcoxon signed-rank + FDR-BH + Cohen dz
#    · Band power canal×banda (T1→T2)
#    · Métricas de red por banda (strength/degree)
#    · Figuras de publicación (heatmaps, red sig., paired means, topo)
#
#  Salida (compatible Fase 14 dashboard):
#    results/longitudinal/{EC|EO}/
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/run_longitudinal_analysis.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 25-07-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/<este-script>.jl …
#
#  Prerrequisito
#  ─────────────
#    results/subjects/… con wPLI T1 y T2
#    data/bids/longitudinal_pairs.csv  (include_longitudinal=true)
#    results/qc/qc_decision_table.csv  (recomendado)
#
#  Config
#  ──────
#    Argumento posicional opcional; defecto: config/pipeline.toml
#
#  Uso
#  ───
#    julia --project=. scripts/run_longitudinal_analysis.jl
#    julia --project=. scripts/run_longitudinal_analysis.jl config/pipeline.toml
#
#  Salida
#  ──────
#    results/longitudinal/{EC|EO}/
#      longitudinal_*.csv · paired_subjects.csv · longitudinal_summary.json
#      tables/spectral/ · tables/network/ · figures/ · config_snapshot.toml

using CSV, DataFrames, Statistics, LinearAlgebra, Dates, TOML, Printf, CairoMakie

include(joinpath(@__DIR__, "..", "src", "viz", "GroupVizCommon.jl"))
using .GroupVizCommon

const PROJ     = dirname(@__DIR__)
const CONFIG_P = length(ARGS) > 0 ? ARGS[1] :
                 joinpath(PROJ, "config", "pipeline.toml")

const QC_ALLOWED = Set(["include", "include_with_warning"])
const N_PAIRED_DESIGN = 30   # Fig. 3.1 — pares EM T1+T2 comparables

cfg_raw   = TOML.parsefile(CONFIG_P)
paths_raw = get(cfg_raw, "paths", Dict{String,Any}())
res_root  = let r = get(paths_raw, "results", "results")
                isabspath(r) ? r : joinpath(PROJ, r)
            end
bids_root = let b = get(paths_raw, "bids_root", "data/bids")
                p = isabspath(b) ? b : joinpath(PROJ, b)
                # Tolerar data/BIDS vs data/bids en macOS/Linux
                isdir(p) ? p : (isdir(joinpath(PROJ, "data", "BIDS")) ?
                    joinpath(PROJ, "data", "BIDS") : p)
            end
bands_cfg = get(cfg_raw, "bands", Dict(
    "DELTA"     => [0.5, 4.0], "THETA" => [4.0, 8.0], "ALPHA" => [7.8, 11.7],
    "BETA_LOW"  => [12.0, 15.0], "BETA_MID" => [15.0, 18.0],
    "BETA_HIGH" => [18.0, 30.0], "GAMMA" => [30.0, 50.0],
))
BANDS = sort(collect(keys(bands_cfg)))

conn_cfg  = get(cfg_raw, "connectivity", Dict{String,Any}())
graph_cfg = get(cfg_raw, "graph", Dict{String,Any}())
WPLI_METHOD = String(get(conn_cfg, "wpli_method", "hilbert"))
USE_DWPLI   = Bool(get(conn_cfg, "use_dwpli", false))
GRAPH_DENS  = Float64(get(graph_cfg, "density", 0.1))
GRAPH_METH  = String(get(graph_cfg, "threshold_method", "proportional"))
N_CH_MONT   = Int(get(get(cfg_raw, "montage", Dict()), "n_channels_analysis", 31))

# ─── Helpers básicos ──────────────────────────────────────────

function norm_cond(c::AbstractString)::String
    lc = lowercase(String(c))
    lc in ("ec", "eyesclosed") && return "eyesclosed"
    lc in ("eo", "eyesopen")   && return "eyesopen"
    return lc
end

function cond_code(c::AbstractString)::String
    nc = norm_cond(c)
    nc == "eyesclosed" && return "EC"
    nc == "eyesopen"   && return "EO"
    return uppercase(String(c))
end

function is_t1_session(s::AbstractString)::Bool
    uppercase(String(s)) in ("T1", "BASELINE", "BL", "V1", "VISIT1", "S1", "PRE")
end

"""True si el ID parece control (MC*, C*, HC*), no paciente EM."""
function is_control_id(sid::AbstractString)::Bool
    s = uppercase(String(sid))
    startswith(s, "MC") && return true
    startswith(s, "HC") && return true
    startswith(s, "C") && !startswith(s, "M") && return true
    return false
end

function is_ms_id(sid::AbstractString)::Bool
    !is_control_id(sid)
end

function export_dir(subj::AbstractString, sess::AbstractString, cond::AbstractString)::String
    joinpath(res_root, "subjects", "sub-$(subj)", "ses-$(sess)", norm_cond(cond))
end

function bh_qvalues(p::Vector{Float64})::Vector{Float64}
    m = length(p); m == 0 && return Float64[]
    ord = sortperm(p); rnk = invperm(ord)
    q   = p .* m ./ rnk
    qs  = q[ord]
    for i in (m - 1):-1:1; qs[i] = min(qs[i], qs[i + 1]); end
    qo = zeros(m); qo[ord] = qs
    return min.(qo, 1.0)
end

function _erf_approx(x::Float64)::Float64
    t = 1.0 / (1.0 + 0.3275911 * abs(x))
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 +
               t * (-1.453152027 + t * 1.061405429))))
    sign(x) * (1.0 - poly * exp(-x * x))
end

function _norm_cdf(z::Float64)::Float64
    0.5 * (1.0 + _erf_approx(z / sqrt(2.0)))
end

function _assign_ranks(v::Vector{Float64})::Vector{Float64}
    n = length(v); order = sortperm(v); ranks = zeros(Float64, n)
    i = 1
    while i <= n
        j = i
        while j < n && v[order[j + 1]] == v[order[i]]; j += 1; end
        r_avg = mean(Float64(i):Float64(j))
        for k in i:j; ranks[order[k]] = r_avg; end
        i = j + 1
    end
    return ranks
end

"""Wilcoxon signed-rank → (p, effect_r)."""
function wilcoxon_p(before::Vector{Float64}, after::Vector{Float64})
    length(before) == length(after) || return (1.0, 0.0)
    diffs = filter(!=(0.0), after .- before)
    nd = length(diffs)
    nd < 2 && return (1.0, 0.0)
    ranks = _assign_ranks(abs.(diffs))
    W_plus  = sum(ranks[diffs .> 0])
    W_minus = sum(ranks[diffs .< 0])
    W = min(W_plus, W_minus)
    μW = nd * (nd + 1) / 4
    σW = sqrt(nd * (nd + 1) * (2nd + 1) / 24)
    σW < 1e-12 && return (1.0, 0.0)
    z = (W - μW) / σW
    p = clamp(2.0 * (1.0 - _norm_cdf(abs(z))), 0.0, 1.0)
    r = abs(z) / sqrt(nd)
    return (p, r)
end

"""Paired t → (p, cohen_dz)."""
function paired_t(before::Vector{Float64}, after::Vector{Float64})
    length(before) == length(after) || return (1.0, 0.0)
    diffs = after .- before
    n = length(diffs)
    n < 2 && return (1.0, 0.0)
    μ = mean(diffs); s = std(diffs)
    s < 1e-12 && return (1.0, 0.0)
    t = μ * sqrt(n) / s
    p = clamp(2.0 * (1.0 - _norm_cdf(abs(t))), 0.0, 1.0)
    return (p, μ / s)
end

function cohen_dz(before::Vector{Float64}, after::Vector{Float64})::Float64
    length(before) == length(after) || return 0.0
    diffs = after .- before
    n = length(diffs); n < 2 && return 0.0
    s = std(diffs); s < 1e-12 && return 0.0
    return mean(diffs) / s
end

function upper_mean(W::Matrix{Float64})::Float64
    n = min(size(W)...); n < 2 && return NaN
    return mean(W[i, j] for i in 1:n for j in (i + 1):n)
end

function realign_matrix(ch::Vector{String}, W::Matrix{Float64},
                        common_ch::Vector{String})::Matrix{Float64}
    length(ch) == size(W, 1) == size(W, 2) ||
        error("Matriz inconsistente: $(length(ch)) ch vs $(size(W))")
    idx = [findfirst(==(c), ch) for c in common_ch]
    any(isnothing, idx) && error("Canales no encontrados")
    return W[idx, idx]
end

function save_mat_csv(path::String, W::Matrix{Float64}, ch::Vector{String})
    df = DataFrame(hcat(ch, W), vcat(["channel"], ch))
    CSV.write(path, df)
end

# ─── Carga de datos por sujeto ────────────────────────────────

function load_wpli(subj_id, sess_id, cond, band)
    path = joinpath(export_dir(subj_id, sess_id, cond),
                    "tables", "connectivity", "wpli_$(band).csv")
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    isempty(df) && return nothing
    row_ch = string.(df[!, 1])
    col_syms = names(df)[2:end]
    col_ch = string.(col_syms)
    common = unique(filter(c -> c in Set(col_ch), row_ch))
    n = length(common); n < 2 && return nothing
    row_idx = [findfirst(==(c), row_ch) for c in common]
    any(isnothing, row_idx) && return nothing
    W = Matrix{Float64}(undef, n, n)
    for (j, c) in enumerate(common)
        col_sym = col_syms[findfirst(==(c), col_ch)]
        hasproperty(df, col_sym) || return nothing
        W[:, j] = Float64.(df[row_idx, col_sym])
    end
    return (common, W)
end

"""Load band_power_summary.csv → Dict(channel => Dict(band => power))."""
function load_band_power(subj_id, sess_id, cond)
    path = joinpath(export_dir(subj_id, sess_id, cond),
                    "tables", "band_power_summary.csv")
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    isempty(df) && return nothing
    ch_col = names(df)[1]
    out = Dict{String, Dict{String,Float64}}()
    for row in eachrow(df)
        ch = string(row[ch_col])
        out[ch] = Dict{String,Float64}()
        for b in BANDS
            if hasproperty(row, Symbol(b))
                out[ch][b] = Float64(row[Symbol(b)])
            end
        end
    end
    return out
end

function load_seg_epochs(subj_id, sess_id, cond)
    path = joinpath(export_dir(subj_id, sess_id, cond),
                    "json", "segmentation_summary.json")
    n_valid = missing; n_total = missing; pct = missing
    if isfile(path)
        txt = read(path, String)
        m = match(r"\"n_valid\"\s*:\s*(\d+)", txt)
        m !== nothing && (n_valid = parse(Int, m.captures[1]))
        m = match(r"\"n_total\"\s*:\s*(\d+)", txt)
        m !== nothing && (n_total = parse(Int, m.captures[1]))
        m = match(r"\"retention_pct\"\s*:\s*([0-9.]+)", txt)
        m !== nothing && (pct = parse(Float64, m.captures[1]))
    end
    return (n_valid, n_total, pct)
end

# ─── QC ───────────────────────────────────────────────────────

function load_qc_table()::Dict{Tuple{String,String,String}, String}
    path = joinpath(res_root, "qc", "qc_decision_table.csv")
    out = Dict{Tuple{String,String,String}, String}()
    isfile(path) || return out
    df = CSV.read(path, DataFrame)
    rename!(df, Dict(n => Symbol(lowercase(string(n))) for n in names(df)))
    sid_col = hasproperty(df, :subject_id) ? :subject_id :
              (hasproperty(df, :bids_id) ? :bids_id : nothing)
    sid_col === nothing && return out
    for row in eachrow(df)
        sid  = string(row[sid_col])
        sess = string(hasproperty(row, :session_id) ? row.session_id : row.session)
        cond = cond_code(string(row.condition))
        dec  = lowercase(string(row.final_decision))
        out[(sid, sess, cond)] = dec
    end
    return out
end

function qc_decision(qc::Dict, sid::String, sess::String, cond::String)::String
    cc = cond_code(cond)
    haskey(qc, (sid, sess, cc)) && return qc[(sid, sess, cc)]
    # Alias M5 ↔ M05 (subject_id vs bids_id)
    if startswith(sid, "M") && !startswith(sid, "MC")
        alt = occursin(r"^M\d$", sid) ? "M0" * sid[2:end] :
              (occursin(r"^M0\d$", sid) ? "M" * sid[3:end] : "")
        !isempty(alt) && haskey(qc, (alt, sess, cc)) && return qc[(alt, sess, cc)]
    end
    return "missing"
end

function qc_ok(dec::String)::Bool
    dec == "missing" && return true  # sin tabla QC → no bloquear
    return dec in QC_ALLOWED
end

# ─── Network helpers (por banda, desde W) ─────────────────────

function _prop_threshold(W::Matrix{Float64}, density::Float64)::Float64
    vals = sort([W[i, j] for i in 1:size(W, 1) for j in (i + 1):size(W, 1)
                 if W[i, j] > 0.0]; rev=true)
    isempty(vals) && return 0.0
    n_keep = max(1, round(Int, density * length(vals)))
    return vals[min(n_keep, length(vals))]
end

function nodal_strength_degree(W::Matrix{Float64}; density::Float64=GRAPH_DENS)
    n = size(W, 1)
    strength = [sum(W[i, j] for j in 1:n if j != i) for i in 1:n]
    thr = _prop_threshold(W, density)
    degree = [count(j -> j != i && W[i, j] > thr, 1:n) for i in 1:n]
    smax = maximum(strength); smax < 1e-12 && (smax = 1.0)
    norm_s = strength ./ smax
    return strength, degree, norm_s, thr
end

# ─── Electrodes / topo ────────────────────────────────────────

function load_electrode_xy(subj_id, sess_id)
    elec_dir = joinpath(bids_root, "electrodes")
    path = joinpath(elec_dir, "sub-$(subj_id)_ses-$(sess_id)_electrodes.tsv")
    isfile(path) || return Dict{String,Tuple{Float64,Float64}}()
    df = CSV.read(path, DataFrame; delim='\t')
    rename!(df, Dict(n => Symbol(lowercase(string(n))) for n in names(df)))
    out = Dict{String,Tuple{Float64,Float64}}()
    for row in eachrow(df)
        hasproperty(row, :name) || continue
        x = hasproperty(row, :x) ? Float64(row.x) : 0.0
        y = hasproperty(row, :y) ? Float64(row.y) : 0.0
        out[string(row.name)] = (x, y)
    end
    return out
end

# ─── Figuras (GroupVizCommon: RdBu_r, topo, bandas fisiológicas) ─

_save_heatmap(path, W, ch, title; cmap=:viridis, diverging=false, colorrange=nothing) =
    save_heatmap(path, W, ch, title; cmap=cmap, diverging=diverging, colorrange=colorrange,
                 colorbar_label=diverging ? "Δ wPLI" : "wPLI")

_save_sig_network(path, ch, sig_rows, title) =
    save_topo_network(path, ch, sig_rows, title)

_save_paired_means(path, subject_means_df) =
    save_paired_means(path, subject_means_df; band_order=BAND_ORDER)

_save_topo_delta(path, ch, delta_vals, xy, title) =
    save_topo_delta(path, ch, delta_vals, xy, title; colorbar_label="Δ power (T2−T1)")

function _save_explore_or_sig_network(figs_dir, band, cond, common_ch, sig_rows, stats_df)
    if !isempty(sig_rows)
        save_topo_network(joinpath(figs_dir, "sig_network_$(band).png"),
                          common_ch, sig_rows, "Edges sig. FDR — $band ($cond)")
    elseif nrow(stats_df) > 0 && hasproperty(stats_df, :effect_d)
        top = first(sort(stats_df, :effect_d; by=abs, rev=true), min(20, nrow(stats_df)))
        rows = [NamedTuple(r) for r in eachrow(top)]
        save_topo_network(joinpath(figs_dir, "explore_network_topN_$(band).png"),
                          common_ch, rows,
                          "Top-20 |dz| (exploratorio) — $band ($cond)")
    end
    return nothing
end

# ─── Pares longitudinales ─────────────────────────────────────

struct SubjPair
    subject_id::String
    session_t1::String
    session_t2::String
    has_ec::Bool
    has_eo::Bool
end

function load_pairs()::Vector{SubjPair}
    # Buscar CSV en bids_root (y variante mayúscula)
    candidates = [
        joinpath(bids_root, "longitudinal_pairs.csv"),
        joinpath(PROJ, "data", "bids", "longitudinal_pairs.csv"),
        joinpath(PROJ, "data", "BIDS", "longitudinal_pairs.csv"),
    ]
    pairs_path = ""
    for c in candidates
        if isfile(c); pairs_path = c; break; end
    end

    all_pairs = SubjPair[]
    n_skipped_not_incl = 0
    n_skipped_ctrl = 0

    if !isempty(pairs_path)
        pdf = CSV.read(pairs_path, DataFrame)
        rename!(pdf, Dict(n => Symbol(lowercase(string(n))) for n in names(pdf)))
        has_new = hasproperty(pdf, :has_t1_ec) || hasproperty(pdf, :include_longitudinal)
        has_bids = hasproperty(pdf, :bids_id)
        has_sess = hasproperty(pdf, :session_t1) && hasproperty(pdf, :session_t2)
        has_incl_col = hasproperty(pdf, :include_longitudinal)

        for row in eachrow(pdf)
            sid = has_bids ? string(row.bids_id) : string(row.subject_id)
            # Diseño: solo EM; controles fuera
            if !is_ms_id(sid)
                n_skipped_ctrl += 1
                continue
            end
            # Diseño: solo pares comparables (include_longitudinal=true)
            if has_incl_col && !_as_bool(row.include_longitudinal)
                n_skipped_not_incl += 1
                continue
            end
            sess_t1 = has_sess ? string(row.session_t1) : "T1"
            sess_t2 = has_sess ? string(row.session_t2) : "T2"
            if has_new
                ht1ec = hasproperty(row, :has_t1_ec) && _as_bool(row.has_t1_ec)
                ht2ec = hasproperty(row, :has_t2_ec) && _as_bool(row.has_t2_ec)
                ht1eo = hasproperty(row, :has_t1_eo) && _as_bool(row.has_t1_eo)
                ht2eo = hasproperty(row, :has_t2_eo) && _as_bool(row.has_t2_eo)
                # Sin T2 en ninguna condición → no comparable
                if !(ht1ec && ht2ec) && !(ht1eo && ht2eo)
                    n_skipped_not_incl += 1
                    continue
                end
                push!(all_pairs, SubjPair(sid, sess_t1, sess_t2, ht1ec && ht2ec, ht1eo && ht2eo))
            else
                push!(all_pairs, SubjPair(sid, sess_t1, sess_t2, true, true))
            end
        end
        println("📋 Pares leídos desde: $pairs_path")
        println("   Comparables (EM, include_longitudinal): $(length(all_pairs))  |  N diseño=$N_PAIRED_DESIGN")
        n_skipped_not_incl > 0 && println("   Omitidos (sin T2 / include=false): $n_skipped_not_incl  — no se analizan")
        n_skipped_ctrl > 0 && println("   Omitidos (controles): $n_skipped_ctrl")
    else
        subj_root = joinpath(res_root, "subjects")
        if isdir(subj_root)
            for sd in readdir(subj_root)
                startswith(sd, "sub-") || continue
                subj_id = sd[5:end]
                is_ms_id(subj_id) || continue
                sess_dir = joinpath(subj_root, sd)
                sessions = filter(s -> startswith(s, "ses-"), readdir(sess_dir))
                t1s = filter(s -> is_t1_session(s[5:end]), sessions)
                t2s = filter(s -> !is_t1_session(s[5:end]), sessions)
                (isempty(t1s) || isempty(t2s)) && continue
                for t1s_ in t1s, t2s_ in t2s
                    has_ec = isdir(joinpath(sess_dir, t1s_, "eyesclosed")) &&
                             isdir(joinpath(sess_dir, t2s_, "eyesclosed"))
                    has_eo = isdir(joinpath(sess_dir, t1s_, "eyesopen")) &&
                             isdir(joinpath(sess_dir, t2s_, "eyesopen"))
                    (has_ec || has_eo) || continue
                    push!(all_pairs, SubjPair(subj_id, t1s_[5:end], t2s_[5:end], has_ec, has_eo))
                end
            end
        end
        println("🔍 Pares auto-detectados (EM con T1+T2): $(length(all_pairs))  |  N diseño=$N_PAIRED_DESIGN")
    end
    return all_pairs
end

_as_bool(x) = x isa Bool ? x : (lowercase(string(x)) == "true")

# ─── Inclusión por condición ──────────────────────────────────

function evaluate_inclusion(sp::SubjPair, cond::String, qc::Dict)
    cc = cond_code(cond)
    eligible = cc == "EC" ? sp.has_ec : sp.has_eo
    if !eligible
        return (false, "Sin par T1/T2 para $cc en longitudinal_pairs",
                "n/a", "n/a", missing, missing, missing, missing)
    end

    qc_t1 = qc_decision(qc, sp.subject_id, sp.session_t1, cc)
    qc_t2 = qc_decision(qc, sp.subject_id, sp.session_t2, cc)
    ev_t1 = load_seg_epochs(sp.subject_id, sp.session_t1, cc)
    ev_t2 = load_seg_epochs(sp.subject_id, sp.session_t2, cc)

    # ¿Hay al menos una banda con wPLI en ambas sesiones?
    n_bands = 0
    for b in BANDS
        r1 = load_wpli(sp.subject_id, sp.session_t1, cc, b)
        r2 = load_wpli(sp.subject_id, sp.session_t2, cc, b)
        (r1 !== nothing && r2 !== nothing) && (n_bands += 1)
    end

    reasons = String[]
    n_bands == 0 && push!(reasons, "Sin wPLI pareado en disco")
    !qc_ok(qc_t1) && push!(reasons, "QC T1=$(qc_t1)")
    !qc_ok(qc_t2) && push!(reasons, "QC T2=$(qc_t2)")

    included = isempty(reasons)
    reason = included ? "" : join(reasons, "; ")
    return (included, reason, qc_t1, qc_t2, ev_t1[1], ev_t2[1], ev_t1[3], ev_t2[3])
end

# ═══════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════

println("=" ^ 62)
println(" NeuroMIND — Análisis longitudinal T1→T2 (solo EM, N diseño=$N_PAIRED_DESIGN)")
println(" $(now())")
println(" Diseño: pares completos; EC y EO en paralelo (sin pooling)")
println(" wPLI=$WPLI_METHOD  dwPLI=$USE_DWPLI  test=wilcoxon  density=$GRAPH_DENS")
println("=" ^ 62)

all_pairs = load_pairs()
if isempty(all_pairs)
    @warn "No se encontraron pares T1/T2."
    exit(1)
end

qc_table = load_qc_table()
println("  QC decisions cargadas: $(length(qc_table))")
println()

for cond in ["EC", "EO"]
    println("── Condición: $cond  (contraste longitudinal T1→T2) " * "─"^20)
    out_dir = joinpath(res_root, "longitudinal", cond)
    figs_dir = joinpath(out_dir, "figures")
    spec_dir = joinpath(out_dir, "tables", "spectral")
    net_dir  = joinpath(out_dir, "tables", "network")
    mkpath(figs_dir); mkpath(spec_dir); mkpath(net_dir)

    # Snapshot de config
    open(joinpath(out_dir, "config_snapshot.toml"), "w") do io
        println(io, "# NeuroMIND longitudinal snapshot — $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
        println(io, "design = \"longitudinal_MS_T1_T2\"")
        println(io, "cohort = \"MS_paired_only\"")
        println(io, "n_paired_design = $N_PAIRED_DESIGN")
        println(io, "condition = \"$cond\"")
        println(io, "note = \"EC and EO run in parallel; same T1→T2 contrast; no pooling\"")
        println(io, "test = \"wilcoxon_signed_rank\"")
        println(io, "fdr = \"bh\"")
        println(io, "wpli_method = \"$WPLI_METHOD\"")
        println(io, "use_dwpli = $USE_DWPLI")
        println(io, "graph_density = $GRAPH_DENS")
        println(io, "graph_threshold_method = \"$GRAPH_METH\"")
        println(io, "n_channels_analysis = $N_CH_MONT")
        println(io, "qc_allowed = [\"include\", \"include_with_warning\"]")
        println(io, "source_config = \"$CONFIG_P\"")
    end

    paired_info = NamedTuple[]
    included_pairs = SubjPair[]

    for sp in all_pairs
        incl, reason, qc1, qc2, ep1, ep2, pct1, pct2 = evaluate_inclusion(sp, cond, qc_table)
        # n_bands_ok se rellena al cargar
        n_ok = 0
        if incl
            for b in BANDS
                r1 = load_wpli(sp.subject_id, sp.session_t1, cond, b)
                r2 = load_wpli(sp.subject_id, sp.session_t2, cond, b)
                (r1 !== nothing && r2 !== nothing) && (n_ok += 1)
            end
            push!(included_pairs, sp)
        end
        push!(paired_info, (
            subject_id      = sp.subject_id,
            session_t1      = sp.session_t1,
            session_t2      = sp.session_t2,
            n_bands_ok      = n_ok,
            included        = incl,
            excluded_reason = reason,
            qc_t1           = qc1,
            qc_t2           = qc2,
            n_epochs_t1     = ep1 === missing ? "" : string(ep1),
            n_epochs_t2     = ep2 === missing ? "" : string(ep2),
            retention_pct_t1 = pct1 === missing ? "" : string(round(pct1, digits=1)),
            retention_pct_t2 = pct2 === missing ? "" : string(round(pct2, digits=1)),
        ))
    end

    n_paired = count(r -> r.included, paired_info)
    n_excl   = count(r -> !r.included, paired_info)
    println("  Comparación longitudinal $cond | Pares incluidos: $n_paired / diseño $N_PAIRED_DESIGN | Excluidos (QC/datos): $n_excl")
    CSV.write(joinpath(out_dir, "paired_subjects.csv"), DataFrame(paired_info))

    if n_paired < 1
        println("  ⚠  Análisis omitido: sin pares completos T1/T2 para $cond")
        open(joinpath(out_dir, "longitudinal_summary.json"), "w") do io
            write(io, """{"n_paired":0,"n_paired_design":$N_PAIRED_DESIGN,"n_candidates":$(length(all_pairs)),"n_t1":0,"n_t2":0,"n_excluded":$n_excl,"n_excluded_qc":0,"n_excluded_data":$n_excl,"n_total_sig":0,"n_bands":0,"cond":"$cond","test":"wilcoxon","cohort":"MS_paired_only","timestamp":"$(Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"))"}""")
        end
        continue
    end

    # ── Acumuladores wPLI ──
    t1_data = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
    t2_data = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
    paired_data = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}, Vector{String}, Matrix{Float64}}}}()
    subject_means = NamedTuple[]

    for sp in included_pairs
        for band in BANDS
            r1 = load_wpli(sp.subject_id, sp.session_t1, cond, band)
            r2 = load_wpli(sp.subject_id, sp.session_t2, cond, band)
            (r1 === nothing || r2 === nothing) && continue
            (ch1, W1) = r1; (ch2, W2) = r2
            haskey(t1_data, band) || (t1_data[band] = [])
            haskey(t2_data, band) || (t2_data[band] = [])
            haskey(paired_data, band) || (paired_data[band] = [])
            push!(t1_data[band], (ch1, W1))
            push!(t2_data[band], (ch2, W2))
            push!(paired_data[band], (ch1, W1, ch2, W2))
            for (W, tp) in ((W1, "T1"), (W2, "T2"))
                μ = upper_mean(W)
                isnan(μ) && continue
                push!(subject_means, (
                    subject_id = sp.subject_id,
                    timepoint  = tp,
                    band       = band,
                    cond       = cond,
                    mean_wpli  = round(μ, digits=5),
                ))
            end
        end
    end

    !isempty(subject_means) &&
        CSV.write(joinpath(out_dir, "subject_band_means.csv"), DataFrame(subject_means))

    band_stats_rows = NamedTuple[]
    all_sig_for_figs = Dict{String, Vector{Any}}()

    # ── wPLI edge-wise ──
    for band in BANDS
        t1_mats = get(t1_data, band, Tuple{Vector{String},Matrix{Float64}}[])
        t2_mats = get(t2_data, band, Tuple{Vector{String},Matrix{Float64}}[])
        pd_mats = get(paired_data, band, Tuple{Vector{String},Matrix{Float64},Vector{String},Matrix{Float64}}[])
        (isempty(t1_mats) || isempty(t2_mats)) && continue

        all_ch_sets = [Set(ch) for (ch, _) in vcat(t1_mats, t2_mats)]
        common_set  = reduce(intersect, all_ch_sets)
        ch_ref = t1_mats[1][1]
        common_ch = filter(c -> c in common_set, ch_ref)
        n = length(common_ch)
        n < 2 && continue

        get_W(ch, W) = realign_matrix(ch, W, common_ch)
        Wt1 = mean(get_W(ch, W) for (ch, W) in t1_mats)
        Wt2 = mean(get_W(ch, W) for (ch, W) in t2_mats)
        Wdiff = Wt2 .- Wt1

        save_mat_csv(joinpath(out_dir, "longitudinal_connectivity_t1_$(band).csv"), Wt1, common_ch)
        save_mat_csv(joinpath(out_dir, "longitudinal_connectivity_t2_$(band).csv"), Wt2, common_ch)
        save_mat_csv(joinpath(out_dir, "longitudinal_difference_$(band).csv"), Wdiff, common_ch)

        upper_idx = [(i, j) for i in 1:n for j in (i + 1):n]
        p_vec  = ones(length(upper_idx))
        pp_vec = ones(length(upper_idx))
        d_vec  = zeros(length(upper_idx))
        n_vec  = zeros(Int, length(upper_idx))

        for (k, (i, j)) in enumerate(upper_idx)
            v1 = Float64[]; v2 = Float64[]
            for (ch1, W1, ch2, W2) in pd_mats
                push!(v1, get_W(ch1, W1)[i, j])
                push!(v2, get_W(ch2, W2)[i, j])
            end
            length(v1) < 2 && continue
            n_vec[k] = length(v1)
            p_vec[k], _ = wilcoxon_p(v1, v2)
            pp_vec[k], _ = paired_t(v1, v2)
            d_vec[k] = cohen_dz(v1, v2)
        end
        q_vec = bh_qvalues(p_vec)

        stat_rows = [(
            ch_a         = common_ch[i],
            ch_b         = common_ch[j],
            t1_mean      = round(Wt1[i, j], digits=5),
            t2_mean      = round(Wt2[i, j], digits=5),
            diff         = round(Wdiff[i, j], digits=5),
            p_value      = round(p_vec[k], digits=5),
            q_value      = round(q_vec[k], digits=5),
            p_parametric = round(pp_vec[k], digits=5),
            effect_d     = round(d_vec[k], digits=4),
            n            = n_vec[k],
        ) for (k, (i, j)) in enumerate(upper_idx)]

        CSV.write(joinpath(out_dir, "longitudinal_statistics_$(band).csv"), DataFrame(stat_rows))

        sig_rows = filter(r -> r.q_value < 0.05, stat_rows)
        sig_df = isempty(sig_rows) ?
            DataFrame(ch_a=String[], ch_b=String[], t1_mean=Float64[],
                      t2_mean=Float64[], diff=Float64[],
                      p_value=Float64[], q_value=Float64[],
                      p_parametric=Float64[], effect_d=Float64[], n=Int[]) :
            DataFrame(sig_rows)
        CSV.write(joinpath(out_dir, "significant_longitudinal_edges_$(band).csv"), sig_df)
        all_sig_for_figs[band] = collect(sig_rows)

        n_sig = length(sig_rows)
        mean_t1 = mean(Wt1[i, j] for (i, j) in upper_idx)
        mean_t2 = mean(Wt2[i, j] for (i, j) in upper_idx)
        @printf("  %-10s  %3d ch  %4d pares  %3d sig (q<0.05)  T1=%.3f  T2=%.3f\n",
                band, n, length(upper_idx), n_sig, mean_t1, mean_t2)

        push!(band_stats_rows, (
            band       = band,
            n_channels = n,
            n_pairs    = length(upper_idx),
            n_sig      = n_sig,
            pct_sig    = round(100.0 * n_sig / max(1, length(upper_idx)), digits=2),
            t1_mean    = round(mean_t1, digits=5),
            t2_mean    = round(mean_t2, digits=5),
            diff_mean  = round(mean_t2 - mean_t1, digits=5),
            mean_p     = round(mean(p_vec), digits=5),
            mean_d     = round(mean(abs, d_vec), digits=4),
            n_subjects = length(pd_mats),
        ))

        # Figuras heatmaps (escala emparejada T1/T2 + Δ) y red FDR/exploratoria
        try
            lim_ab = max(maximum(Wt1), maximum(Wt2), 1e-12)
            _save_heatmap(joinpath(figs_dir, "heatmap_t1_$(band).png"), Wt1, common_ch,
                          "wPLI T1 — $band ($cond)"; colorrange=(0.0, lim_ab))
            _save_heatmap(joinpath(figs_dir, "heatmap_t2_$(band).png"), Wt2, common_ch,
                          "wPLI T2 — $band ($cond)"; colorrange=(0.0, lim_ab))
            _save_heatmap(joinpath(figs_dir, "heatmap_delta_$(band).png"), Wdiff, common_ch,
                          "Δ wPLI (T2−T1) — $band ($cond)"; diverging=true)
            save_heatmap_triplet(joinpath(figs_dir, "heatmap_triplet_$(band).png"),
                Wt1, Wt2, Wdiff, common_ch,
                ("T1 — $band ($cond)", "T2 — $band ($cond)", "Δ (T2−T1) — $band ($cond)"))
            stats_df = DataFrame(stat_rows)
            _save_explore_or_sig_network(figs_dir, band, cond, common_ch, sig_rows, stats_df)
        catch e
            @warn "Figuras wPLI $band fallidas: $e"
        end

        # ── Network por banda ──
        try
            s1, d1, ns1, thr1 = nodal_strength_degree(Wt1)
            s2, d2, ns2, thr2 = nodal_strength_degree(Wt2)
            CSV.write(joinpath(net_dir, "network_metrics_t1_$(band).csv"),
                DataFrame(channel=common_ch, strength=round.(s1, digits=5),
                          degree=d1, norm_strength=round.(ns1, digits=5)))
            CSV.write(joinpath(net_dir, "network_metrics_t2_$(band).csv"),
                DataFrame(channel=common_ch, strength=round.(s2, digits=5),
                          degree=d2, norm_strength=round.(ns2, digits=5)))
            CSV.write(joinpath(net_dir, "network_metrics_delta_$(band).csv"),
                DataFrame(channel=common_ch,
                          delta_strength=round.(s2 .- s1, digits=5),
                          delta_degree=d2 .- d1,
                          delta_norm_strength=round.(ns2 .- ns1, digits=5)))
        catch e
            @warn "Network $band fallido: $e"
        end
    end

    !isempty(band_stats_rows) &&
        CSV.write(joinpath(out_dir, "band_statistics_longitudinal.csv"), DataFrame(band_stats_rows))

    # Network global statistics (mean strength por sujeto×banda×tiempo → Wilcoxon)
    net_global = NamedTuple[]
    for band in BANDS
        pd_mats = get(paired_data, band, Tuple{Vector{String},Matrix{Float64},Vector{String},Matrix{Float64}}[])
        length(pd_mats) < 2 && continue
        # canales comunes
        all_ch_sets = [intersect(Set(ch1), Set(ch2)) for (ch1, _, ch2, _) in pd_mats]
        isempty(all_ch_sets) && continue
        common_set = reduce(intersect, all_ch_sets)
        ch_ref = pd_mats[1][1]
        common_ch = filter(c -> c in common_set, ch_ref)
        length(common_ch) < 2 && continue
        ms1 = Float64[]; ms2 = Float64[]
        for (ch1, W1, ch2, W2) in pd_mats
            Wa = realign_matrix(ch1, W1, common_ch)
            Wb = realign_matrix(ch2, W2, common_ch)
            s1, _, _, _ = nodal_strength_degree(Wa)
            s2, _, _, _ = nodal_strength_degree(Wb)
            push!(ms1, mean(s1)); push!(ms2, mean(s2))
        end
        pw, _ = wilcoxon_p(ms1, ms2)
        dz = cohen_dz(ms1, ms2)
        push!(net_global, (
            band = band,
            metric = "mean_strength",
            t1_mean = round(mean(ms1), digits=5),
            t2_mean = round(mean(ms2), digits=5),
            diff = round(mean(ms2) - mean(ms1), digits=5),
            p_value = round(pw, digits=5),
            effect_d = round(dz, digits=4),
            n = length(ms1),
        ))
    end
    if !isempty(net_global)
        ng_df = DataFrame(net_global)
        ng_df[!, :q_value] = round.(bh_qvalues(Float64.(ng_df.p_value)), digits=5)
        CSV.write(joinpath(net_dir, "network_global_statistics.csv"), ng_df)
    end

    # ── Band power longitudinal (canal×banda → vectores T1/T2) ──
    bp_store = Dict{Tuple{String,String}, Tuple{Vector{Float64},Vector{Float64}}}()
    for sp in included_pairs
        bp1 = load_band_power(sp.subject_id, sp.session_t1, cond)
        bp2 = load_band_power(sp.subject_id, sp.session_t2, cond)
        (bp1 === nothing || bp2 === nothing) && continue
        for ch in intersect(keys(bp1), keys(bp2)), band in BANDS
            haskey(bp1[ch], band) && haskey(bp2[ch], band) || continue
            key = (ch, band)
            haskey(bp_store, key) || (bp_store[key] = (Float64[], Float64[]))
            push!(bp_store[key][1], bp1[ch][band])
            push!(bp_store[key][2], bp2[ch][band])
        end
    end

    bp_stat_rows = NamedTuple[]
    for ((ch, band), (v1, v2)) in bp_store
        length(v1) < 2 && continue
        pw, _ = wilcoxon_p(v1, v2)
        dz = cohen_dz(v1, v2)
        push!(bp_stat_rows, (
            channel  = ch,
            band     = band,
            t1_mean  = round(mean(v1), digits=5),
            t2_mean  = round(mean(v2), digits=5),
            diff     = round(mean(v2) - mean(v1), digits=5),
            p_value  = round(pw, digits=5),
            effect_d = round(dz, digits=4),
            n        = length(v1),
        ))
    end

    if !isempty(bp_stat_rows)
        # FDR por banda
        bp_df = DataFrame(bp_stat_rows)
        q_all = fill(1.0, nrow(bp_df))
        for band in BANDS
            idx = findall(i -> string(bp_df.band[i]) == band, 1:nrow(bp_df))
            isempty(idx) && continue
            q_all[idx] = bh_qvalues(Float64.(bp_df.p_value[idx]))
        end
        bp_df[!, :q_value] = round.(q_all, digits=5)
        CSV.write(joinpath(spec_dir, "band_power_delta_statistics.csv"), bp_df)
        sig_bp_df = filter(r -> r.q_value < 0.05, bp_df)
        CSV.write(joinpath(spec_dir, "significant_band_power_changes.csv"), sig_bp_df)

        # Topo Δ por banda (media grupal del Δ)
        xy = Dict{String,Tuple{Float64,Float64}}()
        if !isempty(included_pairs)
            xy = load_electrode_xy(included_pairs[1].subject_id, included_pairs[1].session_t1)
        end
        for band in BANDS
            sub = filter(r -> string(r.band) == band, eachrow(bp_df))
            isempty(sub) && continue
            chs = [string(r.channel) for r in sub]
            dvs = [Float64(r.diff) for r in sub]
            try
                _save_topo_delta(joinpath(figs_dir, "topo_delta_bandpower_$(band).png"),
                                 chs, dvs, xy, "Δ band power — $band ($cond)")
            catch e
                @warn "Topo band power $band fallido: $e"
            end
        end
    end

    # Paired mean wPLI figure
    if !isempty(subject_means)
        try
            _save_paired_means(joinpath(figs_dir, "paired_mean_wpli_by_band.png"),
                               DataFrame(subject_means))
        catch e
            @warn "Figura paired means fallida: $e"
        end
    end

    n_total_sig = isempty(band_stats_rows) ? 0 : sum(r.n_sig for r in band_stats_rows)
    best_band = honest_best_band(band_stats_rows)

    n_qc_loss = count(r -> !r.included && occursin("QC", r.excluded_reason), paired_info)
    n_data_loss = n_excl - n_qc_loss

    open(joinpath(out_dir, "longitudinal_summary.json"), "w") do io
        pairs_json = join(
            ["  \"$(r.band)_n_sig\": $(r.n_sig), \"$(r.band)_t1\": $(r.t1_mean), \"$(r.band)_t2\": $(r.t2_mean)"
             for r in band_stats_rows], ",\n"
        )
        write(io, """{
  "n_paired": $n_paired,
  "n_paired_design": $N_PAIRED_DESIGN,
  "n_candidates": $(length(all_pairs)),
  "n_t1": $n_paired,
  "n_t2": $n_paired,
  "n_excluded": $n_excl,
  "n_excluded_qc": $n_qc_loss,
  "n_excluded_data": $n_data_loss,
  "n_total_sig": $n_total_sig,
  "n_bands": $(length(band_stats_rows)),
  "best_band": "$best_band",
  "cond": "$cond",
  "cohort": "MS_paired_only",
  "design": "longitudinal_MS_T1_T2",
  "test": "wilcoxon_signed_rank",
  "fdr": "bh",
  "wpli_method": "$WPLI_METHOD",
  "use_dwpli": $USE_DWPLI,
  "graph_density": $GRAPH_DENS,
  "timestamp": "$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))",
  $pairs_json
}""")
    end
    println("  ✅ Guardado en: $out_dir\n")
end

println("✅ Análisis longitudinal completado (EC y EO en paralelo).")
println("   Abre el dashboard → Fase 14 Evaluación Longitudinal")
