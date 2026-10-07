# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Longitudinal
#  Análisis T1→T2 (solo EM) + generación de todas las figuras, en
#  una única pasada: CSV → figura PNG, con datos reales recién
#  calculados (mean_strength por sujeto se captura durante el propio
#  cálculo de network_global_statistics — sin reconstrucción aparte).
# ═══════════════════════════════════════════════════════════════
#
#  Sustituye: scripts/run_longitudinal_analysis.jl (lógica) +
#  src/longitudinal/{LongitudinalFigures,LongitudinalManuscript}.jl +
#  src/visualization/{GroupVizCommon,PublicationCommon,PublicationTheme,SummaryFigures}.jl
#  (solo la porción que usaba el lado longitudinal / la síntesis).
#
#  Diseño experimental (Fig. 3.1):
#    · Solo pacientes EM con par completo T1+T2 (N diseño = 30)
#    · Controles NO entran; eyesclosed y eyesopen en paralelo, sin pooling
#
#  Salida:
#    results/longitudinal/{eyesclosed|eyesopen}/
#      config_snapshot.toml · longitudinal_summary.json
#      tables/   — todos los CSV + figures_manifest.tsv
#      figures/  — todas las figuras PNG (exploratorias +, solo en
#                  eyesclosed, forest/Δ pareado/matrices)
#    results/summary/ — figura de síntesis transversal↔longitudinal,
#      generada por el segundo de los dos análisis que se ejecute
#      (comprueba si el otro dominio ya tiene resultados en disco).
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/longitudinal/Longitudinal.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      27-07-2026
#  Modificado  28-07-2026
# ───────────────────────────────────────────────────────────────

module Longitudinal

using CSV, DataFrames, Statistics, LinearAlgebra, Dates, TOML, Printf, CairoMakie, Random

export run

# ═══════════════════════════════════════════════════════════════
#  Constantes compartidas (duplicadas con Transversal.jl a propósito
#  — ver nota de diseño en el plan: datos de referencia, no lógica)
# ═══════════════════════════════════════════════════════════════

const BAND_ORDER = ["DELTA", "THETA", "ALPHA", "BETA_LOW", "BETA_MID", "BETA_HIGH", "GAMMA"]
const DIVERGING_CMAP = Reverse(:RdBu)   # positivo (T2−T1) = rojo

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

const BOOTSTRAP_N    = 5000
const BOOTSTRAP_SEED = 20260725
const FDR_ALPHA       = 0.05
const WILCOXON_EXACT_MAX_N = 30
const WILCOXON_METHOD = "auto_exact_conditional_dp_n_le_$(WILCOXON_EXACT_MAX_N)"
const LONGITUDINAL_SCHEMA_VERSION = 2
const STATISTICS_SOURCE = "NeuroMIND.Longitudinal.production"
const BOOTSTRAP_METHOD = "paired_percentile_bootstrap"
const QUANTILE_METHOD = "Statistics.quantile_linear_alpha_1_beta_1"
const RRB_METHOD = "paired_signed_ranks_(Wplus-Wminus)/(Wplus+Wminus)"
const EFFECT_DZ_METHOD = "mean(T2-T1)/std(T2-T1)"
const EDGE_FDR_SCOPE = "351_or_available_edges_within_each_band"
const GLOBAL_FDR_SCOPE = "available_frequency_bands"
const POWER_FDR_SCOPE = "available_channels_within_each_band"

const N_PAIRED_DESIGN = 30   # Fig. 3.1 — pares EM T1+T2 comparables
const QC_ALLOWED      = Set(["include", "include_with_warning"])

_ch_xy(ch::AbstractString) = get(CH_POS, uppercase(String(ch)), (0.0, 0.0))

# ── Tema único (ex PublicationTheme.jl), aplicado a TODA figura ────
#  Paleta Okabe–Ito (apta daltonismo); CairoMakie acepta strings hex directamente.

const COLOR_POS  = "#D55E00"
const COLOR_NEG  = "#0072B2"
const COLOR_NS   = "#999999"
const COLOR_ZERO = "#333333"

function apply_theme!()
    set_theme!(Theme(
        fontsize = 11, backgroundcolor = :white,
        Axis = (backgroundcolor=:white, spinewidth=0.8, xtickwidth=0.7, ytickwidth=0.7,
                xticksize=3, yticksize=3, xgridvisible=false, ygridvisible=false,
                titlefont=:regular, titlesize=11, xlabelsize=10, ylabelsize=10,
                xticklabelsize=8, yticklabelsize=8),
        Colorbar = (labelsize=9, ticklabelsize=8, spinewidth=0.6, tickwidth=0.6),
        Legend   = (framevisible=false, labelsize=9, titlesize=9, padding=(2,2,2,2), rowgap=2),
        Lines    = (linewidth=1.2,), Scatter = (strokewidth=0.4,),
    ))
    return nothing
end

# size de Figure ≈ puntos tipográficos (72 dpi); la resolución de impresión
# se aplica solo al exportar vía save_png (px_per_unit = PRINT_DPI/72).
const PRINT_DPI    = 300
const MAKIE_PT_DPI = 72

fig_size_mm(w_mm, h_mm) = (
    round(Int, w_mm / 25.4 * MAKIE_PT_DPI),
    round(Int, h_mm / 25.4 * MAKIE_PT_DPI),
)

fmt_p(x) = x < 0.001 ? "<0.001" : string(round(x; digits=3))
fmt_q(x) = x < 0.001 ? "<0.001" : string(round(x; digits=3))
fmt_d(x) = string(round(x; digits=2))

"""Guarda PNG a `dpi` de impresión; sin PDF (ver nota de diseño §Idea 3)."""
function save_png(fig, stem::String; dpi::Real = PRINT_DPI)
    mkpath(dirname(stem))
    path = stem * ".png"
    save(path, fig; px_per_unit = dpi / MAKIE_PT_DPI)
    return path
end

# ═══════════════════════════════════════════════════════════════
#  Estadística (idéntica a la de run_longitudinal_analysis.jl)
# ═══════════════════════════════════════════════════════════════

function norm_cond(c::AbstractString)::String
    lc = lowercase(String(c))
    lc in ("ec", "eyesclosed") && return "eyesclosed"
    lc in ("eo", "eyesopen")   && return "eyesopen"
    return lc
end

cond_label(c::AbstractString)::String = norm_cond(c) == "eyesclosed" ? "EC" : "EO"

function cond_code(c::AbstractString)::String
    nc = norm_cond(c)
    nc == "eyesclosed" && return "EC"
    nc == "eyesopen"   && return "EO"
    return uppercase(String(c))
end

is_t1_session(s::AbstractString)::Bool =
    uppercase(String(s)) in ("T1", "BASELINE", "BL", "V1", "VISIT1", "S1", "PRE")

"""True si el ID parece control (MC*, C*, HC*), no paciente EM."""
function is_control_id(sid::AbstractString)::Bool
    s = uppercase(String(sid))
    startswith(s, "MC") && return true
    startswith(s, "HC") && return true
    startswith(s, "C") && !startswith(s, "M") && return true
    return false
end
is_ms_id(sid::AbstractString)::Bool = !is_control_id(sid)

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
_norm_cdf(z::Float64)::Float64 = 0.5 * (1.0 + _erf_approx(z / sqrt(2.0)))

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

"""p exacto condicional de Wilcoxon para rangos (incluye empates).

Los rangos se multiplican por 2 para convertir los medios rangos en enteros y
la distribución nula de W+ se obtiene por programación dinámica. Esto evita
enumerar las 2^n asignaciones de signo y coincide con el test exacto estándar
cuando no hay empates.
"""
function _wilcoxon_exact_p(ranks::Vector{Float64}, positive::BitVector)::Float64
    ranks2 = round.(Int, 2 .* ranks)
    observed = sum(ranks2[positive])
    total = sum(ranks2)
    counts = zeros(Int128, total + 1)
    counts[1] = 1
    reached = 0
    for r in ranks2
        for s in reached:-1:0
            counts[s + r + 1] += counts[s + 1]
        end
        reached += r
    end
    denom = 2.0^length(ranks2)
    p_lo = Float64(sum(counts[1:(observed + 1)])) / denom
    p_hi = Float64(sum(counts[(observed + 1):end])) / denom
    return min(1.0, 2.0 * min(p_lo, p_hi))
end

"""Wilcoxon signed-rank → (p, effect_r), exacto para n≤30.

`method=:auto` usa el cálculo exacto condicional para las muestras del diseño
longitudinal actual y una aproximación normal con corrección de empates y
continuidad únicamente si n supera `WILCOXON_EXACT_MAX_N`.
Con `return_method=true` devuelve además el método realmente usado.
"""
function wilcoxon_p(before::Vector{Float64}, after::Vector{Float64};
                     method::Symbol=:auto, return_method::Bool=false)
    if length(before) != length(after)
        return return_method ? (1.0, 0.0, "length_mismatch") : (1.0, 0.0)
    end
    diffs = filter(!=(0.0), after .- before)
    nd = length(diffs)
    if nd < 2
        return return_method ? (1.0, 0.0, "insufficient_n") : (1.0, 0.0)
    end
    ranks = _assign_ranks(abs.(diffs))
    W_plus  = sum(ranks[diffs .> 0])
    W_minus = sum(ranks[diffs .< 0])
    use_exact = method == :exact || (method == :auto && nd <= WILCOXON_EXACT_MAX_N)
    if use_exact
        p = _wilcoxon_exact_p(ranks, diffs .> 0)
        r = abs(W_plus - W_minus) / max(W_plus + W_minus, eps())
        used = "exact_conditional_dp"
    else
        μW = nd * (nd + 1) / 4
        tie_map = Dict{Float64,Int}()
        for v in abs.(diffs)
            tie_map[v] = get(tie_map, v, 0) + 1
        end
        tie_counts = values(tie_map)
        tie_term = sum(t * (t + 1) * (2t + 1) for t in tie_counts if t > 1)
        σ2 = (nd * (nd + 1) * (2nd + 1) - tie_term) / 24
        σW = sqrt(max(σ2, 0.0))
        if σW < 1e-12
            return return_method ? (1.0, 0.0, "asymptotic_tie_corrected") : (1.0, 0.0)
        end
        z = max(0.0, abs(W_plus - μW) - 0.5) / σW
        p = clamp(2.0 * (1.0 - _norm_cdf(z)), 0.0, 1.0)
        r = abs(z) / sqrt(nd)
        used = "asymptotic_tie_continuity_corrected"
    end
    return return_method ? (p, r, used) : (p, r)
end

function cohen_dz(before::Vector{Float64}, after::Vector{Float64})::Float64
    length(before) == length(after) || return 0.0
    diffs = after .- before
    n = length(diffs); n < 2 && return 0.0
    s = std(diffs); s < 1e-12 && return 0.0
    return mean(diffs) / s
end

"""Correlación biserial de rangos pareada (Kerby): (P−N)/(P+N)."""
function paired_rank_biserial(before::Vector{Float64}, after::Vector{Float64})::Float64
    length(before) == length(after) || return 0.0
    diffs = filter(!=(0.0), after .- before)
    nd = length(diffs); nd < 1 && return 0.0
    ranks = _assign_ranks(abs.(diffs))
    W_plus  = sum(ranks[diffs .> 0])
    W_minus = sum(ranks[diffs .< 0])
    denom = W_plus + W_minus; denom < 1e-12 && return 0.0
    return clamp((W_plus - W_minus) / denom, -1.0, 1.0)
end

"""IC95% bootstrap de Cohen d_z (after−before), remuestreo de pares."""
function bootstrap_cohen_dz_ci(before::Vector{Float64}, after::Vector{Float64};
                               n_boot::Int=BOOTSTRAP_N, seed::Int=BOOTSTRAP_SEED, alpha=0.05)
    d0 = cohen_dz(before, after)
    n = length(before); n < 2 && return (d0, NaN, NaN)
    rng = MersenneTwister(seed)
    ds = Vector{Float64}(undef, n_boot)
    for i in 1:n_boot
        idx = rand(rng, 1:n, n)
        ds[i] = cohen_dz(before[idx], after[idx])
    end
    lo, hi = quantile(ds, [alpha/2, 1-alpha/2])
    return (d0, lo, hi)
end

function median_iqr(x::Vector{Float64})
    isempty(x) && return (NaN, NaN, NaN)
    qs = quantile(x, [0.25, 0.5, 0.75])
    return (qs[2], qs[1], qs[3])
end

function mean_ci_bootstrap(x::Vector{Float64}; n_boot::Int=BOOTSTRAP_N, seed::Int=BOOTSTRAP_SEED, alpha=0.05)
    isempty(x) && return (NaN, NaN, NaN)
    μ = mean(x); n = length(x); n < 2 && return (μ, NaN, NaN)
    rng = MersenneTwister(seed + 17)
    ms = [mean(x[rand(rng, 1:n, n)]) for _ in 1:n_boot]
    lo, hi = quantile(ms, [alpha/2, 1-alpha/2])
    return (μ, lo, hi)
end

"""
    paired_change_summary(before, after; n_boot, seed)

Resumen longitudinal único de producción para un vector pareado. Centraliza
los descriptivos y tamaños de efecto que consumen CSV, figuras y visores; el
frontend no debe recalcular ninguno de estos valores.
"""
function paired_change_summary(before::Vector{Float64}, after::Vector{Float64};
                               n_boot::Int=BOOTSTRAP_N, seed::Int=BOOTSTRAP_SEED)
    length(before) == length(after) || error("paired_change_summary: longitudes distintas")
    isempty(before) && error("paired_change_summary: muestra vacía")
    diffs = after .- before
    med, q1, q3 = median_iqr(diffs)
    μ1, t1_lo, t1_hi = mean_ci_bootstrap(before; n_boot=n_boot, seed=seed)
    μ2, t2_lo, t2_hi = mean_ci_bootstrap(after; n_boot=n_boot, seed=seed + 101)
    μd, d_lo, d_hi = mean_ci_bootstrap(diffs; n_boot=n_boot, seed=seed + 202)
    dz, dz_lo, dz_hi = bootstrap_cohen_dz_ci(before, after; n_boot=n_boot, seed=seed + 303)
    return (
        t1_mean=μ1, t1_ci_low=t1_lo, t1_ci_high=t1_hi,
        t2_mean=μ2, t2_ci_low=t2_lo, t2_ci_high=t2_hi,
        diff_mean=μd, diff_ci_low=d_lo, diff_ci_high=d_hi,
        median_diff=med, q1_diff=q1, q3_diff=q3,
        effect_dz=dz, effect_dz_ci_low=dz_lo, effect_dz_ci_high=dz_hi,
        effect_rrb=paired_rank_biserial(before, after),
        n=length(diffs), n_positive=count(>(0.0), diffs),
        n_negative=count(<(0.0), diffs), n_zero=count(==(0.0), diffs),
    )
end

function write_statistics_contract(path::String, condition::String)
    contract = Dict{String,Any}(
        "schema_version" => LONGITUDINAL_SCHEMA_VERSION,
        "statistics_source" => STATISTICS_SOURCE,
        "condition" => condition,
        "fdr_scope_edges" => EDGE_FDR_SCOPE,
        "fdr_scope_global" => GLOBAL_FDR_SCOPE,
        "fdr_scope_power" => POWER_FDR_SCOPE,
        "bootstrap_method" => BOOTSTRAP_METHOD,
        "bootstrap_iterations" => BOOTSTRAP_N,
        "bootstrap_seed" => BOOTSTRAP_SEED,
        "quantile_method" => QUANTILE_METHOD,
        "rrb_method" => RRB_METHOD,
        "effect_dz_method" => EFFECT_DZ_METHOD,
        "generated_at" => Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"),
    )
    open(path, "w") do io
        TOML.print(io, contract; sorted=true)
    end
    return path
end

upper_mean(W::Matrix{Float64})::Float64 = begin
    n = min(size(W)...); n < 2 ? NaN : mean(W[i,j] for i in 1:n for j in (i+1):n)
end

function realign_matrix(ch::AbstractVector{<:AbstractString}, W::Matrix{Float64}, common_ch::Vector{String})::Matrix{Float64}
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

function read_mat_csv(path::String)
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    ch = String.(names(df)[2:end])
    M = Matrix{Float64}(df[:, 2:end])
    return (ch, M)
end

function sort_bands_physio(bands)::Vector{String}
    bs = String.(collect(bands))
    order = Dict(b => i for (i, b) in enumerate(BAND_ORDER))
    return sort(bs; by = b -> get(order, b, 1000 + hash(b) % 100))
end

function honest_best_band(rows; n_sig_field=:n_sig, diff_field=:diff_mean)::String
    isempty(rows) && return ""
    n_sigs = [Int(getfield(r, n_sig_field)) for r in rows]
    sum(n_sigs) == 0 && return ""
    scores = [n_sigs[i]*1000.0 + abs(Float64(getfield(rows[i], diff_field))) for i in eachindex(rows)]
    return string(getfield(rows[argmax(scores)], :band))
end

# ═══════════════════════════════════════════════════════════════
#  Carga de datos
# ═══════════════════════════════════════════════════════════════

function export_dir(res_root, subj::AbstractString, sess::AbstractString, cond::AbstractString)::String
    joinpath(res_root, "subjects", "sub-$(subj)", "ses-$(sess)", norm_cond(cond))
end

function load_wpli(res_root, subj_id, sess_id, cond, band)
    path = joinpath(export_dir(res_root, subj_id, sess_id, cond),
                    "tables", "connectivity", "wpli_$(band).csv")
    isfile(path) || return nothing
    try
        df = CSV.read(path, DataFrame)
        isempty(df) && return nothing
        row_ch = String.(df[!, 1])
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
    catch e
        @warn "load_wpli: fallo al leer/parsear CSV, se omite" subj_id sess_id cond band path exception=e
        return nothing
    end
end

function load_band_power(res_root, subj_id, sess_id, cond, bands)
    path = joinpath(export_dir(res_root, subj_id, sess_id, cond), "tables", "band_power_summary.csv")
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    isempty(df) && return nothing
    ch_col = names(df)[1]
    out = Dict{String, Dict{String,Float64}}()
    for row in eachrow(df)
        ch = string(row[ch_col])
        out[ch] = Dict{String,Float64}()
        for b in bands
            hasproperty(row, Symbol(b)) && (out[ch][b] = Float64(row[Symbol(b)]))
        end
    end
    return out
end

function load_seg_epochs(res_root, subj_id, sess_id, cond)
    path = joinpath(export_dir(res_root, subj_id, sess_id, cond), "json", "segmentation_summary.json")
    n_valid = missing; n_total = missing; pct = missing
    if isfile(path)
        txt = read(path, String)
        m = match(r"\"n_valid\"\s*:\s*(\d+)", txt); m !== nothing && (n_valid = parse(Int, m.captures[1]))
        m = match(r"\"n_total\"\s*:\s*(\d+)", txt); m !== nothing && (n_total = parse(Int, m.captures[1]))
        m = match(r"\"retention_pct\"\s*:\s*([0-9.]+)", txt); m !== nothing && (pct = parse(Float64, m.captures[1]))
    end
    return (n_valid, n_total, pct)
end

function load_qc_table(res_root)::Dict{Tuple{String,String,String}, String}
    path = joinpath(res_root, "qc", "qc_decision_table.csv")
    out = Dict{Tuple{String,String,String}, String}()
    isfile(path) || return out
    df = CSV.read(path, DataFrame)
    rename!(df, Dict(n => Symbol(lowercase(string(n))) for n in names(df)))
    sid_col = hasproperty(df, :subject_id) ? :subject_id : (hasproperty(df, :bids_id) ? :bids_id : nothing)
    sid_col === nothing && return out
    for row in eachrow(df)
        sid  = string(row[sid_col])
        sess = string(hasproperty(row, :session_id) ? row.session_id : row.session)
        cond = cond_code(string(row.condition))
        out[(sid, sess, cond)] = lowercase(string(row.final_decision))
    end
    return out
end

function qc_decision(qc::Dict, sid::String, sess::String, cond::String)::String
    cc = cond_code(cond)
    haskey(qc, (sid, sess, cc)) && return qc[(sid, sess, cc)]
    if startswith(sid, "M") && !startswith(sid, "MC")
        alt = occursin(r"^M\d$", sid) ? "M0"*sid[2:end] : (occursin(r"^M0\d$", sid) ? "M"*sid[3:end] : "")
        !isempty(alt) && haskey(qc, (alt, sess, cc)) && return qc[(alt, sess, cc)]
    end
    return "missing"
end

qc_ok(dec::String)::Bool = dec == "missing" || dec in QC_ALLOWED

function _prop_threshold(W::Matrix{Float64}, density::Float64)::Float64
    vals = sort([W[i,j] for i in 1:size(W,1) for j in (i+1):size(W,1) if W[i,j] > 0.0]; rev=true)
    isempty(vals) && return 0.0
    return vals[min(max(1, round(Int, density*length(vals))), length(vals))]
end

function nodal_strength_degree(W::Matrix{Float64}; density::Float64=0.1)
    n = size(W, 1)
    strength = [sum(W[i,j] for j in 1:n if j != i) for i in 1:n]
    thr = _prop_threshold(W, density)
    degree = [count(j -> j != i && W[i,j] > thr, 1:n) for i in 1:n]
    smax = maximum(strength); smax < 1e-12 && (smax = 1.0)
    return strength, degree, strength ./ smax, thr
end

function load_electrode_xy(bids_root, subj_id, sess_id)
    path = joinpath(bids_root, "electrodes", "sub-$(subj_id)_ses-$(sess_id)_electrodes.tsv")
    isfile(path) || return Dict{String,Tuple{Float64,Float64}}()
    df = CSV.read(path, DataFrame; delim='\t')
    rename!(df, Dict(n => Symbol(lowercase(string(n))) for n in names(df)))
    out = Dict{String,Tuple{Float64,Float64}}()
    for row in eachrow(df)
        hasproperty(row, :name) || continue
        out[string(row.name)] = (hasproperty(row,:x) ? Float64(row.x) : 0.0,
                                  hasproperty(row,:y) ? Float64(row.y) : 0.0)
    end
    return out
end

# ═══════════════════════════════════════════════════════════════
#  Primitivas de figura (duplicadas de src/transversal/Transversal.jl
#  a propósito — ver nota de diseño del plan §6)
# ═══════════════════════════════════════════════════════════════

function save_heatmap(path, W, ch, title; cmap=:viridis, diverging=false, colorrange=nothing, colorbar_label="")
    n = length(ch)
    fig = Figure(size=(720, 640), fontsize=11)
    ax = Axis(fig[1,1]; title=title, xlabel="Canal", ylabel="Canal",
              xticks=(1:n, ch), yticks=(1:n, ch), xticklabelrotation=π/3,
              xticklabelsize=7, yticklabelsize=7)
    if diverging
        lim = colorrange === nothing ? maximum(abs, W) : maximum(abs, colorrange)
        lim < 1e-12 && (lim = 1.0)
        hm = heatmap!(ax, W; colormap=DIVERGING_CMAP, colorrange=(-lim, lim))
        Colorbar(fig[1,2], hm; label=isempty(colorbar_label) ? "Δ" : colorbar_label)
    else
        hm = colorrange === nothing ? heatmap!(ax, W; colormap=cmap) : heatmap!(ax, W; colormap=cmap, colorrange=colorrange)
        Colorbar(fig[1,2], hm; label=isempty(colorbar_label) ? "wPLI" : colorbar_label)
    end
    mkpath(dirname(path)); save(path, fig); return path
end

function save_heatmap_triplet(path, Wa, Wb, Wd, ch, titles; label_a="wPLI", label_d="Δ wPLI")
    n = length(ch)
    lim_ab = max(maximum(Wa), maximum(Wb), 1e-12)
    lim_d = maximum(abs, Wd); lim_d < 1e-12 && (lim_d = 1.0)
    fig = Figure(size=(1400, 480), fontsize=10)
    for (col, (mat, ttl, div)) in enumerate([(Wa,titles[1],false), (Wb,titles[2],false), (Wd,titles[3],true)])
        ax = Axis(fig[1,col]; title=ttl, xticks=(1:n,ch), yticks=(1:n,ch),
                  xticklabelrotation=π/3, xticklabelsize=6, yticklabelsize=6)
        if div
            hm = heatmap!(ax, mat; colormap=DIVERGING_CMAP, colorrange=(-lim_d, lim_d))
            Colorbar(fig[2,col], hm; vertical=false, label=label_d, flipaxis=false)
        else
            hm = heatmap!(ax, mat; colormap=:viridis, colorrange=(0.0, lim_ab))
            Colorbar(fig[2,col], hm; vertical=false, label=label_a, flipaxis=false)
        end
    end
    mkpath(dirname(path)); save(path, fig); return path
end

function save_topo_network(path, ch, edge_rows, title; linewidth_base=1.2)
    isempty(edge_rows) && return nothing
    fig = Figure(size=(720, 700), fontsize=11)
    ax = Axis(fig[1,1]; title=title, aspect=DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    θs = range(0, 2π; length=120)
    lines!(ax, 1.05.*cos.(θs), 1.05.*sin.(θs); color=:gray60, linewidth=1)
    xy = Dict(c => _ch_xy(c) for c in ch)
    diffs = [Float64(r.diff) for r in edge_rows]
    dmax = max(maximum(abs, diffs), 1e-12)
    for r in edge_rows
        a, b = string(r.ch_a), string(r.ch_b)
        (haskey(xy,a) && haskey(xy,b)) || continue
        d = Float64(r.diff)
        col = d >= 0 ? (:firebrick) : (:steelblue)
        lw = linewidth_base + 2.5*abs(d)/dmax
        lines!(ax, [xy[a][1],xy[b][1]], [xy[a][2],xy[b][2]]; color=col, linewidth=lw)
    end
    xs = [xy[c][1] for c in ch]; ys = [xy[c][2] for c in ch]
    scatter!(ax, xs, ys; color=:white, strokecolor=:gray20, strokewidth=1.2, markersize=14)
    for c in ch; text!(ax, xy[c][1], xy[c][2]+0.06, text=c; fontsize=8, align=(:center,:bottom)); end
    xlims!(ax, -1.25, 1.25); ylims!(ax, -1.25, 1.25)
    Label(fig[2,1], "Rojo: Δ>0 (↑)   Azul: Δ<0 (↓)   Grosor ∝ |Δ|"; fontsize=10, color=:gray40, tellwidth=false)
    mkpath(dirname(path)); save(path, fig); return path
end

function save_paired_means(path, subject_means_df; band_order=BAND_ORDER)
    isempty(subject_means_df) && return nothing
    raw = unique(string.(subject_means_df.band))
    bands = sort_bands_physio(intersect(band_order, raw)); isempty(bands) && (bands = sort_bands_physio(raw))
    fig = Figure(size=(900, 460), fontsize=11)
    ax = Axis(fig[1,1]; title="Mean wPLI pareado T1 vs T2 — montaje común", xlabel="Banda",
              ylabel="Mean wPLI equivalente", xticks=(1:length(bands), bands))
    scatter!(ax, [NaN], [NaN]; color=:steelblue, markersize=8, label="T1")
    scatter!(ax, [NaN], [NaN]; color=:darkorange, markersize=8, label="T2")
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
            lines!(ax, [bi-0.15,bi+0.15], [tp["T1"],tp["T2"]]; color=(:gray60,0.45), linewidth=1)
            scatter!(ax, [bi-0.15], [tp["T1"]]; color=:steelblue, markersize=8)
            scatter!(ax, [bi+0.15], [tp["T2"]]; color=:darkorange, markersize=8)
        end
    end
    axislegend(ax; position=:rt, orientation=:horizontal)
    Label(fig[2,1], "Cada línea une el mismo participante. T1 y T2 se restringen al montaje común de la condición.",
          fontsize=8, color=:gray35, tellwidth=false, justification=:center)
    mkpath(dirname(path)); save(path, fig); return path
end

function save_topo_delta(path, ch, delta_vals, xy, title; colorbar_label="Δ",
                         q_values=nothing, n_values=nothing, low_threshold=0.7)
    fig = Figure(size=(560, 600), fontsize=11)
    ax = Axis(fig[1,1]; title=title, aspect=DataAspect())
    hidedecorations!(ax); hidespines!(ax)
    θs = range(0, 2π; length=120)
    lines!(ax, 1.05.*cos.(θs), 1.05.*sin.(θs); color=:gray60, linewidth=1)
    xs=Float64[]; ys=Float64[]; vs=Float64[]; labs=String[]
    for (i,c) in enumerate(ch)
        if haskey(xy, c)
            push!(xs, xy[c][1]); push!(ys, xy[c][2])
        else
            p = _ch_xy(c); push!(xs, p[1]); push!(ys, p[2])
        end
        push!(vs, delta_vals[i]); push!(labs, c)
    end
    if !isempty(xs)
        mx = max(maximum(abs,xs), maximum(abs,ys), 1e-9)
        mx > 1.2 && (xs ./= mx; ys ./= mx)
    end
    lim = maximum(abs, vs); lim < 1e-12 && (lim = 1.0)
    sc = scatter!(ax, xs, ys; color=vs, colormap=DIVERGING_CMAP, colorrange=(-lim,lim), markersize=28)
    for (i,lab) in enumerate(labs); text!(ax, xs[i], ys[i]+0.05, text=lab; fontsize=7, align=(:center,:bottom)); end
    has_q = q_values !== nothing && length(q_values) == length(ch)
    has_n = n_values !== nothing && length(n_values) == length(ch)
    n_max = has_n ? maximum(Int.(n_values)) : 0
    n_sig = has_q ? count(q -> Float64(q) < FDR_ALPHA, q_values) : 0
    n_low = has_n ? count(n -> Int(n) < low_threshold*n_max, n_values) : 0
    for i in eachindex(ch)
        sig = has_q && Float64(q_values[i]) < FDR_ALPHA
        low = has_n && Int(n_values[i]) < low_threshold*n_max
        if sig && low
            scatter!(ax, [xs[i]], [ys[i]]; marker=:cross, markersize=12,
                     color=:orange, strokewidth=0)
        elseif sig
            scatter!(ax, [xs[i]], [ys[i]]; marker=:cross, markersize=12,
                     color=:black, strokewidth=0)
        elseif low
            scatter!(ax, [xs[i]], [ys[i]]; marker=:circle, markersize=34,
                     color=:transparent, strokecolor=:orange, strokewidth=1.4)
        end
    end
    Colorbar(fig[1,2], sc; label=colorbar_label)
    if has_q || has_n
        Label(fig[2,1:2],
              "Cruz negra: q<0.05 (FDR-BH por banda, n=$n_sig canales).\n" *
              "Naranja (cruz o círculo): N reducido (<70% del máximo de esta condición: " *
              "n_max=$n_max; n=$n_low canales) — interpretar con cautela.\n" *
              "N exacto por canal en band_power_delta_statistics.csv.",
              fontsize=7, color=:gray35, tellwidth=false, justification=:center)
    end
    mkpath(dirname(path)); save(path, fig); return path
end

function anatomical_order(channels)::Vector{String}
    chans = unique(String.(collect(channels)))
    known = filter(c -> haskey(CH_POS, uppercase(c)), chans)
    unknown = sort(filter(c -> !haskey(CH_POS, uppercase(c)), chans))
    sort!(known; by = c -> (CH_POS[uppercase(c)][2], CH_POS[uppercase(c)][1]))
    return vcat(known, unknown)
end

"""Heatmap canal×banda de Cohen dz para potencia longitudinal T2−T1.

Replica la política visual del transversal: orden anatómico, cruz negra para
FDR y naranja para celdas con N inferior al 70% del máximo de la condición.
"""
function fig_power_effect_heatmap(bp_df::DataFrame, figs_dir, cl::AbstractString)
    isempty(bp_df) && return nothing
    chans = anatomical_order(unique(String.(bp_df.channel)))
    bnds = sort_bands_physio(unique(string.(bp_df.band)))
    nC, nB = length(chans), length(bnds)
    (nC < 1 || nB < 1) && return nothing
    idxC = Dict(c => i for (i,c) in enumerate(chans))
    idxB = Dict(b => i for (i,b) in enumerate(bnds))
    D = fill(NaN, nC, nB)
    Q = fill(1.0, nC, nB)
    N = fill(0, nC, nB)
    for r in eachrow(bp_df)
        ci, bi = idxC[string(r.channel)], idxB[string(r.band)]
        D[ci,bi] = Float64(r.effect_dz)
        Q[ci,bi] = Float64(r.q_value)
        N[ci,bi] = Int(r.n)
    end
    n_max = maximum(N)
    is_low(ci,bi) = N[ci,bi] < 0.7*n_max
    n_sig = count(!isnan(D[i,j]) && Q[i,j] < FDR_ALPHA for i in 1:nC, j in 1:nB)
    n_low = count(!isnan(D[i,j]) && is_low(i,j) for i in 1:nC, j in 1:nB)

    apply_theme!()
    lim = maximum(x -> isnan(x) ? 0.0 : abs(x), D)
    lim < 1e-9 && (lim = 1.0)
    fig = Figure(size=fig_size_mm(150, 28 + 4.3*nC); figure_padding=8)
    ax = Axis(fig[1,1];
              title="Potencia — cambio por canal×banda (Cohen d_z, T2−T1) — $cl",
              xlabel="Banda", ylabel="Canal", xticks=(1:nB,bnds), yticks=(1:nC,chans),
              xticklabelrotation=π/4, xticklabelsize=8, yticklabelsize=7)
    hm = heatmap!(ax, permutedims(D); colormap=DIVERGING_CMAP,
                  colorrange=(-lim,lim), nan_color=:gray90)
    Colorbar(fig[1,2], hm; label="Cohen d_z")
    for ci in 1:nC, bi in 1:nB
        isnan(D[ci,bi]) && continue
        sig = Q[ci,bi] < FDR_ALPHA
        low = is_low(ci,bi)
        if sig && low
            scatter!(ax, [Float64(bi)], [Float64(ci)]; marker=:cross,
                     markersize=6, color=:orange, strokewidth=0)
        elseif sig
            scatter!(ax, [Float64(bi)], [Float64(ci)]; marker=:cross,
                     markersize=6, color=:black, strokewidth=0)
        elseif low
            scatter!(ax, [Float64(bi)], [Float64(ci)]; marker=:circle,
                     markersize=5, color=:transparent, strokecolor=:orange,
                     strokewidth=0.8)
        end
    end
    Label(fig[2,1:2],
          "Cruz negra: q<0.05 (FDR-BH por banda, n=$n_sig celdas).\n" *
          "Naranja (cruz o círculo): N reducido (<70% del máximo de esta condición: " *
          "n_max=$n_max; n=$n_low celdas) — interpretar con cautela.\n" *
          "N exacto por celda en band_power_delta_statistics.csv.",
          fontsize=7, color=:gray35, tellwidth=false, justification=:center)
    return save_png(fig, joinpath(figs_dir, "power_effect_heatmap_channel_band"))
end

function save_wpli_heatmaps(figs_dir, band, cond_lbl, Wt1, Wt2, Wdiff, common_ch)
    lim_ab = max(maximum(Wt1), maximum(Wt2), 1e-12)
    save_heatmap(joinpath(figs_dir, "heatmap_t1_$(band).png"), Wt1, common_ch,
                 "wPLI T1 — $band ($cond_lbl)"; colorrange=(0.0,lim_ab), colorbar_label="wPLI")
    save_heatmap(joinpath(figs_dir, "heatmap_t2_$(band).png"), Wt2, common_ch,
                 "wPLI T2 — $band ($cond_lbl)"; colorrange=(0.0,lim_ab), colorbar_label="wPLI")
    save_heatmap(joinpath(figs_dir, "heatmap_delta_$(band).png"), Wdiff, common_ch,
                 "Δ wPLI (T2−T1) — $band ($cond_lbl)"; diverging=true, colorbar_label="Δ wPLI")
    save_heatmap_triplet(joinpath(figs_dir, "heatmap_triplet_$(band).png"), Wt1, Wt2, Wdiff, common_ch,
        ("T1 — $band ($cond_lbl)", "T2 — $band ($cond_lbl)", "Δ (T2−T1) — $band ($cond_lbl)"))
    return nothing
end

function save_sig_or_explore_network(figs_dir, band, cond_lbl, common_ch, sig_rows, stats_df)
    if !isempty(sig_rows)
        save_topo_network(joinpath(figs_dir, "sig_network_$(band).png"), common_ch, sig_rows,
                          "Edges sig. FDR — $band ($cond_lbl)")
    elseif nrow(stats_df) > 0 && hasproperty(stats_df, :effect_dz)
        top = first(sort(stats_df, :effect_dz; by=abs, rev=true), min(20, nrow(stats_df)))
        save_topo_network(joinpath(figs_dir, "explore_network_topN_$(band).png"), common_ch,
                          [NamedTuple(r) for r in eachrow(top)],
                          "Top-20 |dz| (exploratorio) — $band ($cond_lbl)")
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════
#  Figuras de manuscrito (ex LongitudinalManuscript.jl) — solo
#  eyesclosed, PNG únicamente, mismas `figures/` que las exploratorias
# ═══════════════════════════════════════════════════════════════

function _band_scores(scores::DataFrame, band::AbstractString)
    sub = filter(r -> string(r.band) == band, scores)
    nrow(sub) == 0 && error("Sin scores C para $band")
    return Float64.(sub.mean_strength_t1), Float64.(sub.mean_strength_t2)
end

function _ng_row(ng::DataFrame, band::AbstractString)
    rows = filter(r -> string(r.band) == band && string(r.metric) == "mean_strength", ng)
    nrow(rows) == 1 || error("Fila network_global ausente para $band")
    return rows[1, :]
end

function fig_forest_by_band_dz(tab_dir, figs_dir, manifest)
    ng = CSV.read(joinpath(tab_dir, "network_global_statistics.csv"), DataFrame)
    scores = CSV.read(joinpath(tab_dir, "mean_strength_scores.csv"), DataFrame)
    bands = BAND_ORDER; nB = length(bands)
    n_ch_values = unique(Int.(scores.n_channels))
    n_ch = length(n_ch_values) == 1 ? only(n_ch_values) : minimum(n_ch_values)
    dz=zeros(nB); lo=zeros(nB); hi=zeros(nB); p=zeros(nB); q=zeros(nB); rrb=zeros(nB)
    n=zeros(Int,nB); dmean=zeros(nB)
    for (i,b) in enumerate(bands)
        row = _ng_row(ng, b)
        dz[i] = Float64(row.effect_dz)
        lo[i] = Float64(row.effect_dz_ci_low)
        hi[i] = Float64(row.effect_dz_ci_high)
        p[i] = Float64(row.p_value); q[i] = Float64(row.q_value)
        rrb[i] = Float64(row.effect_rrb)
        n[i] = Int(row.n); dmean[i] = Float64(row.diff)
    end
    y_of = Dict(b => Float64(nB-i+1) for (i,b) in enumerate(bands))
    apply_theme!()
    fig = Figure(size=fig_size_mm(180,100); figure_padding=8)
    ax = Axis(fig[1,1]; xlabel="Cohen d_z (T2 − T1)", yticks=(Float64.(nB:-1:1), bands),
              title="Cambios globales mean_strength — estimando C (primario)", limits=(nothing,nothing,0.4,nB+0.6))
    vlines!(ax, [0.0]; color=COLOR_ZERO, linewidth=0.9, linestyle=:dash)
    any_sig = any(<(FDR_ALPHA), q)
    for i in 1:nB
        yi = y_of[bands[i]]; sig = q[i] < FDR_ALPHA
        col = (!any_sig) ? COLOR_NS : (sig ? (dz[i]>=0 ? COLOR_POS : COLOR_NEG) : COLOR_NS)
        lines!(ax, [lo[i],hi[i]], [yi,yi]; color=col, linewidth=1.5)
        scatter!(ax, [dz[i]], [yi]; color=col, markersize=7, marker=:circle, strokecolor=:black, strokewidth=0.4)
    end
    ax2 = Axis(fig[1,2]; limits=(0,1.05,0.4,nB+0.6))
    hidedecorations!(ax2); hidespines!(ax2)
    for i in 1:nB
        yi = y_of[bands[i]]
        txt = "n=$(n[i])  dS=$(round(dmean[i]; digits=3))  p=$(fmt_p(p[i]))  q=$(fmt_q(q[i]))  r_rb=$(fmt_d(rrb[i]))"
        text!(ax2, 0.0, yi; text=txt, align=(:left,:center), fontsize=7, color=:gray25)
    end
    colsize!(fig.layout,1,Relative(0.48)); colsize!(fig.layout,2,Relative(0.52))
    Label(fig[2,1:2],
          "Wilcoxon exacto condicional; FDR-BH entre las 7 bandas.\n" *
          "Ninguna banda con color de significación si q_global≥0.05. " *
          "IC95% bootstrap de pares (B=$BOOTSTRAP_N, semilla=$BOOTSTRAP_SEED).\n" *
          "Estimando C sobre montaje común de $n_ch canales.",
          fontsize=7, color=:gray35, tellwidth=false, justification=:center)
    path = save_png(fig, joinpath(figs_dir, "longitudinal_global_changes_by_band"))
    push!(manifest, (figure_name="longitudinal_global_changes_by_band", analysis="longitudinal", band="ALL",
          condition="EC", output_file=path, source_tables="network_global_statistics.csv;mean_strength_scores.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="C (mean_strength, intersected $(n_ch)-ch montage)",
          fdr_family="FDR-BH across 7 bands (global)"))
    return path
end

function _panel_delta!(fig, pos, t1, t2, band, row; letter)
    ax = Axis(fig[pos...]; title="$letter  $band", xlabel="Participante (ordenado por Δ)", ylabel="Δ mean_strength (T2−T1)")
    deltas = t2 .- t1
    ord = sortperm(deltas)
    xs = Float64.(1:length(deltas)); ys = deltas[ord]
    barplot!(ax, xs, ys; color=[y>=0 ? (COLOR_POS,0.75) : (COLOR_NEG,0.75) for y in ys], width=0.8)
    hlines!(ax, [0.0]; color=COLOR_ZERO, linewidth=0.9, linestyle=:dash)
    med = Float64(row.median_diff)
    q1 = Float64(row.q1_diff)
    q3 = Float64(row.q3_diff)
    hlines!(ax, [med]; color=:black, linewidth=1.6)
    μ = Float64(row.diff)
    mlo = Float64(row.diff_ci_low)
    mhi = Float64(row.diff_ci_high)
    n = length(deltas)
    dz = Float64(row.effect_dz)
    rrb = Float64(row.effect_rrb)
    ann = "n=$n  mediana=$(round(med;digits=3)) IQR[$(round(q1;digits=3)), $(round(q3;digits=3))]\n" *
          "media=$(round(μ;digits=3)) IC95%[$(round(mlo;digits=3)), $(round(mhi;digits=3))] (secundario)\n" *
          "d_z=$(fmt_d(dz))  r_rb=$(fmt_d(rrb))  p=$(fmt_p(Float64(row.p_value)))  q_global=$(fmt_q(Float64(row.q_value)))"
    return ax, ann
end

function fig_alpha_delta_pareado(tab_dir, figs_dir, manifest)
    ng = CSV.read(joinpath(tab_dir, "network_global_statistics.csv"), DataFrame)
    scores = CSV.read(joinpath(tab_dir, "mean_strength_scores.csv"), DataFrame)
    apply_theme!()
    fig = Figure(size=fig_size_mm(170,100); figure_padding=8)
    for (col, band, letter) in ((1,"ALPHA","A"), (2,"DELTA","B"))
        t1, t2 = _band_scores(scores, band)
        row = _ng_row(ng, band)
        ax, ann = _panel_delta!(fig, (1,col), t1, t2, band, row; letter=letter)
        if band == "DELTA" && Float64(row.q_value) >= FDR_ALPHA
            Label(fig[2,col], ann * "\nDELTA: mayor |Δ| descriptivo global; no supera FDR.", fontsize=6.5, color=:gray30, tellwidth=false)
        else
            Label(fig[2,col], ann, fontsize=6.5, color=:gray30, tellwidth=false)
        end
    end
    Label(fig[3,1:2],
          "Estimando C (mean_strength, montaje común). Wilcoxon exacto condicional; FDR-BH entre 7 bandas.\n" *
          "No afirmar significación si q_global≥0.05. ALPHA: banda de interés; " *
          "DELTA: mayor cambio descriptivo.",
          fontsize=7, color=:gray35, tellwidth=false, justification=:center)
    path = save_png(fig, joinpath(figs_dir, "longitudinal_alpha_delta_individual_changes"))
    push!(manifest, (figure_name="longitudinal_alpha_delta_individual_changes", analysis="longitudinal", band="ALPHA+DELTA",
          condition="EC", output_file=path, source_tables="mean_strength_scores.csv;network_global_statistics.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="C (mean_strength, intersected montage)",
          fdr_family="FDR-BH across 7 bands (global q)"))
    return path
end

function fig_matrices_alpha(tab_dir, figs_dir, manifest)
    r1 = read_mat_csv(joinpath(tab_dir, "longitudinal_connectivity_t1_ALPHA.csv"))
    r2 = read_mat_csv(joinpath(tab_dir, "longitudinal_connectivity_t2_ALPHA.csv"))
    rd = read_mat_csv(joinpath(tab_dir, "longitudinal_difference_ALPHA.csv"))
    (r1 === nothing || r2 === nothing) && error("Matrices ALPHA longitudinales ausentes")
    ch, Wt1, Wt2 = r1[1], r1[2], r2[2]
    Wd = rd === nothing ? (Wt2 .- Wt1) : rd[2]
    n = length(ch)
    stats = CSV.read(joinpath(tab_dir, "longitudinal_statistics_ALPHA.csv"), DataFrame)
    sig = filter(r -> Float64(r.q_value) < FDR_ALPHA, stats)
    n_edges = nrow(stats)
    n_sig = nrow(sig)
    pdf = CSV.read(joinpath(tab_dir, "paired_subjects.csv"), DataFrame)
    n_subj = count(r -> Bool(r.included) || string(r.included) == "true", eachrow(pdf))

    lim_ab = max(maximum(Wt1), maximum(Wt2), 1e-12)
    lim_d = max(maximum(abs, Wd), 1e-12)
    Wd_tri = copy(Wd)
    for i in 1:n, j in 1:n; i <= j && (Wd_tri[i,j] = NaN); end

    apply_theme!()
    fig = Figure(size=fig_size_mm(180,85); figure_padding=6)
    titles = ("A  T1", "B  T2", "C  T2 − T1")
    mats = (Wt1, Wt2, Wd_tri)
    tick_idx = collect(1:2:n)
    n in tick_idx || push!(tick_idx, n)
    for (col, (mat,ttl,div)) in enumerate(zip(mats, titles, (false,false,true)))
        ax = Axis(fig[1,col]; title=ttl, aspect=1,
                  xticks=(tick_idx,ch[tick_idx]), yticks=(tick_idx,ch[tick_idx]),
                  xticklabelrotation=π/2, xticklabelsize=6, yticklabelsize=6)
        if div
            hm = heatmap!(ax, mat; colormap=DIVERGING_CMAP, colorrange=(-lim_d,lim_d), nan_color=:white)
            Colorbar(fig[2,col], hm; vertical=false, label="Δ wPLI", flipaxis=false, height=8, labelsize=7, ticklabelsize=6)
            for r in eachrow(sig)
                ia = findfirst(==(string(r.ch_a)), ch)
                ib = findfirst(==(string(r.ch_b)), ch)
                (ia === nothing || ib === nothing) && continue
                i, j = min(ia,ib), max(ia,ib)
                scatter!(ax, [Float64(i)], [Float64(j)]; marker=:rect, markersize=4,
                         color=:transparent, strokecolor=:black, strokewidth=0.9)
            end
        else
            hm = heatmap!(ax, mat; colormap=:viridis, colorrange=(0.0,lim_ab))
            Colorbar(fig[2,col], hm; vertical=false, label="wPLI", flipaxis=false, height=8, labelsize=7, ticklabelsize=6)
        end
    end
    fdr_note = n_sig == 0 ? "Ninguna arista supera FDR (n_sig=0). No se marcan hallazgos confirmatorios." :
                             "Contornos: aristas FDR (n_sig=$n_sig)."
    Label(fig[3,1:3],
          "Montaje común ($n canales), n=$n_subj pares. A/B comparten escala; " *
          "C muestra el triángulo inferior.\n" *
          "Wilcoxon exacto condicional; FDR-BH entre $n_edges aristas dentro de ALPHA.\n" *
          "$fdr_note",
          fontsize=7, color=:gray35, tellwidth=false, justification=:center)
    path = save_png(fig, joinpath(figs_dir, "longitudinal_alpha_connectivity_matrices"))
    push!(manifest, (figure_name="longitudinal_alpha_connectivity_matrices", analysis="longitudinal", band="ALPHA",
          condition="EC", output_file=path,
          source_tables="longitudinal_connectivity_*_ALPHA.csv;longitudinal_difference_ALPHA.csv;longitudinal_statistics_ALPHA.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="B (group matrices, intersected montage)",
          fdr_family="FDR-BH among $n_edges edges within ALPHA"))
    return path
end

function write_manifest(path, rows)
    isempty(rows) && return nothing
    mkpath(dirname(path))
    CSV.write(path, DataFrame(rows); delim='\t')
    return path
end

# ═══════════════════════════════════════════════════════════════
#  Figura de síntesis transversal↔longitudinal (ex SummaryFigures.jl)
#  — duplicada a propósito en Transversal.jl (nota de diseño §6);
#  siempre lee ambos dominios de disco, independiente de quién dispare.
# ═══════════════════════════════════════════════════════════════

function generate_summary_if_ready(res_root)
    trans_tab = joinpath(res_root, "transversal", "eyesclosed", "tables")
    long_tab  = joinpath(res_root, "longitudinal", "eyesclosed", "tables")
    gpath = joinpath(trans_tab, "global_mean_wpli_statistics.csv")
    spath = joinpath(trans_tab, "subject_band_means.csv")
    if !isfile(gpath) || !isfile(spath)
        println("  ⤷ Síntesis omitida: aún no hay resultados transversales " *
                "(ejecuta después: run_transversal_analysis.jl).")
        return nothing
    end
    gdf = CSV.read(gpath, DataFrame)
    sm  = CSV.read(spath, DataFrame)
    ng  = CSV.read(joinpath(long_tab, "network_global_statistics.csv"), DataFrame)
    scores = CSV.read(joinpath(long_tab, "mean_strength_scores.csv"), DataFrame)

    bands = BAND_ORDER; nB = length(bands)
    dT=zeros(nB); loT=zeros(nB); hiT=zeros(nB); qT=zeros(nB)
    dL=zeros(nB); loL=zeros(nB); hiL=zeros(nB); qL=zeros(nB)
    for (i,b) in enumerate(bands)
        gsub = filter(r -> string(r.band) == b, gdf)
        nrow(gsub) == 1 || error("Fila transversal ausente: $b")
        grow = gsub[1,:]
        ms = Float64[r.mean_wpli for r in eachrow(sm) if string(r.band)==b && string(r.group)=="MS"]
        ct = Float64[r.mean_wpli for r in eachrow(sm) if string(r.band)==b && string(r.group)=="Control"]
        na, nb_ = length(ms), length(ct)
        d0 = (na<2||nb_<2) ? 0.0 : begin
            μa,μb=mean(ms),mean(ct); sa2,sb2=var(ms),var(ct)
            sp = sqrt(((na-1)*sa2+(nb_-1)*sb2)/max(na+nb_-2,1))
            sp > 1e-12 ? (μa-μb)/sp : 0.0
        end
        dT[i] = Float64(grow.effect_d_pooled)
        rng = MersenneTwister(BOOTSTRAP_SEED)
        ds = Vector{Float64}(undef, BOOTSTRAP_N)
        for k in 1:BOOTSTRAP_N
            aa = ms[rand(rng,1:na,na)]; bb = ct[rand(rng,1:nb_,nb_)]
            μa,μb=mean(aa),mean(bb); sa2,sb2=var(aa),var(bb)
            sp = sqrt(((na-1)*sa2+(nb_-1)*sb2)/max(na+nb_-2,1))
            ds[k] = sp > 1e-12 ? (μa-μb)/sp : 0.0
        end
        l0, h0 = quantile(ds, [0.025, 0.975])
        shift = dT[i] - d0; loT[i] = l0 + shift; hiT[i] = h0 + shift
        qT[i] = Float64(grow.q_value)

        lrow = _ng_row(ng, b)
        dL[i] = Float64(lrow.effect_dz)
        loL[i] = Float64(lrow.effect_dz_ci_low)
        hiL[i] = Float64(lrow.effect_dz_ci_high)
        qL[i] = Float64(lrow.q_value)
    end

    y_of = Dict(b => Float64(nB-i+1) for (i,b) in enumerate(bands))
    apply_theme!()
    fig = Figure(size=fig_size_mm(170,90); figure_padding=8)
    axL = Axis(fig[1,1]; title="Transversal · Cohen d (A)\nEM − Control", xlabel="d",
               yticks=(Float64.(nB:-1:1), bands), limits=(nothing,nothing,0.4,nB+0.6))
    vlines!(axL, [0.0]; color=COLOR_ZERO, linewidth=0.9, linestyle=:dash)
    for i in 1:nB
        yi = y_of[bands[i]]; sig = qT[i] < FDR_ALPHA
        col = dT[i] >= 0 ? COLOR_POS : COLOR_NEG
        lines!(axL, [loT[i],hiT[i]], [yi,yi]; color=(col, sig ? 1.0 : 0.45), linewidth=1.6)
        scatter!(axL, [dT[i]], [yi]; color=col, markersize=sig ? 9 : 6,
                 marker=sig ? :diamond : :circle, strokecolor=sig ? :black : :gray40, strokewidth=sig ? 1.0 : 0.4)
    end
    axR = Axis(fig[1,2]; title="Longitudinal · Cohen d_z (C)\nT2 − T1", xlabel="d_z",
               yticks=(Float64.(nB:-1:1), ["" for _ in 1:nB]), limits=(nothing,nothing,0.4,nB+0.6))
    vlines!(axR, [0.0]; color=COLOR_ZERO, linewidth=0.9, linestyle=:dash)
    for i in 1:nB
        yi = y_of[bands[i]]; sig = qL[i] < FDR_ALPHA
        col = dL[i] >= 0 ? COLOR_POS : COLOR_NEG
        lines!(axR, [loL[i],hiL[i]], [yi,yi]; color=(col, sig ? 1.0 : 0.45), linewidth=1.6)
        scatter!(axR, [dL[i]], [yi]; color=col, markersize=sig ? 9 : 6,
                 marker=sig ? :diamond : :circle, strokecolor=sig ? :black : :gray40, strokewidth=sig ? 1.0 : 0.35)
    end
    Label(fig[2,1:2], "Paneles independientes: d (estimando A, montaje nativo) != d_z (estimando C, montaje comun). " *
          "No son el mismo estimando ni admiten prueba directa entre disenos. Rombo: q_global<0.05 (FDR-BH, 7 bandas de cada analisis). " *
          "Mensaje: diferencias transversales EM-Control no implican cambios longitudinales detectables en la cohorte seguida.",
          fontsize=7, color=:gray30, tellwidth=false)

    out_dir = joinpath(res_root, "summary")
    stem = joinpath(out_dir, "summary_transversal_longitudinal_effects")
    path = save_png(fig, stem)
    manifest = [(figure_name="summary_transversal_longitudinal_effects", analysis="summary", band="ALL", condition="EC",
                 output_file=path, source_tables="global_mean_wpli_statistics.csv;network_global_statistics.csv;mean_strength_scores.csv",
                 generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"),
                 estimand="Transversal A (d) | Longitudinal C (d_z) — separate panels",
                 fdr_family="Each panel: FDR-BH across 7 bands of its own analysis")]
    write_manifest(joinpath(out_dir, "figures_manifest.tsv"), manifest)
    println("  Síntesis → $path")
    return path
end

# ═══════════════════════════════════════════════════════════════
#  Pares longitudinales
# ═══════════════════════════════════════════════════════════════

struct SubjPair
    subject_id::String
    session_t1::String
    session_t2::String
    has_ec::Bool
    has_eo::Bool
end

_as_bool(x) = x isa Bool ? x : (lowercase(string(x)) == "true")

function load_pairs(proj, res_root, bids_root)::Vector{SubjPair}
    candidates = [joinpath(bids_root, "longitudinal_pairs.csv"),
                  joinpath(proj, "data", "bids", "longitudinal_pairs.csv"),
                  joinpath(proj, "data", "BIDS", "longitudinal_pairs.csv")]
    pairs_path = ""
    for c in candidates
        isfile(c) && (pairs_path = c; break)
    end
    all_pairs = SubjPair[]
    n_skipped_not_incl = 0; n_skipped_ctrl = 0

    if !isempty(pairs_path)
        pdf = CSV.read(pairs_path, DataFrame)
        rename!(pdf, Dict(n => Symbol(lowercase(string(n))) for n in names(pdf)))
        has_condition_flags = any(hasproperty(pdf, c) for c in
                                  (:has_t1_ec, :has_t2_ec, :has_t1_eo, :has_t2_eo))
        has_bids = hasproperty(pdf, :bids_id)
        has_sess = hasproperty(pdf, :session_t1) && hasproperty(pdf, :session_t2)
        has_incl_col = hasproperty(pdf, :include_longitudinal)
        for row in eachrow(pdf)
            sid = has_bids ? string(row.bids_id) : string(row.subject_id)
            if !is_ms_id(sid); n_skipped_ctrl += 1; continue; end
            sess_t1 = has_sess ? string(row.session_t1) : "T1"
            sess_t2 = has_sess ? string(row.session_t2) : "T2"
            if has_condition_flags
                ht1ec = hasproperty(row,:has_t1_ec) && _as_bool(row.has_t1_ec)
                ht2ec = hasproperty(row,:has_t2_ec) && _as_bool(row.has_t2_ec)
                ht1eo = hasproperty(row,:has_t1_eo) && _as_bool(row.has_t1_eo)
                ht2eo = hasproperty(row,:has_t2_eo) && _as_bool(row.has_t2_eo)
                if !(ht1ec && ht2ec) && !(ht1eo && ht2eo); n_skipped_not_incl += 1; continue; end
                push!(all_pairs, SubjPair(sid, sess_t1, sess_t2, ht1ec && ht2ec, ht1eo && ht2eo))
            else
                if has_incl_col && !_as_bool(row.include_longitudinal)
                    n_skipped_not_incl += 1
                    continue
                end
                push!(all_pairs, SubjPair(sid, sess_t1, sess_t2, true, true))
            end
        end
        println("📋 Pares leídos desde: $pairs_path")
        println("   Comparables (EM, par EC y/o EO): $(length(all_pairs))  |  N diseño=$N_PAIRED_DESIGN")
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
                    has_ec = isdir(joinpath(sess_dir,t1s_,"eyesclosed")) && isdir(joinpath(sess_dir,t2s_,"eyesclosed"))
                    has_eo = isdir(joinpath(sess_dir,t1s_,"eyesopen")) && isdir(joinpath(sess_dir,t2s_,"eyesopen"))
                    (has_ec || has_eo) || continue
                    push!(all_pairs, SubjPair(subj_id, t1s_[5:end], t2s_[5:end], has_ec, has_eo))
                end
            end
        end
        println("🔍 Pares auto-detectados (EM con T1+T2): $(length(all_pairs))  |  N diseño=$N_PAIRED_DESIGN")
    end
    return all_pairs
end

function evaluate_inclusion(res_root, sp::SubjPair, cond::String, qc::Dict)
    cc = cond_code(cond)
    eligible = cc == "EC" ? sp.has_ec : sp.has_eo
    if !eligible
        return (false, "Sin par T1/T2 para $cc en longitudinal_pairs", "n/a", "n/a",
                missing, missing, missing, missing, 0)
    end
    qc_t1 = qc_decision(qc, sp.subject_id, sp.session_t1, cc)
    qc_t2 = qc_decision(qc, sp.subject_id, sp.session_t2, cc)
    ev_t1 = load_seg_epochs(res_root, sp.subject_id, sp.session_t1, cc)
    ev_t2 = load_seg_epochs(res_root, sp.subject_id, sp.session_t2, cc)
    n_bands = 0
    for b in BAND_ORDER
        r1 = load_wpli(res_root, sp.subject_id, sp.session_t1, cc, b)
        r2 = load_wpli(res_root, sp.subject_id, sp.session_t2, cc, b)
        (r1 !== nothing && r2 !== nothing) && (n_bands += 1)
    end
    reasons = String[]
    n_bands == 0 && push!(reasons, "Sin wPLI pareado en disco")
    qc_ok(qc_t1) || push!(reasons, "QC T1=$(qc_t1)")
    qc_ok(qc_t2) || push!(reasons, "QC T2=$(qc_t2)")
    included = isempty(reasons)
    return (included, included ? "" : join(reasons,"; "), qc_t1, qc_t2,
            ev_t1[1], ev_t2[1], ev_t1[3], ev_t2[3], n_bands)
end

# ═══════════════════════════════════════════════════════════════
#  run(config_path) — punto de entrada único
# ═══════════════════════════════════════════════════════════════

function run(config_path::String)
    proj      = dirname(dirname(@__DIR__))
    cfg_raw   = TOML.parsefile(config_path)
    paths_raw = get(cfg_raw, "paths", Dict{String,Any}())
    res_root  = let r = get(paths_raw, "results", "results")
                    isabspath(r) ? r : joinpath(proj, r)
                end
    bids_root = let b = get(paths_raw, "bids_root", "data/bids")
                    p = isabspath(b) ? b : joinpath(proj, b)
                    isdir(p) ? p : (isdir(joinpath(proj,"data","BIDS")) ? joinpath(proj,"data","BIDS") : p)
                end
    bands_cfg = get(cfg_raw, "bands", Dict(
        "DELTA"=>[0.5,4.0], "THETA"=>[4.0,8.0], "ALPHA"=>[7.8,11.7],
        "BETA_LOW"=>[12.0,15.0], "BETA_MID"=>[15.0,18.0], "BETA_HIGH"=>[18.0,30.0], "GAMMA"=>[30.0,50.0]))
    bands = sort(collect(keys(bands_cfg)))

    conn_cfg  = get(cfg_raw, "connectivity", Dict{String,Any}())
    graph_cfg = get(cfg_raw, "graph", Dict{String,Any}())
    wpli_method = String(get(conn_cfg, "wpli_method", "hilbert"))
    use_dwpli   = Bool(get(conn_cfg, "use_dwpli", false))
    graph_dens  = Float64(get(graph_cfg, "density", 0.1))
    graph_meth  = String(get(graph_cfg, "threshold_method", "proportional"))
    n_ch_mont   = Int(get(get(cfg_raw, "montage", Dict()), "n_channels_analysis", 31))

    println("="^62)
    println(" NeuroMIND — Análisis longitudinal T1→T2 (solo EM, N diseño=$N_PAIRED_DESIGN)")
    println(" $(now())")
    println(" Diseño: pares completos; eyesclosed y eyesopen en paralelo (sin pooling)")
    println(" wPLI=$wpli_method  dwPLI=$use_dwpli  test=wilcoxon ($WILCOXON_METHOD)  density=$graph_dens")
    println("="^62)

    all_pairs = load_pairs(proj, res_root, bids_root)
    if isempty(all_pairs)
        @warn "No se encontraron pares T1/T2."
        exit(1)
    end
    qc_table = load_qc_table(res_root)
    println("  QC decisions cargadas: $(length(qc_table))")
    println()

    for cond in ["eyesclosed", "eyesopen"]
        cl = cond_label(cond)
        println("── Condición: $cl  (contraste longitudinal T1→T2) " * "─"^20)
        out_dir  = joinpath(res_root, "longitudinal", cond)
        figs_dir = joinpath(out_dir, "figures")
        tab_dir  = joinpath(out_dir, "tables")
        mkpath(figs_dir); mkpath(tab_dir)

        open(joinpath(out_dir, "config_snapshot.toml"), "w") do io
            println(io, "# NeuroMIND longitudinal snapshot — $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
            println(io, "schema_version = $LONGITUDINAL_SCHEMA_VERSION")
            println(io, "statistics_source = \"$STATISTICS_SOURCE\"")
            println(io, "design = \"longitudinal_MS_T1_T2\"")
            println(io, "cohort = \"MS_paired_only\"")
            println(io, "n_paired_design = $N_PAIRED_DESIGN")
            println(io, "condition = \"$cl\"")
            println(io, "test = \"wilcoxon_signed_rank\"")
            println(io, "wilcoxon_method = \"$WILCOXON_METHOD\"")
            println(io, "fdr = \"bh\"")
            println(io, "fdr_scope_edges = \"$EDGE_FDR_SCOPE\"")
            println(io, "fdr_scope_global = \"$GLOBAL_FDR_SCOPE\"")
            println(io, "fdr_scope_power = \"$POWER_FDR_SCOPE\"")
            println(io, "bootstrap_method = \"$BOOTSTRAP_METHOD\"")
            println(io, "bootstrap_iterations = $BOOTSTRAP_N")
            println(io, "bootstrap_seed = $BOOTSTRAP_SEED")
            println(io, "quantile_method = \"$QUANTILE_METHOD\"")
            println(io, "rrb_method = \"$RRB_METHOD\"")
            println(io, "effect_dz_method = \"$EFFECT_DZ_METHOD\"")
            println(io, "wpli_method = \"$wpli_method\"")
            println(io, "use_dwpli = $use_dwpli")
            println(io, "graph_density = $graph_dens")
            println(io, "graph_threshold_method = \"$graph_meth\"")
            println(io, "n_channels_analysis = $n_ch_mont")
            println(io, "qc_allowed = [\"include\", \"include_with_warning\"]")
            println(io, "source_config = \"$config_path\"")
        end

        paired_info = NamedTuple[]
        included_pairs = SubjPair[]
        condition_pairs = filter(sp -> cl == "EC" ? sp.has_ec : sp.has_eo, all_pairs)
        for sp in condition_pairs
            incl, reason, qc1, qc2, ep1, ep2, pct1, pct2, n_ok =
                evaluate_inclusion(res_root, sp, cond, qc_table)
            if incl
                push!(included_pairs, sp)
            end
            push!(paired_info, (schema_version=LONGITUDINAL_SCHEMA_VERSION,
                  statistics_source=STATISTICS_SOURCE,
                  subject_id=sp.subject_id, session_t1=sp.session_t1, session_t2=sp.session_t2,
                  n_bands_ok=n_ok, included=incl, excluded_reason=reason, qc_t1=qc1, qc_t2=qc2,
                  n_epochs_t1=ep1===missing ? "" : string(ep1), n_epochs_t2=ep2===missing ? "" : string(ep2),
                  retention_pct_t1=pct1===missing ? "" : string(pct1),
                  retention_pct_t2=pct2===missing ? "" : string(pct2)))
        end

        n_paired = count(r -> r.included, paired_info)
        n_excl   = count(r -> !r.included, paired_info)
        println("  Comparación longitudinal $cl | Pares incluidos: $n_paired / diseño $N_PAIRED_DESIGN | Excluidos (QC/datos): $n_excl")
        CSV.write(joinpath(tab_dir, "paired_subjects.csv"), DataFrame(paired_info))

        if n_paired < 1
            println("  ⚠  Análisis omitido: sin pares completos T1/T2 para $cl")
            open(joinpath(out_dir, "longitudinal_summary.json"), "w") do io
                write(io, """{"n_paired":0,"n_paired_design":$N_PAIRED_DESIGN,"n_candidates":$(length(condition_pairs)),"n_t1":0,"n_t2":0,"n_excluded":$n_excl,"n_excluded_qc":0,"n_excluded_data":$n_excl,"n_total_sig":0,"n_bands":0,"cond":"$cl","test":"wilcoxon_signed_rank","wilcoxon_method":"$WILCOXON_METHOD","cohort":"MS_paired_only","timestamp":"$(Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"))"}""")
            end
            continue
        end

        t1_data = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
        t2_data = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
        paired_data = Dict{String, Vector{NamedTuple}}()

        for sp in included_pairs
            for band in bands
                r1 = load_wpli(res_root, sp.subject_id, sp.session_t1, cond, band)
                r2 = load_wpli(res_root, sp.subject_id, sp.session_t2, cond, band)
                (r1 === nothing || r2 === nothing) && continue
                (ch1, W1) = r1; (ch2, W2) = r2
                haskey(t1_data, band) || (t1_data[band] = [])
                haskey(t2_data, band) || (t2_data[band] = [])
                haskey(paired_data, band) || (paired_data[band] = NamedTuple[])
                push!(t1_data[band], (ch1, W1))
                push!(t2_data[band], (ch2, W2))
                push!(paired_data[band], (subject_id=sp.subject_id, ch1=ch1, W1=W1, ch2=ch2, W2=W2))
            end
        end

        band_stats_rows = NamedTuple[]

        for band in bands
            t1_mats = get(t1_data, band, Tuple{Vector{String},Matrix{Float64}}[])
            t2_mats = get(t2_data, band, Tuple{Vector{String},Matrix{Float64}}[])
            pd_mats = get(paired_data, band, NamedTuple[])
            (isempty(t1_mats) || isempty(t2_mats)) && continue
            all_ch_sets = [Set(ch) for (ch,_) in vcat(t1_mats, t2_mats)]
            common_set  = reduce(intersect, all_ch_sets)
            ch_ref = t1_mats[1][1]
            common_ch = filter(c -> c in common_set, ch_ref)
            n = length(common_ch); n < 2 && continue

            get_W(ch, W) = realign_matrix(ch, W, common_ch)
            Wt1 = mean(get_W(ch,W) for (ch,W) in t1_mats)
            Wt2 = mean(get_W(ch,W) for (ch,W) in t2_mats)
            Wdiff = Wt2 .- Wt1

            save_mat_csv(joinpath(tab_dir, "longitudinal_connectivity_t1_$(band).csv"), Wt1, common_ch)
            save_mat_csv(joinpath(tab_dir, "longitudinal_connectivity_t2_$(band).csv"), Wt2, common_ch)
            save_mat_csv(joinpath(tab_dir, "longitudinal_difference_$(band).csv"), Wdiff, common_ch)

            upper_idx = [(i,j) for i in 1:n for j in (i+1):n]
            p_vec = ones(length(upper_idx))
            d_vec = zeros(length(upper_idx)); n_vec = zeros(Int, length(upper_idx))
            method_vec = fill("", length(upper_idx))
            for (k,(i,j)) in enumerate(upper_idx)
                v1 = [get_W(r.ch1,r.W1)[i,j] for r in pd_mats]
                v2 = [get_W(r.ch2,r.W2)[i,j] for r in pd_mats]
                length(v1) < 2 && continue
                n_vec[k] = length(v1)
                p_vec[k], _, method_vec[k] = wilcoxon_p(v1, v2; return_method=true)
                d_vec[k] = cohen_dz(v1, v2)
            end
            q_vec = bh_qvalues(p_vec)
            effect_order = sortperm(collect(eachindex(d_vec));
                by=k -> (-abs(d_vec[k]), common_ch[upper_idx[k][1]], common_ch[upper_idx[k][2]]))
            effect_rank = zeros(Int, length(d_vec))
            for (rank, k) in enumerate(effect_order)
                effect_rank[k] = rank
            end
            stat_rows = [(
                schema_version=LONGITUDINAL_SCHEMA_VERSION,
                statistics_source=STATISTICS_SOURCE,
                fdr_scope=EDGE_FDR_SCOPE,
                fdr_family_size=length(upper_idx),
                ch_a=common_ch[i], ch_b=common_ch[j],
                t1_mean=Wt1[i,j], t2_mean=Wt2[i,j], diff=Wdiff[i,j],
                p_value=p_vec[k], q_value=q_vec[k],
                effect_dz=d_vec[k], effect_rank_abs_dz=effect_rank[k],
                is_nominal=p_vec[k] < FDR_ALPHA, is_fdr=q_vec[k] < FDR_ALPHA,
                n=n_vec[k], p_method=method_vec[k],
            ) for (k,(i,j)) in enumerate(upper_idx)]
            stats_df = DataFrame(stat_rows)
            CSV.write(joinpath(tab_dir, "longitudinal_statistics_$(band).csv"), stats_df)

            sig_df = filter(:is_fdr => identity, stats_df)
            sig_rows = [NamedTuple(r) for r in eachrow(sig_df)]
            CSV.write(joinpath(tab_dir, "significant_longitudinal_edges_$(band).csv"), sig_df)

            n_sig = length(sig_rows)
            n_nominal = count(<(FDR_ALPHA), p_vec)
            n_top20_nominal = count(k -> effect_rank[k] <= 20 && p_vec[k] < FDR_ALPHA,
                                    eachindex(p_vec))
            mean_t1 = mean(Wt1[i,j] for (i,j) in upper_idx)
            mean_t2 = mean(Wt2[i,j] for (i,j) in upper_idx)
            @printf("  %-10s  %3d ch  %4d pares  %3d sig (q<0.05)  T1=%.3f  T2=%.3f\n",
                    band, n, length(upper_idx), n_sig, mean_t1, mean_t2)

            push!(band_stats_rows, (
                  schema_version=LONGITUDINAL_SCHEMA_VERSION,
                  statistics_source=STATISTICS_SOURCE,
                  fdr_scope=EDGE_FDR_SCOPE,
                  band=band, n_channels=n, n_edges=length(upper_idx),
                  n_pairs=length(upper_idx), fdr_family_size=length(upper_idx),
                  n_nominal=n_nominal, n_sig=n_sig,
                  top20_nominal_overlap=n_top20_nominal,
                  pct_sig=100.0*n_sig/max(1,length(upper_idx)), t1_mean=mean_t1,
                  t2_mean=mean_t2, diff_mean=mean_t2-mean_t1,
                  mean_p=mean(p_vec), mean_abs_dz=mean(abs,d_vec),
                  n_subjects=length(pd_mats)))

            try
                save_wpli_heatmaps(figs_dir, band, cl, Wt1, Wt2, Wdiff, common_ch)
                save_sig_or_explore_network(figs_dir, band, cl, common_ch, sig_rows, stats_df)
            catch e
                @warn "Figuras wPLI $band fallidas: $e"
            end

            try
                s1, d1, ns1, _ = nodal_strength_degree(Wt1; density=graph_dens)
                s2, d2, ns2, _ = nodal_strength_degree(Wt2; density=graph_dens)
                CSV.write(joinpath(tab_dir, "network_metrics_t1_$(band).csv"),
                    DataFrame(schema_version=fill(LONGITUDINAL_SCHEMA_VERSION, n),
                              channel=common_ch, strength=s1, degree=d1, norm_strength=ns1))
                CSV.write(joinpath(tab_dir, "network_metrics_t2_$(band).csv"),
                    DataFrame(schema_version=fill(LONGITUDINAL_SCHEMA_VERSION, n),
                              channel=common_ch, strength=s2, degree=d2, norm_strength=ns2))
                CSV.write(joinpath(tab_dir, "network_metrics_delta_$(band).csv"),
                    DataFrame(schema_version=fill(LONGITUDINAL_SCHEMA_VERSION, n),
                              channel=common_ch, delta_strength=s2.-s1,
                              delta_degree=d2.-d1, delta_norm_strength=ns2.-ns1))
            catch e
                @warn "Network $band fallido: $e"
            end
        end

        !isempty(band_stats_rows) && CSV.write(joinpath(tab_dir, "band_statistics_longitudinal.csv"), DataFrame(band_stats_rows))

        # ── Network global + captura de mean_strength por sujeto (sin
        #    reconstrucción aparte: se calcula una única vez, aquí mismo) ──
        net_global = NamedTuple[]
        scores_rows = NamedTuple[]
        for band in bands
            pd = get(paired_data, band, NamedTuple[])
            length(pd) < 2 && continue
            all_ch_sets = [intersect(Set(r.ch1), Set(r.ch2)) for r in pd]
            common_set = reduce(intersect, all_ch_sets)
            ch_ref = pd[1].ch1
            common_ch = filter(c -> c in common_set, ch_ref)
            length(common_ch) < 2 && continue
            n_ch = length(common_ch)
            ms1 = Float64[]; ms2 = Float64[]
            for r in pd
                Wa = realign_matrix(r.ch1, r.W1, common_ch)
                Wb = realign_matrix(r.ch2, r.W2, common_ch)
                s1,_,_,_ = nodal_strength_degree(Wa; density=graph_dens)
                s2,_,_,_ = nodal_strength_degree(Wb; density=graph_dens)
                m1, m2 = mean(s1), mean(s2)
                push!(ms1, m1); push!(ms2, m2)
                push!(scores_rows, (
                      schema_version=LONGITUDINAL_SCHEMA_VERSION,
                      statistics_source=STATISTICS_SOURCE,
                      subject_id=r.subject_id, band=band, cond=cl, n_channels=n_ch,
                      mean_strength_t1=m1, mean_strength_t2=m2,
                      delta_strength=m2-m1,
                      mean_wpli_equiv_t1=m1/(n_ch-1),
                      mean_wpli_equiv_t2=m2/(n_ch-1),
                      delta_wpli_equiv=(m2-m1)/(n_ch-1)))
            end
            pw, _, p_method = wilcoxon_p(ms1, ms2; return_method=true)
            band_idx = something(findfirst(==(band), bands), 0)
            cond_offset = cl == "EC" ? 0 : 10_000
            s = paired_change_summary(ms1, ms2;
                                      seed=BOOTSTRAP_SEED + cond_offset + band_idx)
            denom = n_ch - 1
            push!(net_global, (
                band=band, metric="mean_strength", n_channels=n_ch,
                t1_mean=s.t1_mean, t1_ci_low=s.t1_ci_low, t1_ci_high=s.t1_ci_high,
                t2_mean=s.t2_mean, t2_ci_low=s.t2_ci_low, t2_ci_high=s.t2_ci_high,
                diff=s.diff_mean, diff_ci_low=s.diff_ci_low, diff_ci_high=s.diff_ci_high,
                median_diff=s.median_diff, q1_diff=s.q1_diff, q3_diff=s.q3_diff,
                mean_wpli_t1=s.t1_mean/denom, mean_wpli_t1_ci_low=s.t1_ci_low/denom,
                mean_wpli_t1_ci_high=s.t1_ci_high/denom,
                mean_wpli_t2=s.t2_mean/denom, mean_wpli_t2_ci_low=s.t2_ci_low/denom,
                mean_wpli_t2_ci_high=s.t2_ci_high/denom,
                diff_wpli=s.diff_mean/denom, diff_wpli_ci_low=s.diff_ci_low/denom,
                diff_wpli_ci_high=s.diff_ci_high/denom,
                median_diff_wpli=s.median_diff/denom,
                q1_diff_wpli=s.q1_diff/denom, q3_diff_wpli=s.q3_diff/denom,
                p_value=pw, effect_dz=s.effect_dz,
                effect_dz_ci_low=s.effect_dz_ci_low, effect_dz_ci_high=s.effect_dz_ci_high,
                effect_rrb=s.effect_rrb, n=s.n, n_positive=s.n_positive,
                n_negative=s.n_negative, n_zero=s.n_zero, p_method=p_method,
                bootstrap_seed=BOOTSTRAP_SEED + cond_offset + band_idx,
            ))
        end
        if !isempty(net_global)
            ng_df = DataFrame(net_global)
            ng_df[!, :q_value] = bh_qvalues(Float64.(ng_df.p_value))
            ng_df[!, :fdr_family_size] = fill(nrow(ng_df), nrow(ng_df))
            ng_df[!, :fdr_scope] = fill(GLOBAL_FDR_SCOPE, nrow(ng_df))
            ng_df[!, :schema_version] = fill(LONGITUDINAL_SCHEMA_VERSION, nrow(ng_df))
            ng_df[!, :statistics_source] = fill(STATISTICS_SOURCE, nrow(ng_df))
            ng_df[!, :bootstrap_method] = fill(BOOTSTRAP_METHOD, nrow(ng_df))
            ng_df[!, :bootstrap_iterations] = fill(BOOTSTRAP_N, nrow(ng_df))
            ng_df[!, :quantile_method] = fill(QUANTILE_METHOD, nrow(ng_df))
            ng_df[!, :rrb_method] = fill(RRB_METHOD, nrow(ng_df))
            ng_df[!, :effect_dz_method] = fill(EFFECT_DZ_METHOD, nrow(ng_df))
            ng_df[!, :is_fdr] = Float64.(ng_df.q_value) .< FDR_ALPHA
            ng_df[!, :is_largest_abs_effect] = falses(nrow(ng_df))
            ng_df[argmax(abs.(Float64.(ng_df.effect_dz))), :is_largest_abs_effect] = true
            CSV.write(joinpath(tab_dir, "network_global_statistics.csv"), ng_df)
        end
        !isempty(scores_rows) && CSV.write(joinpath(tab_dir, "mean_strength_scores.csv"), DataFrame(scores_rows))

        # `subject_band_means` y su figura usan el mismo montaje común que el
        # estimando C primario. Así T1 y T2 nunca mezclan conjuntos de canales.
        subject_means = NamedTuple[]
        for r in scores_rows
            push!(subject_means, (subject_id=r.subject_id, timepoint="T1", band=r.band,
                  schema_version=LONGITUDINAL_SCHEMA_VERSION, statistics_source=STATISTICS_SOURCE,
                  cond=r.cond, mean_wpli=Float64(r.mean_wpli_equiv_t1), n_channels=Int(r.n_channels)))
            push!(subject_means, (subject_id=r.subject_id, timepoint="T2", band=r.band,
                  schema_version=LONGITUDINAL_SCHEMA_VERSION, statistics_source=STATISTICS_SOURCE,
                  cond=r.cond, mean_wpli=Float64(r.mean_wpli_equiv_t2), n_channels=Int(r.n_channels)))
        end
        !isempty(subject_means) &&
            CSV.write(joinpath(tab_dir, "subject_band_means.csv"), DataFrame(subject_means))

        bp_store = Dict{Tuple{String,String}, Tuple{Vector{Float64},Vector{Float64}}}()
        for sp in included_pairs
            bp1 = load_band_power(res_root, sp.subject_id, sp.session_t1, cond, bands)
            bp2 = load_band_power(res_root, sp.subject_id, sp.session_t2, cond, bands)
            (bp1 === nothing || bp2 === nothing) && continue
            for ch in intersect(keys(bp1), keys(bp2)), band in bands
                haskey(bp1[ch], band) && haskey(bp2[ch], band) || continue
                key = (ch, band)
                haskey(bp_store, key) || (bp_store[key] = (Float64[], Float64[]))
                push!(bp_store[key][1], bp1[ch][band]); push!(bp_store[key][2], bp2[ch][band])
            end
        end
        bp_stat_rows = NamedTuple[]
        for ((ch,band), (v1,v2)) in bp_store
            length(v1) < 2 && continue
            pw, _, p_method = wilcoxon_p(v1, v2; return_method=true)
            dz = cohen_dz(v1, v2)
            push!(bp_stat_rows, (schema_version=LONGITUDINAL_SCHEMA_VERSION,
                  statistics_source=STATISTICS_SOURCE, fdr_scope=POWER_FDR_SCOPE,
                  channel=ch, band=band, t1_mean=mean(v1), t2_mean=mean(v2),
                  diff=mean(v2)-mean(v1), p_value=pw, effect_dz=dz, n=length(v1),
                  p_method=p_method))
        end
        if !isempty(bp_stat_rows)
            bp_df = DataFrame(bp_stat_rows)
            q_all = fill(1.0, nrow(bp_df))
            family_size = zeros(Int, nrow(bp_df))
            coverage_threshold = zeros(Int, nrow(bp_df))
            for band in bands
                idx = findall(i -> string(bp_df.band[i]) == band, 1:nrow(bp_df))
                isempty(idx) && continue
                q_all[idx] = bh_qvalues(Float64.(bp_df.p_value[idx]))
                family_size[idx] .= length(idx)
                coverage_threshold[idx] .= ceil(Int, 0.7 * maximum(Int.(bp_df.n[idx])))
            end
            bp_df[!, :q_value] = q_all
            bp_df[!, :fdr_family_size] = family_size
            bp_df[!, :coverage_threshold_n] = coverage_threshold
            bp_df[!, :coverage_pct] = 100.0 .* Int.(bp_df.n) ./ max(n_paired, 1)
            bp_df[!, :coverage_low] = Int.(bp_df.n) .< coverage_threshold
            bp_df[!, :is_fdr] = Float64.(bp_df.q_value) .< FDR_ALPHA
            CSV.write(joinpath(tab_dir, "band_power_delta_statistics.csv"), bp_df)
            CSV.write(joinpath(tab_dir, "significant_band_power_changes.csv"), filter(r -> r.q_value < 0.05, bp_df))

            xy = isempty(included_pairs) ? Dict{String,Tuple{Float64,Float64}}() :
                 load_electrode_xy(bids_root, included_pairs[1].subject_id, included_pairs[1].session_t1)
            for band in bands
                sub = filter(r -> string(r.band) == band, eachrow(bp_df))
                isempty(sub) && continue
                chs = [string(r.channel) for r in sub]
                dvs = [Float64(r.diff) for r in sub]
                qvs = [Float64(r.q_value) for r in sub]
                nvs = [Int(r.n) for r in sub]
                try
                    save_topo_delta(joinpath(figs_dir, "topo_delta_bandpower_$(band).png"), chs, dvs, xy,
                                    "Δ band power — $band ($cl)"; colorbar_label="Δ power (T2−T1)",
                                    q_values=qvs, n_values=nvs)
                catch e
                    @warn "Topo band power $band fallido: $e"
                end
            end
            try
                fig_power_effect_heatmap(bp_df, figs_dir, cl)
            catch e
                @warn "Heatmap canal×banda de potencia fallido: $e"
            end
        end

        if !isempty(subject_means)
            try
                save_paired_means(joinpath(figs_dir, "paired_mean_wpli_by_band.png"), DataFrame(subject_means))
            catch e
                @warn "Figura paired means fallida: $e"
            end
        end

        n_total_sig = isempty(band_stats_rows) ? 0 : sum(r.n_sig for r in band_stats_rows)
        best_band = honest_best_band(band_stats_rows)
        n_qc_loss = count(r -> !r.included && occursin("QC", r.excluded_reason), paired_info)
        n_data_loss = n_excl - n_qc_loss

        open(joinpath(out_dir, "longitudinal_summary.json"), "w") do io
            pairs_json = join(["  \"$(r.band)_n_sig\": $(r.n_sig), \"$(r.band)_t1\": $(r.t1_mean), \"$(r.band)_t2\": $(r.t2_mean)"
                                for r in band_stats_rows], ",\n")
            write(io, """{
  "schema_version": $LONGITUDINAL_SCHEMA_VERSION,
  "statistics_source": "$STATISTICS_SOURCE",
  "n_paired": $n_paired,
  "n_paired_design": $N_PAIRED_DESIGN,
  "n_candidates": $(length(condition_pairs)),
  "n_t1": $n_paired,
  "n_t2": $n_paired,
  "n_excluded": $n_excl,
  "n_excluded_qc": $n_qc_loss,
  "n_excluded_data": $n_data_loss,
  "n_total_sig": $n_total_sig,
  "n_bands": $(length(band_stats_rows)),
  "best_band": "$best_band",
  "cond": "$cl",
  "cohort": "MS_paired_only",
  "design": "longitudinal_MS_T1_T2",
  "test": "wilcoxon_signed_rank",
  "wilcoxon_method": "$WILCOXON_METHOD",
  "fdr": "bh",
  "fdr_scope_edges": "$EDGE_FDR_SCOPE",
  "fdr_scope_global": "$GLOBAL_FDR_SCOPE",
  "fdr_scope_power": "$POWER_FDR_SCOPE",
  "bootstrap_method": "$BOOTSTRAP_METHOD",
  "bootstrap_iterations": $BOOTSTRAP_N,
  "bootstrap_seed": $BOOTSTRAP_SEED,
  "quantile_method": "$QUANTILE_METHOD",
  "rrb_method": "$RRB_METHOD",
  "effect_dz_method": "$EFFECT_DZ_METHOD",
  "wpli_method": "$wpli_method",
  "use_dwpli": $use_dwpli,
  "graph_density": $graph_dens,
  "timestamp": "$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))",
  $pairs_json
}""")
        end
        write_statistics_contract(joinpath(out_dir, "statistics_contract.toml"), cl)

        # ── Figuras de manuscrito + síntesis (solo eyesclosed) ───
        if cond == "eyesclosed"
            manifest = NamedTuple[]
            try
                fig_forest_by_band_dz(tab_dir, figs_dir, manifest)
                fig_alpha_delta_pareado(tab_dir, figs_dir, manifest)
                fig_matrices_alpha(tab_dir, figs_dir, manifest)
                write_manifest(joinpath(tab_dir, "figures_manifest.tsv"), manifest)
                println("  Figuras de manuscrito ($(length(manifest))): " * join([r.figure_name for r in manifest], ", "))
            catch e
                @warn "Figuras de manuscrito longitudinal fallidas: $e"
            end
            try
                generate_summary_if_ready(res_root)
            catch e
                @warn "Figura de síntesis fallida: $e"
            end
        end

        println("  ✅ Guardado en: $out_dir\n")
    end

    println("✅ Análisis longitudinal completado (eyesclosed y eyesopen en paralelo).")
    println("   Abre el dashboard → Fase 14 Evaluación Longitudinal")
    return nothing
end

end # module
