# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Transversal
#  Análisis EM vs Control (T1) + generación de todas las figuras,
#  en una única pasada: CSV → figura PNG, con datos reales recién
#  calculados (sin regenerar a partir de una reconstrucción aparte).
# ═══════════════════════════════════════════════════════════════
#
#  Sustituye: scripts/run_transversal_analysis.jl (lógica) +
#  src/transversal/{TransversalFigures,TransversalManuscript}.jl +
#  src/visualization/{GroupVizCommon,PublicationCommon,PublicationTheme}.jl
#  (solo la porción que usaba el lado transversal).
#
#  Diseño experimental (Fig. 3.1):
#    · Caso-control en T1: EM (N diseño=44) vs Control (N diseño=40)
#    · Sesión T2 NO entra (pertenece al longitudinal)
#    · eyesclosed y eyesopen en paralelo: mismo contraste, sin pooling
#
#  Salida:
#    results/transversal/{eyesclosed|eyesopen}/
#      config_snapshot.toml · statistics_contract.toml · transversal_summary.json
#      tables/   — todos los CSV (estadística edge-wise, red, potencia,
#                  inclusión, medias por sujeto, figures_manifest.tsv)
#      figures/  — todas las figuras PNG: exploratorias por banda (heatmaps,
#                  red FDR/Top-N, raincloud por banda, heatmap de efecto
#                  canal×banda de potencia) +, solo en eyesclosed,
#                  forest/raincloud-ALPHA/matrices/red+potencia (manuscrito)
#    results/transversal/combined/
#      tables/, figures/ — interacción grupo×condición (EC−EO), calculada
#      sobre los sujetos con datos válidos en ambas condiciones (join de
#      subject_band_means.csv de eyesclosed y eyesopen, generado una vez
#      terminadas ambas pasadas del bucle principal)
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/transversal/Transversal.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      27-07-2026
#  Modificado  28-07-2026
# ───────────────────────────────────────────────────────────────

module Transversal

using CSV, DataFrames, Statistics, LinearAlgebra, Dates, TOML, Printf, CairoMakie, Random

export run

# ═══════════════════════════════════════════════════════════════
#  Constantes compartidas (duplicadas con Longitudinal.jl a propósito
#  — ver nota de diseño en el plan: datos de referencia, no lógica)
# ═══════════════════════════════════════════════════════════════

const BAND_ORDER = ["DELTA", "THETA", "ALPHA", "BETA_LOW", "BETA_MID", "BETA_HIGH", "GAMMA"]
const DIVERGING_CMAP = Reverse(:RdBu)   # positivo (MS−Control) = rojo

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
const TRANSVERSAL_SCHEMA_VERSION = 2
const STATISTICS_SOURCE = "NeuroMIND.Transversal.production"
const BOOTSTRAP_METHOD = "stratified_group_percentile_bootstrap"
const QUANTILE_METHOD = "Statistics.quantile_linear_alpha_1_beta_1"
const RRB_METHOD = "2*U_MS/(n_ms*n_ctrl)-1"
const EFFECT_D_POOLED_METHOD = "(mean_MS-mean_Control)/pooled_sample_sd"
const MANNWHITNEY_METHOD = "two_sided_normal_approximation_tie_corrected"
const EDGE_FDR_SCOPE = "available_edges_within_each_band"
const GLOBAL_FDR_SCOPE = "available_frequency_bands"
const POWER_FDR_SCOPE = "available_channels_within_each_band"

const N_MS_DESIGN   = 44   # Fig. 3.1
const N_CTRL_DESIGN = 40
const QC_ALLOWED    = Set(["include", "include_with_warning"])

_ch_xy(ch::AbstractString) = get(CH_POS, uppercase(String(ch)), (0.0, 0.0))

# ── Tema único (ex PublicationTheme.jl), aplicado a TODA figura ────
#  Paleta Okabe–Ito (apta daltonismo); CairoMakie acepta strings hex directamente.

const COLOR_POS  = "#D55E00"   # vermillion (Okabe–Ito, apto daltonismo)
const COLOR_NEG  = "#0072B2"   # blue
const COLOR_NS   = "#999999"   # gray
const COLOR_MS   = "#D55E00"
const COLOR_CTRL = "#56B4E9"   # sky
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
#  Estadística (idéntica a la de run_transversal_analysis.jl)
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
is_ms_group(g::AbstractString)::Bool =
    uppercase(String(g)) in ("MS", "EM", "PATIENT", "PATIENTS", "CASE", "CASES")
is_ctrl_group(g::AbstractString)::Bool =
    uppercase(String(g)) in ("CONTROL", "CONTROLS", "HC", "HEALTHY", "CTRL")

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

"""Welch t-test → (p_value, cohen_d)."""
function welch_t(a::Vector{Float64}, b::Vector{Float64})
    na, nb = length(a), length(b)
    (na < 2 || nb < 2) && return (1.0, 0.0)
    μa, μb = mean(a), mean(b)
    sa2, sb2 = var(a), var(b)
    se = sqrt(sa2 / na + sb2 / nb)
    se < 1e-12 && return (1.0, 0.0)
    t  = (μa - μb) / se
    p  = clamp(2.0 * (1.0 - _norm_cdf(abs(t))), 0.0, 1.0)
    sp = sqrt(((na - 1) * sa2 + (nb - 1) * sb2) / max(na + nb - 2, 1))
    d  = sp > 1e-12 ? (μa - μb) / sp : 0.0
    return (p, d)
end

"""Mann–Whitney U (two-sided, normal approx + corrección de empates) → (p, r_rb)."""
function mannwhitney_p(a::Vector{Float64}, b::Vector{Float64})
    na, nb = length(a), length(b)
    (na < 2 || nb < 2) && return (1.0, 0.0)
    allv = vcat(a, b)
    ranks = _assign_ranks(allv)
    Ra = sum(ranks[1:na])
    Ua = Ra - na * (na + 1) / 2
    Ub = na * nb - Ua
    U  = min(Ua, Ub)
    μU = na * nb / 2
    σ2 = na * nb * (na + nb + 1) / 12
    uniq = unique(allv)
    if length(uniq) < length(allv)
        tie_term = 0.0
        for u in uniq
            t = count(==(u), allv)
            t > 1 && (tie_term += t^3 - t)
        end
        σ2 -= na * nb * tie_term / (12 * (na + nb) * (na + nb - 1))
    end
    σU = sqrt(max(σ2, 0.0))
    σU < 1e-12 && return (1.0, 0.0)
    z = (U - μU) / σU
    p = clamp(2.0 * (1.0 - _norm_cdf(abs(z))), 0.0, 1.0)
    r_rb = clamp(2.0 * Ua / (na * nb) - 1.0, -1.0, 1.0)
    return (p, r_rb)
end

"""IC95% bootstrap de Cohen d (a−b), remuestreo estratificado por grupo."""
function bootstrap_cohen_d_ci(a::Vector{Float64}, b::Vector{Float64};
                              n_boot::Int=BOOTSTRAP_N, seed::Int=BOOTSTRAP_SEED, alpha=0.05)
    _, d0 = welch_t(a, b)
    na, nb = length(a), length(b)
    (na < 2 || nb < 2) && return (d0, NaN, NaN)
    rng = MersenneTwister(seed)
    ds = Vector{Float64}(undef, n_boot)
    for i in 1:n_boot
        _, ds[i] = welch_t(a[rand(rng, 1:na, na)], b[rand(rng, 1:nb, nb)])
    end
    lo, hi = quantile(ds, [alpha/2, 1-alpha/2])
    return (d0, lo, hi)
end

"""IC95% bootstrap de la diferencia de medias bruta (mean(a)−mean(b)), misma
unidad que el dato de entrada (p. ej. wPLI) — a diferencia de
`bootstrap_cohen_d_ci`, que devuelve el IC en escala de Cohen d estandarizada."""
function bootstrap_mean_diff_ci(a::Vector{Float64}, b::Vector{Float64};
                                n_boot::Int=BOOTSTRAP_N, seed::Int=BOOTSTRAP_SEED, alpha=0.05)
    na, nb = length(a), length(b)
    diff0 = (na>0 && nb>0) ? mean(a) - mean(b) : NaN
    (na < 2 || nb < 2) && return (diff0, NaN, NaN)
    rng = MersenneTwister(seed)
    ds = Vector{Float64}(undef, n_boot)
    for i in 1:n_boot
        ds[i] = mean(a[rand(rng, 1:na, na)]) - mean(b[rand(rng, 1:nb, nb)])
    end
    lo, hi = quantile(ds, [alpha/2, 1-alpha/2])
    return (diff0, lo, hi)
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

function independent_group_summary(ms::Vector{Float64}, ctrl::Vector{Float64};
                                   n_boot::Int=BOOTSTRAP_N, seed::Int=BOOTSTRAP_SEED)
    (length(ms) >= 2 && length(ctrl) >= 2) ||
        error("independent_group_summary: se requieren al menos dos sujetos por grupo")
    p_mw, rrb = mannwhitney_p(ms, ctrl)
    p_welch, effect_d = welch_t(ms, ctrl)
    _, d_lo, d_hi = bootstrap_cohen_d_ci(ms, ctrl; n_boot=n_boot, seed=seed)
    _, diff_lo, diff_hi = bootstrap_mean_diff_ci(ms, ctrl; n_boot=n_boot, seed=seed + 101)
    ms_mean, ms_lo, ms_hi = mean_ci_bootstrap(ms; n_boot=n_boot, seed=seed + 202)
    ctrl_mean, ctrl_lo, ctrl_hi = mean_ci_bootstrap(ctrl; n_boot=n_boot, seed=seed + 303)
    ms_med, ms_q1, ms_q3 = median_iqr(ms)
    ctrl_med, ctrl_q1, ctrl_q3 = median_iqr(ctrl)
    return (
        ms_mean=ms_mean, ms_ci_low=ms_lo, ms_ci_high=ms_hi,
        ms_sem=std(ms) / sqrt(length(ms)),
        ms_median=ms_med, ms_q1=ms_q1, ms_q3=ms_q3,
        ctrl_mean=ctrl_mean, ctrl_ci_low=ctrl_lo, ctrl_ci_high=ctrl_hi,
        ctrl_sem=std(ctrl) / sqrt(length(ctrl)),
        ctrl_median=ctrl_med, ctrl_q1=ctrl_q1, ctrl_q3=ctrl_q3,
        diff=ms_mean-ctrl_mean, diff_ci_low=diff_lo, diff_ci_high=diff_hi,
        p_mannwhitney=p_mw, p_welch=p_welch,
        effect_d_pooled=effect_d,
        effect_d_pooled_ci_low=d_lo, effect_d_pooled_ci_high=d_hi,
        effect_rrb=rrb, probability_superiority=(rrb + 1.0) / 2.0,
        n_ms=length(ms), n_ctrl=length(ctrl),
    )
end

function write_statistics_contract(path::String, condition::String)
    contract = Dict{String,Any}(
        "schema_version" => TRANSVERSAL_SCHEMA_VERSION,
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
        "effect_d_pooled_method" => EFFECT_D_POOLED_METHOD,
        "mannwhitney_method" => MANNWHITNEY_METHOD,
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

function load_qc_notes(res_root)::Dict{Tuple{String,String,String}, String}
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
        parts = String[]
        notes = hasproperty(row, :notes) ? string(something(row.notes, "")) : ""
        excl  = hasproperty(row, :exclusion_reason) ? string(something(row.exclusion_reason, "")) : ""
        amp   = hasproperty(row, :amplitude_warning) && (row.amplitude_warning === true || string(row.amplitude_warning) == "true")
        nbad  = hasproperty(row, :n_bad_channels_non_fp2) ? something(tryparse(Int, string(row.n_bad_channels_non_fp2)), 0) : 0
        badch = hasproperty(row, :bad_channels_non_fp2) ? string(something(row.bad_channels_non_fp2, "")) : ""
        !isempty(notes) && notes != "missing" && push!(parts, notes)
        !isempty(excl)  && excl  != "missing" && push!(parts, excl)
        amp && !any(contains(p, "amplitude") for p in parts) && push!(parts, "amplitude_warning")
        nbad > 0 && !isempty(badch) && !any(contains(p, badch) for p in parts) && push!(parts, "bad_ch=$badch")
        out[(sid, sess, cond)] = join(unique(parts), "; ")
    end
    return out
end

function qc_note(notes::Dict, sid::String, sess::String, cond::String)::String
    cc = cond_code(cond)
    haskey(notes, (sid, sess, cc)) && return notes[(sid, sess, cc)]
    if startswith(sid, "M") && !startswith(sid, "MC")
        alt = occursin(r"^M\d$", sid) ? "M0"*sid[2:end] : (occursin(r"^M0\d$", sid) ? "M"*sid[3:end] : "")
        !isempty(alt) && haskey(notes, (alt, sess, cc)) && return notes[(alt, sess, cc)]
    end
    return ""
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
#  Primitivas de figura (duplicadas de src/longitudinal/Longitudinal.jl
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

function save_group_means(path, subject_means_df; band_order=BAND_ORDER)
    isempty(subject_means_df) && return nothing
    raw = unique(string.(subject_means_df.band))
    bands = sort_bands_physio(intersect(band_order, raw)); isempty(bands) && (bands = sort_bands_physio(raw))
    fig = Figure(size=(900, 420), fontsize=11)
    ax = Axis(fig[1,1]; title="Mean wPLI por grupo (MS vs Control)", xlabel="Banda",
              ylabel="Mean wPLI (triángulo superior)", xticks=(1:length(bands), bands))
    for (bi, b) in enumerate(bands)
        sub = filter(r -> string(r.band) == b, eachrow(subject_means_df))
        ms_v = [Float64(r.mean_wpli) for r in sub if lowercase(string(r.group)) in ("ms","em","patient")]
        ct_v = [Float64(r.mean_wpli) for r in sub if lowercase(string(r.group)) in ("control","ctrl","hc")]
        for v in ms_v; scatter!(ax, [bi-0.15], [v]; color=(:firebrick,0.55), markersize=7); end
        for v in ct_v; scatter!(ax, [bi+0.15], [v]; color=(:steelblue,0.55), markersize=7); end
        if !isempty(ms_v)
            m = mean(ms_v); s = length(ms_v)>1 ? std(ms_v)/sqrt(length(ms_v)) : 0.0
            lines!(ax, [bi-0.22,bi-0.08], [m,m]; color=:firebrick, linewidth=2.5)
            lines!(ax, [bi-0.15,bi-0.15], [m-s,m+s]; color=:firebrick, linewidth=1.5)
        end
        if !isempty(ct_v)
            m = mean(ct_v); s = length(ct_v)>1 ? std(ct_v)/sqrt(length(ct_v)) : 0.0
            lines!(ax, [bi+0.08,bi+0.22], [m,m]; color=:steelblue, linewidth=2.5)
            lines!(ax, [bi+0.15,bi+0.15], [m-s,m+s]; color=:steelblue, linewidth=1.5)
        end
    end
    mkpath(dirname(path)); save(path, fig); return path
end

function save_topo_delta(path, ch, delta_vals, xy, title; colorbar_label="Δ")
    fig = Figure(size=(560, 520), fontsize=11)
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
    Colorbar(fig[1,2], sc; label=colorbar_label)
    mkpath(dirname(path)); save(path, fig); return path
end

function save_wpli_heatmaps(figs_dir, band, cond_lbl, Wctrl, Wms, Wdiff, common_ch)
    lim_ab = max(maximum(Wms), maximum(Wctrl), 1e-12)
    save_heatmap(joinpath(figs_dir, "heatmap_ms_$(band).png"), Wms, common_ch,
                 "wPLI MS — $band ($cond_lbl)"; colorrange=(0.0,lim_ab), colorbar_label="wPLI")
    save_heatmap(joinpath(figs_dir, "heatmap_control_$(band).png"), Wctrl, common_ch,
                 "wPLI Control — $band ($cond_lbl)"; colorrange=(0.0,lim_ab), colorbar_label="wPLI")
    save_heatmap(joinpath(figs_dir, "heatmap_diff_$(band).png"), Wdiff, common_ch,
                 "Δ wPLI (MS−Control) — $band ($cond_lbl)"; diverging=true, colorbar_label="Δ wPLI")
    save_heatmap_triplet(joinpath(figs_dir, "heatmap_triplet_$(band).png"), Wctrl, Wms, Wdiff, common_ch,
        ("Control — $band ($cond_lbl)", "MS — $band ($cond_lbl)", "Δ (MS−Ctrl) — $band ($cond_lbl)"))
    return nothing
end

function save_sig_or_explore_network(figs_dir, band, cond_lbl, common_ch, sig_rows, stats_df)
    if !isempty(sig_rows)
        save_topo_network(joinpath(figs_dir, "sig_network_$(band).png"), common_ch, sig_rows,
                          "Edges sig. FDR — $band ($cond_lbl)")
    elseif nrow(stats_df) > 0 && hasproperty(stats_df, :effect_d_pooled)
        top = first(sort(stats_df, :effect_rank_abs_d), min(20, nrow(stats_df)))
        save_topo_network(joinpath(figs_dir, "explore_network_topN_$(band).png"), common_ch,
                          [NamedTuple(r) for r in eachrow(top)],
                          "Top-20 |d| (exploratorio) — $band ($cond_lbl)")
    end
    return nothing
end

# ═══════════════════════════════════════════════════════════════
#  Figuras de manuscrito (ex TransversalManuscript.jl) — solo eyesclosed,
#  PNG únicamente, mismas `figures/` que las exploratorias
# ═══════════════════════════════════════════════════════════════

function _band_row(gdf::DataFrame, band::AbstractString)
    rows = filter(r -> string(r.band) == band, gdf)
    nrow(rows) == 1 || error("Fila global ausente para banda $band")
    return rows[1, :]
end

function _subject_band_vectors(sm::DataFrame, band::AbstractString)
    sub = filter(r -> string(r.band) == band, sm)
    ms = Float64[r.mean_wpli for r in eachrow(sub) if string(r.group) == "MS"]
    ct = Float64[r.mean_wpli for r in eachrow(sub) if string(r.group) == "Control"]
    return ms, ct
end

function fig_forest_by_band(tab_dir, figs_dir, manifest)
    gdf = CSV.read(joinpath(tab_dir, "global_mean_wpli_statistics.csv"), DataFrame)
    sm  = CSV.read(joinpath(tab_dir, "subject_band_means.csv"), DataFrame)
    bands = BAND_ORDER; nB = length(bands)
    d=zeros(nB); lo=zeros(nB); hi=zeros(nB); p=zeros(nB); q=zeros(nB); rrb=zeros(nB)
    nms=zeros(Int,nB); nct=zeros(Int,nB)
    for (i,b) in enumerate(bands)
        row = _band_row(gdf, b)
        d[i] = Float64(row.effect_d_pooled)
        lo[i] = Float64(row.effect_d_pooled_ci_low)
        hi[i] = Float64(row.effect_d_pooled_ci_high)
        p[i] = Float64(row.p_mannwhitney); q[i] = Float64(row.q_value)
        rrb[i] = Float64(row.effect_rrb)
        nms[i] = Int(row.n_ms); nct[i] = Int(row.n_ctrl)
    end
    y_of = Dict(b => Float64(nB-i+1) for (i,b) in enumerate(bands))
    apply_theme!()
    fig = Figure(size=fig_size_mm(170,95); figure_padding=8)
    ax = Axis(fig[1,1]; xlabel="Cohen d (EM − Control)", yticks=(Float64.(nB:-1:1), bands),
              title="Efectos globales wPLI por banda — estimando A", limits=(nothing,nothing,0.4,nB+0.6))
    vlines!(ax, [0.0]; color=COLOR_ZERO, linewidth=0.9, linestyle=:dash)
    for i in 1:nB
        yi = y_of[bands[i]]; sig = q[i] < FDR_ALPHA
        col = sig ? (d[i]>=0 ? COLOR_POS : COLOR_NEG) : COLOR_NS
        lines!(ax, [lo[i],hi[i]], [yi,yi]; color=col, linewidth=sig ? 2.2 : 1.4)
        scatter!(ax, [d[i]], [yi]; color=col, markersize=sig ? 9 : 7,
                 marker=sig ? :diamond : :circle, strokecolor=:black, strokewidth=0.4)
    end
    ax2 = Axis(fig[1,2]; limits=(0,1,0.4,nB+0.6), xgridvisible=false, ygridvisible=false)
    hidedecorations!(ax2); hidespines!(ax2)
    for i in 1:nB
        yi = y_of[bands[i]]
        txt = "n=$(nms[i])/$(nct[i])  p=$(fmt_p(p[i]))  q=$(fmt_q(q[i]))  r_rb=$(fmt_d(rrb[i]))"
        text!(ax2, 0.02, yi; text=txt, align=(:left,:center), fontsize=7, color=q[i]<FDR_ALPHA ? :black : :gray40)
    end
    colsize!(fig.layout,1,Relative(0.55)); colsize!(fig.layout,2,Relative(0.45))
    Label(fig[2,1:2], "FDR-BH entre las 7 bandas. Color intenso solo si q<0.05. Barras: IC95% bootstrap (B=$BOOTSTRAP_N, semilla=$BOOTSTRAP_SEED).",
          fontsize=7, color=:gray35, tellwidth=false)
    path = save_png(fig, joinpath(figs_dir, "transversal_global_effects_by_band"))
    push!(manifest, (figure_name="transversal_global_effects_by_band", analysis="transversal", band="ALL",
          condition="EC", output_file=path, source_tables="global_mean_wpli_statistics.csv;subject_band_means.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="A (mean wPLI, native montage)",
          fdr_family="FDR-BH across 7 bands (global)"))
    return path
end

function fig_raincloud_alpha(tab_dir, figs_dir, manifest)
    gdf = CSV.read(joinpath(tab_dir, "global_mean_wpli_statistics.csv"), DataFrame)
    sm  = CSV.read(joinpath(tab_dir, "subject_band_means.csv"), DataFrame)
    row = _band_row(gdf, "ALPHA")
    ms, ct = _subject_band_vectors(sm, "ALPHA")
    d_tab = Float64(row.effect_d_pooled)
    dlo = Float64(row.effect_d_pooled_ci_low)
    dhi = Float64(row.effect_d_pooled_ci_high)
    rrb = Float64(row.effect_rrb)
    psup = Float64(row.probability_superiority)

    apply_theme!()
    fig = Figure(size=fig_size_mm(95,110); figure_padding=10)
    ax = Axis(fig[1,1]; xlabel="Grupo", ylabel="wPLI medio por participante — estimando A",
              title="ALPHA · distribución individual",
              xticks=([1.0,2.0], ["Control\n(n=$(length(ct)))", "EM\n(n=$(length(ms)))"]))
    function _density_poly(vals)
        n = length(vals); n < 2 && return Float64[], Float64[]
        lo, hi = extrema(vals); pad = 0.05*(hi-lo+eps())
        ys = collect(range(lo-pad, hi+pad; length=60))
        bw = 1.06*std(vals)*n^(-1/5); bw < 1e-8 && (bw = 1e-3)
        dens = [mean(exp(-0.5*((y-v)/bw)^2) for v in vals)/(bw*sqrt(2π)) for y in ys]
        dens ./= (maximum(dens)+eps())
        return ys, dens .* 0.38
    end
    ys_c, dens_c = _density_poly(ct); ys_m, dens_m = _density_poly(ms)
    if !isempty(ys_c)
        poly!(ax, Point2f.(1.0 .- dens_c, ys_c); color=(COLOR_CTRL,0.35), strokewidth=0)
        lines!(ax, 1.0 .- dens_c, ys_c; color=COLOR_CTRL, linewidth=1)
    end
    if !isempty(ys_m)
        poly!(ax, Point2f.(2.0 .+ dens_m, ys_m); color=(COLOR_MS,0.35), strokewidth=0)
        lines!(ax, 2.0 .+ dens_m, ys_m; color=COLOR_MS, linewidth=1)
    end
    rng = MersenneTwister(BOOTSTRAP_SEED)
    for (x0, vals, col) in ((1.0,ct,COLOR_CTRL), (2.0,ms,COLOR_MS))
        jit = 0.06 .* (rand(rng, length(vals)) .- 0.5)
        scatter!(ax, fill(x0,length(vals)).+jit, vals; color=(col,0.75), markersize=5, strokewidth=0.3, strokecolor=:black)
        is_ctrl = x0 == 1.0
        med = Float64(is_ctrl ? row.ctrl_median : row.ms_median)
        q1 = Float64(is_ctrl ? row.ctrl_q1 : row.ms_q1)
        q3 = Float64(is_ctrl ? row.ctrl_q3 : row.ms_q3)
        lines!(ax, [x0-0.12,x0+0.12], [med,med]; color=:black, linewidth=2)
        lines!(ax, [x0,x0], [q1,q3]; color=:black, linewidth=1.4)
        μ = Float64(is_ctrl ? row.ctrl_mean : row.ms_mean)
        mlo = Float64(is_ctrl ? row.ctrl_ci_low : row.ms_ci_low)
        mhi = Float64(is_ctrl ? row.ctrl_ci_high : row.ms_ci_high)
        lines!(ax, [x0+0.18,x0+0.18], [mlo,mhi]; color=:gray40, linewidth=1)
        scatter!(ax, [x0+0.18], [μ]; color=:gray40, markersize=5, marker=:rect)
    end
    ann = "p=$(fmt_p(Float64(row.p_mannwhitney)))  q_global=$(fmt_q(Float64(row.q_value)))\n" *
          "d=$(fmt_d(d_tab)) [IC95% $(fmt_d(dlo)), $(fmt_d(dhi))]  r_rb=$(fmt_d(rrb))  P(EM>Ctrl)=$(fmt_d(psup))\n" *
          "Resumen principal: mediana±IQR; media±IC95% bootstrap (secundario)."
    Label(fig[2,1], ann; fontsize=7, color=:gray30, tellwidth=false)
    path = save_png(fig, joinpath(figs_dir, "transversal_alpha_individual_distribution"))
    push!(manifest, (figure_name="transversal_alpha_individual_distribution", analysis="transversal", band="ALPHA",
          condition="EC", output_file=path, source_tables="subject_band_means.csv;global_mean_wpli_statistics.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="A (mean wPLI, native montage)",
          fdr_family="FDR-BH across 7 bands (global q)"))
    return path
end

function fig_matrices_alpha(tab_dir, figs_dir, manifest)
    rc = read_mat_csv(joinpath(tab_dir, "group_connectivity_control_ALPHA.csv"))
    rm = read_mat_csv(joinpath(tab_dir, "group_connectivity_ms_ALPHA.csv"))
    rd = read_mat_csv(joinpath(tab_dir, "group_difference_ALPHA.csv"))
    (rc === nothing || rm === nothing || rd === nothing) &&
        error("Matrices ALPHA transversales ausentes")
    ch, Wc, Wm = rc[1], rc[2], rm[2]
    Wd = rd[2]
    n = length(ch)
    stats = CSV.read(joinpath(tab_dir, "group_statistics_ALPHA.csv"), DataFrame)
    n_edges = nrow(stats)
    sig = filter(r -> Float64(r.q_value) < FDR_ALPHA, stats)
    n_sig = nrow(sig)
    lim_ab = max(maximum(Wc), maximum(Wm), 1e-12)
    lim_d = max(maximum(abs, Wd), 1e-12)
    Wd_tri = copy(Wd)
    for i in 1:n, j in 1:n; i <= j && (Wd_tri[i,j] = NaN); end

    apply_theme!()
    fig = Figure(size=fig_size_mm(180,70); figure_padding=6)
    titles = ("A  Control", "B  EM", "C  EM − Control")
    mats = (Wc, Wm, Wd_tri)
    for (col, (mat,ttl,div)) in enumerate(zip(mats, titles, (false,false,true)))
        ax = Axis(fig[1,col]; title=ttl, aspect=1, xticks=(1:n,ch), yticks=(1:n,ch),
                  xticklabelrotation=π/2, xticklabelsize=5, yticklabelsize=5)
        if div
            hm = heatmap!(ax, mat; colormap=DIVERGING_CMAP, colorrange=(-lim_d,lim_d), nan_color=:white)
            Colorbar(fig[2,col], hm; vertical=false, label="Δ wPLI", flipaxis=false, height=8, labelsize=7, ticklabelsize=6)
            for r in eachrow(sig)
                ia = findfirst(==(string(r.ch_a)), ch); ib = findfirst(==(string(r.ch_b)), ch)
                (ia === nothing || ib === nothing) && continue
                i, j = min(ia,ib), max(ia,ib)
                scatter!(ax, [Float64(i)], [Float64(j)]; marker=:rect, markersize=4, color=:transparent,
                         strokecolor=:black, strokewidth=0.9)
            end
        else
            hm = heatmap!(ax, mat; colormap=:viridis, colorrange=(0.0,lim_ab))
            Colorbar(fig[2,col], hm; vertical=false, label="wPLI", flipaxis=false, height=8, labelsize=7, ticklabelsize=6)
        end
    end
    Label(fig[3,1:3], "Montaje común ($n canales). A/B misma escala. C: triángulo inferior; contorno = aristas FDR. " *
          "FDR-BH entre $n_edges aristas dentro de ALPHA (n_sig=$n_sig).", fontsize=7, color=:gray35, tellwidth=false)
    path = save_png(fig, joinpath(figs_dir, "transversal_alpha_connectivity_matrices"))
    push!(manifest, (figure_name="transversal_alpha_connectivity_matrices", analysis="transversal", band="ALPHA",
          condition="EC", output_file=path, source_tables="group_connectivity_*_ALPHA.csv;group_difference_ALPHA.csv;group_statistics_ALPHA.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="B (group matrices, intersected montage)",
          fdr_family="FDR-BH among $n_edges edges within ALPHA"))
    return path
end

function fig_network_power_alpha(tab_dir, figs_dir, manifest)
    rc = read_mat_csv(joinpath(tab_dir, "group_connectivity_control_ALPHA.csv"))
    rc === nothing && error("matriz control ALPHA ausente")
    ch = rc[1]
    stats = CSV.read(joinpath(tab_dir, "group_statistics_ALPHA.csv"), DataFrame)
    sig = filter(r -> Float64(r.q_value) < FDR_ALPHA, stats)
    n_edges = nrow(stats)
    bp = CSV.read(joinpath(tab_dir, "band_power_group_statistics.csv"), DataFrame)
    bpA = filter(r -> string(r.band) == "ALPHA", bp)
    ndiff_path = joinpath(tab_dir, "network_metrics_diff_ALPHA.csv")
    has_ns = isfile(ndiff_path)
    ndiff = has_ns ? CSV.read(ndiff_path, DataFrame) : DataFrame()

    apply_theme!()
    fig = Figure(size=fig_size_mm(180,75); figure_padding=6)
    θs = range(0, 2π; length=120)

    axA = Axis(fig[1,1]; title="A  Red FDR ALPHA", aspect=DataAspect())
    hidedecorations!(axA); hidespines!(axA)
    lines!(axA, 1.05.*cos.(θs), 1.05.*sin.(θs); color=:gray60, linewidth=0.8)
    xy = Dict(c => get(CH_POS, uppercase(c), (0.0,0.0)) for c in ch)
    if nrow(sig) > 0
        mags = abs.(Float64.(sig.effect_rrb))
        mmax = max(maximum(mags), 1e-12)
        for (k,r) in enumerate(eachrow(sig))
            a, b = string(r.ch_a), string(r.ch_b)
            (haskey(xy,a) && haskey(xy,b)) || continue
            dlt = Float64(r.diff); col = dlt>=0 ? COLOR_POS : COLOR_NEG
            lw = 1.0 + 3.0*mags[k]/mmax
            lines!(axA, [xy[a][1],xy[b][1]], [xy[a][2],xy[b][2]]; color=col, linewidth=lw)
        end
    end
    xs=[xy[c][1] for c in ch]; ys=[xy[c][2] for c in ch]
    scatter!(axA, xs, ys; color=:white, strokecolor=:gray20, strokewidth=1, markersize=10)
    for c in ch; text!(axA, xy[c][1], xy[c][2]+0.07, text=c; fontsize=5.5, align=(:center,:bottom)); end
    xlims!(axA,-1.25,1.25); ylims!(axA,-1.25,1.25)
    lab_w = hasproperty(sig,:effect_rrb) ? "grosor ∝ |rᵣᵦ|" : "grosor ∝ |d|"
    Label(fig[2,1], "Solo aristas q<0.05. Rojo Δ>0, azul Δ<0; $lab_w. Nodos uniformes.", fontsize=6.5, color=:gray35, tellwidth=false)

    axB = Axis(fig[1,2]; title="B  Δ potencia ALPHA (EM−Ctrl)", aspect=DataAspect())
    hidedecorations!(axB); hidespines!(axB)
    lines!(axB, 1.05.*cos.(θs), 1.05.*sin.(θs); color=:gray60, linewidth=0.8)
    chs = String.(bpA.channel); dvs = Float64.(bpA.diff)
    lim = max(maximum(abs,dvs), 1e-12)
    for (c,dv,qv) in zip(chs, dvs, Float64.(bpA.q_value))
        pos = get(CH_POS, uppercase(c), (0.0,0.0))
        col = dv >= 0 ? (COLOR_POS, 0.35+0.65*abs(dv)/lim) : (COLOR_NEG, 0.35+0.65*abs(dv)/lim)
        msz = Float64(qv) < FDR_ALPHA ? 16.0 : 10.0
        scatter!(axB, [pos[1]], [pos[2]]; color=col, markersize=msz,
                 strokecolor=Float64(qv)<FDR_ALPHA ? :black : :gray50, strokewidth=Float64(qv)<FDR_ALPHA ? 1.2 : 0.4)
        text!(axB, pos[1], pos[2]+0.07, text=c; fontsize=5, align=(:center,:bottom))
    end
    xlims!(axB,-1.25,1.25); ylims!(axB,-1.25,1.25)
    n_fdr_p = count(<(FDR_ALPHA), Float64.(bpA.q_value))
    n_ms = Int(first(bpA.n_ms)); n_ct = Int(first(bpA.n_ctrl))
    Label(fig[2,2], "Potencia de banda (unidades del CSV). Contorno negro: canal FDR. " *
          "FDR-BH entre $(nrow(bpA)) canales en ALPHA (n_FDR=$n_fdr_p). n_EM=$n_ms, n_Ctrl=$n_ct.", fontsize=6.5, color=:gray35, tellwidth=false)

    axC = Axis(fig[1,3]; title="C  Δ mean_strength nodal (descr.)", aspect=DataAspect())
    hidedecorations!(axC); hidespines!(axC)
    lines!(axC, 1.05.*cos.(θs), 1.05.*sin.(θs); color=:gray60, linewidth=0.8)
    if has_ns && nrow(ndiff) > 0 && hasproperty(ndiff, :delta_strength)
        dvs2 = Float64.(ndiff.delta_strength); lim2 = max(maximum(abs,dvs2), 1e-12)
        for (c,dv) in zip(String.(ndiff.channel), dvs2)
            pos = get(CH_POS, uppercase(c), (0.0,0.0))
            col = dv >= 0 ? (COLOR_POS, 0.4+0.6*abs(dv)/lim2) : (COLOR_NEG, 0.4+0.6*abs(dv)/lim2)
            scatter!(axC, [pos[1]], [pos[2]]; color=col, markersize=12, strokecolor=:gray40, strokewidth=0.5)
            text!(axC, pos[1], pos[2]+0.07, text=c; fontsize=5, align=(:center,:bottom))
        end
        Label(fig[2,3], "Descriptivo (matrices grupales); no es contraste nodal por participante. No interpreta hubs inferenciales.",
              fontsize=6.5, color=:gray35, tellwidth=false)
    else
        text!(axC, 0.0, 0.0; text="Sin network_metrics_diff", align=(:center,:center), fontsize=8)
        Label(fig[2,3], "Panel omitible si no aporta.", fontsize=6.5, color=:gray35, tellwidth=false)
    end
    xlims!(axC,-1.25,1.25); ylims!(axC,-1.25,1.25)
    Label(fig[3,1:3], "Potencia y conectividad proceden de los mismos registros; no son validaciones independientes. " *
          "FDR aristas: $n_edges pruebas dentro de ALPHA.", fontsize=7, color=:gray30, tellwidth=false)
    path = save_png(fig, joinpath(figs_dir, "transversal_alpha_network_and_power"))
    push!(manifest, (figure_name="transversal_alpha_network_and_power", analysis="transversal", band="ALPHA",
          condition="EC", output_file=path, source_tables="group_statistics_ALPHA.csv;band_power_group_statistics.csv;network_metrics_diff_ALPHA.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="B edges FDR + spectral power (MS−Control)",
          fdr_family="edges: FDR-BH within ALPHA; power: FDR-BH among channels within ALPHA"))
    return path
end

# ═══════════════════════════════════════════════════════════════
#  Figuras exploratorias nuevas — por condición (EC y EO por igual)
# ═══════════════════════════════════════════════════════════════

"""Raincloud de wPLI medio por sujeto, generalizado a cualquier banda
(la versión ALPHA-EC de manuscrito, `fig_raincloud_alpha`, se conserva
intacta con su propio nombre de archivo por compatibilidad con LaTeX)."""
function fig_raincloud_band(tab_dir, figs_dir, manifest, band::AbstractString, cl::AbstractString)
    gdf = CSV.read(joinpath(tab_dir, "global_mean_wpli_statistics.csv"), DataFrame)
    sm  = CSV.read(joinpath(tab_dir, "subject_band_means.csv"), DataFrame)
    row = _band_row(gdf, band)
    ms, ct = _subject_band_vectors(sm, band)
    (length(ms) < 2 || length(ct) < 2) && return nothing
    d_tab = Float64(row.effect_d_pooled)
    dlo = Float64(row.effect_d_pooled_ci_low)
    dhi = Float64(row.effect_d_pooled_ci_high)
    rrb = Float64(row.effect_rrb)
    psup = Float64(row.probability_superiority)

    apply_theme!()
    fig = Figure(size=fig_size_mm(95,110); figure_padding=10)
    ax = Axis(fig[1,1]; xlabel="Grupo", ylabel="wPLI medio por participante",
              title="$band ($cl) · distribución individual",
              xticks=([1.0,2.0], ["Control\n(n=$(length(ct)))", "EM\n(n=$(length(ms)))"]))
    function _density_poly(vals)
        n = length(vals); n < 2 && return Float64[], Float64[]
        lo, hi = extrema(vals); pad = 0.05*(hi-lo+eps())
        ys = collect(range(lo-pad, hi+pad; length=60))
        bw = 1.06*std(vals)*n^(-1/5); bw < 1e-8 && (bw = 1e-3)
        dens = [mean(exp(-0.5*((y-v)/bw)^2) for v in vals)/(bw*sqrt(2π)) for y in ys]
        dens ./= (maximum(dens)+eps())
        return ys, dens .* 0.38
    end
    ys_c, dens_c = _density_poly(ct); ys_m, dens_m = _density_poly(ms)
    if !isempty(ys_c)
        poly!(ax, Point2f.(1.0 .- dens_c, ys_c); color=(COLOR_CTRL,0.35), strokewidth=0)
        lines!(ax, 1.0 .- dens_c, ys_c; color=COLOR_CTRL, linewidth=1)
    end
    if !isempty(ys_m)
        poly!(ax, Point2f.(2.0 .+ dens_m, ys_m); color=(COLOR_MS,0.35), strokewidth=0)
        lines!(ax, 2.0 .+ dens_m, ys_m; color=COLOR_MS, linewidth=1)
    end
    rng = MersenneTwister(BOOTSTRAP_SEED)
    for (x0, vals, col) in ((1.0,ct,COLOR_CTRL), (2.0,ms,COLOR_MS))
        jit = 0.06 .* (rand(rng, length(vals)) .- 0.5)
        scatter!(ax, fill(x0,length(vals)).+jit, vals; color=(col,0.75), markersize=5, strokewidth=0.3, strokecolor=:black)
        is_ctrl = x0 == 1.0
        med = Float64(is_ctrl ? row.ctrl_median : row.ms_median)
        q1 = Float64(is_ctrl ? row.ctrl_q1 : row.ms_q1)
        q3 = Float64(is_ctrl ? row.ctrl_q3 : row.ms_q3)
        lines!(ax, [x0-0.12,x0+0.12], [med,med]; color=:black, linewidth=2)
        lines!(ax, [x0,x0], [q1,q3]; color=:black, linewidth=1.4)
        μ = Float64(is_ctrl ? row.ctrl_mean : row.ms_mean)
        mlo = Float64(is_ctrl ? row.ctrl_ci_low : row.ms_ci_low)
        mhi = Float64(is_ctrl ? row.ctrl_ci_high : row.ms_ci_high)
        lines!(ax, [x0+0.18,x0+0.18], [mlo,mhi]; color=:gray40, linewidth=1)
        scatter!(ax, [x0+0.18], [μ]; color=:gray40, markersize=5, marker=:rect)
    end
    ann = "p=$(fmt_p(Float64(row.p_mannwhitney)))  q_global=$(fmt_q(Float64(row.q_value)))\n" *
          "d=$(fmt_d(d_tab)) [IC95% $(fmt_d(dlo)), $(fmt_d(dhi))]  r_rb=$(fmt_d(rrb))  P(EM>Ctrl)=$(fmt_d(psup))"
    Label(fig[2,1], ann; fontsize=7, color=:gray30, tellwidth=false)
    path = save_png(fig, joinpath(figs_dir, "raincloud_$(band)"))
    push!(manifest, (figure_name="raincloud_$(band)", analysis="transversal", band=band,
          condition=cl, output_file=path, source_tables="subject_band_means.csv;global_mean_wpli_statistics.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="A (mean wPLI, native montage)",
          fdr_family="FDR-BH across 7 bands (global)"))
    return path
end

"""Ordena canales anatómicamente (anterior→posterior, izq→dcha) usando
`CH_POS`; canales sin posición conocida se añaden al final, alfabéticos."""
function anatomical_order(chans::Vector{String})::Vector{String}
    known   = filter(c -> haskey(CH_POS, uppercase(c)), chans)
    unknown = sort(filter(c -> !haskey(CH_POS, uppercase(c)), chans))
    # Ascendente en y: occipital primero (índice bajo → abajo en el heatmap,
    # eje y creciente hacia arriba) y frontal último (índice alto → arriba),
    # para que la fila superior sea frontal, como en la lectura habitual F→C→P→O.
    sort!(known; by = c -> (CH_POS[uppercase(c)][2], CH_POS[uppercase(c)][1]))
    return vcat(known, unknown)
end

"""Heatmap canal×banda de Cohen d (potencia, EM−Control) — resume los
14 topomapas Δ-potencia por banda en una única figura, con contorno
en las celdas que sobreviven FDR-BH (calculado por banda, ya en
`band_power_group_statistics.csv`). Canales en orden anatómico. Las
celdas con N reducido (<70% del n máximo de esta condición) se marcan
aparte — un N bajo puede producir significación FDR espuria."""
function fig_power_effect_heatmap(bp_df::DataFrame, figs_dir, manifest, cl::AbstractString)
    isempty(bp_df) && return nothing
    chans = anatomical_order(unique(String.(bp_df.channel)))
    bnds  = sort_bands_physio(unique(string.(bp_df.band)))
    nC, nB = length(chans), length(bnds)
    (nC < 1 || nB < 1) && return nothing
    idxC = Dict(c => i for (i,c) in enumerate(chans))
    idxB = Dict(b => i for (i,b) in enumerate(bnds))
    D = fill(NaN, nC, nB); Q = fill(1.0, nC, nB)
    Nms = fill(0, nC, nB); Nct = fill(0, nC, nB)
    for r in eachrow(bp_df)
        ci, bi = idxC[string(r.channel)], idxB[string(r.band)]
        D[ci,bi] = Float64(r.effect_d_pooled); Q[ci,bi] = Float64(r.q_value)
        Nms[ci,bi] = Int(r.n_ms); Nct[ci,bi] = Int(r.n_ctrl)
    end
    n_ms_max = maximum(Nms); n_ct_max = maximum(Nct)
    low_thr = 0.7
    is_low(ci,bi) = Nms[ci,bi] < low_thr*n_ms_max || Nct[ci,bi] < low_thr*n_ct_max
    n_sig = count(!isnan(D[i,j]) && Q[i,j] < FDR_ALPHA for i in 1:nC, j in 1:nB)
    n_low = count(!isnan(D[i,j]) && is_low(i,j) for i in 1:nC, j in 1:nB)

    apply_theme!()
    lim = maximum(x -> isnan(x) ? 0.0 : abs(x), D); lim < 1e-9 && (lim = 1.0)
    fig = Figure(size=fig_size_mm(150, 28 + 4.3*nC); figure_padding=8)
    ax = Axis(fig[1,1]; title="Potencia — tamaño de efecto por canal×banda (Cohen d, EM−Control) — $cl",
              xlabel="Banda", ylabel="Canal", xticks=(1:nB, bnds), yticks=(1:nC, chans),
              xticklabelrotation=π/4, xticklabelsize=8, yticklabelsize=7)
    Dt = permutedims(D)
    hm = heatmap!(ax, Dt; colormap=DIVERGING_CMAP, colorrange=(-lim,lim), nan_color=:gray90)
    Colorbar(fig[1,2], hm; label="Cohen d")
    for ci in 1:nC, bi in 1:nB
        isnan(D[ci,bi]) && continue
        sig = Q[ci,bi] < FDR_ALPHA
        low = is_low(ci,bi)
        if sig && low
            scatter!(ax, [Float64(bi)], [Float64(ci)]; marker=:cross, markersize=6, color=:orange, strokewidth=0)
        elseif sig
            scatter!(ax, [Float64(bi)], [Float64(ci)]; marker=:cross, markersize=6, color=:black, strokewidth=0)
        elseif low
            scatter!(ax, [Float64(bi)], [Float64(ci)]; marker=:circle, markersize=5, color=:transparent,
                     strokecolor=:orange, strokewidth=0.8)
        end
    end
    Label(fig[2,1:2], "Cruz negra: q<0.05 (FDR-BH por banda, n=$n_sig celdas).\n" *
          "Naranja (cruz o círculo): N reducido (<70% del máximo de esta condición: " *
          "n_EM=$n_ms_max, n_Ctrl=$n_ct_max; n=$n_low celdas) — interpretar con cautela.\n" *
          "N exacto por celda en band_power_group_statistics.csv.",
          fontsize=7, color=:gray35, tellwidth=false, justification=:center)
    path = save_png(fig, joinpath(figs_dir, "power_effect_heatmap_channel_band"))
    push!(manifest, (figure_name="power_effect_heatmap_channel_band", analysis="transversal", band="ALL",
          condition=cl, output_file=path, source_tables="band_power_group_statistics.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="Band power effect size (MS−Control)",
          fdr_family="FDR-BH among channels within each band"))
    return path
end

# ═══════════════════════════════════════════════════════════════
#  Interacción grupo × condición (EC×EO) — cruza eyesclosed/eyesopen,
#  se genera UNA vez, después del bucle principal, leyendo los CSV ya
#  escritos por ambas condiciones (mismo principio CSV→figura de todo
#  el módulo). Sale a results/transversal/combined/{tables,figures}/,
#  al no pertenecer a ninguna de las dos carpetas de condición.
# ═══════════════════════════════════════════════════════════════

function _fig_interaction_forest(idf::DataFrame, figs_dir, manifest, n_ms_s::Int, n_ct_s::Int)
    bands = sort_bands_physio(idf.band); nB = length(bands)
    y_of = Dict(b => Float64(nB-i+1) for (i,b) in enumerate(bands))
    apply_theme!()
    fig = Figure(size=fig_size_mm(170,95); figure_padding=8)
    ax = Axis(fig[1,1]; xlabel="Diff-of-diff (wPLI): (EM: EC−EO) − (Control: EC−EO)",
              yticks=(Float64.(nB:-1:1), bands),
              title="Interacción grupo×condición (EC/EO) — n=$n_ms_s EM / $n_ct_s Control",
              limits=(nothing,nothing,0.4,nB+0.6))
    vlines!(ax, [0.0]; color=COLOR_ZERO, linewidth=0.9, linestyle=:dash)
    for row in eachrow(idf)
        yi = y_of[string(row.band)]
        sig = Float64(row.q_value) < FDR_ALPHA
        col = sig ? (Float64(row.diff_of_diff)>=0 ? COLOR_POS : COLOR_NEG) : COLOR_NS
        lines!(ax, [Float64(row.ci_lo),Float64(row.ci_hi)], [yi,yi]; color=col, linewidth=sig ? 2.2 : 1.4)
        scatter!(ax, [Float64(row.diff_of_diff)], [yi]; color=col, markersize=sig ? 9 : 7,
                 marker=sig ? :diamond : :circle, strokecolor=:black, strokewidth=0.4)
    end
    ax2 = Axis(fig[1,2]; limits=(0,1,0.4,nB+0.6), xgridvisible=false, ygridvisible=false)
    hidedecorations!(ax2); hidespines!(ax2)
    for row in eachrow(idf)
        yi = y_of[string(row.band)]
        txt = "n=$(row.n_ms)/$(row.n_ctrl)  p=$(fmt_p(Float64(row.p_mannwhitney)))  q=$(fmt_q(Float64(row.q_value)))"
        text!(ax2, 0.02, yi; text=txt, align=(:left,:center), fontsize=7, color=Float64(row.q_value)<FDR_ALPHA ? :black : :gray40)
    end
    colsize!(fig.layout,1,Relative(0.55)); colsize!(fig.layout,2,Relative(0.45))
    Label(fig[2,1:2], "Diferencia de diferencias EC−EO entre grupos (mismo estimando de interacción grupo×condición;\n" *
          "Mann–Whitney sobre los cambios EC−EO por sujeto, no formalmente equivalente a un ANOVA 2×2).\n" *
          "Punto e IC en wPLI bruto (misma escala). FDR-BH entre las 7 bandas. IC95% bootstrap (B=$BOOTSTRAP_N, semilla=$BOOTSTRAP_SEED).",
          fontsize=7, color=:gray35, tellwidth=false, justification=:center)
    path = save_png(fig, joinpath(figs_dir, "interaction_ec_eo_forest"))
    push!(manifest, (figure_name="interaction_ec_eo_forest", analysis="transversal", band="ALL",
          condition="EC_EO", output_file=path, source_tables="interaction_ec_eo_by_band.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="Group×condition interaction (EC−EO diff-of-diff)",
          fdr_family="FDR-BH across 7 bands"))
    return path
end

function _fig_interaction_slopeplots(pdf::DataFrame, figs_dir, manifest, n_ms_s::Int, n_ct_s::Int)
    bands = sort_bands_physio(unique(string.(pdf.band))); nB = length(bands)
    ncols = 4; nrows = cld(nB, ncols)
    apply_theme!()
    fig = Figure(size=fig_size_mm(220, 55.0*nrows + 12); figure_padding=8)
    for (k,b) in enumerate(bands)
        r = div(k-1, ncols) + 1; c = mod(k-1, ncols) + 1
        sub = filter(row -> string(row.band)==b, pdf)
        ax = Axis(fig[r,c]; title=b, titlesize=9, xticks=([1.0,2.0], ["EC","EO"]), xticklabelsize=8,
                  ylabel = c==1 ? "wPLI medio" : "", yticklabelsize=7)
        for row in eachrow(sub)
            col = row.group == "MS" ? (COLOR_MS, 0.35) : (COLOR_CTRL, 0.35)
            lines!(ax, [1.0,2.0], [Float64(row.wpli_ec), Float64(row.wpli_eo)]; color=col, linewidth=1.0)
        end
        for (grp, col) in (("MS",COLOR_MS), ("Control",COLOR_CTRL))
            g = filter(row -> row.group==grp, sub)
            isempty(g) && continue
            mec = mean(Float64.(g.wpli_ec)); meo = mean(Float64.(g.wpli_eo))
            lines!(ax, [1.0,2.0], [mec,meo]; color=col, linewidth=3)
            scatter!(ax, [1.0,2.0], [mec,meo]; color=col, markersize=9, strokecolor=:black, strokewidth=0.5)
        end
    end
    Label(fig[nrows+1, 1:ncols], "Línea fina: un sujeto (rojo=EM, azul=Control). Línea gruesa: media de grupo. " *
          "n=$n_ms_s EM / $n_ct_s Control con datos válidos en EC y EO.", fontsize=7, color=:gray35, tellwidth=false)
    path = save_png(fig, joinpath(figs_dir, "interaction_ec_eo_slopeplots"))
    push!(manifest, (figure_name="interaction_ec_eo_slopeplots", analysis="transversal", band="ALL",
          condition="EC_EO", output_file=path, source_tables="subject_ec_eo_paired.csv",
          generated_at=Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"), estimand="Within-subject EC→EO trajectories by group",
          fdr_family="n/a (descriptive)"))
    return path
end

"""Sensibilidad de montaje común: repite la interacción EC×EO restringiendo
cada sujeto a los canales que tiene en común entre EC y EO (montaje
*individual*, no el agregado por banda que usa el resto del módulo) —
descarta que la interacción sea un artefacto de comparar montajes distintos
entre condiciones. Requiere recargar wPLI crudo (subject_band_means.csv no
guarda qué canales usó cada sujeto)."""
function fig_interaction_common_channels(res_root, bands, tab_dir)
    ec_incl_path = joinpath(res_root, "transversal", "eyesclosed", "tables", "subject_inclusion.csv")
    eo_incl_path = joinpath(res_root, "transversal", "eyesopen",   "tables", "subject_inclusion.csv")
    if !(isfile(ec_incl_path) && isfile(eo_incl_path))
        @warn "Sensibilidad montaje común omitida: falta subject_inclusion.csv"
        return nothing
    end
    ec_i = filter(r -> string(r.included)=="true", CSV.read(ec_incl_path, DataFrame))
    eo_i = filter(r -> string(r.included)=="true", CSV.read(eo_incl_path, DataFrame))
    eo_by_sid = Dict(string(r.subject_id) => (string(r.session_id), string(r.group)) for r in eachrow(eo_i))

    rows = NamedTuple[]
    n_common_all = Int[]
    for r in eachrow(ec_i)
        sid = string(r.subject_id)
        haskey(eo_by_sid, sid) || continue
        sess_eo, grp = eo_by_sid[sid]
        sess_ec = string(r.session_id)
        for band in bands
            rec = load_wpli(res_root, sid, sess_ec, "eyesclosed", band)
            reo = load_wpli(res_root, sid, sess_eo, "eyesopen", band)
            (rec === nothing || reo === nothing) && continue
            ch_ec, W_ec = rec; ch_eo, W_eo = reo
            common = intersect(ch_ec, ch_eo)
            length(common) < 2 && continue
            m_ec = upper_mean(realign_matrix(ch_ec, W_ec, common))
            m_eo = upper_mean(realign_matrix(ch_eo, W_eo, common))
            (isnan(m_ec) || isnan(m_eo)) && continue
            push!(n_common_all, length(common))
            push!(rows, (schema_version=TRANSVERSAL_SCHEMA_VERSION,
                  statistics_source=STATISTICS_SOURCE,
                  subject_id=sid, group=grp, band=band, n_common_channels=length(common),
                  wpli_ec_common=m_ec, wpli_eo_common=m_eo,
                  delta_ec_eo_common=m_ec-m_eo))
        end
    end
    if isempty(rows)
        @warn "Sensibilidad montaje común omitida: sin datos válidos"
        return nothing
    end
    cdf = DataFrame(rows)
    CSV.write(joinpath(tab_dir, "subject_ec_eo_common_channels.csv"), cdf)

    inter_rows = NamedTuple[]
    for band in bands
        sub = filter(row -> string(row.band)==band, cdf)
        ms_d = Float64[r.delta_ec_eo_common for r in eachrow(sub) if r.group=="MS"]
        ct_d = Float64[r.delta_ec_eo_common for r in eachrow(sub) if r.group=="Control"]
        (length(ms_d) < 2 || length(ct_d) < 2) && continue
        pmw, rrb = mannwhitney_p(ms_d, ct_d)
        _, d0 = welch_t(ms_d, ct_d)
        _, lo, hi = bootstrap_mean_diff_ci(ms_d, ct_d)
        push!(inter_rows, (schema_version=TRANSVERSAL_SCHEMA_VERSION,
              statistics_source=STATISTICS_SOURCE, fdr_scope=GLOBAL_FDR_SCOPE,
              band=band, n_ms=length(ms_d), n_ctrl=length(ct_d),
              delta_ec_eo_ms=mean(ms_d), delta_ec_eo_ctrl=mean(ct_d),
              diff_of_diff=mean(ms_d)-mean(ct_d),
              p_mannwhitney=pmw, effect_d_pooled=d0, effect_rrb=rrb,
              ci_lo=lo, ci_hi=hi))
    end
    if !isempty(inter_rows)
        idf = DataFrame(inter_rows)
        idf[!, :q_value] = bh_qvalues(Float64.(idf.p_mannwhitney))
        idf[!, :fdr_family_size] = fill(nrow(idf), nrow(idf))
        idf[!, :is_fdr] = Float64.(idf.q_value) .< FDR_ALPHA
        CSV.write(joinpath(tab_dir, "interaction_ec_eo_by_band_common_channels.csv"), idf)
    end
    med_n = round(Int, median(n_common_all))
    println("  Sensibilidad montaje común: $(length(unique(cdf.subject_id))) sujetos, " *
            "canales por sujeto mediana=$med_n, rango=$(minimum(n_common_all))–$(maximum(n_common_all))")
    return nothing
end

"""Interacción grupo×condición: lee subject_band_means.csv de eyesclosed
y eyesopen (ya escritos por el bucle principal), empareja por sujeto+banda
y compara Δ(EC−EO) entre EM y Control. No requiere sesión T2 — usa los
mismos sujetos T1 con datos válidos en ambas condiciones oculares."""
function fig_interaction_ec_eo(res_root, bands)
    ec_path = joinpath(res_root, "transversal", "eyesclosed", "tables", "subject_band_means.csv")
    eo_path = joinpath(res_root, "transversal", "eyesopen",   "tables", "subject_band_means.csv")
    if !(isfile(ec_path) && isfile(eo_path))
        @warn "Interacción EC×EO omitida: falta subject_band_means.csv en eyesclosed y/o eyesopen"
        return nothing
    end
    ec = CSV.read(ec_path, DataFrame); eo = CSV.read(eo_path, DataFrame)
    out_dir  = joinpath(res_root, "transversal", "combined")
    tab_dir  = joinpath(out_dir, "tables"); figs_dir = joinpath(out_dir, "figures")
    mkpath(tab_dir); mkpath(figs_dir)

    paired = NamedTuple[]
    for band in bands
        ec_b = filter(r -> string(r.band)==band, ec)
        eo_b = filter(r -> string(r.band)==band, eo)
        eo_by_sid = Dict(string(r.subject_id) => r for r in eachrow(eo_b))
        for r in eachrow(ec_b)
            sid = string(r.subject_id)
            haskey(eo_by_sid, sid) || continue
            r2 = eo_by_sid[sid]
            push!(paired, (schema_version=TRANSVERSAL_SCHEMA_VERSION,
                  statistics_source=STATISTICS_SOURCE,
                  subject_id=sid, group=string(r.group), band=band,
                  wpli_ec=Float64(r.mean_wpli), wpli_eo=Float64(r2.mean_wpli),
                  delta_ec_eo=Float64(r.mean_wpli)-Float64(r2.mean_wpli)))
        end
    end
    if isempty(paired)
        @warn "Interacción EC×EO omitida: ningún sujeto con datos válidos en ambas condiciones"
        return nothing
    end
    pdf = DataFrame(paired)
    CSV.write(joinpath(tab_dir, "subject_ec_eo_paired.csv"), pdf)

    inter_rows = NamedTuple[]
    for band in bands
        sub = filter(row -> string(row.band)==band, pdf)
        ms_d = Float64[r.delta_ec_eo for r in eachrow(sub) if r.group=="MS"]
        ct_d = Float64[r.delta_ec_eo for r in eachrow(sub) if r.group=="Control"]
        (length(ms_d) < 2 || length(ct_d) < 2) && continue
        pmw, rrb = mannwhitney_p(ms_d, ct_d)
        _, d0 = welch_t(ms_d, ct_d)
        # IC bruto (misma unidad que diff_of_diff, wPLI) — NO el IC de Cohen d de
        # bootstrap_cohen_d_ci, que está en escala estandarizada y no es comparable
        # con el punto graficado (bug detectado en revisión, corregido aquí).
        _, lo, hi = bootstrap_mean_diff_ci(ms_d, ct_d)
        push!(inter_rows, (schema_version=TRANSVERSAL_SCHEMA_VERSION,
              statistics_source=STATISTICS_SOURCE, fdr_scope=GLOBAL_FDR_SCOPE,
              band=band, n_ms=length(ms_d), n_ctrl=length(ct_d),
              delta_ec_eo_ms=mean(ms_d), delta_ec_eo_ctrl=mean(ct_d),
              diff_of_diff=mean(ms_d)-mean(ct_d),
              p_mannwhitney=pmw, effect_d_pooled=d0, effect_rrb=rrb,
              ci_lo=lo, ci_hi=hi))
    end
    if isempty(inter_rows)
        @warn "Interacción EC×EO omitida: ninguna banda con ≥2 sujetos por grupo"
        return nothing
    end
    idf = DataFrame(inter_rows)
    idf[!, :q_value] = bh_qvalues(Float64.(idf.p_mannwhitney))
    idf[!, :fdr_family_size] = fill(nrow(idf), nrow(idf))
    idf[!, :is_fdr] = Float64.(idf.q_value) .< FDR_ALPHA
    CSV.write(joinpath(tab_dir, "interaction_ec_eo_by_band.csv"), idf)

    n_ms_s = length(unique(string.(filter(r -> r.group=="MS", pdf).subject_id)))
    n_ct_s = length(unique(string.(filter(r -> r.group=="Control", pdf).subject_id)))
    println("  Interacción EC×EO: $(n_ms_s+n_ct_s) sujetos con ambas condiciones (EM=$n_ms_s, Control=$n_ct_s)")

    try
        fig_interaction_common_channels(res_root, bands, tab_dir)
    catch e
        @warn "Sensibilidad montaje común fallida: $e"
    end

    manifest = NamedTuple[]
    try
        _fig_interaction_forest(idf, figs_dir, manifest, n_ms_s, n_ct_s)
        _fig_interaction_slopeplots(pdf, figs_dir, manifest, n_ms_s, n_ct_s)
        write_manifest(joinpath(tab_dir, "figures_manifest.tsv"), manifest)
        println("  Figuras interacción ($(length(manifest))): " * join([r.figure_name for r in manifest], ", "))
    catch e
        @warn "Figuras de interacción EC×EO fallidas: $e"
    end
    println("  ✅ Guardado en: $out_dir\n")
    return nothing
end

function write_manifest(path, rows)
    isempty(rows) && return nothing
    mkpath(dirname(path))
    CSV.write(path, DataFrame(rows); delim='\t')
    return path
end

# ═══════════════════════════════════════════════════════════════
#  Cohorte T1 desde groups.csv
# ═══════════════════════════════════════════════════════════════

struct SubjT1
    subject_id::String
    session_id::String
    group::String   # "MS" | "Control"
end

function load_t1_cohort(bids_root)::Vector{SubjT1}
    groups_path = joinpath(bids_root, "groups.csv")
    if !isfile(groups_path)
        @warn "No se encontró groups.csv en: $groups_path"
        exit(1)
    end
    gdf = CSV.read(groups_path, DataFrame)
    rename!(gdf, Dict(n => Symbol(lowercase(string(n))) for n in names(gdf)))
    if !hasproperty(gdf, :session_id) && hasproperty(gdf, :session); rename!(gdf, :session => :session_id); end
    if !hasproperty(gdf, :session_id) && hasproperty(gdf, :ses); rename!(gdf, :ses => :session_id); end
    if hasproperty(gdf, :bids_id)
        hasproperty(gdf, :subject_id) && select!(gdf, Not(:subject_id))
        rename!(gdf, :bids_id => :subject_id)
    end
    missing_cols = filter(c -> !hasproperty(gdf, c), [:subject_id, :group, :session_id])
    isempty(missing_cols) || error("groups.csv: faltan columnas: $(join(missing_cols, ", "))")
    filter!(row -> is_t1_session(string(row.session_id)), gdf)
    unique!(gdf, [:subject_id, :session_id, :group])
    cohort = SubjT1[]
    for row in eachrow(gdf)
        sid, ses, grp = string(row.subject_id), string(row.session_id), string(row.group)
        if is_ms_group(grp); push!(cohort, SubjT1(sid, ses, "MS"))
        elseif is_ctrl_group(grp); push!(cohort, SubjT1(sid, ses, "Control"))
        else; @warn "Grupo desconocido '$grp' para $sid — omitido"; end
    end
    seen = Set{String}(); uniq = SubjT1[]
    for s in cohort
        s.subject_id in seen && continue
        push!(seen, s.subject_id); push!(uniq, s)
    end
    return uniq
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
    println(" NeuroMIND — Análisis transversal EM vs Control (T1)")
    println(" $(now())")
    println(" Diseño: caso-control T1 (EM=$N_MS_DESIGN vs Ctrl=$N_CTRL_DESIGN)")
    println(" eyesclosed y eyesopen en paralelo (mismo contraste, sin pooling)")
    println(" wPLI=$wpli_method  dwPLI=$use_dwpli  test=mannwhitney  density=$graph_dens")
    println("="^62)

    cohort = load_t1_cohort(bids_root)
    n_ms_cand   = count(s -> s.group == "MS", cohort)
    n_ctrl_cand = count(s -> s.group == "Control", cohort)
    println("📋 Cohorte T1: EM=$n_ms_cand | Control=$n_ctrl_cand  (candidatos; N diseño $N_MS_DESIGN/$N_CTRL_DESIGN)")

    qc_table = load_qc_table(res_root)
    qc_notes = load_qc_notes(res_root)
    println("  QC decisions cargadas: $(length(qc_table)) · notes: $(length(qc_notes))")
    println()

    for cond in ["eyesclosed", "eyesopen"]
        cl = cond_label(cond)   # "EC" | "EO" — solo para títulos/prints
        println("── Condición: $cl  (contraste transversal EM vs Control) " * "─"^12)
        out_dir  = joinpath(res_root, "transversal", cond)
        figs_dir = joinpath(out_dir, "figures")
        tab_dir  = joinpath(out_dir, "tables")
        mkpath(figs_dir); mkpath(tab_dir)
        write_statistics_contract(joinpath(out_dir, "statistics_contract.toml"), cl)
        manifest = NamedTuple[]

        open(joinpath(out_dir, "config_snapshot.toml"), "w") do io
            println(io, "# NeuroMIND transversal snapshot — $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
            println(io, "schema_version = $TRANSVERSAL_SCHEMA_VERSION")
            println(io, "statistics_source = \"$STATISTICS_SOURCE\"")
            println(io, "design = \"case_control_T1\"")
            println(io, "session_policy = \"T1_only\"")
            println(io, "condition = \"$cl\"")
            println(io, "n_ms_design = $N_MS_DESIGN")
            println(io, "n_ctrl_design = $N_CTRL_DESIGN")
            println(io, "test = \"mannwhitney\"")
            println(io, "fdr = \"bh\"")
            println(io, "fdr_scope_edges = \"$EDGE_FDR_SCOPE\"")
            println(io, "fdr_scope_global = \"$GLOBAL_FDR_SCOPE\"")
            println(io, "fdr_scope_power = \"$POWER_FDR_SCOPE\"")
            println(io, "bootstrap_method = \"$BOOTSTRAP_METHOD\"")
            println(io, "bootstrap_iterations = $BOOTSTRAP_N")
            println(io, "bootstrap_seed = $BOOTSTRAP_SEED")
            println(io, "quantile_method = \"$QUANTILE_METHOD\"")
            println(io, "rrb_method = \"$RRB_METHOD\"")
            println(io, "effect_d_pooled_method = \"$EFFECT_D_POOLED_METHOD\"")
            println(io, "mannwhitney_method = \"$MANNWHITNEY_METHOD\"")
            println(io, "wpli_method = \"$wpli_method\"")
            println(io, "use_dwpli = $use_dwpli")
            println(io, "graph_density = $graph_dens")
            println(io, "graph_threshold_method = \"$graph_meth\"")
            println(io, "n_channels_analysis = $n_ch_mont")
            println(io, "qc_allowed = [\"include\", \"include_with_warning\"]")
            println(io, "source_config = \"$config_path\"")
        end

        ms_data   = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
        ctrl_data = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
        subject_means = NamedTuple[]
        inclusion     = NamedTuple[]
        included_subj = SubjT1[]

        for s in cohort
            qc_dec = qc_decision(qc_table, s.subject_id, s.session_id, cond)
            ep = load_seg_epochs(res_root, s.subject_id, s.session_id, cond)
            bands_ok = String[]
            for band in bands
                r = load_wpli(res_root, s.subject_id, s.session_id, cond, band)
                r === nothing || push!(bands_ok, band)
            end
            reasons = String[]
            isempty(bands_ok) && push!(reasons, "Sin datos wPLI ($cl)")
            qc_ok(qc_dec) || push!(reasons, "QC=$(qc_dec)")
            included = isempty(reasons)
            warn_txt = ""
            if included && qc_dec == "include_with_warning"
                warn_txt = qc_note(qc_notes, s.subject_id, s.session_id, cond)
                isempty(warn_txt) && (warn_txt = "include_with_warning")
            end
            push!(inclusion, (schema_version=TRANSVERSAL_SCHEMA_VERSION,
                  statistics_source=STATISTICS_SOURCE,
                  subject_id=s.subject_id, session_id=s.session_id, group=s.group,
                  n_bands_ok=length(bands_ok), included=included,
                  excluded_reason=included ? "" : join(reasons, "; "), qc_decision=qc_dec,
                  warning_reason=warn_txt, n_epochs_valid=ep[1]===missing ? "" : string(ep[1])))
            included || continue
            push!(included_subj, s)
            for band in bands
                r = load_wpli(res_root, s.subject_id, s.session_id, cond, band)
                r === nothing && continue
                (ch, W) = r
                target = s.group == "MS" ? ms_data : ctrl_data
                haskey(target, band) || (target[band] = [])
                push!(target[band], (ch, W))
                μ = upper_mean(W); isnan(μ) && continue
                push!(subject_means, (schema_version=TRANSVERSAL_SCHEMA_VERSION,
                      statistics_source=STATISTICS_SOURCE,
                      subject_id=s.subject_id, group=s.group, band=band, cond=cl,
                      mean_wpli=μ))
            end
        end

        n_ms   = count(r -> r.group=="MS" && r.included, inclusion)
        n_ctrl = count(r -> r.group=="Control" && r.included, inclusion)
        n_excl = count(r -> !r.included, inclusion)
        n_excl_qc = count(r -> !r.included && occursin("QC=", r.excluded_reason), inclusion)
        n_excl_data = n_excl - n_excl_qc
        println("  Comparación transversal $cl | EM: $n_ms / $N_MS_DESIGN | Control: $n_ctrl / $N_CTRL_DESIGN | Excluidos: $n_excl")

        CSV.write(joinpath(tab_dir, "subject_inclusion.csv"), DataFrame(inclusion))
        !isempty(subject_means) && CSV.write(joinpath(tab_dir, "subject_band_means.csv"), DataFrame(subject_means))

        if n_ms < 1 || n_ctrl < 1
            println("  ⚠  Análisis omitido: se necesitan sujetos en ambos grupos para $cl")
            open(joinpath(out_dir, "transversal_summary.json"), "w") do io
                write(io, """{"schema_version":$TRANSVERSAL_SCHEMA_VERSION,"statistics_source":"$STATISTICS_SOURCE","n_ms":$n_ms,"n_ctrl":$n_ctrl,"n_ms_design":$N_MS_DESIGN,"n_ctrl_design":$N_CTRL_DESIGN,"n_included":$(n_ms+n_ctrl),"n_excluded":$n_excl,"n_excluded_qc":$n_excl_qc,"n_excluded_data":$n_excl_data,"n_total_sig":0,"n_bands":0,"cond":"$cl","design":"case_control_T1","status":"insufficient_groups","timestamp":"$(Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"))"}""")
            end
            continue
        end

        band_stats_rows = NamedTuple[]

        for band in bands
            ms_mats   = get(ms_data,   band, Tuple{Vector{String},Matrix{Float64}}[])
            ctrl_mats = get(ctrl_data, band, Tuple{Vector{String},Matrix{Float64}}[])
            (isempty(ms_mats) || isempty(ctrl_mats)) && continue
            all_ch_sets = [Set(ch) for (ch,_) in vcat(ms_mats, ctrl_mats)]
            common_set  = reduce(intersect, all_ch_sets)
            ch_ref      = ctrl_mats[1][1]
            common_ch   = filter(c -> c in common_set, ch_ref)
            n = length(common_ch); n < 2 && continue

            get_W(ch, W) = realign_matrix(ch, W, common_ch)
            Wms   = mean(get_W(ch,W) for (ch,W) in ms_mats)
            Wctrl = mean(get_W(ch,W) for (ch,W) in ctrl_mats)
            Wdiff = Wms .- Wctrl

            save_mat_csv(joinpath(tab_dir, "group_connectivity_control_$(band).csv"), Wctrl, common_ch)
            save_mat_csv(joinpath(tab_dir, "group_connectivity_ms_$(band).csv"), Wms, common_ch)
            save_mat_csv(joinpath(tab_dir, "group_difference_$(band).csv"), Wdiff, common_ch)

            upper_idx = [(i,j) for i in 1:n for j in (i+1):n]
            p_mw = ones(length(upper_idx)); p_w = ones(length(upper_idx))
            d_vec = zeros(length(upper_idx)); rrb_vec = zeros(length(upper_idx))
            for (k,(i,j)) in enumerate(upper_idx)
                va = [get_W(ch,W)[i,j] for (ch,W) in ms_mats]
                vb = [get_W(ch,W)[i,j] for (ch,W) in ctrl_mats]
                p_mw[k], rrb_vec[k] = mannwhitney_p(va, vb)
                p_w[k], d_vec[k] = welch_t(va, vb)
            end
            q_vec = bh_qvalues(p_mw)
            n_ms_b = length(ms_mats); n_ct_b = length(ctrl_mats)
            effect_order = sortperm(collect(eachindex(d_vec));
                by=k -> (-abs(d_vec[k]), common_ch[upper_idx[k][1]], common_ch[upper_idx[k][2]]))
            effect_rank = zeros(Int, length(d_vec))
            for (rank, k) in enumerate(effect_order)
                effect_rank[k] = rank
            end
            stat_rows = [(
                schema_version=TRANSVERSAL_SCHEMA_VERSION,
                statistics_source=STATISTICS_SOURCE,
                fdr_scope=EDGE_FDR_SCOPE,
                fdr_family_size=length(upper_idx),
                ch_a=common_ch[i], ch_b=common_ch[j],
                ctrl_mean=Wctrl[i,j], ms_mean=Wms[i,j], diff=Wdiff[i,j],
                p_value=p_mw[k], p_mannwhitney=p_mw[k], p_welch=p_w[k],
                q_value=q_vec[k], effect_d_pooled=d_vec[k], effect_rrb=rrb_vec[k],
                effect_rank_abs_d=effect_rank[k],
                is_nominal=p_mw[k] < FDR_ALPHA, is_fdr=q_vec[k] < FDR_ALPHA,
                n_ms=n_ms_b, n_ctrl=n_ct_b,
            ) for (k,(i,j)) in enumerate(upper_idx)]
            stats_df = DataFrame(stat_rows)
            CSV.write(joinpath(tab_dir, "group_statistics_$(band).csv"), stats_df)

            sig_df = filter(:is_fdr => identity, stats_df)
            sig_rows = [NamedTuple(r) for r in eachrow(sig_df)]
            CSV.write(joinpath(tab_dir, "significant_edges_$(band).csv"), sig_df)

            n_sig = length(sig_rows)
            n_nominal = count(<(FDR_ALPHA), p_mw)
            n_top20_nominal = count(k -> effect_rank[k] <= 20 && p_mw[k] < FDR_ALPHA,
                                    eachindex(p_mw))
            mean_ctrl = mean(Wctrl[i,j] for (i,j) in upper_idx)
            mean_ms   = mean(Wms[i,j] for (i,j) in upper_idx)
            mean_d    = mean(abs, d_vec)
            @printf("  %-10s  %3d ch  %4d pares  %3d sig (q<0.05)  ctrl=%.3f  MS=%.3f\n",
                    band, n, length(upper_idx), n_sig, mean_ctrl, mean_ms)

            push!(band_stats_rows, (
                  schema_version=TRANSVERSAL_SCHEMA_VERSION,
                  statistics_source=STATISTICS_SOURCE, fdr_scope=EDGE_FDR_SCOPE,
                  band=band, n_channels=n, n_edges=length(upper_idx),
                  n_pairs=length(upper_idx), fdr_family_size=length(upper_idx),
                  n_nominal=n_nominal, n_sig=n_sig,
                  top20_nominal_overlap=n_top20_nominal,
                  pct_sig=100.0*n_sig/max(1,length(upper_idx)),
                  ctrl_mean=mean_ctrl, ms_mean=mean_ms, diff_mean=mean_ms-mean_ctrl,
                  mean_p=mean(p_mw), mean_abs_d_pooled=mean_d,
                  n_ms=n_ms_b, n_ctrl=n_ct_b))

            try
                save_wpli_heatmaps(figs_dir, band, cl, Wctrl, Wms, Wdiff, common_ch)
                save_sig_or_explore_network(figs_dir, band, cl, common_ch, sig_rows, stats_df)
            catch e
                @warn "Figuras wPLI $band fallidas: $e"
            end

            try
                s_ms, d_ms, ns_ms, _ = nodal_strength_degree(Wms; density=graph_dens)
                s_ct, d_ct, ns_ct, _ = nodal_strength_degree(Wctrl; density=graph_dens)
                CSV.write(joinpath(tab_dir, "network_metrics_ms_$(band).csv"),
                    DataFrame(schema_version=fill(TRANSVERSAL_SCHEMA_VERSION, n),
                              channel=common_ch, strength=s_ms, degree=d_ms, norm_strength=ns_ms))
                CSV.write(joinpath(tab_dir, "network_metrics_control_$(band).csv"),
                    DataFrame(schema_version=fill(TRANSVERSAL_SCHEMA_VERSION, n),
                              channel=common_ch, strength=s_ct, degree=d_ct, norm_strength=ns_ct))
                CSV.write(joinpath(tab_dir, "network_metrics_diff_$(band).csv"),
                    DataFrame(schema_version=fill(TRANSVERSAL_SCHEMA_VERSION, n),
                              channel=common_ch, delta_strength=s_ms.-s_ct,
                              delta_degree=d_ms.-d_ct, delta_norm_strength=ns_ms.-ns_ct))
            catch e
                @warn "Network $band fallido: $e"
            end
        end

        !isempty(band_stats_rows) && CSV.write(joinpath(tab_dir, "band_statistics.csv"), DataFrame(band_stats_rows))

        global_rows = NamedTuple[]
        for band in bands
            ms_μ = [r.mean_wpli for r in subject_means if r.band==band && r.group=="MS"]
            ct_μ = [r.mean_wpli for r in subject_means if r.band==band && r.group=="Control"]
            (length(ms_μ) < 2 || length(ct_μ) < 2) && continue
            band_idx = something(findfirst(==(band), bands), 0)
            cond_offset = cl == "EC" ? 0 : 10_000
            seed = BOOTSTRAP_SEED + cond_offset + band_idx
            s = independent_group_summary(Float64.(ms_μ), Float64.(ct_μ); seed=seed)
            push!(global_rows, (
                schema_version=TRANSVERSAL_SCHEMA_VERSION,
                statistics_source=STATISTICS_SOURCE,
                fdr_scope=GLOBAL_FDR_SCOPE,
                band=band,
                ms_mean=s.ms_mean, ms_ci_low=s.ms_ci_low, ms_ci_high=s.ms_ci_high,
                ms_sem=s.ms_sem, ms_median=s.ms_median, ms_q1=s.ms_q1, ms_q3=s.ms_q3,
                ctrl_mean=s.ctrl_mean, ctrl_ci_low=s.ctrl_ci_low, ctrl_ci_high=s.ctrl_ci_high,
                ctrl_sem=s.ctrl_sem, ctrl_median=s.ctrl_median,
                ctrl_q1=s.ctrl_q1, ctrl_q3=s.ctrl_q3,
                diff=s.diff, diff_ci_low=s.diff_ci_low, diff_ci_high=s.diff_ci_high,
                p_mannwhitney=s.p_mannwhitney, p_welch=s.p_welch,
                effect_d_pooled=s.effect_d_pooled,
                effect_d_pooled_ci_low=s.effect_d_pooled_ci_low,
                effect_d_pooled_ci_high=s.effect_d_pooled_ci_high,
                effect_rrb=s.effect_rrb,
                probability_superiority=s.probability_superiority,
                n_ms=s.n_ms, n_ctrl=s.n_ctrl, bootstrap_seed=seed,
            ))
        end
        if !isempty(global_rows)
            gdf = DataFrame(global_rows)
            gdf[!, :q_value] = bh_qvalues(Float64.(gdf.p_mannwhitney))
            gdf[!, :fdr_family_size] = fill(nrow(gdf), nrow(gdf))
            gdf[!, :bootstrap_method] = fill(BOOTSTRAP_METHOD, nrow(gdf))
            gdf[!, :bootstrap_iterations] = fill(BOOTSTRAP_N, nrow(gdf))
            gdf[!, :quantile_method] = fill(QUANTILE_METHOD, nrow(gdf))
            gdf[!, :rrb_method] = fill(RRB_METHOD, nrow(gdf))
            gdf[!, :effect_d_pooled_method] = fill(EFFECT_D_POOLED_METHOD, nrow(gdf))
            gdf[!, :mannwhitney_method] = fill(MANNWHITNEY_METHOD, nrow(gdf))
            gdf[!, :is_fdr] = Float64.(gdf.q_value) .< FDR_ALPHA
            gdf[!, :is_largest_abs_effect] = falses(nrow(gdf))
            gdf[argmax(abs.(Float64.(gdf.effect_d_pooled))), :is_largest_abs_effect] = true
            CSV.write(joinpath(tab_dir, "global_mean_wpli_statistics.csv"), gdf)
            for band in string.(gdf.band)
                try
                    fig_raincloud_band(tab_dir, figs_dir, manifest, band, cl)
                catch e
                    @warn "Raincloud $band ($cl) fallido: $e"
                end
            end
        end

        net_global = NamedTuple[]
        for band in bands
            ms_mats   = get(ms_data,   band, Tuple{Vector{String},Matrix{Float64}}[])
            ctrl_mats = get(ctrl_data, band, Tuple{Vector{String},Matrix{Float64}}[])
            (length(ms_mats) < 2 || length(ctrl_mats) < 2) && continue
            all_ch_sets = [Set(ch) for (ch,_) in vcat(ms_mats, ctrl_mats)]
            common_set = reduce(intersect, all_ch_sets)
            ch_ref = ctrl_mats[1][1]
            common_ch = filter(c -> c in common_set, ch_ref)
            length(common_ch) < 2 && continue
            ms_s = Float64[]; ct_s = Float64[]
            for (ch,W) in ms_mats
                s,_,_,_ = nodal_strength_degree(realign_matrix(ch,W,common_ch); density=graph_dens)
                push!(ms_s, mean(s))
            end
            for (ch,W) in ctrl_mats
                s,_,_,_ = nodal_strength_degree(realign_matrix(ch,W,common_ch); density=graph_dens)
                push!(ct_s, mean(s))
            end
            band_idx = something(findfirst(==(band), bands), 0)
            cond_offset = cl == "EC" ? 0 : 10_000
            seed = BOOTSTRAP_SEED + 100_000 + cond_offset + band_idx
            s = independent_group_summary(ms_s, ct_s; seed=seed)
            push!(net_global, (
                schema_version=TRANSVERSAL_SCHEMA_VERSION,
                statistics_source=STATISTICS_SOURCE, fdr_scope=GLOBAL_FDR_SCOPE,
                band=band, metric="mean_strength", n_channels=length(common_ch),
                ms_mean=s.ms_mean, ms_ci_low=s.ms_ci_low, ms_ci_high=s.ms_ci_high,
                ctrl_mean=s.ctrl_mean, ctrl_ci_low=s.ctrl_ci_low, ctrl_ci_high=s.ctrl_ci_high,
                diff=s.diff, diff_ci_low=s.diff_ci_low, diff_ci_high=s.diff_ci_high,
                p_value=s.p_mannwhitney, p_welch=s.p_welch,
                effect_d_pooled=s.effect_d_pooled,
                effect_d_pooled_ci_low=s.effect_d_pooled_ci_low,
                effect_d_pooled_ci_high=s.effect_d_pooled_ci_high,
                effect_rrb=s.effect_rrb,
                probability_superiority=s.probability_superiority,
                n_ms=s.n_ms, n_ctrl=s.n_ctrl, bootstrap_seed=seed,
            ))
        end
        if !isempty(net_global)
            ng_df = DataFrame(net_global)
            ng_df[!, :q_value] = bh_qvalues(Float64.(ng_df.p_value))
            ng_df[!, :fdr_family_size] = fill(nrow(ng_df), nrow(ng_df))
            ng_df[!, :bootstrap_method] = fill(BOOTSTRAP_METHOD, nrow(ng_df))
            ng_df[!, :bootstrap_iterations] = fill(BOOTSTRAP_N, nrow(ng_df))
            ng_df[!, :rrb_method] = fill(RRB_METHOD, nrow(ng_df))
            ng_df[!, :effect_d_pooled_method] = fill(EFFECT_D_POOLED_METHOD, nrow(ng_df))
            ng_df[!, :mannwhitney_method] = fill(MANNWHITNEY_METHOD, nrow(ng_df))
            ng_df[!, :is_fdr] = Float64.(ng_df.q_value) .< FDR_ALPHA
            CSV.write(joinpath(tab_dir, "network_global_statistics.csv"), ng_df)
        end

        bp_store = Dict{Tuple{String,String}, Tuple{Vector{Float64},Vector{Float64}}}()
        for s in included_subj
            bp = load_band_power(res_root, s.subject_id, s.session_id, cond, bands)
            bp === nothing && continue
            for ch in keys(bp), band in bands
                haskey(bp[ch], band) || continue
                key = (ch, band)
                haskey(bp_store, key) || (bp_store[key] = (Float64[], Float64[]))
                (s.group == "MS" ? push!(bp_store[key][1], bp[ch][band]) : push!(bp_store[key][2], bp[ch][band]))
            end
        end
        bp_stat_rows = NamedTuple[]
        for ((ch,band), (vms,vct)) in bp_store
            (length(vms) < 2 || length(vct) < 2) && continue
            pmw, _ = mannwhitney_p(vms, vct)
            _, d = welch_t(vms, vct)
            push!(bp_stat_rows, (
                  schema_version=TRANSVERSAL_SCHEMA_VERSION,
                  statistics_source=STATISTICS_SOURCE, fdr_scope=POWER_FDR_SCOPE,
                  channel=ch, band=band, ms_mean=mean(vms), ctrl_mean=mean(vct),
                  diff=mean(vms)-mean(vct), p_value=pmw, effect_d_pooled=d,
                  n_ms=length(vms), n_ctrl=length(vct)))
        end
        if !isempty(bp_stat_rows)
            bp_df = DataFrame(bp_stat_rows)
            q_all = fill(1.0, nrow(bp_df))
            family_size = zeros(Int, nrow(bp_df))
            coverage_low = falses(nrow(bp_df))
            coverage_ms_pct = zeros(Float64, nrow(bp_df))
            coverage_ctrl_pct = zeros(Float64, nrow(bp_df))
            for band in bands
                idx = findall(i -> string(bp_df.band[i]) == band, 1:nrow(bp_df))
                isempty(idx) && continue
                q_all[idx] = bh_qvalues(Float64.(bp_df.p_value[idx]))
                family_size[idx] .= length(idx)
                max_ms = maximum(Int.(bp_df.n_ms[idx]))
                max_ctrl = maximum(Int.(bp_df.n_ctrl[idx]))
                coverage_ms_pct[idx] = 100.0 .* Int.(bp_df.n_ms[idx]) ./ max(max_ms, 1)
                coverage_ctrl_pct[idx] = 100.0 .* Int.(bp_df.n_ctrl[idx]) ./ max(max_ctrl, 1)
                coverage_low[idx] = (Int.(bp_df.n_ms[idx]) .< ceil(Int, 0.7 * max_ms)) .|
                                    (Int.(bp_df.n_ctrl[idx]) .< ceil(Int, 0.7 * max_ctrl))
            end
            bp_df[!, :q_value] = q_all
            bp_df[!, :fdr_family_size] = family_size
            bp_df[!, :coverage_ms_pct] = coverage_ms_pct
            bp_df[!, :coverage_ctrl_pct] = coverage_ctrl_pct
            bp_df[!, :coverage_low] = coverage_low
            bp_df[!, :is_fdr] = Float64.(bp_df.q_value) .< FDR_ALPHA
            CSV.write(joinpath(tab_dir, "band_power_group_statistics.csv"), bp_df)
            CSV.write(joinpath(tab_dir, "significant_band_power_differences.csv"), filter(r -> r.q_value < 0.05, bp_df))

            try
                fig_power_effect_heatmap(bp_df, figs_dir, manifest, cl)
            catch e
                @warn "Heatmap efecto potencia canal×banda ($cl) fallido: $e"
            end

            xy = isempty(included_subj) ? Dict{String,Tuple{Float64,Float64}}() :
                 load_electrode_xy(bids_root, included_subj[1].subject_id, included_subj[1].session_id)
            for band in bands
                sub = filter(r -> string(r.band) == band, eachrow(bp_df))
                isempty(sub) && continue
                chs = [string(r.channel) for r in sub]; dvs = [Float64(r.diff) for r in sub]
                try
                    save_topo_delta(joinpath(figs_dir, "topo_diff_bandpower_$(band).png"), chs, dvs, xy,
                                    "Δ band power (MS−Control) — $band ($cl)"; colorbar_label="Δ power (MS−Control)")
                catch e
                    @warn "Topo band power $band fallido: $e"
                end
            end
        end

        if !isempty(subject_means)
            try
                save_group_means(joinpath(figs_dir, "group_mean_wpli_by_band.png"), DataFrame(subject_means))
            catch e
                @warn "Figura group means fallida: $e"
            end
        end

        n_total_sig = isempty(band_stats_rows) ? 0 : sum(r.n_sig for r in band_stats_rows)
        best_band = honest_best_band(band_stats_rows)

        open(joinpath(out_dir, "transversal_summary.json"), "w") do io
            pairs_json = join(["  \"$(r.band)_n_sig\": $(r.n_sig), \"$(r.band)_ctrl\": $(r.ctrl_mean), \"$(r.band)_ms\": $(r.ms_mean)"
                                for r in band_stats_rows], ",\n")
            write(io, """{
  "schema_version": $TRANSVERSAL_SCHEMA_VERSION,
  "statistics_source": "$STATISTICS_SOURCE",
  "n_ms": $n_ms,
  "n_ctrl": $n_ctrl,
  "n_ms_design": $N_MS_DESIGN,
  "n_ctrl_design": $N_CTRL_DESIGN,
  "n_included": $(n_ms + n_ctrl),
  "n_excluded": $n_excl,
  "n_excluded_qc": $n_excl_qc,
  "n_excluded_data": $n_excl_data,
  "n_total_sig": $n_total_sig,
  "n_bands": $(length(band_stats_rows)),
  "best_band": "$best_band",
  "cond": "$cl",
  "design": "case_control_T1",
  "session_policy": "T1_only",
  "test": "mannwhitney",
  "mannwhitney_method": "$MANNWHITNEY_METHOD",
  "fdr": "bh",
  "fdr_scope_edges": "$EDGE_FDR_SCOPE",
  "fdr_scope_global": "$GLOBAL_FDR_SCOPE",
  "fdr_scope_power": "$POWER_FDR_SCOPE",
  "bootstrap_method": "$BOOTSTRAP_METHOD",
  "bootstrap_iterations": $BOOTSTRAP_N,
  "bootstrap_seed": $BOOTSTRAP_SEED,
  "quantile_method": "$QUANTILE_METHOD",
  "rrb_method": "$RRB_METHOD",
  "effect_d_pooled_method": "$EFFECT_D_POOLED_METHOD",
  "wpli_method": "$wpli_method",
  "use_dwpli": $use_dwpli,
  "timestamp": "$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))",
  $pairs_json
}""")
        end
        # ── Figuras de manuscrito (solo eyesclosed) ──────────────
        if cond == "eyesclosed"
            try
                fig_forest_by_band(tab_dir, figs_dir, manifest)
                fig_raincloud_alpha(tab_dir, figs_dir, manifest)
                fig_matrices_alpha(tab_dir, figs_dir, manifest)
                fig_network_power_alpha(tab_dir, figs_dir, manifest)
            catch e
                @warn "Figuras de manuscrito transversal fallidas: $e"
            end
        end

        if !isempty(manifest)
            write_manifest(joinpath(tab_dir, "figures_manifest.tsv"), manifest)
            println("  Figuras generadas ($(length(manifest))): " * join([r.figure_name for r in manifest], ", "))
        end

        println("  ✅ Guardado en: $out_dir\n")
    end

    println("── Interacción grupo×condición (EC×EO) " * "─"^12)
    try
        fig_interaction_ec_eo(res_root, bands)
    catch e
        @warn "Interacción EC×EO fallida: $e"
    end

    println("✅ Análisis transversal completado (eyesclosed y eyesopen en paralelo).")
    println("   Abre el dashboard → Fase 13 Evaluación Transversal")
    return nothing
end

end # module
