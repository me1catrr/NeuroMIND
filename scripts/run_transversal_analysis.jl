# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Análisis transversal profesional (EM vs Control)
# ═══════════════════════════════════════════════════════════════
#
#  Diseño experimental (Fig. 3.1):
#    · Caso-control en T1: EM (N diseño=44) vs Control (N diseño=40)
#    · Sesión T2 NO entra (pertenece al longitudinal)
#    · EC y EO en paralelo: mismo contraste EM vs Control, sin pooling
#
#  Lee derivados del pipeline y produce estadísticas grupales:
#    · Inclusión T1 + QC (qc_decision_table)
#    · wPLI edge-wise: Mann–Whitney + Welch (ref.) + FDR-BH + Cohen d
#    · Band power canal×banda (MS vs Control)
#    · Métricas de red por banda (strength/degree)
#    · Figuras de publicación
#
#  Salida (compatible Fase 13 dashboard):
#    results/transversal/{EC|EO}/
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/run_transversal_analysis.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 25-07-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/<este-script>.jl …
#
#  Prerrequisito
#  ─────────────
#    results/subjects/…          (pipeline por sujeto, sesión T1)
#    data/bids/groups.csv        columnas: subject_id/bids_id, group, session
#    results/qc/qc_decision_table.csv  (recomendado)
#
#  Uso
#  ───
#    julia --project=. scripts/run_transversal_analysis.jl
#    julia --project=. scripts/run_transversal_analysis.jl config/pipeline.toml
#
#  Salida
#  ──────
#    results/transversal/{EC|EO}/
#      group_*.csv · significant_edges_*.csv · band_statistics.csv
#      subject_inclusion.csv · transversal_summary.json · config_snapshot.toml
#      tables/spectral/ · tables/network/ · figures/

using CSV, DataFrames, Statistics, LinearAlgebra, Dates, TOML, Printf, CairoMakie

include(joinpath(@__DIR__, "..", "src", "viz", "GroupVizCommon.jl"))
using .GroupVizCommon

const PROJ     = dirname(@__DIR__)
const CONFIG_P = length(ARGS) > 0 ? ARGS[1] :
                 joinpath(PROJ, "config", "pipeline.toml")

const QC_ALLOWED     = Set(["include", "include_with_warning"])
const N_MS_DESIGN    = 44   # Fig. 3.1
const N_CTRL_DESIGN  = 40

cfg_raw   = TOML.parsefile(CONFIG_P)
paths_raw = get(cfg_raw, "paths", Dict{String,Any}())
res_root  = let r = get(paths_raw, "results", "results")
                isabspath(r) ? r : joinpath(PROJ, r)
            end
bids_root = let b = get(paths_raw, "bids_root", "data/bids")
                p = isabspath(b) ? b : joinpath(PROJ, b)
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

function is_ms_group(g::AbstractString)::Bool
    uppercase(String(g)) in ("MS", "EM", "PATIENT", "PATIENTS", "CASE", "CASES")
end

function is_ctrl_group(g::AbstractString)::Bool
    uppercase(String(g)) in ("CONTROL", "CONTROLS", "HC", "HEALTHY", "CTRL")
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

"""Mann–Whitney U (two-sided, normal approx with tie correction) → (p, effect_r)."""
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
    # Tie correction
    r_all = ranks
    σ2 = na * nb * (na + nb + 1) / 12
    # simple tie term
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
    r = abs(z) / sqrt(na + nb)
    return (p, r)
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

# ─── Carga de datos ───────────────────────────────────────────

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
    if startswith(sid, "M") && !startswith(sid, "MC")
        alt = occursin(r"^M\d$", sid) ? "M0" * sid[2:end] :
              (occursin(r"^M0\d$", sid) ? "M" * sid[3:end] : "")
        !isempty(alt) && haskey(qc, (alt, sess, cc)) && return qc[(alt, sess, cc)]
    end
    return "missing"
end

function qc_ok(dec::String)::Bool
    dec == "missing" && return true
    return dec in QC_ALLOWED
end

# ─── Network ──────────────────────────────────────────────────

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

# ─── Figuras (GroupVizCommon) ─────────────────────────────────

_save_heatmap(path, W, ch, title; cmap=:viridis, diverging=false, colorrange=nothing) =
    save_heatmap(path, W, ch, title; cmap=cmap, diverging=diverging, colorrange=colorrange,
                 colorbar_label=diverging ? "Δ wPLI" : "wPLI")

_save_sig_network(path, ch, sig_rows, title) =
    save_topo_network(path, ch, sig_rows, title)

_save_group_means(path, subject_means_df) =
    save_group_means(path, subject_means_df; band_order=BAND_ORDER)

_save_topo_delta(path, ch, delta_vals, xy, title) =
    save_topo_delta(path, ch, delta_vals, xy, title; colorbar_label="Δ power (MS−Control)")

function _save_explore_or_sig_network(figs_dir, band, cond, common_ch, sig_rows, stats_df)
    if !isempty(sig_rows)
        save_topo_network(joinpath(figs_dir, "sig_network_$(band).png"),
                          common_ch, sig_rows, "Edges sig. FDR — $band ($cond)")
    elseif nrow(stats_df) > 0 && hasproperty(stats_df, :effect_d)
        top = first(sort(stats_df, :effect_d; by=abs, rev=true), min(20, nrow(stats_df)))
        rows = [NamedTuple(r) for r in eachrow(top)]
        save_topo_network(joinpath(figs_dir, "explore_network_topN_$(band).png"),
                          common_ch, rows,
                          "Top-20 |d| (exploratorio) — $band ($cond)")
    end
    return nothing
end

# ─── Cohorte T1 desde groups.csv ──────────────────────────────

struct SubjT1
    subject_id::String
    session_id::String
    group::String   # "MS" | "Control"
end

function load_t1_cohort()::Vector{SubjT1}
    groups_path = joinpath(bids_root, "groups.csv")
    if !isfile(groups_path)
        @warn """
        No se encontró groups.csv en:
          $groups_path
        Crea el archivo con columnas: subject_id/bids_id, group, session
        """
        exit(1)
    end

    gdf = CSV.read(groups_path, DataFrame)
    rename!(gdf, Dict(n => Symbol(lowercase(string(n))) for n in names(gdf)))

    if !hasproperty(gdf, :session_id) && hasproperty(gdf, :session)
        rename!(gdf, :session => :session_id)
    end
    if !hasproperty(gdf, :session_id) && hasproperty(gdf, :ses)
        rename!(gdf, :ses => :session_id)
    end
    if hasproperty(gdf, :bids_id)
        if hasproperty(gdf, :subject_id)
            select!(gdf, Not(:subject_id))
        end
        rename!(gdf, :bids_id => :subject_id)
    end

    required_cols = [:subject_id, :group, :session_id]
    missing_cols  = filter(c -> !hasproperty(gdf, c), required_cols)
    if !isempty(missing_cols)
        error("groups.csv: faltan columnas: $(join(missing_cols, ", "))")
    end

    # Solo T1 (diseño caso-control)
    filter!(row -> is_t1_session(string(row.session_id)), gdf)
    unique!(gdf, [:subject_id, :session_id, :group])

    cohort = SubjT1[]
    for row in eachrow(gdf)
        sid = string(row.subject_id)
        ses = string(row.session_id)
        grp = string(row.group)
        if is_ms_group(grp)
            push!(cohort, SubjT1(sid, ses, "MS"))
        elseif is_ctrl_group(grp)
            push!(cohort, SubjT1(sid, ses, "Control"))
        else
            @warn "Grupo desconocido '$grp' para $sid — omitido"
        end
    end
    # Una fila por sujeto (si hay duplicados de grupo, conservar primero)
    seen = Set{String}()
    uniq = SubjT1[]
    for s in cohort
        s.subject_id in seen && continue
        push!(seen, s.subject_id)
        push!(uniq, s)
    end
    return uniq
end

# ═══════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════

println("=" ^ 62)
println(" NeuroMIND — Análisis transversal EM vs Control (T1)")
println(" $(now())")
println(" Diseño: caso-control T1 (EM=$N_MS_DESIGN vs Ctrl=$N_CTRL_DESIGN)")
println(" EC y EO en paralelo (mismo contraste, sin pooling)")
println(" wPLI=$WPLI_METHOD  dwPLI=$USE_DWPLI  test=mannwhitney  density=$GRAPH_DENS")
println("=" ^ 62)

cohort = load_t1_cohort()
n_ms_cand   = count(s -> s.group == "MS", cohort)
n_ctrl_cand = count(s -> s.group == "Control", cohort)
println("📋 Cohorte T1: EM=$n_ms_cand | Control=$n_ctrl_cand  (candidatos; N diseño $N_MS_DESIGN/$N_CTRL_DESIGN)")

qc_table = load_qc_table()
println("  QC decisions cargadas: $(length(qc_table))")
println()

for cond in ["EC", "EO"]
    println("── Condición: $cond  (contraste transversal EM vs Control) " * "─"^12)
    out_dir  = joinpath(res_root, "transversal", cond)
    figs_dir = joinpath(out_dir, "figures")
    spec_dir = joinpath(out_dir, "tables", "spectral")
    net_dir  = joinpath(out_dir, "tables", "network")
    mkpath(figs_dir); mkpath(spec_dir); mkpath(net_dir)

    open(joinpath(out_dir, "config_snapshot.toml"), "w") do io
        println(io, "# NeuroMIND transversal snapshot — $(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))")
        println(io, "design = \"case_control_T1\"")
        println(io, "session_policy = \"T1_only\"")
        println(io, "condition = \"$cond\"")
        println(io, "n_ms_design = $N_MS_DESIGN")
        println(io, "n_ctrl_design = $N_CTRL_DESIGN")
        println(io, "note = \"EC and EO run in parallel; same EM vs Control contrast; no pooling\"")
        println(io, "test = \"mannwhitney\"")
        println(io, "fdr = \"bh\"")
        println(io, "wpli_method = \"$WPLI_METHOD\"")
        println(io, "use_dwpli = $USE_DWPLI")
        println(io, "graph_density = $GRAPH_DENS")
        println(io, "graph_threshold_method = \"$GRAPH_METH\"")
        println(io, "n_channels_analysis = $N_CH_MONT")
        println(io, "qc_allowed = [\"include\", \"include_with_warning\"]")
        println(io, "source_config = \"$CONFIG_P\"")
    end

    ms_data   = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
    ctrl_data = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
    subject_means = NamedTuple[]
    inclusion     = NamedTuple[]
    included_subj = SubjT1[]

    for s in cohort
        qc_dec = qc_decision(qc_table, s.subject_id, s.session_id, cond)
        ep = load_seg_epochs(s.subject_id, s.session_id, cond)

        bands_ok = String[]
        for band in BANDS
            r = load_wpli(s.subject_id, s.session_id, cond, band)
            r === nothing && continue
            push!(bands_ok, band)
        end

        reasons = String[]
        isempty(bands_ok) && push!(reasons, "Sin datos wPLI ($cond)")
        !qc_ok(qc_dec) && push!(reasons, "QC=$(qc_dec)")

        included = isempty(reasons)
        push!(inclusion, (
            subject_id      = s.subject_id,
            session_id      = s.session_id,
            group           = s.group,
            n_bands_ok      = length(bands_ok),
            included        = included,
            excluded_reason = included ? "" : join(reasons, "; "),
            qc_decision     = qc_dec,
            n_epochs_valid  = ep[1] === missing ? "" : string(ep[1]),
        ))

        included || continue
        push!(included_subj, s)

        for band in BANDS
            r = load_wpli(s.subject_id, s.session_id, cond, band)
            r === nothing && continue
            (ch, W) = r
            target = s.group == "MS" ? ms_data : ctrl_data
            haskey(target, band) || (target[band] = [])
            push!(target[band], (ch, W))
            μ = upper_mean(W)
            isnan(μ) && continue
            push!(subject_means, (
                subject_id = s.subject_id,
                group      = s.group,
                band       = band,
                cond       = cond,
                mean_wpli  = round(μ, digits=5),
            ))
        end
    end

    n_ms   = count(r -> r.group == "MS"      && r.included, inclusion)
    n_ctrl = count(r -> r.group == "Control" && r.included, inclusion)
    n_excl = count(r -> !r.included, inclusion)
    n_excl_qc = count(r -> !r.included && occursin("QC=", r.excluded_reason), inclusion)
    n_excl_data = n_excl - n_excl_qc

    println("  Comparación transversal $cond | EM: $n_ms / $N_MS_DESIGN | Control: $n_ctrl / $N_CTRL_DESIGN | Excluidos: $n_excl")

    CSV.write(joinpath(out_dir, "subject_inclusion.csv"), DataFrame(inclusion))
    !isempty(subject_means) &&
        CSV.write(joinpath(out_dir, "subject_band_means.csv"), DataFrame(subject_means))

    if n_ms < 1 || n_ctrl < 1
        println("  ⚠  Análisis omitido: se necesitan sujetos en ambos grupos para $cond")
        open(joinpath(out_dir, "transversal_summary.json"), "w") do io
            write(io, """{"n_ms":$n_ms,"n_ctrl":$n_ctrl,"n_ms_design":$N_MS_DESIGN,"n_ctrl_design":$N_CTRL_DESIGN,"n_included":$(n_ms+n_ctrl),"n_excluded":$n_excl,"n_excluded_qc":$n_excl_qc,"n_excluded_data":$n_excl_data,"n_total_sig":0,"n_bands":0,"cond":"$cond","design":"case_control_T1","timestamp":"$(Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"))"}""")
        end
        continue
    end

    band_stats_rows = NamedTuple[]
    all_sig_for_figs = Dict{String, Vector{Any}}()

    for band in BANDS
        ms_mats   = get(ms_data,   band, Tuple{Vector{String}, Matrix{Float64}}[])
        ctrl_mats = get(ctrl_data, band, Tuple{Vector{String}, Matrix{Float64}}[])
        (isempty(ms_mats) || isempty(ctrl_mats)) && continue

        all_ch_sets = [Set(ch) for (ch, _) in vcat(ms_mats, ctrl_mats)]
        common_set  = reduce(intersect, all_ch_sets)
        ch_ref      = ctrl_mats[1][1]
        common_ch   = filter(c -> c in common_set, ch_ref)
        n           = length(common_ch)
        n < 2 && continue

        get_W(ch, W) = realign_matrix(ch, W, common_ch)
        Wms   = mean(get_W(ch, W) for (ch, W) in ms_mats)
        Wctrl = mean(get_W(ch, W) for (ch, W) in ctrl_mats)
        Wdiff = Wms .- Wctrl

        save_mat_csv(joinpath(out_dir, "group_connectivity_control_$(band).csv"), Wctrl, common_ch)
        save_mat_csv(joinpath(out_dir, "group_connectivity_ms_$(band).csv"),      Wms,   common_ch)
        save_mat_csv(joinpath(out_dir, "group_difference_$(band).csv"),            Wdiff, common_ch)

        upper_idx = [(i, j) for i in 1:n for j in (i + 1):n]
        p_mw  = ones(length(upper_idx))
        p_w   = ones(length(upper_idx))
        d_vec = zeros(length(upper_idx))

        for (k, (i, j)) in enumerate(upper_idx)
            va = [get_W(ch, W)[i, j] for (ch, W) in ms_mats]
            vb = [get_W(ch, W)[i, j] for (ch, W) in ctrl_mats]
            p_mw[k], _ = mannwhitney_p(va, vb)
            p_w[k], d_vec[k] = welch_t(va, vb)
        end
        q_vec = bh_qvalues(p_mw)

        n_ms_b = length(ms_mats); n_ct_b = length(ctrl_mats)
        stat_rows = [(
            ch_a           = common_ch[i],
            ch_b           = common_ch[j],
            ctrl_mean      = round(Wctrl[i, j], digits=5),
            ms_mean        = round(Wms[i, j],   digits=5),
            diff           = round(Wdiff[i, j], digits=5),
            p_value        = round(p_mw[k],     digits=5),
            p_mannwhitney  = round(p_mw[k],     digits=5),
            p_welch        = round(p_w[k],      digits=5),
            q_value        = round(q_vec[k],    digits=5),
            effect_d       = round(d_vec[k],    digits=4),
            n_ms           = n_ms_b,
            n_ctrl         = n_ct_b,
        ) for (k, (i, j)) in enumerate(upper_idx)]

        CSV.write(joinpath(out_dir, "group_statistics_$(band).csv"), DataFrame(stat_rows))

        sig_rows = filter(r -> r.q_value < 0.05, stat_rows)
        sig_df = isempty(sig_rows) ?
            DataFrame(ch_a=String[], ch_b=String[], ctrl_mean=Float64[],
                      ms_mean=Float64[], diff=Float64[],
                      p_value=Float64[], p_mannwhitney=Float64[], p_welch=Float64[],
                      q_value=Float64[], effect_d=Float64[], n_ms=Int[], n_ctrl=Int[]) :
            DataFrame(sig_rows)
        CSV.write(joinpath(out_dir, "significant_edges_$(band).csv"), sig_df)
        all_sig_for_figs[band] = collect(sig_rows)

        # Top edges by |d| (exploratory)
        top_ord = sortperm([abs(r.effect_d) for r in stat_rows]; rev=true)
        top_n = min(20, length(stat_rows))
        top_rows = [stat_rows[top_ord[k]] for k in 1:top_n]
        CSV.write(joinpath(out_dir, "top_edges_by_effect_$(band).csv"), DataFrame(top_rows))

        n_sig = length(sig_rows)
        mean_ctrl = mean(Wctrl[i, j] for (i, j) in upper_idx)
        mean_ms   = mean(Wms[i, j]   for (i, j) in upper_idx)
        mean_d    = mean(abs, d_vec)

        @printf("  %-10s  %3d ch  %4d pares  %3d sig (q<0.05)  ctrl=%.3f  MS=%.3f\n",
                band, n, length(upper_idx), n_sig, mean_ctrl, mean_ms)

        push!(band_stats_rows, (
            band       = band,
            n_channels = n,
            n_pairs    = length(upper_idx),
            n_sig      = n_sig,
            pct_sig    = round(100.0 * n_sig / max(1, length(upper_idx)), digits=2),
            ctrl_mean  = round(mean_ctrl, digits=5),
            ms_mean    = round(mean_ms,   digits=5),
            diff_mean  = round(mean_ms - mean_ctrl, digits=5),
            mean_p     = round(mean(p_mw), digits=5),
            mean_d     = round(mean_d,     digits=4),
            n_ms       = n_ms_b,
            n_ctrl     = n_ct_b,
        ))

        try
            lim_ab = max(maximum(Wms), maximum(Wctrl), 1e-12)
            _save_heatmap(joinpath(figs_dir, "heatmap_ms_$(band).png"), Wms, common_ch,
                          "wPLI MS — $band ($cond)"; colorrange=(0.0, lim_ab))
            _save_heatmap(joinpath(figs_dir, "heatmap_control_$(band).png"), Wctrl, common_ch,
                          "wPLI Control — $band ($cond)"; colorrange=(0.0, lim_ab))
            _save_heatmap(joinpath(figs_dir, "heatmap_diff_$(band).png"), Wdiff, common_ch,
                          "Δ wPLI (MS−Control) — $band ($cond)"; diverging=true)
            save_heatmap_triplet(joinpath(figs_dir, "heatmap_triplet_$(band).png"),
                Wctrl, Wms, Wdiff, common_ch,
                ("Control — $band ($cond)", "MS — $band ($cond)", "Δ (MS−Ctrl) — $band ($cond)"))
            _save_explore_or_sig_network(figs_dir, band, cond, common_ch, sig_rows, DataFrame(stat_rows))
        catch e
            @warn "Figuras wPLI $band fallidas: $e"
        end

        try
            s_ms, d_ms, ns_ms, _ = nodal_strength_degree(Wms)
            s_ct, d_ct, ns_ct, _ = nodal_strength_degree(Wctrl)
            CSV.write(joinpath(net_dir, "network_metrics_ms_$(band).csv"),
                DataFrame(channel=common_ch, strength=round.(s_ms, digits=5),
                          degree=d_ms, norm_strength=round.(ns_ms, digits=5)))
            CSV.write(joinpath(net_dir, "network_metrics_control_$(band).csv"),
                DataFrame(channel=common_ch, strength=round.(s_ct, digits=5),
                          degree=d_ct, norm_strength=round.(ns_ct, digits=5)))
            CSV.write(joinpath(net_dir, "network_metrics_diff_$(band).csv"),
                DataFrame(channel=common_ch,
                          delta_strength=round.(s_ms .- s_ct, digits=5),
                          delta_degree=d_ms .- d_ct,
                          delta_norm_strength=round.(ns_ms .- ns_ct, digits=5)))
        catch e
            @warn "Network $band fallido: $e"
        end
    end

    !isempty(band_stats_rows) &&
        CSV.write(joinpath(out_dir, "band_statistics.csv"), DataFrame(band_stats_rows))

    # Global mean wPLI statistics (one test per band)
    global_rows = NamedTuple[]
    for band in BANDS
        ms_μ = [r.mean_wpli for r in subject_means if r.band == band && r.group == "MS"]
        ct_μ = [r.mean_wpli for r in subject_means if r.band == band && r.group == "Control"]
        (length(ms_μ) < 2 || length(ct_μ) < 2) && continue
        pmw, _ = mannwhitney_p(Float64.(ms_μ), Float64.(ct_μ))
        pw, d = welch_t(Float64.(ms_μ), Float64.(ct_μ))
        push!(global_rows, (
            band = band,
            ms_mean = round(mean(ms_μ), digits=5),
            ctrl_mean = round(mean(ct_μ), digits=5),
            diff = round(mean(ms_μ) - mean(ct_μ), digits=5),
            p_mannwhitney = round(pmw, digits=5),
            p_welch = round(pw, digits=5),
            effect_d = round(d, digits=4),
            n_ms = length(ms_μ),
            n_ctrl = length(ct_μ),
        ))
    end
    if !isempty(global_rows)
        gdf = DataFrame(global_rows)
        gdf[!, :q_value] = round.(bh_qvalues(Float64.(gdf.p_mannwhitney)), digits=5)
        CSV.write(joinpath(out_dir, "global_mean_wpli_statistics.csv"), gdf)
    end

    # Network global
    net_global = NamedTuple[]
    for band in BANDS
        ms_mats   = get(ms_data,   band, Tuple{Vector{String},Matrix{Float64}}[])
        ctrl_mats = get(ctrl_data, band, Tuple{Vector{String},Matrix{Float64}}[])
        (length(ms_mats) < 2 || length(ctrl_mats) < 2) && continue
        all_ch_sets = [Set(ch) for (ch, _) in vcat(ms_mats, ctrl_mats)]
        common_set = reduce(intersect, all_ch_sets)
        ch_ref = ctrl_mats[1][1]
        common_ch = filter(c -> c in common_set, ch_ref)
        length(common_ch) < 2 && continue
        ms_s = Float64[]; ct_s = Float64[]
        for (ch, W) in ms_mats
            Wa = realign_matrix(ch, W, common_ch)
            s, _, _, _ = nodal_strength_degree(Wa)
            push!(ms_s, mean(s))
        end
        for (ch, W) in ctrl_mats
            Wa = realign_matrix(ch, W, common_ch)
            s, _, _, _ = nodal_strength_degree(Wa)
            push!(ct_s, mean(s))
        end
        pmw, _ = mannwhitney_p(ms_s, ct_s)
        _, d = welch_t(ms_s, ct_s)
        push!(net_global, (
            band = band,
            metric = "mean_strength",
            ms_mean = round(mean(ms_s), digits=5),
            ctrl_mean = round(mean(ct_s), digits=5),
            diff = round(mean(ms_s) - mean(ct_s), digits=5),
            p_value = round(pmw, digits=5),
            effect_d = round(d, digits=4),
            n_ms = length(ms_s),
            n_ctrl = length(ct_s),
        ))
    end
    if !isempty(net_global)
        ng_df = DataFrame(net_global)
        ng_df[!, :q_value] = round.(bh_qvalues(Float64.(ng_df.p_value)), digits=5)
        CSV.write(joinpath(net_dir, "network_global_statistics.csv"), ng_df)
    end

    # Band power group stats
    bp_store = Dict{Tuple{String,String}, Tuple{Vector{Float64},Vector{Float64}}}()
    for s in included_subj
        bp = load_band_power(s.subject_id, s.session_id, cond)
        bp === nothing && continue
        for ch in keys(bp), band in BANDS
            haskey(bp[ch], band) || continue
            key = (ch, band)
            haskey(bp_store, key) || (bp_store[key] = (Float64[], Float64[]))
            if s.group == "MS"
                push!(bp_store[key][1], bp[ch][band])
            else
                push!(bp_store[key][2], bp[ch][band])
            end
        end
    end

    bp_stat_rows = NamedTuple[]
    for ((ch, band), (vms, vct)) in bp_store
        (length(vms) < 2 || length(vct) < 2) && continue
        pmw, _ = mannwhitney_p(vms, vct)
        _, d = welch_t(vms, vct)
        push!(bp_stat_rows, (
            channel   = ch,
            band      = band,
            ms_mean   = round(mean(vms), digits=5),
            ctrl_mean = round(mean(vct), digits=5),
            diff      = round(mean(vms) - mean(vct), digits=5),
            p_value   = round(pmw, digits=5),
            effect_d  = round(d, digits=4),
            n_ms      = length(vms),
            n_ctrl    = length(vct),
        ))
    end

    if !isempty(bp_stat_rows)
        bp_df = DataFrame(bp_stat_rows)
        q_all = fill(1.0, nrow(bp_df))
        for band in BANDS
            idx = findall(i -> string(bp_df.band[i]) == band, 1:nrow(bp_df))
            isempty(idx) && continue
            q_all[idx] = bh_qvalues(Float64.(bp_df.p_value[idx]))
        end
        bp_df[!, :q_value] = round.(q_all, digits=5)
        CSV.write(joinpath(spec_dir, "band_power_group_statistics.csv"), bp_df)
        sig_bp = filter(r -> r.q_value < 0.05, bp_df)
        CSV.write(joinpath(spec_dir, "significant_band_power_differences.csv"), sig_bp)

        xy = Dict{String,Tuple{Float64,Float64}}()
        if !isempty(included_subj)
            xy = load_electrode_xy(included_subj[1].subject_id, included_subj[1].session_id)
        end
        for band in BANDS
            sub = filter(r -> string(r.band) == band, eachrow(bp_df))
            isempty(sub) && continue
            chs = [string(r.channel) for r in sub]
            dvs = [Float64(r.diff) for r in sub]
            try
                _save_topo_delta(joinpath(figs_dir, "topo_diff_bandpower_$(band).png"),
                                 chs, dvs, xy, "Δ band power (MS−Control) — $band ($cond)")
            catch e
                @warn "Topo band power $band fallido: $e"
            end
        end
    end

    if !isempty(subject_means)
        try
            _save_group_means(joinpath(figs_dir, "group_mean_wpli_by_band.png"),
                              DataFrame(subject_means))
        catch e
            @warn "Figura group means fallida: $e"
        end
    end

    n_total_sig = isempty(band_stats_rows) ? 0 : sum(r.n_sig for r in band_stats_rows)
    best_band = honest_best_band(band_stats_rows)

    open(joinpath(out_dir, "transversal_summary.json"), "w") do io
        pairs_json = join(
            ["  \"$(r.band)_n_sig\": $(r.n_sig), \"$(r.band)_ctrl\": $(r.ctrl_mean), \"$(r.band)_ms\": $(r.ms_mean)"
             for r in band_stats_rows],
            ",\n"
        )
        write(io, """{
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
  "cond": "$cond",
  "design": "case_control_T1",
  "session_policy": "T1_only",
  "test": "mannwhitney",
  "fdr": "bh",
  "wpli_method": "$WPLI_METHOD",
  "use_dwpli": $USE_DWPLI,
  "timestamp": "$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))",
  $pairs_json
}""")
    end

    println("  ✅ Guardado en: $out_dir\n")
end

println("✅ Análisis transversal completado (EC y EO en paralelo).")
println("   Abre el dashboard → Fase 13 Evaluación Transversal")
