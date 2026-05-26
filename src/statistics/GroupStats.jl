# NeuroMIND/src/statistics/GroupStats.jl
# Estadística grupal: Mann-Whitney, Wilcoxon, Spearman.
# Estrategia del estudio BRAIN (BRAIN.pdf pág. 33):
#   - Comparación grupal: Mann-Whitney U (no paramétrico)
#   - Seguimiento T1→T2: Wilcoxon signed-rank (muestras relacionadas)
#   - Correlaciones clínicas: Spearman (EDSS, fatiga, cognición, carga lesional)

"""
    mann_whitney_test(a, b; alpha=0.05) -> StatResult

Test Mann-Whitney U para comparar dos grupos independientes (MS vs Control).
"""
function mann_whitney_test(a::Vector{Float64}, b::Vector{Float64}; alpha::Float64=0.05)::StatResult
    na, nb = length(a), length(b)
    (na < 2 || nb < 2) && return _empty_stat("mann_whitney")

    # U statistic via rank sum
    all_vals = vcat(a, b)
    labels   = vcat(fill(1, na), fill(2, nb))
    order    = sortperm(all_vals)
    ranks    = _assign_ranks(all_vals[order])

    R1   = sum(ranks[labels[order] .== 1])
    U1   = R1 - na * (na + 1) / 2
    U2   = na * nb - U1
    U    = min(U1, U2)

    # Aproximación normal para n grandes
    μU   = na * nb / 2
    σU   = sqrt(na * nb * (na + nb + 1) / 12)
    z    = (U - μU) / σU
    p    = 2 * (1 - _norm_cdf(abs(z)))

    # Effect size r = z / sqrt(N)
    r    = abs(z) / sqrt(na + nb)

    StatResult("mann_whitney", U, p, p,  # p_adjusted se fija desde fuera con FDR
               r, p <= alpha,
               mean(a), mean(b), na, nb)
end

"""
    wilcoxon_signed_rank(before, after; alpha=0.05) -> StatResult

Test Wilcoxon signed-rank para muestras relacionadas (T1 vs T2 del mismo paciente).
"""
function wilcoxon_signed_rank(before::Vector{Float64}, after::Vector{Float64}; alpha::Float64=0.05)::StatResult
    length(before) == length(after) || error("Vectores de distinto tamaño")
    n     = length(before)
    diffs = after .- before
    diffs = filter(!=(0.0), diffs)  # eliminar ceros
    nd    = length(diffs)
    nd < 2 && return _empty_stat("wilcoxon")

    ranks  = _assign_ranks(abs.(diffs))
    W_plus  = sum(ranks[diffs .> 0])
    W_minus = sum(ranks[diffs .< 0])
    W       = min(W_plus, W_minus)

    μW  = nd * (nd + 1) / 4
    σW  = sqrt(nd * (nd + 1) * (2nd + 1) / 24)
    z   = (W - μW) / σW
    p   = 2 * (1 - _norm_cdf(abs(z)))
    r   = abs(z) / sqrt(nd)

    StatResult("wilcoxon", W, p, p, r, p <= alpha,
               mean(before), mean(after), n, n)
end

"""
    spearman_correlation(x, y; alpha=0.05) -> StatResult

Correlación de Spearman entre una métrica EEG y una variable clínica (EDSS, etc.).
"""
function spearman_correlation(x::Vector{Float64}, y::Vector{Float64}; alpha::Float64=0.05)::StatResult
    n = length(x)
    (n != length(y) || n < 3) && return _empty_stat("spearman")

    rx   = _assign_ranks(x)
    ry   = _assign_ranks(y)
    ρ    = cor(rx, ry)
    t    = ρ * sqrt((n - 2) / (1 - ρ^2 + 1e-12))
    p    = 2 * (1 - _t_cdf(abs(t), n - 2))

    StatResult("spearman", ρ, p, p, ρ, p <= alpha,
               mean(x), mean(y), n, n)
end

"""
    compare_groups_stats(subjects, session_id, condition, band, cfg) -> DataFrame

Tabla de estadísticos por par de canales (MS vs Control).
"""
function compare_groups_stats(
    subjects::Vector{Subject},
    session_id::String,
    condition::String,
    band::String,
    cfg::PipelineConfig
)::DataFrame

    ms_subjs   = filter(s -> s.group == "MS",      subjects)
    ctrl_subjs = filter(s -> s.group == "Control",  subjects)
    alpha      = get(cfg.statistics, "alpha", 0.05)

    ms_mats   = _collect_wpli(ms_subjs,   session_id, condition, band)
    ctrl_mats = _collect_wpli(ctrl_subjs, session_id, condition, band)

    (isempty(ms_mats) || isempty(ctrl_mats)) &&
        return DataFrame(from=String[], to=String[], U=Float64[], p=Float64[],
                         p_fdr=Float64[], effect_r=Float64[], significant=Bool[],
                         mean_ms=Float64[], mean_ctrl=Float64[])

    n_ch  = size(ms_mats[1], 1)
    ch    = ms_subjs[1].sessions[session_id].connectivity[condition].channel_names

    rows = NamedTuple[]
    for i in 1:n_ch, j in (i+1):n_ch
        a = [m[i,j] for m in ms_mats]
        b = [m[i,j] for m in ctrl_mats]
        sr = mann_whitney_test(a, b; alpha)
        push!(rows, (from=ch[i], to=ch[j], U=sr.statistic, p=sr.p_value,
                     p_fdr=sr.p_adjusted, effect_r=sr.effect_size,
                     significant=sr.significant, mean_ms=sr.group_a_mean,
                     mean_ctrl=sr.group_b_mean))
    end

    df = DataFrame(rows)

    # Aplicar FDR sobre todos los p-valores del triángulo superior
    p_vec  = df.p
    α_fdr  = get(cfg.statistics, "fdr_q", 0.05)
    thr    = fdr_correction(p_vec; alpha=α_fdr)
    df.p_fdr      .= p_vec
    df.significant .= p_vec .<= thr

    return sort(df, :p)
end

# ─── Helpers privados ─────────────────────────────────────────

function _assign_ranks(v::Vector{Float64})::Vector{Float64}
    n     = length(v)
    order = sortperm(v)
    ranks = zeros(Float64, n)
    i     = 1
    while i <= n
        j = i
        while j < n && v[order[j+1]] == v[order[i]]
            j += 1
        end
        r_avg = mean(i:j)
        for k in i:j
            ranks[order[k]] = r_avg
        end
        i = j + 1
    end
    return ranks
end

# Abramowitz & Stegun 7.1.26 — max |ε| ≤ 1.5×10⁻⁷ (no SpecialFunctions needed)
function _erf_approx(x::Float64)::Float64
    t = 1.0 / (1.0 + 0.3275911 * abs(x))
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 +
               t * (-1.453152027 + t * 1.061405429))))
    sign(x) * (1.0 - poly * exp(-x * x))
end

function _norm_cdf(z::Float64)::Float64
    0.5 * (1.0 + _erf_approx(z / sqrt(2.0)))
end

function _t_cdf(t::Float64, df::Real)::Float64
    df > 30 && return _norm_cdf(t)
    x = df / (df + t^2)
    clamp(1.0 - 0.5 * x^(df/2), 0.0, 1.0)
end

function _empty_stat(name::String)::StatResult
    StatResult(name, NaN, 1.0, 1.0, 0.0, false, NaN, NaN, 0, 0)
end

function _collect_wpli(subjects, session_id, condition, band)
    [s.sessions[session_id].connectivity[condition].matrices[band]
     for s in subjects
     if haskey(s.sessions, session_id) &&
        haskey(s.sessions[session_id].connectivity, condition) &&
        haskey(s.sessions[session_id].connectivity[condition].matrices, band)]
end
