# NeuroMIND/src/longitudinal/LongitudinalAnalysis.jl
# Análisis temporal y de grupo: longitudinal + cross-seccional.

"""
    compute_longitudinal(subj::Subject, condition::String, band::String) -> LongitudinalAnalysis

Recopila la conectividad de todas las visitas de un sujeto y construye la trayectoria temporal.
"""
function compute_longitudinal(
    subj::Subject,
    condition::String,
    band::String
)::LongitudinalAnalysis

    visits   = sort(collect(keys(subj.sessions)))
    matrices = Matrix{Float64}[]
    metrics  = GraphMetrics[]
    ch_names = String[]

    for sess_id in visits
        sess = subj.sessions[sess_id]
        conn = get(sess.connectivity, condition, nothing)
        conn === nothing && continue
        haskey(conn.matrices, band) || continue

        push!(matrices, conn.matrices[band])
        isempty(ch_names) && (ch_names = conn.channel_names)

        gm = get(sess.graph_metrics, "$(condition)_$(band)", nothing)
        gm !== nothing && push!(metrics, gm)
    end

    return LongitudinalAnalysis(
        subj.id, condition, band,
        visits, matrices, metrics, ch_names
    )
end

"""
    group_mean_connectivity(subjects::Vector{Subject},
                            session_id::String, condition::String, band::String) -> Matrix{Float64}

Calcula la media de conectividad de un grupo de sujetos para una sesión/condición/banda.
"""
function group_mean_connectivity(
    subjects::Vector{Subject},
    session_id::String,
    condition::String,
    band::String
)::Matrix{Float64}

    mats = Matrix{Float64}[]
    for subj in subjects
        sess = get(subj.sessions, session_id, nothing)
        sess === nothing && continue
        conn = get(sess.connectivity, condition, nothing)
        conn === nothing && continue
        haskey(conn.matrices, band) || continue
        push!(mats, conn.matrices[band])
    end

    isempty(mats) && error("No hay datos de conectividad disponibles para los sujetos dados")
    return mean(cat(mats...; dims=3), dims=3)[:, :, 1]
end

"""
    compare_groups(subjects_a, subjects_b, session_id, condition, band) -> GroupAnalysis

Compara dos grupos calculando la diferencia de medias y un estadístico t por edge.
"""
function compare_groups(
    subjects_a::Vector{Subject},
    subjects_b::Vector{Subject},
    session_id::String,
    condition::String,
    band::String;
    group_a_name::String = "A",
    group_b_name::String = "B"
)::GroupAnalysis

    mean_a = group_mean_connectivity(subjects_a, session_id, condition, band)
    mean_b = group_mean_connectivity(subjects_b, session_id, condition, band)

    n_ch = size(mean_a, 1)

    # Recopilar matrices individuales para estadística
    mats_a = _collect_matrices(subjects_a, session_id, condition, band)
    mats_b = _collect_matrices(subjects_b, session_id, condition, band)

    p_values  = ones(Float64, n_ch, n_ch)
    effect    = zeros(Float64, n_ch, n_ch)

    if !isempty(mats_a) && !isempty(mats_b)
        for i in 1:n_ch, j in (i+1):n_ch
            va = [m[i,j] for m in mats_a]
            vb = [m[i,j] for m in mats_b]
            p, d = _welch_t_test(va, vb)
            p_values[i,j] = p_values[j,i] = p
            effect[i,j]   = effect[j,i]   = d
        end
    end

    ch_names = isempty(mats_a) ? String[] :
               subjects_a[1].sessions[session_id].connectivity[condition].channel_names

    return GroupAnalysis(
        session_id, condition, band,
        group_a_name, group_b_name,
        mean_a, mean_b,
        p_values, effect,
        ch_names,
        length(subjects_a), length(subjects_b)
    )
end

# ─── Helpers privados ─────────────────────────────────────────

function _collect_matrices(subjects, session_id, condition, band)
    [subj.sessions[session_id].connectivity[condition].matrices[band]
     for subj in subjects
     if haskey(subj.sessions, session_id) &&
        haskey(subj.sessions[session_id].connectivity, condition) &&
        haskey(subj.sessions[session_id].connectivity[condition].matrices, band)]
end

function _welch_t_test(a::Vector{Float64}, b::Vector{Float64})
    na, nb = length(a), length(b)
    (na < 2 || nb < 2) && return (1.0, 0.0)

    μa, μb = mean(a), mean(b)
    sa2, sb2 = var(a), var(b)

    se  = sqrt(sa2/na + sb2/nb)
    se < 1e-12 && return (1.0, 0.0)

    t   = (μa - μb) / se
    df  = (sa2/na + sb2/nb)^2 / ((sa2/na)^2/(na-1) + (sb2/nb)^2/(nb-1))
    p   = 2.0 * (1.0 - _t_cdf(abs(t), df))  # aproximación

    # Cohen's d
    sp  = sqrt(((na-1)*sa2 + (nb-1)*sb2) / (na+nb-2))
    d   = sp > 1e-12 ? (μa - μb) / sp : 0.0

    return clamp(p, 0.0, 1.0), d
end

# Normal CDF sin SpecialFunctions — aproximación polinomial A&S 26.2.17, error máx 7.5e-8
function _norm_cdf_approx(z::Float64)::Float64
    z < 0.0 && return 1.0 - _norm_cdf_approx(-z)
    t = 1.0 / (1.0 + 0.2316419 * z)
    poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    return 1.0 - exp(-0.5 * z * z) * poly / sqrt(2.0 * π)
end

# Aproximación normal para df > 30; para df pequeños es menos preciso
function _t_cdf(t::Float64, df::Float64)::Float64
    df > 30.0 && return _norm_cdf_approx(t)
    # Aproximación via regularized incomplete beta
    x = df / (df + t^2)
    return 1.0 - 0.5 * _incomplete_beta_approx(x, df/2, 0.5)
end

function _incomplete_beta_approx(x, a, b)
    # Continuada fraction (single step approximation — solo para orientación)
    clamp(x^a * (1-x)^b, 0.0, 1.0)
end
