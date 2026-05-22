# NeuroMIND/scripts/run_longitudinal_analysis.jl
# Análisis longitudinal: T1 → T2 intra-sujeto (pacientes EM).
# Lee resultados individuales (wpli_*.csv) ya generados por el pipeline
# y produce archivos longitudinales en results/group/longitudinal/{EC|EO}/
#
# Detecta automáticamente sujetos con sesiones T1 y T2 en results/subjects/
# o lee un archivo opcional data/BIDS/longitudinal_pairs.csv con columnas:
#   subject_id, session_t1, session_t2
# (si el archivo existe se usa; si no, se auto-detecta T1/T2 por nombre de sesión)
#
# Uso:
#   julia --project=. scripts/run_longitudinal_analysis.jl
#   julia --project=. scripts/run_longitudinal_analysis.jl config/single_subject.toml

using CSV, DataFrames, Statistics, LinearAlgebra, Dates, TOML, Printf

const PROJ      = dirname(@__DIR__)
const CONFIG_P  = length(ARGS) > 0 ? ARGS[1] :
                  joinpath(PROJ, "config", "single_subject.toml")

cfg_raw   = TOML.parsefile(CONFIG_P)
paths_raw = get(cfg_raw, "paths", Dict{String,Any}())
res_root  = let r = get(paths_raw, "results", "results")
                isabspath(r) ? r : joinpath(PROJ, r)
            end
bids_root = let b = get(paths_raw, "bids_root", "data/BIDS")
                isabspath(b) ? b : joinpath(PROJ, b)
            end
bands_cfg = get(cfg_raw, "bands", Dict(
    "DELTA"    => [0.5, 4.0], "THETA" => [4.0, 8.0], "ALPHA"    => [7.8, 11.7],
    "BETA_LOW" => [12.0, 15.0], "BETA_MID" => [15.0, 18.0],
    "BETA_HIGH"=> [18.0, 30.0], "GAMMA"    => [30.0, 50.0],
))
BANDS = sort(collect(keys(bands_cfg)))

# ─── Helpers ──────────────────────────────────────────────────

function norm_cond(c::String)::String
    lc = lowercase(c)
    lc in ("ec","eyesclosed") && return "eyesclosed"
    lc in ("eo","eyesopen")   && return "eyesopen"
    return lc
end

function is_t1_session(s::String)::Bool
    uppercase(s) in ("T1","BASELINE","BL","V1","VISIT1","S1","PRE")
end

function load_wpli(subj_id::String, sess_id::String, cond::String, band::String)
    path = joinpath(res_root, "subjects",
                    "sub-$(subj_id)", "ses-$(sess_id)",
                    norm_cond(cond), "wpli_$(band).csv")
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    isempty(df) && return nothing
    ch = string.(df[!, 1])
    n  = length(ch)
    W  = Matrix{Float64}(undef, n, n)
    for (j, c) in enumerate(ch)
        col_sym = Symbol(c)
        hasproperty(df, col_sym) || return nothing
        W[:, j] = Float64.(df[!, col_sym])
    end
    return (ch, W)
end

function bh_qvalues(p::Vector{Float64})::Vector{Float64}
    m = length(p); m == 0 && return Float64[]
    ord = sortperm(p); rnk = invperm(ord)
    q   = p .* m ./ rnk
    qs  = q[ord]
    for i in (m-1):-1:1; qs[i] = min(qs[i], qs[i+1]); end
    qo = zeros(m); qo[ord] = qs
    return min.(qo, 1.0)
end

"""Paired t-test: returns (p_value, cohen_dz)."""
function paired_t(before::Vector{Float64}, after::Vector{Float64})
    length(before) == length(after) || return (1.0, 0.0)
    diffs = after .- before
    n     = length(diffs)
    n < 2 && return (1.0, 0.0)
    μ = mean(diffs)
    s = std(diffs)
    s < 1e-12 && return (1.0, 0.0)
    t = μ * sqrt(n) / s
    p = clamp(2.0 * (1.0 - 0.5*(1.0 + erf(abs(t)/sqrt(2.0)))), 0.0, 1.0)
    d = μ / s   # Cohen's dz
    return (p, d)
end

function realign_matrix(ch::Vector{String}, W::Matrix{Float64},
                         common_ch::Vector{String})::Matrix{Float64}
    idx = [findfirst(==(c), ch) for c in common_ch]
    any(isnothing, idx) && error("Canales no encontrados")
    W[idx, idx]
end

function save_mat_csv(path::String, W::Matrix{Float64}, ch::Vector{String})
    df = DataFrame(hcat(ch, W), vcat(["channel"], ch))
    CSV.write(path, df)
end

# ─── Detectar / cargar pares longitudinales ───────────────────

pairs_path = joinpath(bids_root, "longitudinal_pairs.csv")

struct SubjPair
    subject_id::String
    session_t1::String
    session_t2::String
end

all_pairs = SubjPair[]

if isfile(pairs_path)
    # Usar archivo explícito
    pdf = CSV.read(pairs_path, DataFrame)
    rename!(pdf, Dict(n => Symbol(lowercase(string(n))) for n in names(pdf)))
    for row in eachrow(pdf)
        push!(all_pairs, SubjPair(
            string(row.subject_id),
            string(row.session_t1),
            string(row.session_t2),
        ))
    end
    println("📋 Pares longitudinales leídos desde: $pairs_path")
else
    # Auto-detectar en results/subjects/
    subj_root = joinpath(res_root, "subjects")
    if isdir(subj_root)
        for sd in readdir(subj_root)
            !startswith(sd, "sub-") && continue
            subj_id   = sd[5:end]
            sess_dir  = joinpath(subj_root, sd)
            sessions  = filter(s -> startswith(s, "ses-"),
                               readdir(sess_dir))
            t1_sesses = filter(s -> is_t1_session(s[5:end]), sessions)
            t2_sesses = filter(s -> !is_t1_session(s[5:end]), sessions)
            for t1s in t1_sesses, t2s in t2_sesses
                push!(all_pairs, SubjPair(subj_id, t1s[5:end], t2s[5:end]))
            end
        end
    end
    println("🔍 Pares auto-detectados en results/subjects/: $(length(all_pairs))")
    println("   (Crea data/BIDS/longitudinal_pairs.csv para control explícito)")
end

if isempty(all_pairs)
    @warn """
    No se encontraron pares T1/T2.

    Opciones:
    1. Crea data/BIDS/longitudinal_pairs.csv con columnas:
       subject_id, session_t1, session_t2
    2. Asegúrate de que las sesiones se llaman T1 (baseline) y T2 (seguimiento)
    """
    exit(1)
end

println("   Pares encontrados: $(length(all_pairs))")
println()

# ─── Proceso por condición ────────────────────────────────────

for cond in ["EC", "EO"]
    println("── Condición: $cond " * "─"^40)
    out_dir = joinpath(res_root, "group", "longitudinal", cond)
    mkpath(out_dir)

    t1_data   = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
    t2_data   = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
    # Para paired t-test guardamos por sujeto y banda
    paired_data = Dict{String, Vector{Tuple{Matrix{Float64}, Matrix{Float64}}}}()

    subject_means = NamedTuple[]
    paired_info   = NamedTuple[]

    for sp in all_pairs
        t1_ok   = Dict{String,Bool}()
        t2_ok   = Dict{String,Bool}()
        pair_ok = 0

        for band in BANDS
            r1 = load_wpli(sp.subject_id, sp.session_t1, cond, band)
            r2 = load_wpli(sp.subject_id, sp.session_t2, cond, band)
            r1 === nothing && (t1_ok[band] = false; continue)
            r2 === nothing && (t2_ok[band] = false; continue)
            t1_ok[band] = true; t2_ok[band] = true; pair_ok += 1

            (ch1, W1) = r1; (ch2, W2) = r2
            haskey(t1_data, band) || (t1_data[band] = [])
            haskey(t2_data, band) || (t2_data[band] = [])
            haskey(paired_data, band) || (paired_data[band] = [])
            push!(t1_data[band], (ch1, W1))
            push!(t2_data[band], (ch2, W2))
            push!(paired_data[band], (W1, W2))

            n = length(ch1)
            up = [(i, j) for i in 1:n for j in (i+1):n]
            for (W, tp) in [(W1,"T1"),(W2,"T2")]
                isempty(up) && continue
                μ = mean(W[i,j] for (i,j) in up)
                push!(subject_means, (
                    subject_id = sp.subject_id,
                    timepoint  = tp,
                    band       = band,
                    cond       = cond,
                    mean_wpli  = round(μ, digits=5),
                ))
            end
        end

        included = pair_ok > 0
        push!(paired_info, (
            subject_id      = sp.subject_id,
            session_t1      = sp.session_t1,
            session_t2      = sp.session_t2,
            n_bands_ok      = pair_ok,
            included        = included,
            excluded_reason = included ? "" : "Sin datos wPLI pareados ($cond)",
        ))
    end

    n_paired = count(r -> r.included, paired_info)
    n_loss   = count(r -> !r.included, paired_info)
    println("  Pares completos: $n_paired | Pérdidas: $n_loss")

    CSV.write(joinpath(out_dir, "paired_subjects.csv"), DataFrame(paired_info))
    !isempty(subject_means) &&
        CSV.write(joinpath(out_dir, "subject_band_means.csv"), DataFrame(subject_means))

    if n_paired < 1
        println("  ⚠  Análisis omitido: sin pares completos T1/T2")
        open(joinpath(out_dir, "longitudinal_summary.json"), "w") do io
            write(io, """{"n_paired":0,"n_t1":$(length(all_pairs)),"n_t2":0,"n_loss":$n_loss,"n_total_sig":0,"n_bands":0,"cond":"$cond","timestamp":"$(Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"))"}""")
        end
        continue
    end

    band_stats_rows = NamedTuple[]

    for band in BANDS
        t1_mats  = get(t1_data,     band, Tuple{Vector{String},Matrix{Float64}}[])
        t2_mats  = get(t2_data,     band, Tuple{Vector{String},Matrix{Float64}}[])
        pd_mats  = get(paired_data, band, Tuple{Matrix{Float64},Matrix{Float64}}[])
        (isempty(t1_mats) || isempty(t2_mats)) && continue

        # Alinear al conjunto común de canales
        all_ch_sets = [Set(ch) for (ch, _) in vcat(t1_mats, t2_mats)]
        common_set  = reduce(intersect, all_ch_sets)
        ch_ref      = t1_mats[1][1]
        common_ch   = filter(c -> c in common_set, ch_ref)
        n           = length(common_ch)
        n < 2 && continue

        get_W(ch, W) = realign_matrix(ch, W, common_ch)

        Wt1   = mean(get_W(ch, W) for (ch, W) in t1_mats)
        Wt2   = mean(get_W(ch, W) for (ch, W) in t2_mats)
        Wdiff = Wt2 .- Wt1

        save_mat_csv(joinpath(out_dir, "longitudinal_connectivity_t1_$(band).csv"), Wt1, common_ch)
        save_mat_csv(joinpath(out_dir, "longitudinal_connectivity_t2_$(band).csv"), Wt2, common_ch)
        save_mat_csv(joinpath(out_dir, "longitudinal_difference_$(band).csv"),      Wdiff, common_ch)

        upper_idx = [(i, j) for i in 1:n for j in (i+1):n]
        p_vec = zeros(length(upper_idx))
        d_vec = zeros(length(upper_idx))

        for (k, (i, j)) in enumerate(upper_idx)
            # Extraer la conectividad del par en T1 y T2 para cada sujeto
            v1 = [get_W(ch1, W1)[i, j] for (ch1, W1) in t1_mats]
            v2 = [get_W(ch2, W2)[i, j] for (ch2, W2) in t2_mats]
            # Sólo usamos pares completos
            n_min = min(length(v1), length(v2))
            n_min < 2 && continue
            p_vec[k], d_vec[k] = paired_t(v1[1:n_min], v2[1:n_min])
        end
        q_vec = bh_qvalues(p_vec)

        stat_rows = [(
            ch_a     = common_ch[i],
            ch_b     = common_ch[j],
            t1_mean  = round(Wt1[i, j],   digits=5),
            t2_mean  = round(Wt2[i, j],   digits=5),
            diff     = round(Wdiff[i, j], digits=5),
            p_value  = round(p_vec[k],    digits=5),
            q_value  = round(q_vec[k],    digits=5),
            effect_d = round(d_vec[k],    digits=4),
        ) for (k, (i, j)) in enumerate(upper_idx)]

        CSV.write(joinpath(out_dir, "longitudinal_statistics_$(band).csv"), DataFrame(stat_rows))

        sig_rows = filter(r -> r.q_value < 0.05, stat_rows)
        sig_df   = isempty(sig_rows) ?
            DataFrame(ch_a=String[], ch_b=String[], t1_mean=Float64[],
                      t2_mean=Float64[], diff=Float64[],
                      p_value=Float64[], q_value=Float64[], effect_d=Float64[]) :
            DataFrame(sig_rows)
        CSV.write(joinpath(out_dir, "significant_longitudinal_edges_$(band).csv"), sig_df)

        n_sig      = length(sig_rows)
        mean_t1    = mean(Wt1[i, j] for (i, j) in upper_idx)
        mean_t2    = mean(Wt2[i, j] for (i, j) in upper_idx)
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
        ))
    end

    !isempty(band_stats_rows) &&
        CSV.write(joinpath(out_dir, "band_statistics_longitudinal.csv"), DataFrame(band_stats_rows))

    n_total_sig = isempty(band_stats_rows) ? 0 : sum(r.n_sig for r in band_stats_rows)
    best_band   = isempty(band_stats_rows) ? "" :
                  band_stats_rows[argmax(abs(r.diff_mean) for r in band_stats_rows)].band

    open(joinpath(out_dir, "longitudinal_summary.json"), "w") do io
        pairs_json = join(
            ["  \"$(r.band)_n_sig\": $(r.n_sig), \"$(r.band)_t1\": $(r.t1_mean), \"$(r.band)_t2\": $(r.t2_mean)"
             for r in band_stats_rows], ",\n"
        )
        write(io, """{
  "n_paired": $n_paired,
  "n_t1": $(length(all_pairs)),
  "n_t2": $n_paired,
  "n_loss": $n_loss,
  "n_total_sig": $n_total_sig,
  "n_bands": $(length(band_stats_rows)),
  "best_band": "$best_band",
  "cond": "$cond",
  "timestamp": "$(Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"))",
  $pairs_json
}""")
    end
    println("  ✅ Guardado en: $out_dir\n")
end

println("✅ Análisis longitudinal completado.")
println("   Abre el dashboard → Fase 14 Evaluación Longitudinal")
