# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Análisis transversal (MS vs control)
# ═══════════════════════════════════════════════════════════════
#
#  Lee wpli_*.csv ya generados por el pipeline y produce
#  estadísticas grupales en results/transversal/{EC|EO}/.
#
# ───────────────────────────────────────────────────────────────
#  Fichero    scripts/run_transversal_analysis.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 22-07-2026
# ───────────────────────────────────────────────────────────────
#
#  Invocación: julia --project=. scripts/<este-script>.jl …
#  Sin shebang: #!/usr/bin/env julia no activaría --project=.
#
#  Prerrequisito
#  ─────────────
#    results/subjects/…          (pipeline por sujeto)
#    data/bids/groups.csv        columnas: subject_id, group, session_id
#      group ∈ {"MS","EM","Patient"} | {"Control","HC"}
#
#  Config
#  ──────
#    Argumento posicional opcional; defecto: config/pipeline.toml
#    (lee [paths] y [bands])
#
#  Uso
#  ───
#    julia --project=. scripts/run_transversal_analysis.jl
#    julia --project=. scripts/run_transversal_analysis.jl config/pipeline.toml
#
#  Salida
#  ──────
#    results/transversal/{eyesclosed|eyesopen}/

using CSV, DataFrames, Statistics, LinearAlgebra, Dates, TOML, Printf

const PROJ      = dirname(@__DIR__)
const CONFIG_P  = length(ARGS) > 0 ? ARGS[1] :
                  joinpath(PROJ, "config", "pipeline.toml")

cfg_raw   = TOML.parsefile(CONFIG_P)
paths_raw = get(cfg_raw, "paths", Dict{String,Any}())
res_root  = let r = get(paths_raw, "results", "results")
                isabspath(r) ? r : joinpath(PROJ, r)
            end
bids_root = let b = get(paths_raw, "bids_root", "data/BIDS")
                isabspath(b) ? b : joinpath(PROJ, b)
            end
bands_cfg = get(cfg_raw, "bands", Dict(
    "DELTA"    => [0.5, 4.0],
    "THETA"    => [4.0, 8.0],
    "ALPHA"    => [7.8, 11.7],
    "BETA_LOW" => [12.0, 15.0],
    "BETA_MID" => [15.0, 18.0],
    "BETA_HIGH"=> [18.0, 30.0],
    "GAMMA"    => [30.0, 50.0],
))
BANDS = sort(collect(keys(bands_cfg)))

# ─── Helpers ──────────────────────────────────────────────────

function norm_cond(c::AbstractString)::String
    lc = lowercase(String(c))
    lc in ("ec", "eyesclosed") && return "eyesclosed"
    lc in ("eo", "eyesopen")   && return "eyesopen"
    return lc
end

function is_ms_group(g::AbstractString)::Bool
    uppercase(String(g)) in ("MS", "EM", "PATIENT", "PATIENTS", "CASE", "CASES")
end

"""Load wPLI matrix from individual subject results. Returns (channel_names, matrix) or nothing."""
function load_wpli(subj_id::AbstractString, sess_id::AbstractString,
                   cond::AbstractString, band::AbstractString)
    path = joinpath(res_root, "subjects",
                    "sub-$(subj_id)", "ses-$(sess_id)",
                    norm_cond(cond), "tables", "connectivity", "wpli_$(band).csv")
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

"""Compute BH q-values from p-values vector."""
function bh_qvalues(p::Vector{Float64})::Vector{Float64}
    m = length(p); m == 0 && return Float64[]
    ord = sortperm(p); rnk = invperm(ord)
    q   = p .* m ./ rnk
    qs  = q[ord]
    for i in (m - 1):-1:1; qs[i] = min(qs[i], qs[i + 1]); end
    qo = zeros(m); qo[ord] = qs
    return min.(qo, 1.0)
end

"""Welch t-test: returns (p_value, cohen_d)."""
function welch_t(a::Vector{Float64}, b::Vector{Float64})
    na, nb = length(a), length(b)
    (na < 2 || nb < 2) && return (1.0, 0.0)
    μa, μb = mean(a), mean(b)
    sa2, sb2 = var(a), var(b)
    se = sqrt(sa2 / na + sb2 / nb)
    se < 1e-12 && return (1.0, 0.0)
    t  = (μa - μb) / se
    # Normal approximation (good for large N)
    p  = clamp(2.0 * (1.0 - _norm_cdf(abs(t))), 0.0, 1.0)
    sp = sqrt(((na - 1) * sa2 + (nb - 1) * sb2) / max(na + nb - 2, 1))
    d  = sp > 1e-12 ? (μa - μb) / sp : 0.0
    return (p, d)
end
# Normal CDF sin SpecialFunctions — aproximación polinomial A&S 26.2.17, error máx 7.5e-8
function _norm_cdf(z::Float64)::Float64
    z < 0.0 && return 1.0 - _norm_cdf(-z)
    t = 1.0 / (1.0 + 0.2316419 * z)
    poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
    return 1.0 - exp(-0.5 * z * z) * poly / sqrt(2.0 * π)
end

"""Realign wPLI matrix to a common channel set."""
function realign_matrix(ch::Vector{String}, W::Matrix{Float64},
                         common_ch::Vector{String})::Matrix{Float64}
    idx = [findfirst(==(c), ch) for c in common_ch]
    any(isnothing, idx) && error("Canales no encontrados")
    return W[idx, idx]
end

"""Save n×n matrix as CSV with channel column."""
function save_mat_csv(path::String, W::Matrix{Float64}, ch::Vector{String})
    df = DataFrame(hcat(ch, W), vcat(["channel"], ch))
    CSV.write(path, df)
end

# ─── Load groups.csv ──────────────────────────────────────────

groups_path = joinpath(bids_root, "groups.csv")
if !isfile(groups_path)
    @warn """
    No se encontró groups.csv en:
      $groups_path

    Crea el archivo con columnas: subject_id, group, session_id
    Ejemplo (ver groups.example.csv):
      subject_id,group,session_id
      M05,MS,T2
      M10,Control,T2
    """
    exit(1)
end

gdf = CSV.read(groups_path, DataFrame)
rename!(gdf, Dict(n => Symbol(lowercase(string(n))) for n in names(gdf)))

# Aliases tolerantes para nombres de columnas comunes en groups.csv
if !hasproperty(gdf, :session_id) && hasproperty(gdf, :session)
    rename!(gdf, :session => :session_id)
end
if !hasproperty(gdf, :session_id) && hasproperty(gdf, :ses)
    rename!(gdf, :ses => :session_id)
end

# Si existe `bids_id` (sub-M05 con cero a la izquierda), usarla como
# subject_id ya que los resultados del pipeline se guardan bajo sub-{bids_id}.
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

# Deduplicar filas (algunos groups.csv repiten sujeto por condición EC/EO,
# pero el análisis ya itera ambas internamente).
unique!(gdf, [:subject_id, :session_id, :group])

println("🧬 NeuroMIND — Análisis Transversal")
println("   Config: $CONFIG_P")
println("   Grupos: $groups_path  ($(nrow(gdf)) sujetos)")
println()

# ─── Proceso por condición ─────────────────────────────────────

for cond in ["EC", "EO"]
    println("── Condición: $cond " * "─"^40)
    out_dir = joinpath(res_root, "transversal", cond)
    mkpath(out_dir)

    # ── Cargar matrices por banda y grupo ──────────────────────
    ms_data   = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
    ctrl_data = Dict{String, Vector{Tuple{Vector{String}, Matrix{Float64}}}}()
    subject_means = NamedTuple[]   # para distribución per-subject
    inclusion     = NamedTuple[]

    for row in eachrow(gdf)
        sid = string(row.subject_id)
        ses = string(row.session_id)
        grp = string(row.group)
        is_ms = is_ms_group(grp)

        bands_ok = String[]
        for band in BANDS
            r = load_wpli(sid, ses, cond, band)
            r === nothing && continue
            (ch, W) = r
            target = is_ms ? ms_data : ctrl_data
            haskey(target, band) || (target[band] = [])
            push!(target[band], (ch, W))
            push!(bands_ok, band)

            # Per-subject mean connectivity (upper triangle)
            n  = length(ch)
            up = [(i, j) for i in 1:n for j in (i+1):n]
            isempty(up) && continue
            μ  = mean(W[i, j] for (i, j) in up)
            push!(subject_means, (
                subject_id = sid,
                group      = is_ms ? "MS" : "Control",
                band       = band,
                cond       = cond,
                mean_wpli  = round(μ, digits=5),
            ))
        end

        included = !isempty(bands_ok)
        push!(inclusion, (
            subject_id      = sid,
            session_id      = ses,
            group           = is_ms ? "MS" : "Control",
            n_bands_ok      = length(bands_ok),
            included        = included,
            excluded_reason = included ? "" : "Sin datos wPLI ($cond)",
        ))
    end

    n_ms   = count(r -> r.group == "MS"      && r.included, inclusion)
    n_ctrl = count(r -> r.group == "Control" && r.included, inclusion)
    n_excl = count(r -> !r.included, inclusion)
    println("  EM: $n_ms sujetos | Control: $n_ctrl sujetos | Excluidos: $n_excl")

    # ── Guardar subject_inclusion.csv ─────────────────────────
    CSV.write(joinpath(out_dir, "subject_inclusion.csv"), DataFrame(inclusion))
    !isempty(subject_means) &&
        CSV.write(joinpath(out_dir, "subject_band_means.csv"), DataFrame(subject_means))

    if n_ms < 1 || n_ctrl < 1
        println("  ⚠  Análisis omitido: se necesitan sujetos en ambos grupos")
        open(joinpath(out_dir, "transversal_summary.json"), "w") do io
            write(io, """{"n_ms":$n_ms,"n_ctrl":$n_ctrl,"n_included":$(n_ms+n_ctrl),"n_excluded":$n_excl,"n_total_sig":0,"n_bands":0,"cond":"$cond","timestamp":"$(Dates.format(now(),"yyyy-mm-ddTHH:MM:SS"))"}""")
        end
        continue
    end

    band_stats_rows = NamedTuple[]

    for band in BANDS
        ms_mats   = get(ms_data,   band, Tuple{Vector{String}, Matrix{Float64}}[])
        ctrl_mats = get(ctrl_data, band, Tuple{Vector{String}, Matrix{Float64}}[])
        (isempty(ms_mats) || isempty(ctrl_mats)) && continue

        # ── Alinear al conjunto común de canales ───────────────
        all_ch_sets = [Set(ch) for (ch, _) in vcat(ms_mats, ctrl_mats)]
        common_set  = reduce(intersect, all_ch_sets)
        ch_ref      = ctrl_mats[1][1]     # referencia de orden
        common_ch   = filter(c -> c in common_set, ch_ref)
        n           = length(common_ch)
        n < 2 && continue

        get_W(ch, W) = realign_matrix(ch, W, common_ch)

        # ── Matrices medias ────────────────────────────────────
        Wms   = mean(get_W(ch, W) for (ch, W) in ms_mats)
        Wctrl = mean(get_W(ch, W) for (ch, W) in ctrl_mats)
        Wdiff = Wms .- Wctrl

        save_mat_csv(joinpath(out_dir, "group_connectivity_control_$(band).csv"), Wctrl, common_ch)
        save_mat_csv(joinpath(out_dir, "group_connectivity_ms_$(band).csv"),      Wms,   common_ch)
        save_mat_csv(joinpath(out_dir, "group_difference_$(band).csv"),            Wdiff, common_ch)

        # ── Estadística por par ────────────────────────────────
        upper_idx = [(i, j) for i in 1:n for j in (i+1):n]
        p_vec = zeros(length(upper_idx))
        d_vec = zeros(length(upper_idx))

        for (k, (i, j)) in enumerate(upper_idx)
            va = [get_W(ch, W)[i, j] for (ch, W) in ms_mats]
            vb = [get_W(ch, W)[i, j] for (ch, W) in ctrl_mats]
            p_vec[k], d_vec[k] = welch_t(va, vb)
        end
        q_vec = bh_qvalues(p_vec)

        stat_rows = [(
            ch_a      = common_ch[i],
            ch_b      = common_ch[j],
            ctrl_mean = round(Wctrl[i, j], digits=5),
            ms_mean   = round(Wms[i, j],   digits=5),
            diff      = round(Wdiff[i, j], digits=5),
            p_value   = round(p_vec[k],    digits=5),
            q_value   = round(q_vec[k],    digits=5),
            effect_d  = round(d_vec[k],    digits=4),
        ) for (k, (i, j)) in enumerate(upper_idx)]

        CSV.write(joinpath(out_dir, "group_statistics_$(band).csv"), DataFrame(stat_rows))

        sig_rows = filter(r -> r.q_value < 0.05, stat_rows)
        sig_df   = isempty(sig_rows) ?
            DataFrame(ch_a=String[], ch_b=String[], ctrl_mean=Float64[],
                      ms_mean=Float64[], diff=Float64[],
                      p_value=Float64[], q_value=Float64[], effect_d=Float64[]) :
            DataFrame(sig_rows)
        CSV.write(joinpath(out_dir, "significant_edges_$(band).csv"), sig_df)

        n_sig     = length(sig_rows)
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
            mean_p     = round(mean(p_vec), digits=5),
            mean_d     = round(mean_d,      digits=4),
        ))
    end

    !isempty(band_stats_rows) &&
        CSV.write(joinpath(out_dir, "band_statistics.csv"), DataFrame(band_stats_rows))

    # ── transversal_summary.json ───────────────────────────────
    n_total_sig = isempty(band_stats_rows) ? 0 : sum(r.n_sig for r in band_stats_rows)
    best_band   = isempty(band_stats_rows) ? "" :
                  band_stats_rows[argmax(abs(r.diff_mean) for r in band_stats_rows)].band

    open(joinpath(out_dir, "transversal_summary.json"), "w") do io
        pairs_json = join(
            ["  \"$(r.band)_n_sig\": $(r.n_sig), \"$(r.band)_ctrl\": $(r.ctrl_mean), \"$(r.band)_ms\": $(r.ms_mean)"
             for r in band_stats_rows],
            ",\n"
        )
        write(io, """{
  "n_ms": $n_ms,
  "n_ctrl": $n_ctrl,
  "n_included": $(n_ms + n_ctrl),
  "n_excluded": $n_excl,
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

println("✅ Análisis transversal completado.")
println("   Abre el dashboard → Fase 13 Evaluación Transversal")
