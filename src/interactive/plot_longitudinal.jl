#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor interactivo longitudinal (T1 → T2)
# ═══════════════════════════════════════════════════════════════
#
#  Explora results/longitudinal/{eyesclosed|eyesopen}/ tras
#  run_longitudinal_analysis.jl (CLI acepta EC|EO como atajo).
#  · Resumen de producción · cambios pareados + IC bootstrap · heatmaps T1/T2/Δ
#  · red FDR / top-|dz| · volcano · potencia Δ · métricas de red
#
#  Uso:
#    julia --project=. src/interactive/plot_longitudinal.jl
#    julia --project=. src/interactive/plot_longitudinal.jl EO
#    # → http://127.0.0.1:8780/
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/plot_longitudinal.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      25-07-2026
#  Modificado  28-07-2026
# ───────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "viewer_support.jl"))
using .ViewerSupport
using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics, TOML

const PROJ = normpath(joinpath(@__DIR__, "..", ".."))
const HOST = "127.0.0.1"
const PORT = 8780
const COMMON_JS = joinpath(@__DIR__, "viewer_common.js")
const EXPECTED_SCHEMA_VERSION = 2
const INCOMPATIBLE_MESSAGE =
    "Resultados incompatibles con el visor actual. Regenere el análisis longitudinal."

struct IncompatibleResultsError <: Exception
    detail::String
end
Base.showerror(io::IO, e::IncompatibleResultsError) =
    print(io, INCOMPATIBLE_MESSAGE, " ", e.detail)

# ── Tipos ──────────────────────────────────────────────────────

mutable struct CondStore
    cond::String
    dir::String
    bands::Vector{String}
    channels::Dict{String,Vector{String}}          # band → ch
    mat_t1::Dict{String,Matrix{Float64}}
    mat_t2::Dict{String,Matrix{Float64}}
    mat_diff::Dict{String,Matrix{Float64}}
    stats::Dict{String,DataFrame}                  # edge stats
    band_stats::DataFrame
    subject_means::DataFrame
    paired::DataFrame
    summary::Dict{String,Any}
    spectral::DataFrame
    net_global::DataFrame
    net_t1::Dict{String,DataFrame}
    net_t2::Dict{String,DataFrame}
    net_delta::Dict{String,DataFrame}
end

mutable struct LongViewer
    results_root::String
    stores::Dict{String,CondStore}   # "EC"|"EO"
    errors::Dict{String,String}
    default_cond::String
end

# ── Loaders ────────────────────────────────────────────────────

"""Carpeta en disco ("eyesclosed"/"eyesopen") para el código corto "EC"/"EO"."""
_cond_dir(cond::String)::String = cond == "EC" ? "eyesclosed" : "eyesopen"

function _require_columns(df::DataFrame, required::Vector{Symbol}, artifact::String)
    missing_cols = filter(c -> !hasproperty(df, c), required)
    isempty(missing_cols) || throw(IncompatibleResultsError(
        "$artifact: faltan campos obligatorios $(join(string.(missing_cols), ", "))."))
    nrow(df) > 0 || throw(IncompatibleResultsError("$artifact está vacío."))
    if :schema_version in required
        all(Int.(df.schema_version) .== EXPECTED_SCHEMA_VERSION) ||
            throw(IncompatibleResultsError("$artifact usa una versión de esquema no compatible."))
    end
    return df
end

function _required_csv(path::String, required::Vector{Symbol})::DataFrame
    isfile(path) || throw(IncompatibleResultsError("Falta el artefacto obligatorio $(basename(path))."))
    df = try
        CSV.read(path, DataFrame)
    catch e
        throw(IncompatibleResultsError("No se puede leer $(basename(path)): $(sprint(showerror, e))."))
    end
    return _require_columns(df, required, basename(path))
end

function _validate_provenance(df::DataFrame, contract::Dict{String,Any},
                              artifact::String; fdr_key::Union{Nothing,String}=nothing)
    expected_source = string(contract["statistics_source"])
    all(string.(df.statistics_source) .== expected_source) ||
        throw(IncompatibleResultsError("$artifact no coincide con statistics_source del contrato."))
    if fdr_key !== nothing
        expected_scope = string(contract[fdr_key])
        all(string.(df.fdr_scope) .== expected_scope) ||
            throw(IncompatibleResultsError("$artifact no coincide con $fdr_key del contrato."))
    end
    return df
end

function _required_matrix(path::String)
    result = try
        read_mat_csv(path)
    catch e
        throw(IncompatibleResultsError("No se puede leer $(basename(path)): $(sprint(showerror, e))."))
    end
    result === nothing &&
        throw(IncompatibleResultsError("Falta o es inválido el artefacto obligatorio $(basename(path))."))
    return result
end

function validate_condition_contract(dir::String)::Dict{String,Any}
    path = joinpath(dir, "statistics_contract.toml")
    isfile(path) || throw(IncompatibleResultsError("Falta statistics_contract.toml."))
    contract = try
        TOML.parsefile(path)
    catch e
        throw(IncompatibleResultsError("statistics_contract.toml no es legible: $(sprint(showerror, e))."))
    end
    required = [
        "schema_version", "statistics_source",
        "fdr_scope_edges", "fdr_scope_global", "fdr_scope_power",
        "bootstrap_method", "bootstrap_iterations", "bootstrap_seed",
        "quantile_method", "rrb_method", "effect_dz_method",
    ]
    missing_fields = filter(k -> !haskey(contract, k), required)
    isempty(missing_fields) || throw(IncompatibleResultsError(
        "statistics_contract.toml: faltan campos obligatorios $(join(missing_fields, ", "))."))
    Int(contract["schema_version"]) == EXPECTED_SCHEMA_VERSION ||
        throw(IncompatibleResultsError(
            "schema_version=$(contract["schema_version"]); se requiere $EXPECTED_SCHEMA_VERSION."))
    return contract
end

function load_condition(dir::String, cond::String)::CondStore
    contract = validate_condition_contract(dir)
    tab_dir = joinpath(dir, "tables")
    band_stats = _required_csv(joinpath(tab_dir, "band_statistics_longitudinal.csv"), [
        :schema_version, :statistics_source, :fdr_scope, :band, :n_channels,
        :n_edges, :fdr_family_size, :n_nominal, :n_sig, :t1_mean, :t2_mean,
        :top20_nominal_overlap, :diff_mean, :mean_abs_dz, :n_subjects,
    ])
    _validate_provenance(band_stats, contract, "band_statistics_longitudinal.csv";
                         fdr_key="fdr_scope_edges")
    bands = [b for b in BAND_ORDER if b in String.(band_stats.band)]
    bands == BAND_ORDER || throw(IncompatibleResultsError(
        "band_statistics_longitudinal.csv debe contener exactamente las 7 bandas de producción."))
    channels = Dict{String,Vector{String}}()
    mat_t1 = Dict{String,Matrix{Float64}}()
    mat_t2 = Dict{String,Matrix{Float64}}()
    mat_diff = Dict{String,Matrix{Float64}}()
    stats = Dict{String,DataFrame}()
    net_t1 = Dict{String,DataFrame}()
    net_t2 = Dict{String,DataFrame}()
    net_delta = Dict{String,DataFrame}()

    for b in bands
        r1 = _required_matrix(joinpath(tab_dir, "longitudinal_connectivity_t1_$(b).csv"))
        r2 = _required_matrix(joinpath(tab_dir, "longitudinal_connectivity_t2_$(b).csv"))
        rd = _required_matrix(joinpath(tab_dir, "longitudinal_difference_$(b).csv"))
        r1[1] == r2[1] == rd[1] || throw(IncompatibleResultsError(
            "Las matrices T1, T2 y Δ de $b no comparten el mismo montaje."))
        channels[b] = r1[1]
        mat_t1[b] = r1[2]
        mat_t2[b] = r2[2]
        mat_diff[b] = rd[2]
        sdf = _required_csv(joinpath(tab_dir, "longitudinal_statistics_$(b).csv"), [
            :schema_version, :statistics_source, :fdr_scope, :fdr_family_size,
            :ch_a, :ch_b, :t1_mean, :t2_mean, :diff, :p_value, :q_value,
            :effect_dz, :effect_rank_abs_dz, :is_nominal, :is_fdr, :n, :p_method,
        ])
        _validate_provenance(sdf, contract, "longitudinal_statistics_$(b).csv";
                             fdr_key="fdr_scope_edges")
        all(Int.(sdf.fdr_family_size) .== nrow(sdf)) ||
            throw(IncompatibleResultsError(
                "longitudinal_statistics_$(b).csv declara un tamaño de familia FDR incoherente."))
        for c in (:ch_a, :ch_b)
            sdf[!, c] = String.(sdf[!, c])
        end
        stats[b] = sdf
        net_t1[b] = _required_csv(joinpath(tab_dir, "network_metrics_t1_$(b).csv"),
            [:schema_version, :channel, :strength, :degree, :norm_strength])
        net_t2[b] = _required_csv(joinpath(tab_dir, "network_metrics_t2_$(b).csv"),
            [:schema_version, :channel, :strength, :degree, :norm_strength])
        net_delta[b] = _required_csv(joinpath(tab_dir, "network_metrics_delta_$(b).csv"),
            [:schema_version, :channel, :delta_strength, :delta_degree, :delta_norm_strength])
    end

    subject_means = _required_csv(joinpath(tab_dir, "subject_band_means.csv"),
        [:schema_version, :statistics_source, :subject_id, :timepoint, :band, :cond,
         :mean_wpli, :n_channels])
    _validate_provenance(subject_means, contract, "subject_band_means.csv")
    paired = _required_csv(joinpath(tab_dir, "paired_subjects.csv"),
        [:schema_version, :statistics_source, :subject_id, :session_t1, :session_t2,
         :n_bands_ok, :included, :excluded_reason, :qc_t1, :qc_t2])
    _validate_provenance(paired, contract, "paired_subjects.csv")
    spectral = _required_csv(joinpath(tab_dir, "band_power_delta_statistics.csv"), [
        :schema_version, :statistics_source, :fdr_scope, :fdr_family_size,
        :channel, :band, :t1_mean, :t2_mean, :diff, :p_value, :q_value,
        :effect_dz, :n, :coverage_threshold_n, :coverage_pct, :coverage_low, :is_fdr,
    ])
    _validate_provenance(spectral, contract, "band_power_delta_statistics.csv";
                         fdr_key="fdr_scope_power")
    for b in bands
        sub = filter(r -> string(r.band) == b, spectral)
        nrow(sub) > 0 || throw(IncompatibleResultsError(
            "band_power_delta_statistics.csv no contiene la banda $b."))
        all(Int.(sub.fdr_family_size) .== nrow(sub)) ||
            throw(IncompatibleResultsError(
                "band_power_delta_statistics.csv declara una familia FDR incoherente para $b."))
    end
    net_global = _required_csv(joinpath(tab_dir, "network_global_statistics.csv"), [
        :schema_version, :statistics_source, :fdr_scope, :fdr_family_size,
        :band, :metric, :n_channels,
        :t1_mean, :t1_ci_low, :t1_ci_high, :t2_mean, :t2_ci_low, :t2_ci_high,
        :diff, :diff_ci_low, :diff_ci_high, :median_diff, :q1_diff, :q3_diff,
        :mean_wpli_t1, :mean_wpli_t2, :diff_wpli,
        :median_diff_wpli, :q1_diff_wpli, :q3_diff_wpli,
        :p_value, :q_value, :effect_dz, :effect_dz_ci_low, :effect_dz_ci_high,
        :effect_rrb, :n, :p_method, :is_fdr, :is_largest_abs_effect,
        :bootstrap_method, :bootstrap_iterations, :bootstrap_seed,
        :quantile_method, :rrb_method, :effect_dz_method,
    ])
    _validate_provenance(net_global, contract, "network_global_statistics.csv";
                         fdr_key="fdr_scope_global")
    nrow(net_global) == length(BAND_ORDER) || throw(IncompatibleResultsError(
        "network_global_statistics.csv debe contener una fila por cada una de las 7 bandas."))
    all(Int.(net_global.fdr_family_size) .== length(BAND_ORDER)) ||
        throw(IncompatibleResultsError(
            "network_global_statistics.csv debe declarar una familia FDR de 7 bandas."))
    summary = parse_summary_json(joinpath(dir, "longitudinal_summary.json"))
    Int(get(summary, "schema_version", 0)) == EXPECTED_SCHEMA_VERSION ||
        throw(IncompatibleResultsError("longitudinal_summary.json no declara schema_version compatible."))
    string(get(summary, "statistics_source", "")) == string(contract["statistics_source"]) ||
        throw(IncompatibleResultsError("La procedencia estadística no coincide entre contrato y resumen."))

    CondStore(
        cond, dir, bands, channels, mat_t1, mat_t2, mat_diff, stats,
        band_stats, subject_means, paired, summary, spectral, net_global,
        net_t1, net_t2, net_delta,
    )
end

function load_viewer(; results_root::Union{Nothing,String}=nothing,
                       default_cond::String="EC")::LongViewer
    root = resolve_results_root(PROJ, results_root)
    stores = Dict{String,CondStore}()
    errors = Dict{String,String}()
    for cond in ("EC", "EO")
        d = joinpath(root, "longitudinal", _cond_dir(cond))
        isdir(d) || continue
        try
            stores[cond] = load_condition(d, cond)
        catch e
            errors[cond] = sprint(showerror, e)
        end
    end
    if isempty(stores) && isempty(errors)
        errors["ALL"] = INCOMPATIBLE_MESSAGE *
            " No hay artefactos en longitudinal/{eyesclosed|eyesopen}."
    end
    dc = haskey(stores, default_cond) ? default_cond :
         (!isempty(stores) ? first(sort(collect(keys(stores)))) : default_cond)
    return LongViewer(root, stores, errors, dc)
end

# ── HTTP helpers ───────────────────────────────────────────────

function _send(sock, status::Int, body::String; content_type::String="text/html; charset=utf-8")
    bytes = codeunits(body)
    header = "HTTP/1.1 $status\r\nContent-Type: $content_type\r\n" *
             "Content-Length: $(length(bytes))\r\nConnection: close\r\n" *
             "Access-Control-Allow-Origin: *\r\n\r\n"
    write(sock, header)
    write(sock, bytes)
end

function _send_bytes(sock, status::Int, bytes::Vector{UInt8}; content_type::String)
    header = "HTTP/1.1 $status\r\nContent-Type: $content_type\r\n" *
             "Content-Length: $(length(bytes))\r\nConnection: close\r\n" *
             "Access-Control-Allow-Origin: *\r\n" *
             "Cache-Control: no-cache\r\n\r\n"
    write(sock, header)
    write(sock, bytes)
end

function _sort_band_stats(df::DataFrame)::DataFrame
    nrow(df) == 0 && return df
    hasproperty(df, :band) || return df
    order = Dict(b => i for (i, b) in enumerate(BAND_ORDER))
    return sort(df, :band; by = b -> get(order, string(b), 1000 + hash(string(b)) % 100))
end

function _honest_best_from_summary(summ::Dict{String,Any})::String
    return string(get(summ, "best_band", ""))
end

function state_json(viewer::LongViewer, cond::String, band::String;
                    edge_mode::String="top_dz", topn::Int=20)::String
    haskey(viewer.stores, cond) || return "{\"ok\":false,\"error\":\"condición desconocida\"}"
    st = viewer.stores[cond]
    band in st.bands || (band = first(st.bands))
    ch = st.channels[band]
    n = length(ch)
    ch_json = "[" * join(["\"$(json_escape(c))\"" for c in ch], ",") * "]"

    edf = get(st.stats, band, DataFrame())
    edges_out = DataFrame()
    if nrow(edf) > 0
        if edge_mode == "fdr"
            edges_out = filter(:is_fdr => identity, edf)
            edges_out = sort(edges_out, :effect_rank_abs_dz)
        elseif edge_mode == "p05"
            edges_out = filter(:is_nominal => identity, edf)
            edges_out = sort(edges_out, :effect_rank_abs_dz)
        else
            edges_out = filter(r -> Int(r.effect_rank_abs_dz) <= topn, edf)
            edges_out = sort(edges_out, :effect_rank_abs_dz)
        end
    end
    edge_cols = Symbol[:ch_a, :ch_b, :t1_mean, :t2_mean, :diff, :p_value, :q_value,
                       :effect_dz, :effect_rank_abs_dz, :is_nominal, :is_fdr, :n]
    edge_cols_full = edge_cols

    sm = st.subject_means
    sm_band = nrow(sm) > 0 && hasproperty(sm, :band) ?
              filter(r -> string(r.band) == band, sm) : DataFrame()
    sm_cols = Symbol[:subject_id, :timepoint, :band, :mean_wpli]

    sp = st.spectral
    sp_band = nrow(sp) > 0 && hasproperty(sp, :band) ?
              filter(r -> string(r.band) == band, sp) : DataFrame()
    sp_cols = Symbol[:channel, :band, :t1_mean, :t2_mean, :diff, :p_value, :effect_dz,
                     :q_value, :n, :fdr_family_size, :coverage_threshold_n,
                     :coverage_pct, :coverage_low, :is_fdr]

    nt1 = get(st.net_t1, band, DataFrame())
    nt2 = get(st.net_t2, band, DataFrame())
    nd  = get(st.net_delta, band, DataFrame())
    for df in (nt1, nt2, nd)
        nrow(df) > 0 && hasproperty(df, :channel) && (df.channel = String.(df.channel))
    end

    summ = st.summary
    band_stats_sorted = _sort_band_stats(st.band_stats)
    best_band = _honest_best_from_summary(summ)
    n_total_sig = try
        Int(round(Float64(get(summ, "n_total_sig", 0))))
    catch
        0
    end
    n_design = try
        Int(round(Float64(get(summ, "n_paired_design", 30))))
    catch
        30
    end
    n_cand = try
        Int(round(Float64(get(summ, "n_candidates", nrow(st.paired)))))
    catch
        nrow(st.paired)
    end
    n_paired = try
        Int(round(Float64(get(summ, "n_paired", 0))))
    catch
        0
    end
    n_excl_qc = try
        Int(round(Float64(get(summ, "n_excluded_qc", get(summ, "n_qc_loss", 0)))))
    catch
        0
    end
    n_excl_data = try
        Int(round(Float64(get(summ, "n_excluded_data", get(summ, "n_data_loss", 0)))))
    catch
        0
    end
    n_not_candidate = max(0, n_design - n_cand)

    n_warn = 0
    warn_ids = String[]
    if nrow(st.paired) > 0 && hasproperty(st.paired, :included)
        for r in eachrow(st.paired)
            incl = r.included === true || string(r.included) == "true"
            incl || continue
            qc1 = hasproperty(r, :qc_t1) ? string(something(r.qc_t1, "")) : ""
            qc2 = hasproperty(r, :qc_t2) ? string(something(r.qc_t2, "")) : ""
            if occursin("warning", lowercase(qc1)) || occursin("warning", lowercase(qc2))
                n_warn += 1
                hasproperty(r, :subject_id) && push!(warn_ids, string(r.subject_id))
            end
        end
    end
    warn_ids_json = "[" * join(["\"$(json_escape(x))\"" for x in warn_ids], ",") * "]"

    n_spec_ch = 0
    n_spec_min = 0
    n_spec_max = 0
    if nrow(sp_band) > 0 && hasproperty(sp_band, :channel)
        n_spec_ch = length(unique(String.(sp_band.channel)))
        if hasproperty(sp_band, :n)
            ns = Int.(coalesce.(sp_band.n, 0))
            n_spec_min = minimum(ns)
            n_spec_max = maximum(ns)
        end
    end

    kpi = """{
      "n_paired": $n_paired,
      "n_paired_design": $n_design,
      "n_candidates": $n_cand,
      "n_not_candidate": $n_not_candidate,
      "n_excluded": $(get(summ, "n_excluded", n_excl_qc + n_excl_data)),
      "n_excluded_qc": $n_excl_qc,
      "n_excluded_data": $n_excl_data,
      "n_include_with_warning": $n_warn,
      "warning_subjects": $warn_ids_json,
      "n_total_sig": $n_total_sig,
      "best_band": "$(json_escape(best_band))",
      "test": "$(json_escape(string(get(summ, "test", "wilcoxon"))))",
      "primary_global": "C",
      "wpli_method": "$(json_escape(string(get(summ, "wpli_method", ""))))",
      "timestamp": "$(json_escape(string(get(summ, "timestamp", ""))))",
      "cohort_note": "EC and EO analysed independently (parallel, no pooling). Candidates need T1∩T2 in at least one condition; N included can differ by condition.",
      "fdr_note": "Edge FDR-BH within each band over n_pairs edges (e.g. 351); not a joint BH over 7×n_pairs edge×band tests. Global network q-values: FDR-BH across the 7 bands.",
      "effect_note": "Primary global contrast = mean_strength on the intersected montage. mean_wPLI is its exact normalized expression: mean_strength=(n_channels-1)*mean_wPLI. effect_rrb is the Wilcoxon-aligned paired rank-biserial effect; effect_dz is complementary Cohen dz.",
      "spectral_note": "Absolute band power (μV²) from PSD integral; no log/dB; diff = T2−T1; FDR-BH within each band across channels; effect_dz = Cohen dz; n may vary per channel if pair incomplete.",
      "qc_warning_note": "include_with_warning subjects pass QC thresholds but carry amplitude/channel alerts; reasons in qc_decision_table.csv (not a separate warning_reason column here)."
    }"""

    bs_cols = Symbol[:band, :n_channels, :n_pairs, :fdr_family_size, :n_nominal,
                     :n_sig, :top20_nominal_overlap, :pct_sig, :t1_mean, :t2_mean,
                     :diff_mean, :mean_p, :mean_abs_dz, :n_subjects]
    paired_use = filter(c -> hasproperty(st.paired, c),
        [:subject_id, :session_t1, :session_t2, :n_bands_ok, :included,
         :excluded_reason, :qc_t1, :qc_t2])
    netg_cols = Symbol[
        :band, :metric, :n_channels,
        :t1_mean, :t1_ci_low, :t1_ci_high, :t2_mean, :t2_ci_low, :t2_ci_high,
        :diff, :diff_ci_low, :diff_ci_high, :median_diff, :q1_diff, :q3_diff,
        :mean_wpli_t1, :mean_wpli_t1_ci_low, :mean_wpli_t1_ci_high,
        :mean_wpli_t2, :mean_wpli_t2_ci_low, :mean_wpli_t2_ci_high,
        :diff_wpli, :diff_wpli_ci_low, :diff_wpli_ci_high,
        :median_diff_wpli, :q1_diff_wpli, :q3_diff_wpli,
        :p_value, :q_value, :effect_dz, :effect_dz_ci_low, :effect_dz_ci_high,
        :effect_rrb, :n, :n_positive, :n_negative, :n_zero, :p_method,
        :fdr_family_size, :fdr_scope, :is_fdr, :is_largest_abs_effect,
        :bootstrap_method, :bootstrap_iterations, :bootstrap_seed,
        :quantile_method, :rrb_method, :effect_dz_method,
    ]

    bs_band = filter(r -> string(r.band) == band, band_stats_sorted)
    nrow(bs_band) == 1 || throw(IncompatibleResultsError(
        "Falta el resumen persistido de la banda $band."))
    bs_row = bs_band[1, :]
    n_sig_band = Int(bs_row.n_sig)
    n_p_uncorr = Int(bs_row.n_nominal)
    top_all = filter(r -> Int(r.effect_rank_abs_dz) <= topn, edf)
    n_top_nominal = topn == 20 ? Int(bs_row.top20_nominal_overlap) :
                    count(r -> Bool(r.is_nominal), eachrow(top_all))
    n_top_actual = nrow(top_all)
    M1 = st.mat_t1[band]
    M2 = st.mat_t2[band]
    shared_lim_t12 = max(maximum(abs, M1), maximum(abs, M2), 1e-12)

    sp_use = filter(c -> nrow(sp_band) == 0 || hasproperty(sp_band, c), sp_cols)
    ng_use = filter(c -> nrow(st.net_global) == 0 || hasproperty(st.net_global, c), netg_cols)
    nt1_use = filter(c -> nrow(nt1) == 0 || hasproperty(nt1, c), [:channel, :strength, :degree, :norm_strength])
    nt2_use = filter(c -> nrow(nt2) == 0 || hasproperty(nt2, c), [:channel, :strength, :degree, :norm_strength])
    nd_use  = filter(c -> nrow(nd) == 0 || hasproperty(nd, c),
                     [:channel, :delta_strength, :delta_degree, :delta_norm_strength])
    edge_use = filter(c -> nrow(edges_out) == 0 || hasproperty(edges_out, c), edge_cols)
    estat_use = filter(c -> nrow(edf) == 0 || hasproperty(edf, c), edge_cols_full)
    bs_use = filter(c -> nrow(band_stats_sorted) == 0 || hasproperty(band_stats_sorted, c), bs_cols)
    n_global_fdr = count(identity, Bool.(st.net_global.is_fdr))
    descriptive_band = ""
    if nrow(st.net_global) > 0 && hasproperty(st.net_global, :is_largest_abs_effect)
        idx = findfirst(==(true), Bool.(st.net_global.is_largest_abs_effect))
        idx === nothing || (descriptive_band = string(st.net_global.band[idx]))
    end

    return """{
  "ok": true,
  "cond": "$cond",
  "band": "$band",
  "bands": [$(join(["\"$b\"" for b in st.bands], ","))],
  "conditions": [$(join(["\"$c\"" for c in sort(collect(keys(viewer.stores)))], ","))],
  "n": $n,
  "n_channels": $n,
  "channels": $ch_json,
  "positions": $(pos_json()),
  "matrix_t1": $(mat_flat(M1)),
  "matrix_t2": $(mat_flat(M2)),
  "matrix_diff": $(mat_flat(st.mat_diff[band])),
  "shared_lim_t12": $(repr(shared_lim_t12)),
  "n_sig_band": $n_sig_band,
  "n_p_uncorr": $n_p_uncorr,
  "edge_counts": {
    "total": $(Int(bs_row.n_edges)),
    "nominal": $n_p_uncorr,
    "fdr": $n_sig_band,
    "top_n": $n_top_actual,
    "top_nominal_overlap": $n_top_nominal,
    "fdr_family_size": $(Int(bs_row.fdr_family_size))
  },
  "global_counts": {
    "total_bands": $(nrow(st.net_global)),
    "fdr_bands": $n_global_fdr,
    "fdr_family_size": $(Int(first(st.net_global.fdr_family_size))),
    "largest_descriptive_band": "$(json_escape(descriptive_band))"
  },
  "n_spectral_channels": $n_spec_ch,
  "n_spectral_min": $n_spec_min,
  "n_spectral_max": $n_spec_max,
  "edge_mode": "$(json_escape(edge_mode))",
  "edges": $(df_rows_json(edges_out, edge_use)),
  "edge_stats": $(df_rows_json(edf, estat_use)),
  "kpi": $kpi,
  "band_stats": $(df_rows_json(band_stats_sorted, bs_use)),
  "subject_means": $(df_rows_json(sm_band, sm_cols)),
  "paired": $(df_rows_json(st.paired, paired_use)),
  "spectral": $(df_rows_json(sp_band, sp_use)),
  "net_global": $(df_rows_json(st.net_global, ng_use)),
  "net_t1": $(df_rows_json(nt1, nt1_use)),
  "net_t2": $(df_rows_json(nt2, nt2_use)),
  "net_delta": $(df_rows_json(nd, nd_use))
}"""
end

# ── PNG / CSV export ───────────────────────────────────────────

function save_diff_heatmap(st::CondStore, band::String)::String
    M = st.mat_diff[band]
    ch = st.channels[band]
    out = joinpath(st.dir, "figures", "viewer_delta_$(band).png")
    save_heatmap(out, M, ch, "Δ wPLI (T2−T1) — $band ($(st.cond))";
                 diverging=true, colorbar_label="Δ wPLI")
    return out
end

function save_edges_csv(st::CondStore, band::String, edge_mode::String, topn::Int)::String
    edf = get(st.stats, band, DataFrame())
    nrow(edf) == 0 && error("Sin estadísticas para $band")
    if edge_mode == "fdr"
        out_df = sort(filter(:is_fdr => identity, edf), :effect_rank_abs_dz)
    elseif edge_mode == "p05"
        out_df = sort(filter(:is_nominal => identity, edf), :effect_rank_abs_dz)
    else
        out_df = sort(filter(r -> Int(r.effect_rank_abs_dz) <= topn, edf),
                      :effect_rank_abs_dz)
    end
    out = joinpath(st.dir, "figures", "viewer_edges_$(edge_mode)_$(band).csv")
    mkpath(dirname(out))
    CSV.write(out, out_df)
    return out
end

# ── HTML SPA ───────────────────────────────────────────────────

function html_page()::String
    # JS uses string concat (+) — no template literals — to avoid Julia $ escaping.
    return """<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>NeuroMIND — Longitudinal T1→T2</title>
<style>
:root {
  --bg: #f1f5f9; --card: #fff; --ink: #0f172a; --muted: #64748b;
  --line: #e2e8f0; --accent: #0f766e; --ms: #b91c1c; --t1: #1d4ed8; --t2: #c2410c;
  --ok: #15803d; --warn: #a16207; --fdr: #b91c1c; --explore: #c2410c; --soft: #f8fafc;
  --shadow: 0 8px 28px rgba(15, 23, 42, 0.06);
}
* { box-sizing: border-box; }
body { margin:0; font-family: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
  background: var(--bg); color: var(--ink); }
header { background: linear-gradient(120deg, #0f172a 0%, #134e4a 60%, #0f766e 100%);
  color:#fff; padding: 1.15rem max(1.25rem, calc((100vw - 1560px) / 2)); }
header h1 { margin:0; font-size:1.35rem; font-weight:650; letter-spacing:-0.025em; }
header p { margin:0.3rem 0 0; opacity:0.82; font-size:0.86rem; }
.toolbar { display:flex; flex-wrap:wrap; gap:0.6rem; align-items:center;
  padding:0.7rem max(1.25rem, calc((100vw - 1560px) / 2)); background:var(--card);
  border-bottom:1px solid var(--line); position:sticky; top:0; z-index:20;
  box-shadow:0 2px 12px rgba(15,23,42,0.05); }
.toolbar label { font-size:0.72rem; color:var(--muted); margin-right:0.2rem;
  text-transform:uppercase; letter-spacing:0.04em; font-weight:600; }
.toolbar select, .toolbar button, .seg button, .analysis-controls select,
.analysis-controls button, .analysis-controls input {
  font: inherit; font-size:0.85rem; padding:0.35rem 0.65rem; border-radius:6px;
  border:1px solid var(--line); background:#fff; cursor:pointer; }
.seg { display:inline-flex; border:1px solid var(--line); border-radius:8px; overflow:hidden; }
.seg button { border:0; border-radius:0; background:#f8fafc; }
.seg button.active { background:var(--accent); color:#fff; }
.badge { font-size:0.75rem; padding:0.25rem 0.55rem; border-radius:999px;
  background:#ecfdf5; color:#065f46; border:1px solid #a7f3d0; }
.badge.mode-fdr { background:#fef2f2; color:#991b1b; border-color:#fecaca; }
.badge.mode-explore { background:#fff7ed; color:#9a3412; border-color:#fed7aa; }
.badge.mode-nom { background:#f8fafc; color:#475569; border-color:#cbd5e1; }
#topn-wrap { display:none; align-items:center; gap:0.25rem; }
.tabs { display:flex; gap:0.35rem; padding:0.85rem max(1.25rem, calc((100vw - 1560px) / 2)) 0;
  overflow-x:auto; }
.tabs button { border:1px solid transparent; background:transparent; padding:0.5rem 0.85rem;
  font:inherit; font-size:0.85rem; color:var(--muted); cursor:pointer; border-radius:8px 8px 0 0; }
.tabs button.active { background:var(--card); color:var(--ink); border-color:var(--line);
  border-bottom-color:var(--card); font-weight:600; }
main { max-width:1600px; margin:0 auto; padding:0 1.25rem 2.5rem; }
.panel { display:none; background:var(--card); border:1px solid var(--line);
  border-radius:0 12px 12px 12px; padding:1.1rem; min-height:400px; box-shadow:var(--shadow); }
.panel.active { display:block; }
.section-head { display:flex; align-items:flex-start; justify-content:space-between; gap:1rem;
  margin:0 0 1rem; }
.section-head h2 { margin:0; font-size:1.05rem; letter-spacing:-0.01em; }
.section-head p { margin:0.25rem 0 0; color:var(--muted); font-size:0.84rem; }
.grid { display:grid; gap:1rem; }
.grid-3 { grid-template-columns: repeat(3, 1fr); }
.grid-2 { grid-template-columns: 1fr 1fr; }
@media (max-width: 1100px) { .grid-3, .grid-2 { grid-template-columns: 1fr; } }
.card { border:1px solid var(--line); border-radius:10px; padding:0.85rem 0.95rem; background:var(--soft); }
.card h3 { margin:0 0 0.5rem; font-size:0.92rem; }
.result-hero { display:grid; grid-template-columns:minmax(260px,1.25fr) minmax(260px,0.75fr);
  gap:1rem; border:1px solid #a7f3d0; background:linear-gradient(135deg,#f0fdfa,#fff);
  border-radius:12px; padding:1rem 1.1rem; margin-bottom:1rem; }
.result-hero.no-sig { border-color:#cbd5e1; background:linear-gradient(135deg,#f8fafc,#fff); }
.result-hero .eyebrow { color:var(--accent); font-size:0.7rem; font-weight:700;
  letter-spacing:0.07em; text-transform:uppercase; }
.result-hero h2 { margin:0.2rem 0 0.35rem; font-size:1.25rem; letter-spacing:-0.02em; }
.result-hero p { margin:0; color:#475569; font-size:0.88rem; line-height:1.5; }
.result-facts { display:grid; grid-template-columns:1fr 1fr; gap:0.55rem; align-content:center; }
.result-fact { background:rgba(255,255,255,0.78); border:1px solid var(--line);
  border-radius:8px; padding:0.55rem 0.65rem; }
.result-fact .v { font-weight:700; font-size:0.95rem; }
.result-fact .l { color:var(--muted); font-size:0.68rem; text-transform:uppercase;
  letter-spacing:0.04em; margin-top:0.12rem; }
.analysis-controls { display:flex; flex-wrap:wrap; align-items:end; gap:0.65rem;
  padding:0.75rem; margin-bottom:1rem; background:#fff; border:1px solid var(--line);
  border-radius:9px; }
.analysis-controls label { display:flex; flex-direction:column; gap:0.25rem; color:var(--muted);
  font-size:0.7rem; font-weight:650; text-transform:uppercase; letter-spacing:0.04em; }
.analysis-controls .spacer { flex:1; }
.analysis-controls button { color:#334155; }
.analysis-controls button.primary { background:var(--accent); border-color:var(--accent); color:#fff; }
.kpi-row { display:flex; flex-wrap:wrap; gap:0.65rem; margin-bottom:1rem; }
.kpi { flex:1; min-width:130px; padding:0.65rem 0.8rem; border-radius:10px; background:#fff; border:1px solid var(--line); }
.kpi .v { font-size:1.25rem; font-weight:700; }
.kpi .l { font-size:0.68rem; color:var(--muted); text-transform:uppercase; letter-spacing:0.04em; }
.flow { display:flex; flex-wrap:wrap; gap:0.4rem; align-items:stretch; margin-bottom:1rem; }
.flow-step { flex:1; min-width:130px; background:#fff; border:1px solid var(--line);
  border-radius:10px; padding:0.65rem 0.75rem; }
.flow-step .step-l { font-size:0.68rem; text-transform:uppercase; letter-spacing:0.05em;
  color:var(--muted); font-weight:650; }
.flow-step .step-v { font-size:1.05rem; font-weight:700; margin-top:0.15rem; }
.flow-step .step-s { font-size:0.72rem; color:var(--muted); margin-top:0.12rem; }
.flow-arrow { display:flex; align-items:center; color:#94a3b8; font-size:1.05rem; }
.flow-step.ok { border-color:#86efac; background:linear-gradient(180deg,#fff,#f0fdf4); }
.flow-step.warn { border-color:#fcd34d; background:linear-gradient(180deg,#fff,#fffbeb); }
.canvas-wrap { width:100%; }
.canvas-wrap canvas { display:block; background:#fff; border-radius:8px;
  border:1px solid var(--line); cursor:crosshair; }
table { width:100%; border-collapse:collapse; font-size:0.78rem; }
th, td { padding:0.35rem 0.45rem; border-bottom:1px solid var(--line); text-align:left; }
th { color:var(--muted); font-weight:600; position:sticky; top:0; background:#fff; }
.result-table tbody tr.selected { background:#ecfdf5; }
.result-table td:first-child { font-weight:650; }
.sig-chip { display:inline-block; min-width:3.2rem; text-align:center; padding:0.12rem 0.35rem;
  border-radius:999px; font-size:0.68rem; font-weight:650; }
.sig-chip.yes { color:#991b1b; background:#fee2e2; }
.sig-chip.no { color:#475569; background:#e2e8f0; }
.scroll { max-height:320px; overflow:auto; }
.muted { color:var(--muted); font-size:0.85rem; }
.tip { position:fixed; pointer-events:none; background:#0f172a; color:#fff; font-size:0.75rem;
  padding:0.4rem 0.55rem; border-radius:6px; z-index:50; display:none; max-width:320px; white-space:pre-line; }
.status { font-size:0.8rem; color:var(--muted); margin-left:auto; }
.status.ok { color:var(--ok); } .status.err { color:var(--ms); }
.note { font-size:0.85rem; line-height:1.45; color:#334155; background:#f8fafc;
  border:1px solid #cbd5e1; border-radius:8px; padding:0.75rem; margin-top:0.75rem; }
.banner { display:none; margin:0 0 1rem; padding:0.65rem 0.9rem; border-radius:8px;
  background:#fffbeb; border:1px solid #fcd34d; color:#92400e; font-size:0.85rem; }
.banner.show { display:block; }
.banner.explore { background:#fff7ed; border-color:#fdba74; color:#9a3412; }
.empty-cta { text-align:center; padding:2rem 1rem; color:var(--muted); }
.compat-error { display:none; max-width:900px; margin:2rem auto; padding:1.2rem 1.35rem;
  background:#fff7ed; border:1px solid #fdba74; border-radius:12px; color:#7c2d12; }
.compat-error h2 { margin:0 0 0.45rem; font-size:1.1rem; }
.compat-error p { margin:0.25rem 0; line-height:1.5; }
.empty-cta button { margin-top:0.75rem; font:inherit; padding:0.45rem 0.9rem; border-radius:6px;
  border:1px solid var(--line); background:var(--accent); color:#fff; cursor:pointer; }
.legend-row { display:flex; flex-wrap:wrap; gap:0.75rem; margin-top:0.5rem; font-size:0.78rem; color:#475569; }
.legend-row span { display:inline-flex; align-items:center; gap:0.35rem; }
.swatch { width:18px; height:3px; border-radius:2px; display:inline-block; }
.swatch.fdr { background:var(--fdr); height:4px; }
.swatch.explore { background:var(--explore); height:3px; opacity:0.75;
  background-image: repeating-linear-gradient(90deg, var(--explore) 0 4px, transparent 4px 7px); }
.swatch.up { background:#b91c1c; }
.swatch.down { background:#1d4ed8; }
.pill { display:inline-block; font-size:0.7rem; padding:0.15rem 0.45rem; border-radius:999px;
  background:#ecfeff; color:#0f766e; border:1px solid #99f6e4; font-weight:600; margin-left:0.35rem; }
.net-title { display:flex; align-items:baseline; gap:0.5rem; flex-wrap:wrap; margin-bottom:0.75rem; }
.net-title h3 { margin:0; }
details { margin-top:1rem; border:1px solid var(--line); border-radius:10px; background:var(--soft); }
summary { cursor:pointer; padding:0.75rem 0.9rem; font-weight:650; font-size:0.88rem; }
details > .details-body { padding:0 0.9rem 0.9rem; }
@media (max-width: 760px) {
  .result-hero { grid-template-columns:1fr; }
  .toolbar .status { flex-basis:100%; margin-left:0; }
  .flow-arrow { display:none; }
  .analysis-controls .spacer { display:none; }
  .analysis-controls button { flex:1; }
}
</style>
<script src="/static/viewer_common.js"></script>
</head>
<body>
<header>
  <h1>NeuroMIND — Evaluación longitudinal T1 → T2</h1>
  <p>EM · EC/EO por separado · Wilcoxon + FDR-BH · Δ = T2−T1</p>
</header>
<div class="compat-error" id="compat-error">
  <h2>Resultados incompatibles con el visor actual</h2>
  <p id="compat-message">Regenera el análisis longitudinal con la versión de producción correspondiente.</p>
  <p><code>julia --project=. scripts/run_longitudinal_analysis.jl</code></p>
</div>
<div class="toolbar">
  <div class="seg" id="cond-seg">
    <button type="button" data-cond="EC" class="active" title="Ojos cerrados">EC · ojos cerrados</button>
    <button type="button" data-cond="EO" title="Ojos abiertos">EO · ojos abiertos</button>
  </div>
  <label>Banda <select id="band"></select></label>
  <span class="badge" id="badge-n">—</span>
  <span class="status" id="status">Cargando…</span>
</div>
<div class="tabs" id="tabs">
  <button type="button" data-tab="overview" class="active">Resumen</button>
  <button type="button" data-tab="global">Cambio T1 → T2</button>
  <button type="button" data-tab="conn">Conectividad detallada</button>
  <button type="button" data-tab="spec">Potencia y nodos</button>
</div>
<main>
  <section class="panel active" id="panel-overview">
    <div id="result-summary" class="result-hero no-sig"></div>
    <div class="flow" id="cohort-flow"></div>
    <div class="kpi-row" id="kpis"></div>
    <div class="grid grid-2">
      <div class="card">
        <h3>Resultado por banda <span class="pill">contraste primario: mean strength</span></h3>
        <div class="scroll"><table id="tbl-band-summary" class="result-table"><thead></thead><tbody></tbody></table></div>
      </div>
      <div class="card">
        <h3>Δ mean wPLI por banda <span class="pill" title="Normalización exacta de mean strength por el número de conexiones posibles de cada nodo.">expresión normalizada</span></h3>
        <div class="canvas-wrap"><canvas id="cv-diffbar" data-h="260"></canvas></div>
        <p class="muted" style="margin-bottom:0">Rojo = aumento en T2; azul = disminución.</p>
      </div>
    </div>
    <div style="display:none"><canvas id="cv-nsig" data-h="220"></canvas></div>
    <div class="note" id="interp"></div>
    <details>
      <summary>Ver inclusión y control de calidad de los pares</summary>
      <div class="details-body"><div class="scroll"><table id="tbl-paired"><thead></thead><tbody></tbody></table></div></div>
    </details>
  </section>
  <section class="panel" id="panel-global">
    <div class="section-head"><div><h2>Evolución longitudinal</h2>
      <p>Primero el cambio individual; después el contraste global sobre el montaje común.</p></div></div>
    <div class="grid grid-2">
      <div class="card"><h3>Mean wPLI por participante <span class="pill">montaje común</span></h3>
        <p class="muted" style="margin:0 0 0.45rem">Cada línea une T1 y T2 para el mismo participante; usa el mismo montaje común del contraste primario.</p>
        <div class="canvas-wrap"><canvas id="cv-spaghetti" data-h="300"></canvas></div></div>
      <div class="card"><h3>Distribución del cambio individual Δᵢ</h3>
        <p class="muted" style="margin:0 0 0.45rem">Δᵢ = T2−T1; puntos individuales y resumen descriptivo.</p>
        <div class="canvas-wrap"><canvas id="cv-delta" data-h="220"></canvas></div>
        <p class="muted" id="delta-stats"></p></div>
    </div>
    <div class="grid grid-2" style="margin-top:1rem">
      <div class="card"><h3>Δ mean wPLI en todas las bandas <span class="pill">mean strength / (canales−1)</span></h3>
      <div class="canvas-wrap"><canvas id="cv-bandmeans" data-h="260"></canvas></div>
      <p class="muted">Rojo ↑ T2, azul ↓ T2. Es la escala normalizada del contraste primario, no una evidencia independiente.</p></div>
      <div class="card"><h3>Contraste global de red <span class="pill">mean strength · primario</span></h3>
        <p class="muted" style="margin:0 0 0.4rem">Wilcoxon sobre mean strength, con FDR-BH entre las siete bandas.</p>
        <div class="scroll"><table id="tbl-netg"><thead></thead><tbody></tbody></table></div></div>
    </div>
  </section>
  <section class="panel" id="panel-conn">
    <div class="section-head"><div><h2>Conectividad por pares de canales</h2>
      <p>Las matrices muestran todo el montaje; el filtro solo afecta al grafo y a su tabla.</p></div></div>
    <div class="analysis-controls">
      <label>Aristas del grafo
        <select id="edge-mode">
          <option value="top_dz" selected>Top-N por |dz| · exploratorio</option>
          <option value="p05">p &lt; 0.05 · nominal</option>
          <option value="fdr">q &lt; 0.05 · FDR</option>
        </select>
      </label>
      <span id="topn-wrap"><label>N aristas
        <input id="topn" type="number" min="5" max="100" value="20" style="width:5rem"/>
      </label></span>
      <span class="badge" id="badge-mode">—</span>
      <span class="spacer"></span>
      <button type="button" id="btn-png">Guardar matriz Δ (PNG)</button>
      <button type="button" id="btn-csv" class="primary">Guardar aristas (CSV)</button>
    </div>
    <div id="fdr-banner" class="banner"></div>
    <div class="grid grid-3">
      <div class="card"><h3>T1</h3><div class="canvas-wrap"><canvas id="cv-hm-t1"></canvas></div></div>
      <div class="card"><h3>T2</h3><div class="canvas-wrap"><canvas id="cv-hm-t2"></canvas></div></div>
      <div class="card"><h3 id="hm-d-title">Δ (T2−T1) · FDR outline · triángulo inferior</h3>
        <div class="canvas-wrap"><canvas id="cv-hm-d"></canvas></div></div>
    </div>
    <p class="muted" id="hm-tip-hint">Hover: T1, T2, Δ, p, q y dz. Contorno negro = q&lt;0.05.</p>
    <div class="net-title" style="margin-top:1.25rem">
      <h3 id="net-heading">Red exploratoria</h3>
      <span class="badge" id="net-mode-badge">—</span>
    </div>
    <div id="net-empty" class="empty-cta" style="display:none">
      <p>No se detectaron aristas que superaran FDR en esta banda.</p>
      <p class="muted">Modo confirmatorio vacío — cambia a Top-N |dz| o a p&lt;0.05 nominal para exploración.</p>
      <button type="button" id="btn-to-topn">Usar Top-N |dz|</button>
    </div>
    <div id="net-content" class="grid grid-2">
      <div class="card"><h3>Grafo (edges filtrados)</h3>
        <div class="canvas-wrap"><canvas id="cv-graph"></canvas></div>
        <div class="legend-row">
          <span><i class="swatch up"></i> Δ&gt;0 (aumento T2)</span>
          <span><i class="swatch down"></i> Δ&lt;0 (disminución T2)</span>
          <span><i class="swatch fdr"></i> FDR q&lt;0.05</span>
          <span><i class="swatch explore"></i> No FDR / exploratorio</span>
          <span>Grosor ∝ |dz|</span>
        </div>
      </div>
      <div class="card"><h3>Tabla de edges</h3>
        <div class="scroll"><table id="tbl-edges"><thead></thead><tbody></tbody></table></div>
        <p class="muted" id="rows-edges" style="margin-bottom:0"></p>
      </div>
    </div>
    <div class="card" style="margin-top:1rem"><h3>Mapa de efectos — Cohen dz vs −log10(p) <span class="pill">color = estado FDR</span></h3>
      <div class="canvas-wrap"><canvas id="cv-volcano" data-h="320"></canvas></div>
      <p class="muted">Rojo = q&lt;0.05 · ámbar = p&lt;0.05 sin FDR · gris = no significativo.</p>
    </div>
  </section>
  <section class="panel" id="panel-spec">
    <div class="section-head"><div><h2>Potencia espectral y métricas nodales</h2>
      <p>Resultados complementarios de la banda seleccionada.</p></div></div>
    <div class="grid grid-2">
      <div class="card"><h3>Topo Δ band power (T2−T1)</h3>
        <div class="canvas-wrap"><canvas id="cv-topo"></canvas></div>
        <p class="muted" id="spec-ch-note"></p>
        <p class="muted">Potencia absoluta (μV²). Opacidad = cobertura; contorno grueso = q&lt;0.05.</p></div>
      <div class="card"><h3>Canal × banda</h3>
        <div class="scroll"><table id="tbl-spec"><thead></thead><tbody></tbody></table></div>
        <p class="muted" id="rows-spec" style="margin-bottom:0"></p>
      </div>
    </div>
    <div class="grid grid-2" style="margin-top:1rem">
      <div class="card"><h3>Strength nodal T1 vs T2 <span class="pill">descriptivo</span></h3>
        <p class="muted" style="margin:0 0 0.4rem">Cada punto es un canal; color = Δ strength.</p>
        <div class="canvas-wrap"><canvas id="cv-strength" data-h="360"></canvas></div></div>
      <div class="card"><h3>Cómo interpretar esta vista</h3>
        <p class="muted">El topograma resume el cambio de potencia absoluta por canal. El gráfico nodal compara la strength de cada electrodo entre T1 y T2.</p>
        <div class="note">Estas vistas son complementarias. La inferencia longitudinal primaria se encuentra en <b>Cambio T1 → T2</b>, usa mean strength y aplica FDR-BH entre bandas.</div></div>
    </div>
  </section>
</main>
<div class="tip" id="tip"></div>
<script>
(function () {
  if (!window.NMV) {
    document.getElementById('status').textContent = 'Error: NMV no cargado';
    document.getElementById('status').className = 'status err';
    return;
  }
  var gid = function (id) { return document.getElementById(id); };
  var S = null;
  var cond = 'EC';
  var tab = 'overview';
  var edgeLookup = null;

  function setStatus(msg, cls) {
    var el = gid('status');
    el.textContent = msg;
    el.className = 'status' + (cls ? ' ' + cls : '');
  }

  function edgeModeMeta(mode, topn) {
    if (mode === 'fdr') {
      return { title: 'Red FDR', short: 'FDR q<0.05', cls: 'mode-fdr', explore: false };
    }
    if (mode === 'p05') {
      return { title: 'Red nominal (p<0.05)', short: 'Nominal p<0.05', cls: 'mode-nom', explore: true };
    }
    var n = topn || 20;
    return { title: 'Red exploratoria Top-' + n + ' por |dz|', short: 'Top-' + n + ' |dz|', cls: 'mode-explore', explore: true };
  }

  function syncEdgeModeUI() {
    var mode = gid('edge-mode').value;
    var topn = Number(gid('topn').value) || 20;
    var meta = edgeModeMeta(mode, topn);
    gid('topn-wrap').style.display = (mode === 'top_dz') ? 'inline-flex' : 'none';
    gid('net-heading').textContent = meta.title;
    var bm = gid('badge-mode');
    bm.textContent = meta.short;
    bm.className = 'badge ' + meta.cls;
    var nmb = gid('net-mode-badge');
    nmb.textContent = meta.short;
    nmb.className = 'badge ' + meta.cls;
    return meta;
  }

  function updateBanner(meta) {
    var ban = gid('fdr-banner');
    var counts = (S && S.edge_counts) || {};
    var nSig = Number(counts.fdr) || 0;
    var nUnc = Number(counts.nominal) || 0;
    var nTotal = Number(counts.total) || 0;
    var nTop = Number(counts.top_n) || 0;
    var nOverlap = Number(counts.top_nominal_overlap) || 0;
    if (nSig === 0) {
      ban.className = 'banner show';
      ban.innerHTML = '<b>' + nUnc + '/' + nTotal + '</b> aristas con p&lt;0,05 sin corregir y ' +
        '<b>0/' + nTotal + '</b> con q&lt;0,05. De forma independiente, el grafo Top-N muestra ' +
        nTop + ' aristas por |dz|; ' + nOverlap + ' de ellas pertenecen al conjunto nominal.';
    } else if (meta && meta.explore) {
      ban.className = 'banner show explore';
      ban.textContent = 'Modo exploratorio: ' + meta.short + '. En esta banda hay ' +
        nSig + '/' + nTotal + ' aristas con q<0,05.';
    } else {
      ban.className = 'banner';
      ban.textContent = '';
    }
  }

  async function refresh() {
    setStatus('Cargando…');
    var meta = syncEdgeModeUI();
    var mode = gid('edge-mode').value;
    var topn = Number(gid('topn').value) || 20;
    var band = gid('band').value || 'ALPHA';
    var url = '/api/state?cond=' + encodeURIComponent(cond) +
      '&band=' + encodeURIComponent(band) +
      '&edge_mode=' + encodeURIComponent(mode) +
      '&topn=' + encodeURIComponent(String(topn));
    var res = await fetch(url);
    S = await res.json();
    if (!S.ok) { setStatus(S.error || 'Error', 'err'); return; }
    edgeLookup = NMV.edgeLookupFactory(S.edge_stats || []);
    var sel = gid('band');
    var cur = sel.value;
    sel.innerHTML = (S.bands || []).map(function (b) {
      return '<option value="' + b + '">' + b + '</option>';
    }).join('');
    if ((S.bands || []).indexOf(cur) >= 0) sel.value = cur;
    else if (S.band) sel.value = S.band;
    document.querySelectorAll('#cond-seg button').forEach(function (b) {
      b.classList.toggle('active', b.dataset.cond === cond);
      b.style.display = (S.conditions || []).indexOf(b.dataset.cond) >= 0 ? '' : 'none';
    });
    var kpi = S.kpi || {};
    var nCh = S.n_channels || S.n || '—';
    gid('badge-n').textContent = (kpi.n_paired != null ? kpi.n_paired : '—') +
      ' pares incluidos / ' + (kpi.n_candidates != null ? kpi.n_candidates : '—') +
      ' candidatos · ' + nCh + ' canales comunes';
    meta = syncEdgeModeUI();
    updateBanner(meta);
    renderAll();
    var ec = S.edge_counts || {};
    setStatus(cond + ' · ' + S.band + ' · ' + (ec.fdr || 0) +
      ' aristas con q<0,05 · ' + meta.short, 'ok');
  }

  function renderAll() {
    if (!S) return;
    renderResultSummary();
    renderFlow();
    renderKPIs();
    renderBars();
    renderBandSummary();
    renderInterp();
    renderPaired();
    renderSpaghetti();
    renderDeltaDist();
    renderBandMeans();
    renderHeatmaps();
    renderNet();
    renderVolcano();
    renderTopo();
    renderSpec();
    renderStrength();
    renderNetGlobal();
  }

  function formatNumber(v, digits) {
    var n = Number(v);
    return isFinite(n) ? n.toFixed(digits == null ? 3 : digits) : '—';
  }

  function renameHeaders(tableId, labels) {
    gid(tableId).querySelectorAll('thead th').forEach(function (th) {
      th.textContent = labels[th.textContent] || th.textContent;
    });
  }

  function globalRows() {
    return (S.net_global || []).filter(function (r) {
      return !r.metric || String(r.metric) === 'mean_strength';
    });
  }

  function selectedGlobal() {
    return globalRows().find(function (r) {
      return String(r.band) === String(S.band);
    }) || {};
  }

  function renderResultSummary() {
    var k = S.kpi || {};
    var rows = globalRows();
    var globalSigCount = Number((S.global_counts || {}).fdr_bands) || 0;
    var selected = selectedGlobal();
    var descriptiveBand = (S.global_counts || {}).largest_descriptive_band || '';
    var descriptive = rows.find(function (r) {
      return String(r.band) === String(descriptiveBand);
    }) || {};
    var direction = Number(selected.diff) > 0 ? 'aumento' :
      Number(selected.diff) < 0 ? 'disminución' : 'sin cambio';
    var primaryText = globalSigCount
      ? globalSigCount + ' banda' + (globalSigCount === 1 ? '' : 's') +
        ' supera' + (globalSigCount === 1 ? '' : 'n') + ' FDR en el contraste global.'
      : 'No se detectan cambios globales estadísticamente significativos tras FDR-BH.';
    var edgeText = (Number(k.n_total_sig) || 0) > 0
      ? k.n_total_sig + ' aristas superan FDR en el conjunto de bandas.'
      : 'Ninguna arista supera FDR en el conjunto de bandas.';
    var descriptiveText = descriptive.band
      ? ' <b>' + descriptive.band + '</b> presenta el mayor efecto descriptivo global' +
        ' (r_rb=' + formatNumber(descriptive.effect_rrb, 3) +
        '; q=' + formatNumber(descriptive.q_value, 3) + '), sin denominarlo tendencia significativa.'
      : '';
    var hero = gid('result-summary');
    hero.className = 'result-hero' + (globalSigCount ? '' : ' no-sig');
    hero.innerHTML =
      '<div><div class="eyebrow">Conclusión principal · ' + cond + '</div>' +
      '<h2>' + primaryText + '</h2>' +
      '<p>' + edgeText + ' La banda seleccionada, <b>' + S.band +
      '</b>, muestra ' + (direction === 'disminución' ? 'una disminución descriptiva' :
        direction === 'aumento' ? 'un aumento descriptivo' : 'ausencia de cambio descriptivo') +
      '.' + descriptiveText + ' Los modos Top-N y p&lt;0,05 se presentan únicamente como exploratorios.</p></div>' +
      '<div class="result-facts">' +
        '<div class="result-fact"><div class="v">' + formatNumber(selected.diff_wpli, 4) +
          '</div><div class="l">Δ mean wPLI normalizado</div></div>' +
        '<div class="result-fact"><div class="v">' + formatNumber(selected.effect_rrb, 3) +
          '</div><div class="l">r_rb · efecto primario</div></div>' +
        '<div class="result-fact"><div class="v">' + formatNumber(selected.q_value, 3) +
          '</div><div class="l">q global · ' + S.band + '</div></div>' +
        '<div class="result-fact"><div class="v">' + (k.n_paired != null ? k.n_paired : '—') +
          '</div><div class="l">Pares incluidos</div></div>' +
      '</div>';
  }

  function renderBandSummary() {
    var rows = NMV.sortBandsPhysio(globalRows());
    var table = gid('tbl-band-summary');
    table.querySelector('thead').innerHTML =
      '<tr><th>Banda</th><th>Δ wPLI</th><th>Δ strength</th><th>r_rb</th><th>dz</th><th>q global</th><th>FDR</th></tr>';
    table.querySelector('tbody').innerHTML = rows.map(function (r) {
      var sig = r.is_fdr === true || String(r.is_fdr) === 'true';
      var cls = String(r.band) === String(S.band) ? ' class="selected"' : '';
      return '<tr' + cls + '><td>' + r.band + '</td>' +
        '<td>' + formatNumber(r.diff_wpli, 4) + '</td>' +
        '<td>' + formatNumber(r.diff, 4) + '</td>' +
        '<td>' + formatNumber(r.effect_rrb, 3) + '</td>' +
        '<td>' + formatNumber(r.effect_dz, 3) + '</td>' +
        '<td>' + formatNumber(r.q_value, 3) + '</td>' +
        '<td><span class="sig-chip ' + (sig ? 'yes' : 'no') + '">' +
          (sig ? 'Sí' : 'No') + '</span></td></tr>';
    }).join('');
  }

  function renderFlow() {
    var k = S.kpi || {};
    var nNot = k.n_not_candidate != null ? k.n_not_candidate : 0;
    var nWarn = k.n_include_with_warning != null ? k.n_include_with_warning : 0;
    var steps = [
      { l: 'Diseño', v: k.n_paired_design != null ? k.n_paired_design : 30,
        s: 'N protocolo (Fig. 3.1)', cls: '' },
      { l: 'Candidatos', v: k.n_candidates != null ? k.n_candidates : '—',
        s: 'No en roster analítico: ' + nNot, cls: 'warn' },
      { l: 'Excl. QC', v: k.n_excluded_qc != null ? k.n_excluded_qc : '—',
        s: 'Datos: ' + (k.n_excluded_data != null ? k.n_excluded_data : 0), cls: 'warn' },
      { l: 'Incluidos', v: k.n_paired != null ? k.n_paired : '—',
        s: 'con warning: ' + nWarn + ' · ' + cond, cls: 'ok' }
    ];
    var html = '';
    steps.forEach(function (st, i) {
      if (i) html += '<div class="flow-arrow">→</div>';
      html += '<div class="flow-step ' + st.cls + '"><div class="step-l">' + st.l +
        '</div><div class="step-v">' + st.v + '</div><div class="step-s">' + st.s + '</div></div>';
    });
    html += '<p class="muted" style="flex-basis:100%;margin:0.35rem 0 0">' +
      'EC y EO se analizan por separado. Contraste primario: mean strength sobre montaje común.</p>';
    gid('cohort-flow').innerHTML = html;
  }

  function renderKPIs() {
    var k = S.kpi || {};
    var nGlobal = Number((S.global_counts || {}).fdr_bands) || 0;
    var items = [
      ['Pares incluidos', k.n_paired],
      ['Canales comunes', S.n_channels || S.n],
      ['Bandas FDR · global', nGlobal + ' / 7'],
      ['Aristas FDR · total', k.n_total_sig != null ? k.n_total_sig : '—']
    ];
    gid('kpis').innerHTML = items.map(function (it) {
      return '<div class="kpi"><div class="v">' + (it[1] != null ? it[1] : '—') +
        '</div><div class="l">' + it[0] + '</div></div>';
    }).join('');
  }

  function renderInterp() {
    var mode = gid('edge-mode').value;
    var kind = (mode === 'fdr') ? 'confirmatoria (q<0.05)' :
      (mode === 'p05') ? 'nominal (p<0.05)' : 'exploratoria (Top-N |dz|)';
    gid('interp').innerHTML =
      '<b>Lectura:</b> el contraste primario es la mean strength sobre el montaje común ' +
      '(Wilcoxon; FDR-BH entre 7 bandas). Mean wPLI es su normalización exacta: ' +
      'mean strength = (canales−1) × mean wPLI. Las aristas aplican BH dentro de cada banda. ' +
      'El filtro del grafo está en modo <b>' + kind + '</b>.';
  }

  function renderBars() {
    var rows = NMV.sortBandsPhysio(globalRows());
    var labs = rows.map(function (r) { return r.band; });
    NMV.barChart(gid('cv-diffbar'), labs, rows.map(function (r) { return Number(r.diff_wpli) || 0; }),
      function (v) { return v >= 0 ? '#c2410c' : '#1d4ed8'; }, { height: 260, signed: true });
  }

  function renderPaired() {
    NMV.fillTable(gid('tbl-paired'),
      ['subject_id', 'included', 'qc_t1', 'qc_t2', 'n_bands_ok', 'excluded_reason'],
      S.paired || [], {
        excluded_reason: function (v) {
          if (v == null || v === '' || v === 'missing') return '—';
          return v;
        },
        included: function (v) {
          return (v === true || v === 'true') ? 'sí' : 'no';
        },
        qc_t1: function (v) {
          return (v == null || v === '') ? '—' : v;
        },
        qc_t2: function (v) {
          return (v == null || v === '') ? '—' : v;
        }
      });
    var labels = {
      subject_id: 'Participante', included: 'Incluido', qc_t1: 'QC T1', qc_t2: 'QC T2',
      n_bands_ok: 'Bandas válidas', excluded_reason: 'Motivo de exclusión'
    };
    gid('tbl-paired').querySelectorAll('thead th').forEach(function (th) {
      th.textContent = labels[th.textContent] || th.textContent;
    });
  }

  function pairedDeltas() {
    var rows = S.subject_means || [];
    var by = {};
    rows.forEach(function (r) {
      var sid = String(r.subject_id);
      by[sid] = by[sid] || {};
      by[sid][String(r.timepoint)] = Number(r.mean_wpli);
    });
    return Object.keys(by).map(function (sid) {
      return { sid: sid, T1: by[sid].T1, T2: by[sid].T2,
        d: (by[sid].T2 != null && by[sid].T1 != null) ? by[sid].T2 - by[sid].T1 : NaN };
    }).filter(function (p) {
      return p.T1 != null && p.T2 != null && !isNaN(p.T1) && !isNaN(p.T2);
    });
  }

  function renderSpaghetti() {
    var cv = gid('cv-spaghetti');
    var rsz = NMV.resizeCanvas(cv, { height: 300, fallbackW: 700 });
    var ctx = rsz.ctx, w = rsz.w, h = rsz.h;
    var pairs = pairedDeltas();
    if (!pairs.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin subject_band_means', 20, 40);
      return;
    }
    var vals = [];
    pairs.forEach(function (p) { vals.push(p.T1, p.T2); });
    var gs = selectedGlobal();
    [
      gs.mean_wpli_t1_ci_low, gs.mean_wpli_t1_ci_high,
      gs.mean_wpli_t2_ci_low, gs.mean_wpli_t2_ci_high
    ].forEach(function (v) { if (v != null && isFinite(Number(v))) vals.push(Number(v)); });
    var ymin = Math.min.apply(null, vals), ymax = Math.max.apply(null, vals);
    var pad = { l: 62, r: 30, t: 28, b: 44 };
    var yscale = function (v) {
      return pad.t + (1 - (v - ymin) / Math.max(ymax - ymin, 1e-9)) * (h - pad.t - pad.b);
    };
    var x1 = w * 0.32, x2 = w * 0.72;
    // Y grid + numeric scale
    ctx.strokeStyle = '#f1f5f9';
    ctx.fillStyle = '#94a3b8';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'right';
    for (var ti = 0; ti <= 4; ti++) {
      var yv = ymin + (ymax - ymin) * ti / 4;
      var yy = yscale(yv);
      ctx.beginPath();
      ctx.moveTo(pad.l, yy);
      ctx.lineTo(w - pad.r, yy);
      ctx.stroke();
      ctx.fillText(yv.toFixed(3), pad.l - 6, yy + 3);
    }
    ctx.strokeStyle = '#e2e8f0';
    ctx.beginPath();
    ctx.moveTo(x1, pad.t); ctx.lineTo(x1, h - pad.b);
    ctx.moveTo(x2, pad.t); ctx.lineTo(x2, h - pad.b);
    ctx.stroke();
    ctx.fillStyle = '#64748b';
    ctx.font = '12px sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText('T1', x1, h - 14);
    ctx.fillText('T2', x2, h - 14);
    ctx.save();
    ctx.translate(14, (pad.t + h - pad.b) / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.font = '11px sans-serif';
    ctx.fillStyle = '#475569';
    ctx.fillText('Mean wPLI por participante', 0, 0);
    ctx.restore();
    pairs.forEach(function (p) {
      ctx.strokeStyle = 'rgba(100,116,139,0.35)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(x1, yscale(p.T1));
      ctx.lineTo(x2, yscale(p.T2));
      ctx.stroke();
      ctx.fillStyle = '#1d4ed8';
      ctx.beginPath(); ctx.arc(x1, yscale(p.T1), 3.5, 0, Math.PI * 2); ctx.fill();
      ctx.fillStyle = '#c2410c';
      ctx.beginPath(); ctx.arc(x2, yscale(p.T2), 3.5, 0, Math.PI * 2); ctx.fill();
    });
    var ms1 = { m: Number(gs.mean_wpli_t1), lo: Number(gs.mean_wpli_t1_ci_low),
      hi: Number(gs.mean_wpli_t1_ci_high) };
    var ms2 = { m: Number(gs.mean_wpli_t2), lo: Number(gs.mean_wpli_t2_ci_low),
      hi: Number(gs.mean_wpli_t2_ci_high) };
    ctx.strokeStyle = '#0f766e';
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    ctx.moveTo(x1, yscale(ms1.m));
    ctx.lineTo(x2, yscale(ms2.m));
    ctx.stroke();
    [[x1, ms1], [x2, ms2]].forEach(function (pair) {
      var x = pair[0], ms = pair[1];
      if (!isFinite(ms.lo) || !isFinite(ms.hi)) return;
      ctx.strokeStyle = '#0f766e';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(x, yscale(ms.lo));
      ctx.lineTo(x, yscale(ms.hi));
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(x - 6, yscale(ms.lo));
      ctx.lineTo(x + 6, yscale(ms.lo));
      ctx.moveTo(x - 6, yscale(ms.hi));
      ctx.lineTo(x + 6, yscale(ms.hi));
      ctx.stroke();
    });
    ctx.fillStyle = '#0f766e';
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText('Media e IC95% bootstrap de producción (N=' + pairs.length + ')', pad.l, 16);
  }

  function renderDeltaDist() {
    var cv = gid('cv-delta');
    var rsz = NMV.resizeCanvas(cv, { height: 220, fallbackW: 700 });
    var ctx = rsz.ctx, w = rsz.w, h = rsz.h;
    var pairs = pairedDeltas();
    var ds = pairs.map(function (p) { return p.d; }).filter(function (v) { return !isNaN(v); });
    if (!ds.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin Δᵢ', 20, 40);
      gid('delta-stats').textContent = '';
      return;
    }
    var gs = selectedGlobal();
    var med = Number(gs.median_diff_wpli);
    var q1 = Number(gs.q1_diff_wpli);
    var q3 = Number(gs.q3_diff_wpli);
    var meanD = Number(gs.diff_wpli);
    var ciLo = Number(gs.diff_wpli_ci_low);
    var ciHi = Number(gs.diff_wpli_ci_high);
    var absMax = Math.max.apply(null, ds.map(Math.abs).concat(
      [Math.abs(q1), Math.abs(q3), Math.abs(ciLo), Math.abs(ciHi), 1e-6]));
    var pad = { l: 50, r: 24, t: 20, b: 36 };
    var X = function (v) {
      return pad.l + (v + absMax) / (2 * absMax) * (w - pad.l - pad.r);
    };
    var y0 = h / 2;
    ctx.strokeStyle = '#e2e8f0';
    ctx.beginPath();
    ctx.moveTo(pad.l, y0); ctx.lineTo(w - pad.r, y0); ctx.stroke();
    ctx.strokeStyle = '#94a3b8';
    ctx.setLineDash([3, 3]);
    ctx.beginPath();
    ctx.moveTo(X(0), pad.t); ctx.lineTo(X(0), h - pad.b); ctx.stroke();
    ctx.setLineDash([]);
    // IQR box
    ctx.fillStyle = 'rgba(15,118,110,0.12)';
    ctx.fillRect(X(q1), y0 - 28, X(q3) - X(q1), 56);
    ctx.strokeStyle = '#0f766e';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(X(med), y0 - 28); ctx.lineTo(X(med), y0 + 28); ctx.stroke();
    ds.forEach(function (d, i) {
      var jitter = ((i % 7) - 3) * 3.5;
      ctx.beginPath();
      ctx.arc(X(d), y0 + jitter, 3.2, 0, Math.PI * 2);
      ctx.fillStyle = d >= 0 ? 'rgba(194,65,12,0.7)' : 'rgba(29,78,216,0.7)';
      ctx.fill();
    });
    // Media e IC95% bootstrap calculados por Julia.
    ctx.strokeStyle = '#0f172a';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(X(ciLo), y0 - 40);
    ctx.lineTo(X(ciHi), y0 - 40);
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(X(meanD), y0 - 40, 3.5, 0, Math.PI * 2);
    ctx.fillStyle = '#0f172a';
    ctx.fill();
    ctx.fillStyle = '#64748b';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    [-absMax, 0, absMax].forEach(function (xv) {
      ctx.fillText(xv.toFixed(3), X(xv), h - 12);
    });
    ctx.textAlign = 'left';
    ctx.fillText('Δᵢ = T2−T1', pad.l, 14);
    gid('delta-stats').textContent =
      'Producción Julia · N=' + (gs.n != null ? gs.n : ds.length) +
      ' · mediana=' + formatNumber(med, 4) +
      ' · IQR=[' + formatNumber(q1, 4) + ', ' + formatNumber(q3, 4) + ']' +
      ' · media=' + formatNumber(meanD, 4) +
      ' · IC95% bootstrap=[' + formatNumber(ciLo, 4) + ', ' + formatNumber(ciHi, 4) + ']' +
      ' · r_rb=' + formatNumber(gs.effect_rrb, 3) +
      ' · dz=' + formatNumber(gs.effect_dz, 3) +
      ' · p=' + formatNumber(gs.p_value, 4) + ' · q=' + formatNumber(gs.q_value, 4) +
      ' · signos +/−/0=' + (gs.n_positive || 0) + '/' + (gs.n_negative || 0) +
      '/' + (gs.n_zero || 0);
  }

  function renderBandMeans() {
    var rows = NMV.sortBandsPhysio(globalRows());
    // Neutral intensity (direction only) — no global FDR survived
    NMV.barChart(gid('cv-bandmeans'),
      rows.map(function (r) { return r.band; }),
      rows.map(function (r) { return Number(r.diff_wpli) || 0; }),
      function (v) { return v >= 0 ? 'rgba(194,65,12,0.55)' : 'rgba(29,78,216,0.55)'; },
      { height: 260, signed: true });
  }

  function renderHeatmaps() {
    var lim = S.shared_lim_t12 || null;
    var fdrSet = NMV.buildFdrSet(S.edge_stats || []);
    var base = { n: S.n, channels: S.channels, lim: lim };
    NMV.drawHeatmap(gid('cv-hm-t1'), S.matrix_t1, Object.assign({}, base, { diverging: false }));
    NMV.drawHeatmap(gid('cv-hm-t2'), S.matrix_t2, Object.assign({}, base, { diverging: false }));
    NMV.drawHeatmap(gid('cv-hm-d'), S.matrix_diff, {
      n: S.n, channels: S.channels, diverging: true, fdrSet: fdrSet, triangle: 'lower'
    });
    var title = gid('hm-d-title');
    var counts = S.edge_counts || {};
    if ((counts.fdr || 0) === 0) {
      title.textContent = 'Δ completo (T2−T1) · 0/' + (counts.total || 0) +
        ' aristas con q<0,05 · triángulo inferior';
    } else {
      title.textContent = 'Δ (T2−T1) · contorno FDR · triángulo inferior';
    }
  }

  function renderNet() {
    var mode = gid('edge-mode').value;
    var edges = S.edges || [];
    var empty = mode === 'fdr' && edges.length === 0;
    gid('net-empty').style.display = empty ? 'block' : 'none';
    gid('net-content').style.display = empty ? 'none' : '';
    if (empty) return;
    NMV.drawGraph(gid('cv-graph'), S, edges, {
      thicknessKey: 'effect_dz',
      legendDir: 'rojo/ámbar Δ>0 · azul/slate Δ<0',
      legendThick: 'grosor ∝ |dz|'
    });
    var fmt = {
      t1_mean: function (v) { return Number(v).toFixed(4); },
      t2_mean: function (v) { return Number(v).toFixed(4); },
      diff: function (v) { return Number(v).toFixed(4); },
      p_value: function (v) { return Number(v).toFixed(4); },
      q_value: function (v) { return Number(v).toFixed(4); },
      effect_dz: function (v) { return Number(v).toFixed(3); }
    };
    NMV.fillTable(gid('tbl-edges'),
      ['ch_a', 'ch_b', 't1_mean', 't2_mean', 'diff', 'p_value', 'q_value', 'effect_dz', 'n'],
      edges, fmt);
    renameHeaders('tbl-edges', {
      ch_a: 'Canal A', ch_b: 'Canal B', t1_mean: 'wPLI T1', t2_mean: 'wPLI T2',
      diff: 'Δ wPLI', p_value: 'p nominal', q_value: 'q FDR', effect_dz: 'dz', n: 'N'
    });
    gid('rows-edges').textContent = edges.length + ' filas mostradas · filtro: ' +
      edgeModeMeta(mode, Number(gid('topn').value) || 20).short;
  }

  function renderVolcano() {
    NMV.drawVolcano(gid('cv-volcano'), S.edge_stats || [], { yKey: 'p', xLabel: 'Cohen dz' });
  }

  function renderTopo() {
    NMV.drawTopo(gid('cv-topo'), S.spectral || [], S.positions || {}, { showN: true });
  }

  function renderSpec() {
    var fmt = {
      diff: function (v) { return Number(v).toFixed(4); },
      p_value: function (v) { return Number(v).toFixed(4); },
      q_value: function (v) { return Number(v).toFixed(4); },
      effect_dz: function (v) { return Number(v).toFixed(3); }
    };
    var rows = (S.spectral || []).slice().sort(function (a, b) {
      return Math.abs(Number(b.diff) || 0) - Math.abs(Number(a.diff) || 0);
    });
    NMV.fillTable(gid('tbl-spec'),
      ['channel', 'diff', 'p_value', 'q_value', 'effect_dz', 'n', 'coverage_pct'],
      rows.slice(0, 80), fmt);
    renameHeaders('tbl-spec', {
      channel: 'Canal', diff: 'Δ potencia (μV²)', p_value: 'p nominal',
      q_value: 'q FDR', effect_dz: 'dz', n: 'N', coverage_pct: 'Cobertura %'
    });
    gid('rows-spec').textContent = rows.length + ' canales mostrados para ' + S.band + '.';
    var note = gid('spec-ch-note');
    var nmin = S.n_spectral_min, nmax = S.n_spectral_max;
    var nch = S.n_spectral_channels || (rows.length ? rows.length : '—');
    var family = rows.length && rows[0].fdr_family_size != null ?
      rows[0].fdr_family_size : nch;
    var threshold = rows.length && rows[0].coverage_threshold_n != null ?
      rows[0].coverage_threshold_n : '—';
    var msg = 'Canales con potencia: ' + nch + ' (wPLI usa ' + (S.n_channels || S.n) +
      '). FDR-BH entre ' + family + ' canales dentro de ' + S.band +
      '; cobertura baja si N<' + threshold + '.';
    if (nmin != null && nmax != null && nmin !== nmax) {
      msg += ' n por canal varía ' + nmin + '–' + nmax +
        ' (pares incompletos canal-específicos; no todos los ' +
        ((S.kpi && S.kpi.n_paired) || '?') + ' pares aportan cada canal).';
    } else if (nmax) {
      msg += ' n=' + nmax + ' por canal.';
    }
    note.textContent = msg;
  }

  function renderStrength() {
    NMV.drawStrength(gid('cv-strength'), S.net_t1 || [], S.net_t2 || [], 'T1', 'T2');
  }

  function renderNetGlobal() {
    var fmt = {
      diff_wpli: function (v) { return Number(v).toFixed(4); },
      diff: function (v) { return Number(v).toFixed(4); },
      p_value: function (v) { return Number(v).toFixed(4); },
      q_value: function (v) { return Number(v).toFixed(4); },
      effect_dz: function (v) { return Number(v).toFixed(3); },
      effect_rrb: function (v) { return Number(v).toFixed(3); }
    };
    NMV.fillTable(gid('tbl-netg'),
      ['band', 'diff_wpli', 'diff', 'effect_rrb', 'effect_dz', 'p_value', 'q_value', 'n'],
      S.net_global || [], fmt);
    renameHeaders('tbl-netg', {
      band: 'Banda', diff_wpli: 'Δ mean wPLI', diff: 'Δ strength',
      effect_rrb: 'r_rb', effect_dz: 'dz', p_value: 'p Wilcoxon',
      q_value: 'q (BH · 7 bandas)', n: 'N'
    });
  }

  ['cv-hm-t1', 'cv-hm-t2', 'cv-hm-d'].forEach(function (id) {
    NMV.bindHeatmap(gid(id), function (a, b) {
      return edgeLookup ? edgeLookup(a, b) : null;
    });
  });

  document.querySelectorAll('#cond-seg button').forEach(function (b) {
    b.addEventListener('click', function () {
      cond = b.dataset.cond;
      refresh();
    });
  });
  gid('band').addEventListener('change', refresh);
  gid('edge-mode').addEventListener('change', function () { syncEdgeModeUI(); refresh(); });
  gid('topn').addEventListener('change', refresh);
  gid('btn-to-topn').addEventListener('click', function () {
    gid('edge-mode').value = 'top_dz';
    syncEdgeModeUI();
    refresh();
  });

  document.querySelectorAll('#tabs button').forEach(function (b) {
    b.addEventListener('click', function () {
      tab = b.dataset.tab;
      document.querySelectorAll('#tabs button').forEach(function (x) {
        x.classList.toggle('active', x === b);
      });
      document.querySelectorAll('.panel').forEach(function (p) {
        p.classList.toggle('active', p.id === 'panel-' + tab);
      });
      setTimeout(renderAll, 30);
    });
  });

  gid('btn-png').addEventListener('click', async function () {
    setStatus('Exportando PNG…');
    try {
      var r = await fetch('/api/export_png?cond=' + encodeURIComponent(cond) +
        '&band=' + encodeURIComponent(gid('band').value));
      var j = await r.json();
      if (!j.ok) throw new Error(j.error || 'fail');
      setStatus('PNG: ' + j.file, 'ok');
    } catch (e) {
      setStatus(String(e.message || e), 'err');
    }
  });
  gid('btn-csv').addEventListener('click', async function () {
    setStatus('Exportando CSV…');
    try {
      var mode = gid('edge-mode').value;
      var topn = gid('topn').value;
      var r = await fetch('/api/export_csv?cond=' + encodeURIComponent(cond) +
        '&band=' + encodeURIComponent(gid('band').value) +
        '&edge_mode=' + encodeURIComponent(mode) +
        '&topn=' + encodeURIComponent(topn));
      var j = await r.json();
      if (!j.ok) throw new Error(j.error || 'fail');
      setStatus('CSV: ' + j.file, 'ok');
    } catch (e) {
      setStatus(String(e.message || e), 'err');
    }
  });

  window.addEventListener('resize', function () {
    clearTimeout(window._rt);
    window._rt = setTimeout(renderAll, 150);
  });

  (async function init() {
    var meta = await (await fetch('/api/meta')).json();
    if (!meta.ok) {
      gid('compat-error').style.display = 'block';
      gid('compat-message').textContent = meta.error ||
        'Regenera el análisis longitudinal con la versión de producción correspondiente.';
      document.querySelector('.toolbar').style.display = 'none';
      gid('tabs').style.display = 'none';
      document.querySelector('main').style.display = 'none';
      return;
    }
    cond = meta.default_cond || 'EC';
    await refresh();
  })();
})();
</script>
</body>
</html>
"""
end

# ── HTTP ───────────────────────────────────────────────────────

function _parse_qs(q::AbstractString)::Dict{String,String}
    out = Dict{String,String}()
    isempty(q) && return out
    for part in split(q, '&')
        kv = split(part, '='; limit=2)
        length(kv) == 2 || continue
        key = String(kv[1])
        val = String(replace(kv[2], '+' => ' '))
        val = replace(val, r"%([0-9A-Fa-f]{2})" =>
            m -> string(Char(parse(Int, m.captures[1]; base=16))))
        out[key] = val
    end
    return out
end

function handle_request(sock, viewer::LongViewer)
    try
        line = readline(sock)
        isempty(line) && return
        m = match(r"^(GET|POST)\s+(\S+)", line)
        m === nothing && return _send(sock, 400, "Bad Request")
        method, pathfull = m.captures[1], m.captures[2]
        path, qs = let sp = split(pathfull, '?'; limit=2)
            length(sp) == 1 ? (sp[1], Dict{String,String}()) : (sp[1], _parse_qs(sp[2]))
        end
        while true
            h = readline(sock)
            (isempty(h) || h == "\r") && break
        end

        if path == "/" || path == "/index.html"
            return _send(sock, 200, html_page())
        elseif path == "/static/viewer_common.js"
            isfile(COMMON_JS) || return _send(sock, 404, "JS not found")
            return _send_bytes(sock, 200, read(COMMON_JS);
                               content_type="application/javascript; charset=utf-8")
        elseif path == "/api/meta"
            ok = !isempty(viewer.stores) && isempty(viewer.errors)
            errors_json = "{" * join([
                "\"$(json_escape(k))\":\"$(json_escape(v))\"" for (k, v) in
                sort(collect(viewer.errors); by=first)
            ], ",") * "}"
            first_error = isempty(viewer.errors) ? "" : first(values(viewer.errors))
            body = "{\"ok\":$(ok ? "true" : "false")," *
                   "\"default_cond\":\"$(viewer.default_cond)\"," *
                   "\"conditions\":[$(join(["\"$c\"" for c in sort(collect(keys(viewer.stores)))], ","))]," *
                   "\"errors\":$errors_json," *
                   "\"error\":\"$(json_escape(first_error))\"}"
            return _send(sock, 200, body; content_type="application/json")
        elseif path == "/api/state"
            cond = get(qs, "cond", viewer.default_cond)
            band = get(qs, "band", "ALPHA")
            mode = get(qs, "edge_mode", "top_dz")
            topn = something(tryparse(Int, get(qs, "topn", "20")), 20)
            return _send(sock, 200, state_json(viewer, cond, band; edge_mode=mode, topn=topn);
                         content_type="application/json")
        elseif path == "/api/export_png"
            cond = get(qs, "cond", viewer.default_cond)
            band = get(qs, "band", "ALPHA")
            st = viewer.stores[cond]
            out = save_diff_heatmap(st, band)
            body = "{\"ok\":true,\"path\":\"$(json_escape(out))\",\"file\":\"$(json_escape(basename(out)))\"}"
            return _send(sock, 200, body; content_type="application/json")
        elseif path == "/api/export_csv"
            cond = get(qs, "cond", viewer.default_cond)
            band = get(qs, "band", "ALPHA")
            mode = get(qs, "edge_mode", "top_dz")
            topn = something(tryparse(Int, get(qs, "topn", "20")), 20)
            st = viewer.stores[cond]
            out = save_edges_csv(st, band, mode, topn)
            body = "{\"ok\":true,\"path\":\"$(json_escape(out))\",\"file\":\"$(json_escape(basename(out)))\"}"
            return _send(sock, 200, body; content_type="application/json")
        else
            return _send(sock, 404, "Not Found")
        end
    catch e
        msg = sprint(showerror, e)
        try
            _send(sock, 500, "{\"ok\":false,\"error\":\"$(json_escape(msg))\"}";
                  content_type="application/json")
        catch
        end
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

"""
    launch_longitudinal_viewer(; results_root=nothing, port=8780, cond="EC", open=true)

Carga `results/longitudinal/{eyesclosed|eyesopen}/` y sirve el visor interactivo (CLI acepta `EC|EO` como alias).
"""
function launch_longitudinal_viewer(; results_root=nothing, port::Int=PORT,
                                      cond::String="EC", open::Bool=true)
    viewer = load_viewer(; results_root=results_root, default_cond=uppercase(cond))
    println("NeuroMIND longitudinal viewer")
    println("  results: $(viewer.results_root)/longitudinal/")
    for (c, st) in viewer.stores
        println("  $c: $(length(st.bands)) bandas · pares=$(get(st.summary, "n_paired", "?"))")
    end
    for (c, msg) in viewer.errors
        println("  $c: ⚠ $msg")
    end
    server, p = _listen_available(HOST, port)
    url = "http://$HOST:$p/"
    println()
    println("UI → $url")
    p != port && println("(puerto $port ocupado; usando $p)")
    println("Ctrl+C para detener.")
    open && open_browser(url)
    try
        while true
            sock = accept(server)
            @async begin
                try
                    handle_request(sock, viewer)
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

function main()
    cond = "EC"
    for a in ARGS
        ua = uppercase(a)
        ua in ("EC", "EO") && (cond = ua)
    end
    root = nothing
    for a in ARGS
        if isdir(a) || endswith(a, "results")
            root = a
        end
    end
    launch_longitudinal_viewer(; results_root=root, cond=cond)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
