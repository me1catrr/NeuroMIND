#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor interactivo transversal (MS vs Control)
# ═══════════════════════════════════════════════════════════════
#
#  Explora results/transversal/{eyesclosed|eyesopen}/ tras
#  run_transversal_analysis.jl (CLI acepta EC|EO como atajo).
#  · Overview KPIs · strip MS vs Control · heatmaps Ctrl/MS/Δ
#  · red FDR / top-|d| · volcano · potencia Δ · métricas de red
#
#  Uso:
#    julia --project=. src/interactive/plot_transversal.jl
#    julia --project=. src/interactive/plot_transversal.jl EO
#    # → http://127.0.0.1:8781/
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/interactive/plot_transversal.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      25-07-2026
#  Modificado  28-07-2026
# ───────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "viewer_support.jl"))
using .ViewerSupport
using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics, TOML

const PROJ = normpath(joinpath(@__DIR__, "..", ".."))
const HOST = "127.0.0.1"
const PORT = 8781
const COMMON_JS = joinpath(@__DIR__, "viewer_common.js")
const EXPECTED_SCHEMA_VERSION = 2
const INCOMPATIBLE_MESSAGE =
    "Resultados incompatibles con el visor actual. Regenere el análisis transversal."

struct IncompatibleTransversalResultsError <: Exception
    detail::String
end
Base.showerror(io::IO, e::IncompatibleTransversalResultsError) =
    print(io, INCOMPATIBLE_MESSAGE, " ", e.detail)

# ── Tipos ──────────────────────────────────────────────────────

mutable struct CondStore
    cond::String
    dir::String
    bands::Vector{String}
    channels::Dict{String,Vector{String}}          # band → ch
    mat_ctrl::Dict{String,Matrix{Float64}}
    mat_ms::Dict{String,Matrix{Float64}}
    mat_diff::Dict{String,Matrix{Float64}}
    stats::Dict{String,DataFrame}                  # edge stats
    band_stats::DataFrame
    global_band::DataFrame                         # global_mean_wpli_statistics
    subject_means::DataFrame
    inclusion::DataFrame
    summary::Dict{String,Any}
    spectral::DataFrame
    net_global::DataFrame
    net_ctrl::Dict{String,DataFrame}
    net_ms::Dict{String,DataFrame}
    net_diff::Dict{String,DataFrame}
end

mutable struct TransViewer
    results_root::String
    stores::Dict{String,CondStore}   # "EC"|"EO"
    errors::Dict{String,String}
    default_cond::String
end

# ── IO ─────────────────────────────────────────────────────────

"""Carpeta en disco ("eyesclosed"/"eyesopen") para el código corto "EC"/"EO"."""
_cond_dir(cond::String)::String = cond == "EC" ? "eyesclosed" : "eyesopen"

function _require_columns(df::DataFrame, required::Vector{Symbol}, artifact::String)
    missing_cols = filter(c -> !hasproperty(df, c), required)
    isempty(missing_cols) || throw(IncompatibleTransversalResultsError(
        "$artifact: faltan campos obligatorios $(join(string.(missing_cols), ", "))."))
    nrow(df) > 0 || throw(IncompatibleTransversalResultsError("$artifact está vacío."))
    if :schema_version in required
        all(Int.(df.schema_version) .== EXPECTED_SCHEMA_VERSION) ||
            throw(IncompatibleTransversalResultsError(
                "$artifact usa una versión de esquema no compatible."))
    end
    return df
end

function _required_csv(path::String, required::Vector{Symbol})::DataFrame
    isfile(path) || throw(IncompatibleTransversalResultsError(
        "Falta el artefacto obligatorio $(basename(path))."))
    df = try
        CSV.read(path, DataFrame)
    catch e
        throw(IncompatibleTransversalResultsError(
            "No se puede leer $(basename(path)): $(sprint(showerror, e))."))
    end
    return _require_columns(df, required, basename(path))
end

function _required_matrix(path::String)
    result = try
        read_mat_csv(path)
    catch e
        throw(IncompatibleTransversalResultsError(
            "No se puede leer $(basename(path)): $(sprint(showerror, e))."))
    end
    result === nothing && throw(IncompatibleTransversalResultsError(
        "Falta o es inválido el artefacto obligatorio $(basename(path))."))
    return result
end

function _validate_provenance(df::DataFrame, contract::Dict{String,Any},
                              artifact::String; fdr_key::Union{Nothing,String}=nothing)
    expected_source = string(contract["statistics_source"])
    all(string.(df.statistics_source) .== expected_source) ||
        throw(IncompatibleTransversalResultsError(
            "$artifact no coincide con statistics_source del contrato."))
    if fdr_key !== nothing
        expected_scope = string(contract[fdr_key])
        all(string.(df.fdr_scope) .== expected_scope) ||
            throw(IncompatibleTransversalResultsError(
                "$artifact no coincide con $fdr_key del contrato."))
    end
    return df
end

function validate_condition_contract(dir::String)::Dict{String,Any}
    path = joinpath(dir, "statistics_contract.toml")
    isfile(path) || throw(IncompatibleTransversalResultsError(
        "Falta statistics_contract.toml."))
    contract = try
        TOML.parsefile(path)
    catch e
        throw(IncompatibleTransversalResultsError(
            "statistics_contract.toml no es legible: $(sprint(showerror, e))."))
    end
    required = [
        "schema_version", "statistics_source",
        "fdr_scope_edges", "fdr_scope_global", "fdr_scope_power",
        "bootstrap_method", "bootstrap_iterations", "bootstrap_seed",
        "quantile_method", "rrb_method", "effect_d_pooled_method",
        "mannwhitney_method",
    ]
    missing_fields = filter(k -> !haskey(contract, k), required)
    isempty(missing_fields) || throw(IncompatibleTransversalResultsError(
        "statistics_contract.toml: faltan campos obligatorios $(join(missing_fields, ", "))."))
    Int(contract["schema_version"]) == EXPECTED_SCHEMA_VERSION ||
        throw(IncompatibleTransversalResultsError(
            "schema_version=$(contract["schema_version"]); se requiere $EXPECTED_SCHEMA_VERSION."))
    return contract
end

function load_condition(dir::String, cond::String)::CondStore
    contract = validate_condition_contract(dir)
    tab_dir = joinpath(dir, "tables")
    band_stats = _required_csv(joinpath(tab_dir, "band_statistics.csv"), [
        :schema_version, :statistics_source, :fdr_scope, :band, :n_channels,
        :n_edges, :n_pairs, :fdr_family_size, :n_nominal, :n_sig,
        :top20_nominal_overlap, :ctrl_mean, :ms_mean, :diff_mean,
        :mean_abs_d_pooled, :n_ms, :n_ctrl,
    ])
    _validate_provenance(band_stats, contract, "band_statistics.csv";
                         fdr_key="fdr_scope_edges")
    bands = [b for b in BAND_ORDER if b in String.(band_stats.band)]
    bands == BAND_ORDER || throw(IncompatibleTransversalResultsError(
        "band_statistics.csv debe contener exactamente las 7 bandas de producción."))
    channels = Dict{String,Vector{String}}()
    mat_ctrl = Dict{String,Matrix{Float64}}()
    mat_ms = Dict{String,Matrix{Float64}}()
    mat_diff = Dict{String,Matrix{Float64}}()
    stats = Dict{String,DataFrame}()
    net_ctrl = Dict{String,DataFrame}()
    net_ms = Dict{String,DataFrame}()
    net_diff = Dict{String,DataFrame}()

    for b in bands
        rc = _required_matrix(joinpath(tab_dir, "group_connectivity_control_$(b).csv"))
        rm = _required_matrix(joinpath(tab_dir, "group_connectivity_ms_$(b).csv"))
        rd = _required_matrix(joinpath(tab_dir, "group_difference_$(b).csv"))
        rc[1] == rm[1] == rd[1] || throw(IncompatibleTransversalResultsError(
            "Las matrices Control, MS y Δ de $b no comparten el mismo montaje."))
        channels[b] = rc[1]
        mat_ctrl[b] = rc[2]
        mat_ms[b] = rm[2]
        mat_diff[b] = rd[2]
        sdf = _required_csv(joinpath(tab_dir, "group_statistics_$(b).csv"), [
            :schema_version, :statistics_source, :fdr_scope, :fdr_family_size,
            :ch_a, :ch_b, :ctrl_mean, :ms_mean, :diff,
            :p_value, :p_mannwhitney, :p_welch, :q_value,
            :effect_d_pooled, :effect_rrb, :effect_rank_abs_d,
            :is_nominal, :is_fdr, :n_ms, :n_ctrl,
        ])
        _validate_provenance(sdf, contract, "group_statistics_$(b).csv";
                             fdr_key="fdr_scope_edges")
        all(Int.(sdf.fdr_family_size) .== nrow(sdf)) ||
            throw(IncompatibleTransversalResultsError(
                "group_statistics_$(b).csv declara un tamaño de familia FDR incoherente."))
        for c in (:ch_a, :ch_b)
            sdf[!, c] = String.(sdf[!, c])
        end
        stats[b] = sdf
        net_ctrl[b] = _required_csv(joinpath(tab_dir, "network_metrics_control_$(b).csv"),
            [:schema_version, :channel, :strength, :degree, :norm_strength])
        net_ms[b] = _required_csv(joinpath(tab_dir, "network_metrics_ms_$(b).csv"),
            [:schema_version, :channel, :strength, :degree, :norm_strength])
        net_diff[b] = _required_csv(joinpath(tab_dir, "network_metrics_diff_$(b).csv"),
            [:schema_version, :channel, :delta_strength, :delta_degree, :delta_norm_strength])
    end

    global_band = _required_csv(joinpath(tab_dir, "global_mean_wpli_statistics.csv"), [
        :schema_version, :statistics_source, :fdr_scope, :fdr_family_size,
        :band, :ms_mean, :ms_ci_low, :ms_ci_high, :ms_sem,
        :ms_median, :ms_q1, :ms_q3,
        :ctrl_mean, :ctrl_ci_low, :ctrl_ci_high, :ctrl_sem,
        :ctrl_median, :ctrl_q1, :ctrl_q3,
        :diff, :diff_ci_low, :diff_ci_high,
        :p_mannwhitney, :p_welch, :q_value,
        :effect_d_pooled, :effect_d_pooled_ci_low, :effect_d_pooled_ci_high,
        :effect_rrb, :probability_superiority, :n_ms, :n_ctrl,
        :bootstrap_method, :bootstrap_iterations, :bootstrap_seed,
        :quantile_method, :rrb_method, :effect_d_pooled_method,
        :mannwhitney_method, :is_fdr, :is_largest_abs_effect,
    ])
    _validate_provenance(global_band, contract, "global_mean_wpli_statistics.csv";
                         fdr_key="fdr_scope_global")
    nrow(global_band) == length(BAND_ORDER) ||
        throw(IncompatibleTransversalResultsError(
            "global_mean_wpli_statistics.csv debe contener las 7 bandas."))
    all(Int.(global_band.fdr_family_size) .== length(BAND_ORDER)) ||
        throw(IncompatibleTransversalResultsError(
            "global_mean_wpli_statistics.csv debe declarar una familia FDR de 7 bandas."))

    subject_means = _required_csv(joinpath(tab_dir, "subject_band_means.csv"),
        [:schema_version, :statistics_source, :subject_id, :group, :band, :cond, :mean_wpli])
    _validate_provenance(subject_means, contract, "subject_band_means.csv")
    inclusion = _required_csv(joinpath(tab_dir, "subject_inclusion.csv"),
        [:schema_version, :statistics_source, :subject_id, :session_id, :group,
         :n_bands_ok, :included, :excluded_reason, :qc_decision])
    _validate_provenance(inclusion, contract, "subject_inclusion.csv")
    spectral = _required_csv(joinpath(tab_dir, "band_power_group_statistics.csv"), [
        :schema_version, :statistics_source, :fdr_scope, :fdr_family_size,
        :channel, :band, :ms_mean, :ctrl_mean, :diff, :p_value, :q_value,
        :effect_d_pooled, :n_ms, :n_ctrl,
        :coverage_ms_pct, :coverage_ctrl_pct, :coverage_low, :is_fdr,
    ])
    _validate_provenance(spectral, contract, "band_power_group_statistics.csv";
                         fdr_key="fdr_scope_power")
    for b in bands
        sub = filter(r -> string(r.band) == b, spectral)
        nrow(sub) > 0 || throw(IncompatibleTransversalResultsError(
            "band_power_group_statistics.csv no contiene la banda $b."))
        all(Int.(sub.fdr_family_size) .== nrow(sub)) ||
            throw(IncompatibleTransversalResultsError(
                "band_power_group_statistics.csv declara una familia FDR incoherente para $b."))
    end
    net_global = _required_csv(joinpath(tab_dir, "network_global_statistics.csv"), [
        :schema_version, :statistics_source, :fdr_scope, :fdr_family_size,
        :band, :metric, :n_channels, :ms_mean, :ctrl_mean, :diff,
        :p_value, :p_welch, :q_value,
        :effect_d_pooled, :effect_d_pooled_ci_low, :effect_d_pooled_ci_high,
        :effect_rrb, :probability_superiority, :n_ms, :n_ctrl,
        :bootstrap_method, :bootstrap_iterations, :bootstrap_seed,
        :rrb_method, :effect_d_pooled_method, :mannwhitney_method, :is_fdr,
    ])
    _validate_provenance(net_global, contract, "network_global_statistics.csv";
                         fdr_key="fdr_scope_global")
    nrow(net_global) == length(BAND_ORDER) ||
        throw(IncompatibleTransversalResultsError(
            "network_global_statistics.csv debe contener las 7 bandas."))
    summary = parse_summary_json(joinpath(dir, "transversal_summary.json"))
    Int(get(summary, "schema_version", 0)) == EXPECTED_SCHEMA_VERSION ||
        throw(IncompatibleTransversalResultsError(
            "transversal_summary.json no declara schema_version compatible."))
    string(get(summary, "statistics_source", "")) == string(contract["statistics_source"]) ||
        throw(IncompatibleTransversalResultsError(
            "La procedencia estadística no coincide entre contrato y resumen."))

    CondStore(
        cond, dir, bands, channels, mat_ctrl, mat_ms, mat_diff, stats,
        band_stats, global_band, subject_means, inclusion, summary, spectral, net_global,
        net_ctrl, net_ms, net_diff,
    )
end

function load_viewer(; results_root::Union{Nothing,String}=nothing,
                       default_cond::String="EC")::TransViewer
    root = resolve_results_root(PROJ, results_root)
    stores = Dict{String,CondStore}()
    errors = Dict{String,String}()
    for cond in ("EC", "EO")
        d = joinpath(root, "transversal", _cond_dir(cond))
        isdir(d) || continue
        try
            stores[cond] = load_condition(d, cond)
        catch e
            errors[cond] = sprint(showerror, e)
        end
    end
    if isempty(stores) && isempty(errors)
        errors["ALL"] = INCOMPATIBLE_MESSAGE *
            " No hay artefactos en transversal/{eyesclosed|eyesopen}."
    end
    dc = haskey(stores, default_cond) ? default_cond :
         (!isempty(stores) ? first(sort(collect(keys(stores)))) : default_cond)
    return TransViewer(root, stores, errors, dc)
end

# ── JSON / HTTP helpers ────────────────────────────────────────

function _send(sock, status::Int, body::String; content_type::String="text/html; charset=utf-8")
    bytes = codeunits(body)
    header = "HTTP/1.1 $status\r\nContent-Type: $content_type\r\n" *
             "Content-Length: $(length(bytes))\r\nConnection: close\r\n" *
             "Access-Control-Allow-Origin: *\r\n\r\n"
    write(sock, header)
    write(sock, bytes)
end

function _sort_band_df(df::DataFrame)::DataFrame
    nrow(df) == 0 && return df
    hasproperty(df, :band) || return df
    order = Dict(b => i for (i, b) in enumerate(BAND_ORDER))
    return sort(df, :band; by = b -> get(order, string(b), 1000))
end

function state_json(viewer::TransViewer, cond::String, band::String;
                    edge_mode::String="fdr", topn::Int=20)::String
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
            edges_out = sort(filter(:is_fdr => identity, edf), :effect_rank_abs_d)
        elseif edge_mode == "p05"
            edges_out = sort(filter(:is_nominal => identity, edf), :effect_rank_abs_d)
        else
            edges_out = sort(
                filter(r -> Int(r.effect_rank_abs_d) <= topn, edf),
                :effect_rank_abs_d)
        end
    end
    edge_cols = Symbol[:ch_a, :ch_b, :ctrl_mean, :ms_mean, :diff, :p_value,
                       :p_mannwhitney, :p_welch, :q_value,
                       :effect_d_pooled, :effect_rrb, :effect_rank_abs_d,
                       :is_nominal, :is_fdr, :n_ms, :n_ctrl]

    sm = st.subject_means
    sm_band = nrow(sm) > 0 && hasproperty(sm, :band) ?
              filter(r -> string(r.band) == band, sm) : DataFrame()
    sm_cols = Symbol[:subject_id, :group, :band, :cond, :mean_wpli]

    sp = st.spectral
    sp_band = nrow(sp) > 0 && hasproperty(sp, :band) ?
              filter(r -> string(r.band) == band, sp) : DataFrame()
    sp_cols = Symbol[:channel, :band, :ms_mean, :ctrl_mean, :diff, :p_value,
                     :effect_d_pooled, :q_value, :n_ms, :n_ctrl,
                     :fdr_family_size, :coverage_ms_pct, :coverage_ctrl_pct,
                     :coverage_low, :is_fdr]

    nc = get(st.net_ctrl, band, DataFrame())
    nm = get(st.net_ms, band, DataFrame())
    nd = get(st.net_diff, band, DataFrame())
    for df in (nc, nm, nd)
        nrow(df) > 0 && hasproperty(df, :channel) && (df.channel = String.(df.channel))
    end

    summ = st.summary
    n_total_sig = Int(get(summ, "n_total_sig", 0))
    best_raw = string(get(summ, "best_band", ""))
    best_band = n_total_sig == 0 ? "" : best_raw
    best_n_sig = 0
    if !isempty(best_band) && nrow(st.band_stats) > 0 && hasproperty(st.band_stats, :band)
        br = filter(r -> string(r.band) == best_band, st.band_stats)
        best_n_sig = nrow(br) > 0 ? Int(br[1, :n_sig]) : 0
    end

    n_ms_cand = 0
    n_ctrl_cand = 0
    n_excl_qc_ms = 0
    n_excl_qc_ctrl = 0
    n_warn_ms = 0
    n_warn_ctrl = 0
    if nrow(st.inclusion) > 0 && hasproperty(st.inclusion, :group)
        for r in eachrow(st.inclusion)
            g = lowercase(string(r.group))
            is_ms = g in ("ms", "em", "patient")
            is_ct = g in ("control", "ctrl", "hc")
            if is_ms
                n_ms_cand += 1
            elseif is_ct
                n_ctrl_cand += 1
            end
            incl = hasproperty(r, :included) ? r.included : true
            included = incl === true || string(incl) == "true"
            qc = hasproperty(r, :qc_decision) ? string(r.qc_decision) : ""
            if included && qc == "include_with_warning"
                is_ms && (n_warn_ms += 1)
                is_ct && (n_warn_ctrl += 1)
            end
            if !included
                reason = hasproperty(r, :excluded_reason) ? string(something(r.excluded_reason, "")) : ""
                if occursin("QC=", reason) || occursin("QC=", uppercase(reason)) ||
                   (qc in ("exclude", "manual_review"))
                    is_ms && (n_excl_qc_ms += 1)
                    is_ct && (n_excl_qc_ctrl += 1)
                end
            end
        end
    end
    n_ms_incl = Int(get(summ, "n_ms", 0))
    n_ctrl_incl = Int(get(summ, "n_ctrl", 0))
    n_ms_design = Int(get(summ, "n_ms_design", 44))
    n_ctrl_design = Int(get(summ, "n_ctrl_design", 40))
    n_ms_no_t1 = max(0, n_ms_design - n_ms_cand)
    n_ctrl_no_t1 = max(0, n_ctrl_design - n_ctrl_cand)

    kpi = """{
      "n_ms": $n_ms_incl,
      "n_ctrl": $n_ctrl_incl,
      "n_ms_design": $n_ms_design,
      "n_ctrl_design": $n_ctrl_design,
      "n_ms_candidates": $n_ms_cand,
      "n_ctrl_candidates": $n_ctrl_cand,
      "n_ms_no_t1": $n_ms_no_t1,
      "n_ctrl_no_t1": $n_ctrl_no_t1,
      "n_excluded_qc_ms": $n_excl_qc_ms,
      "n_excluded_qc_ctrl": $n_excl_qc_ctrl,
      "n_warn_ms": $n_warn_ms,
      "n_warn_ctrl": $n_warn_ctrl,
      "n_included": $(Int(get(summ, "n_included", 0))),
      "n_excluded": $(Int(get(summ, "n_excluded", 0))),
      "n_excluded_qc": $(Int(get(summ, "n_excluded_qc", 0))),
      "n_excluded_data": $(Int(get(summ, "n_excluded_data", 0))),
      "n_total_sig": $n_total_sig,
      "best_band": "$(json_escape(best_band))",
      "best_band_n_sig": $best_n_sig,
      "test": "$(json_escape(string(get(summ, "test", "mannwhitney"))))",
      "design": "$(json_escape(string(get(summ, "design", "case_control_T1"))))",
      "wpli_method": "$(json_escape(string(get(summ, "wpli_method", ""))))",
      "timestamp": "$(json_escape(string(get(summ, "timestamp", ""))))",
      "fdr_note": "BH within each band (n_pairs edges); n_total_sig is the sum of edge×band discoveries, not a joint 7×n_pairs FDR",
      "spectral_note": "Absolute band power (μV²) from PSD integral; no log/dB; FDR-BH within each band across channels; effect_d_pooled = Cohen d with pooled SD; diff = MS−Control"
    }"""

    bs_cols = Symbol[:band, :n_channels, :n_edges, :n_pairs, :fdr_family_size,
                     :n_nominal, :n_sig, :top20_nominal_overlap, :pct_sig,
                     :ctrl_mean, :ms_mean, :diff_mean, :mean_p,
                     :mean_abs_d_pooled, :n_ms, :n_ctrl]
    gb_cols = Symbol[
        :band, :ms_mean, :ms_ci_low, :ms_ci_high, :ms_sem,
        :ms_median, :ms_q1, :ms_q3,
        :ctrl_mean, :ctrl_ci_low, :ctrl_ci_high, :ctrl_sem,
        :ctrl_median, :ctrl_q1, :ctrl_q3,
        :diff, :diff_ci_low, :diff_ci_high,
        :p_mannwhitney, :p_welch, :q_value,
        :effect_d_pooled, :effect_d_pooled_ci_low, :effect_d_pooled_ci_high,
        :effect_rrb, :probability_superiority, :n_ms, :n_ctrl,
        :fdr_family_size, :is_fdr, :is_largest_abs_effect,
        :bootstrap_method, :bootstrap_iterations, :bootstrap_seed,
        :quantile_method, :rrb_method, :effect_d_pooled_method, :mannwhitney_method,
    ]
    incl_use = filter(c -> hasproperty(st.inclusion, c),
        [:subject_id, :session_id, :group, :n_bands_ok, :included,
         :excluded_reason, :qc_decision, :warning_reason, :n_epochs_valid])
    netg_cols = Symbol[:band, :metric, :n_channels, :ms_mean, :ctrl_mean, :diff,
                       :diff_ci_low, :diff_ci_high, :p_value, :p_welch, :q_value,
                       :effect_d_pooled, :effect_d_pooled_ci_low,
                       :effect_d_pooled_ci_high, :effect_rrb,
                       :probability_superiority, :n_ms, :n_ctrl,
                       :fdr_family_size, :is_fdr]

    bs_sorted = _sort_band_df(st.band_stats)
    bs_band = filter(r -> string(r.band) == band, bs_sorted)
    nrow(bs_band) == 1 || throw(IncompatibleTransversalResultsError(
        "Falta el resumen persistido de la banda $band."))
    bs_row = bs_band[1, :]
    n_sig_band = Int(bs_row.n_sig)
    n_p_uncorr = Int(bs_row.n_nominal)
    top_all = filter(r -> Int(r.effect_rank_abs_d) <= topn, edf)
    n_top_nominal = topn == 20 ? Int(bs_row.top20_nominal_overlap) :
                    count(r -> Bool(r.is_nominal), eachrow(top_all))
    Mc = st.mat_ctrl[band]
    Mm = st.mat_ms[band]
    shared_lim = max(maximum(abs, Mc), maximum(abs, Mm), 1e-12)

    gb_sorted = _sort_band_df(st.global_band)
    n_global_fdr = count(identity, Bool.(st.global_band.is_fdr))
    descriptive_band = ""
    idx_desc = findfirst(==(true), Bool.(st.global_band.is_largest_abs_effect))
    idx_desc === nothing || (descriptive_band = string(st.global_band.band[idx_desc]))

    sp_use = filter(c -> nrow(sp_band) == 0 || hasproperty(sp_band, c), sp_cols)
    ng_use = filter(c -> nrow(st.net_global) == 0 || hasproperty(st.net_global, c), netg_cols)
    nc_use = filter(c -> nrow(nc) == 0 || hasproperty(nc, c), [:channel, :strength, :degree, :norm_strength])
    nm_use = filter(c -> nrow(nm) == 0 || hasproperty(nm, c), [:channel, :strength, :degree, :norm_strength])
    nd_use = filter(c -> nrow(nd) == 0 || hasproperty(nd, c),
                    [:channel, :delta_strength, :delta_degree, :delta_norm_strength])
    gb_use = filter(c -> nrow(gb_sorted) == 0 || hasproperty(gb_sorted, c), gb_cols)
    bs_use = filter(c -> nrow(bs_sorted) == 0 || hasproperty(bs_sorted, c), bs_cols)
    edge_use = filter(c -> nrow(edf) == 0 || hasproperty(edf, c), edge_cols)
    edges_use_out = filter(c -> nrow(edges_out) == 0 || hasproperty(edges_out, c), edge_cols)

    n_spec_ch = 0
    if nrow(sp_band) > 0 && hasproperty(sp_band, :channel)
        n_spec_ch = length(unique(String.(sp_band.channel)))
    end

    return """{
  "ok": true,
  "cond": "$cond",
  "band": "$band",
  "bands": [$(join(["\"$b\"" for b in st.bands], ","))],
  "conditions": [$(join(["\"$c\"" for c in sort(collect(keys(viewer.stores)))], ","))],
  "n": $n,
  "n_channels": $n,
  "n_spectral_channels": $n_spec_ch,
  "channels": $ch_json,
  "positions": $(pos_json()),
  "matrix_ctrl": $(mat_flat(st.mat_ctrl[band])),
  "matrix_ms": $(mat_flat(st.mat_ms[band])),
  "matrix_diff": $(mat_flat(st.mat_diff[band])),
  "shared_lim_ctrl_ms": $(repr(shared_lim)),
  "n_sig_band": $n_sig_band,
  "n_p_uncorr": $n_p_uncorr,
  "edge_counts": {
    "total": $(Int(bs_row.n_edges)),
    "nominal": $n_p_uncorr,
    "fdr": $n_sig_band,
    "top_n": $(nrow(top_all)),
    "top_nominal_overlap": $n_top_nominal,
    "fdr_family_size": $(Int(bs_row.fdr_family_size))
  },
  "global_counts": {
    "total_bands": $(nrow(st.global_band)),
    "fdr_bands": $n_global_fdr,
    "fdr_family_size": $(Int(first(st.global_band.fdr_family_size))),
    "largest_descriptive_band": "$(json_escape(descriptive_band))"
  },
  "edges": $(df_rows_json(edges_out, edges_use_out)),
  "edge_stats": $(df_rows_json(edf, edge_use)),
  "kpi": $kpi,
  "band_stats": $(df_rows_json(bs_sorted, bs_use)),
  "global_band": $(df_rows_json(gb_sorted, gb_use)),
  "subject_means": $(df_rows_json(sm_band, sm_cols)),
  "inclusion": $(df_rows_json(st.inclusion, incl_use)),
  "spectral": $(df_rows_json(sp_band, sp_use)),
  "net_global": $(df_rows_json(st.net_global, ng_use)),
  "net_ctrl": $(df_rows_json(nc, nc_use)),
  "net_ms": $(df_rows_json(nm, nm_use)),
  "net_diff": $(df_rows_json(nd, nd_use))
}"""
end

# ── PNG export ─────────────────────────────────────────────────

function save_diff_heatmap(st::CondStore, band::String)::String
    M = st.mat_diff[band]
    ch = st.channels[band]
    out = joinpath(st.dir, "figures", "viewer_delta_$(band).png")
    save_heatmap(out, M, ch, "Δ wPLI (MS−Control) — $band ($(st.cond))";
                 diverging=true, colorbar_label="Δ wPLI")
    return out
end

function save_edges_csv(st::CondStore, band::String, edge_mode::String, topn::Int)::String
    edf = get(st.stats, band, DataFrame())
    nrow(edf) == 0 && error("Sin estadísticas para $band")
    if edge_mode == "fdr"
        out_df = sort(filter(:is_fdr => identity, edf), :effect_rank_abs_d)
    elseif edge_mode == "p05"
        out_df = sort(filter(:is_nominal => identity, edf), :effect_rank_abs_d)
    else
        out_df = sort(filter(r -> Int(r.effect_rank_abs_d) <= topn, edf),
                      :effect_rank_abs_d)
    end
    out = joinpath(st.dir, "figures", "viewer_edges_$(edge_mode)_$(band).csv")
    mkpath(dirname(out))
    CSV.write(out, out_df)
    return out
end

# ── HTML SPA ───────────────────────────────────────────────────

function html_page()::String
    # JS uses string concat (no template literals) to avoid Julia `$` interpolation.
    return """<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>NeuroMIND — Transversal MS vs Control</title>
<style>
:root {
  --bg: #eef2f6; --card: #ffffff; --ink: #0f172a; --muted: #64748b;
  --line: #e2e8f0; --accent: #1e3a5f; --accent2: #0f766e;
  --ms: #b91c1c; --ctrl: #1d4ed8; --ok: #15803d; --warn: #a16207;
  --fdr: #b91c1c; --explore: #c2410c; --soft: #f8fafc;
}
* { box-sizing: border-box; }
body { margin:0; font-family: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
  background:
    radial-gradient(1200px 500px at 10% -10%, #dbeafe 0%, transparent 55%),
    radial-gradient(900px 400px at 100% 0%, #ccfbf1 0%, transparent 50%),
    var(--bg);
  color: var(--ink); }
header { background: linear-gradient(125deg, #0b1220 0%, #1e3a5f 52%, #0f766e 120%);
  color:#fff; padding: 1.1rem 1.35rem 0.95rem;
  box-shadow: 0 8px 24px rgba(15,23,42,0.18); }
header h1 { margin:0; font-size:1.3rem; font-weight:650; letter-spacing:-0.02em; }
header p { margin:0.3rem 0 0; opacity:0.88; font-size:0.86rem; }
.toolbar { display:flex; flex-wrap:wrap; gap:0.65rem; align-items:center;
  padding:0.8rem 1.25rem; background:rgba(255,255,255,0.92); backdrop-filter: blur(8px);
  border-bottom:1px solid var(--line); position:sticky; top:0; z-index:20; }
.toolbar label { font-size:0.72rem; color:var(--muted); margin-right:0.15rem;
  text-transform:uppercase; letter-spacing:0.04em; font-weight:600; }
.toolbar select, .toolbar button, .seg button {
  font: inherit; font-size:0.85rem; padding:0.38rem 0.7rem; border-radius:8px;
  border:1px solid var(--line); background:#fff; cursor:pointer; transition: all .15s ease; }
.toolbar select:hover, .toolbar button:hover { border-color:#94a3b8; }
.toolbar button:active { transform: translateY(1px); }
.seg { display:inline-flex; border:1px solid var(--line); border-radius:9px; overflow:hidden; }
.seg button { border:0; border-radius:0; background:#f8fafc; }
.seg button.active { background:var(--accent); color:#fff; }
.badge { font-size:0.75rem; padding:0.28rem 0.6rem; border-radius:999px;
  background:#eff6ff; color:#1e3a5f; border:1px solid #bfdbfe; font-weight:500; }
.badge.mode-fdr { background:#fef2f2; color:#991b1b; border-color:#fecaca; }
.badge.mode-explore { background:#fff7ed; color:#9a3412; border-color:#fed7aa; }
.badge.mode-nom { background:#f8fafc; color:#475569; border-color:#cbd5e1; }
#topn-wrap { display:none; align-items:center; gap:0.25rem; }
.tabs { display:flex; gap:0.25rem; padding:0.55rem 1.25rem 0; overflow-x:auto; }
.tabs button { border:1px solid transparent; background:transparent; padding:0.55rem 0.9rem;
  font:inherit; font-size:0.85rem; color:var(--muted); cursor:pointer; border-radius:10px 10px 0 0;
  transition: color .15s, background .15s; }
.tabs button:hover { color:var(--ink); }
.tabs button.active { background:var(--card); color:var(--ink); border-color:var(--line);
  border-bottom-color:var(--card); font-weight:650; box-shadow: 0 -2px 0 var(--accent2) inset; }
main { padding:0 1.25rem 2.5rem; }
.panel { display:none; background:var(--card); border:1px solid var(--line);
  border-radius:0 12px 12px 12px; padding:1.1rem; min-height:420px;
  box-shadow: 0 10px 30px rgba(15,23,42,0.04); }
.panel.active { display:block; animation: fadeIn .25s ease; }
@keyframes fadeIn { from { opacity:0; transform:translateY(4px); } to { opacity:1; transform:none; } }
.grid { display:grid; gap:1rem; }
.grid-3 { grid-template-columns: repeat(3, 1fr); }
.grid-2 { grid-template-columns: 1fr 1fr; }
@media (max-width: 1100px) { .grid-3, .grid-2 { grid-template-columns: 1fr; } }
.card { border:1px solid var(--line); border-radius:12px; padding:0.85rem 0.9rem; background:var(--soft); }
.card h3 { margin:0 0 0.55rem; font-size:0.92rem; letter-spacing:-0.01em; }
.kpi-row { display:flex; flex-wrap:wrap; gap:0.7rem; margin-bottom:1rem; }
.kpi { min-width:108px; padding:0.7rem 0.85rem; border-radius:12px; background:#fff;
  border:1px solid var(--line); box-shadow: 0 1px 2px rgba(15,23,42,0.03); }
.kpi .v { font-size:1.32rem; font-weight:700; letter-spacing:-0.02em; }
.kpi .l { font-size:0.68rem; color:var(--muted); text-transform:uppercase; letter-spacing:0.04em; margin-top:0.15rem; }
.flow { display:flex; flex-wrap:wrap; gap:0.45rem; align-items:stretch; margin-bottom:1rem; }
.flow-step { flex:1; min-width:140px; background:#fff; border:1px solid var(--line);
  border-radius:12px; padding:0.7rem 0.8rem; position:relative; }
.flow-step .step-l { font-size:0.68rem; text-transform:uppercase; letter-spacing:0.05em;
  color:var(--muted); font-weight:650; }
.flow-step .step-v { font-size:1.05rem; font-weight:700; margin-top:0.2rem; }
.flow-step .step-s { font-size:0.72rem; color:var(--muted); margin-top:0.15rem; }
.flow-arrow { display:flex; align-items:center; color:#94a3b8; font-size:1.1rem; padding:0 0.1rem; }
.flow-step.ok { border-color:#86efac; background:linear-gradient(180deg,#fff,#f0fdf4); }
.flow-step.warn { border-color:#fcd34d; background:linear-gradient(180deg,#fff,#fffbeb); }
.canvas-wrap { width:100%; }
.canvas-wrap canvas { display:block; background:#fff; border-radius:10px;
  border:1px solid var(--line); cursor:crosshair; }
table { width:100%; border-collapse:collapse; font-size:0.78rem; }
th, td { padding:0.38rem 0.48rem; border-bottom:1px solid var(--line); text-align:left; }
th { color:var(--muted); font-weight:650; position:sticky; top:0; background:#fff; }
.scroll { max-height:320px; overflow:auto; border-radius:8px; }
.muted { color:var(--muted); font-size:0.85rem; }
.tip { position:fixed; pointer-events:none; background:#0f172a; color:#fff; font-size:0.75rem;
  padding:0.4rem 0.6rem; border-radius:8px; z-index:50; display:none; max-width:300px;
  white-space:pre-line; box-shadow: 0 8px 20px rgba(0,0,0,0.25); }
.status { font-size:0.8rem; color:var(--muted); margin-left:auto; }
.status.ok { color:var(--ok); } .status.err { color:var(--ms); }
.note { font-size:0.85rem; line-height:1.5; color:#334155; background:linear-gradient(180deg,#f8fafc,#f1f5f9);
  border:1px solid #cbd5e1; border-radius:10px; padding:0.85rem; margin-top:0.85rem; }
.banner { display:none; margin:0.75rem 1.25rem 0; padding:0.7rem 0.95rem; border-radius:10px;
  background:#fffbeb; border:1px solid #fcd34d; color:#92400e; font-size:0.85rem; }
.banner.show { display:block; }
.banner.explore { background:#fff7ed; border-color:#fdba74; color:#9a3412; }
.empty-cta { padding:1.5rem; text-align:center; color:var(--muted); font-size:0.9rem; }
.legend-row { display:flex; flex-wrap:wrap; gap:0.75rem; margin-top:0.55rem; font-size:0.78rem; color:#475569; }
.legend-row span { display:inline-flex; align-items:center; gap:0.35rem; }
.swatch { width:18px; height:3px; border-radius:2px; display:inline-block; }
.swatch.fdr { background:var(--fdr); height:4px; }
.swatch.explore { background:var(--explore); height:3px; opacity:0.75;
  background-image: repeating-linear-gradient(90deg, var(--explore) 0 4px, transparent 4px 7px); }
.pill { display:inline-block; font-size:0.7rem; padding:0.15rem 0.45rem; border-radius:999px;
  background:#ecfeff; color:#0f766e; border:1px solid #99f6e4; font-weight:600; margin-left:0.35rem; }
.net-title { display:flex; align-items:baseline; gap:0.5rem; flex-wrap:wrap; }
.net-title h3 { margin:0; }
.fdr-families { width:100%; border-collapse:collapse; font-size:0.78rem; margin:0.4rem 0 0; }
.fdr-families th, .fdr-families td { padding:0.35rem 0.5rem; border-bottom:1px solid var(--line); text-align:left; }
.fdr-families th { color:var(--muted); font-weight:650; background:#fff; }
.fdr-families code { font-size:0.72rem; background:#f1f5f9; padding:0.1rem 0.3rem; border-radius:4px; }
.compat-error { display:none; max-width:900px; margin:2rem auto; padding:1.2rem 1.35rem;
  background:#fff7ed; border:1px solid #fdba74; border-radius:12px; color:#7c2d12; }
.compat-error h2 { margin:0 0 0.45rem; font-size:1.1rem; }
.compat-error p { margin:0.25rem 0; line-height:1.5; }
</style>
</head>
<body>
<header>
  <h1>NeuroMIND — Evaluación transversal MS vs Control</h1>
  <p>Caso-control T1 · EC y EO en paralelo · Mann–Whitney + FDR-BH por banda</p>
</header>
<div class="compat-error" id="compat-error">
  <h2>Resultados incompatibles con el visor actual</h2>
  <p id="compat-message">Regenera el análisis transversal con la versión de producción correspondiente.</p>
  <p><code>julia --project=. scripts/run_transversal_analysis.jl</code></p>
</div>
<div class="toolbar">
  <div class="seg" id="cond-seg">
    <button type="button" data-cond="EC" class="active">EC</button>
    <button type="button" data-cond="EO">EO</button>
  </div>
  <label>Banda <select id="band"></select></label>
  <label>Edges
    <select id="edge-mode">
      <option value="fdr" selected>Solo q&lt;0.05 (FDR)</option>
      <option value="p05">Solo p&lt;0.05 (nominal)</option>
      <option value="top_d">Top-N |d| (exploratorio)</option>
    </select>
  </label>
  <span id="topn-wrap"><label>N <input id="topn" type="number" min="5" max="100" value="20" style="width:4rem;padding:0.3rem;border:1px solid var(--line);border-radius:6px"/></label></span>
  <span class="badge" id="badge-n">—</span>
  <span class="badge" id="badge-mode">—</span>
  <button type="button" id="btn-png">Export PNG Δ</button>
  <button type="button" id="btn-csv">Export CSV edges</button>
  <span class="status" id="status">Cargando…</span>
</div>
<div class="banner" id="fdr-banner">Sin aristas FDR (q&lt;0.05) en esta banda. Vista exploratoria disponible vía Top-N |d|. Confirmatorio solo con FDR.</div>
<div class="tabs" id="tabs">
  <button type="button" data-tab="overview" class="active">Overview</button>
  <button type="button" data-tab="global">Comparación global</button>
  <button type="button" data-tab="conn">Conectividad</button>
  <button type="button" data-tab="net" id="tab-net">Red FDR</button>
  <button type="button" data-tab="volcano">Volcano</button>
  <button type="button" data-tab="spec">Espectros</button>
  <button type="button" data-tab="hubs">Métricas de red</button>
</div>
<main>
  <section class="panel active" id="panel-overview">
    <div class="flow" id="cohort-flow"></div>
    <div class="kpi-row" id="kpis"></div>
    <div class="grid grid-2">
      <div class="card"><h3>Edges FDR (q&lt;0.05) por banda <span class="pill">BH dentro de cada banda</span></h3>
        <div class="canvas-wrap"><canvas id="cv-nsig" data-h="230"></canvas></div></div>
      <div class="card"><h3>Δ wPLI matriz grupal (MS−Ctrl) · ch intersectados <span class="pill" title="Media del triángulo superior de la matriz grupal (promedio de matrices realineadas al montaje común intersectado). Estimando B.">estimando B</span></h3>
        <div class="canvas-wrap"><canvas id="cv-diffbar" data-h="230"></canvas></div></div>
    </div>
    <div class="note" id="interp"></div>
    <div class="card" style="margin-top:1rem">
      <h3>Familias FDR-BH (independientes)</h3>
      <table class="fdr-families">
        <thead><tr><th>Análisis</th><th>Familia</th><th>Etiqueta</th></tr></thead>
        <tbody>
          <tr><td>Aristas wPLI</td><td id="fdr-fam-edges">n_pares aristas dentro de cada banda</td><td><code>q_edge</code></td></tr>
          <tr><td>Potencia espectral</td><td id="fdr-fam-spec">n canales dentro de cada banda</td><td><code>q_power</code></td></tr>
          <tr><td>Comparación global / red</td><td>7 bandas</td><td><code>q_global</code></td></tr>
        </tbody>
      </table>
      <p class="muted" style="margin:0.5rem 0 0">No confundir <code>q_global</code> (7 tests) con <code>q_edge</code> (p. ej. 276 aristas). «Resultados arista×banda FDR» suma descubrimientos por banda, no es una FDR conjunta.</p>
    </div>
    <div class="card" style="margin-top:1rem"><h3>Inclusión de sujetos</h3>
      <div class="scroll"><table id="tbl-incl"><thead></thead><tbody></tbody></table></div>
    </div>
  </section>
  <section class="panel" id="panel-global">
    <div class="card"><h3>Mean wPLI por sujeto — montaje nativo <span class="pill" title="Media del triángulo superior de la matriz wPLI de cada participante en su montaje nativo (antes del hard-intersect de canales). Mann–Whitney sobre esos valores.">estimando A</span></h3>
      <p class="muted" style="margin:0 0 0.5rem">Cada punto = un participante. Barras = media ± SEM. También mediana e IQR. Efectos no paramétricos: correlación biserial de rangos (r_rb) y P(MS&gt;Ctrl).</p>
      <div class="canvas-wrap"><canvas id="cv-strip" data-h="340"></canvas></div>
      <p class="muted" id="mw-effect-note"></p></div>
    <div class="card" style="margin-top:1rem"><h3>Cohen d (MS − Control) · wPLI global por banda</h3>
      <div class="canvas-wrap"><canvas id="cv-bandmeans" data-h="280"></canvas></div>
      <p class="muted" id="cohen-q-note"></p>
      <p class="muted" style="margin-top:0.35rem"><code>q_global</code>: FDR-BH entre las 7 bandas. Barras de error = IC 95 % bootstrap de Cohen d. Color intenso = q_global&lt;0.05.</p></div>
  </section>
  <section class="panel" id="panel-conn">
    <div class="grid grid-3">
      <div class="card"><h3>Control</h3>
        <div class="canvas-wrap"><canvas id="cv-hm-ctrl"></canvas></div></div>
      <div class="card"><h3>MS</h3>
        <div class="canvas-wrap"><canvas id="cv-hm-ms"></canvas></div></div>
      <div class="card"><h3>Δ (MS−Control) · FDR outline · triángulo inferior</h3>
        <div class="canvas-wrap"><canvas id="cv-hm-d"></canvas></div></div>
    </div>
    <p class="muted" id="hm-tip-hint">Hover: medias, Δ, p, q_edge, d, n. Contorno negro = q_edge&lt;0.05 (una vez por arista). Familia FDR: n_pares aristas dentro de la banda. Ctrl/MS matriz completa; Δ solo triángulo inferior.</p>
  </section>
  <section class="panel" id="panel-net">
    <div class="net-title" style="margin-bottom:0.75rem">
      <h3 id="net-heading">Red FDR</h3>
      <span class="badge" id="net-mode-badge">—</span>
    </div>
    <div class="grid grid-2">
      <div class="card"><h3>Grafo (edges filtrados)</h3>
        <div class="canvas-wrap"><canvas id="cv-graph"></canvas></div>
        <div class="legend-row">
          <span><i class="swatch fdr"></i> FDR q&lt;0.05 (confirmatorio)</span>
          <span id="net-legend-explore"><i class="swatch explore"></i> No FDR / exploratorio</span>
          <span>Grosor ∝ |Δ wPLI| · rojo/ámbar = MS&gt;Ctrl</span>
        </div>
        <div class="empty-cta" id="net-empty" style="display:none">Sin edges FDR. Cambia a <b>Top-N |d|</b> para exploración o a <b>p&lt;0.05</b> nominal.</div>
      </div>
      <div class="card"><h3>Tabla de edges</h3>
        <div class="scroll"><table id="tbl-edges"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
  </section>
  <section class="panel" id="panel-volcano">
    <div class="card"><h3>Volcano: efecto vs −log10(p) <span class="pill">color = FDR persistida</span></h3>
      <label class="muted" style="display:inline-flex;align-items:center;gap:0.35rem;margin:0 0 0.5rem">
        Eje X
        <select id="volcano-x" style="font:inherit;font-size:0.85rem;padding:0.25rem 0.45rem;border-radius:6px;border:1px solid var(--line)">
          <option value="rrb" selected>r_rb (rank-biserial, Mann–Whitney)</option>
          <option value="d">Cohen d pooled (complementario)</option>
        </select>
      </label>
      <div class="canvas-wrap"><canvas id="cv-volcano" data-h="320"></canvas></div>
      <p class="muted">Rojo = <code>is_fdr</code> de producción · ámbar = <code>is_nominal</code> sin FDR · gris = no significativo. Línea = p=0.05. Familia BH: n_pares aristas dentro de la banda. r_rb &gt; 0 ⇒ MS &gt; Ctrl.</p>
    </div>
  </section>
  <section class="panel" id="panel-spec">
    <div class="grid grid-2">
      <div class="card"><h3>Topo Δ band power (MS−Control)</h3>
        <div class="canvas-wrap"><canvas id="cv-topo"></canvas></div>
        <p class="muted" id="spec-ch-note"></p>
        <p class="muted">Potencia <b>absoluta</b> de banda (μV²), integral del PSD; sin log/dB. diff = MS−Control. effect_d_pooled = Cohen d con SD combinada. <code>q_power</code>: FDR-BH dentro de cada banda entre canales. Contorno grueso = q_power&lt;0.05.</p></div>
      <div class="card"><h3>Canal × banda</h3>
        <p class="muted" id="spec-n-note" style="margin:0 0 0.4rem"></p>
        <div class="scroll"><table id="tbl-spec"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
  </section>
  <section class="panel" id="panel-hubs">
    <div class="grid grid-2">
      <div class="card"><h3>Strength nodal Control vs MS <span class="pill" title="Descriptivo: cada punto es un canal de la matriz grupal media, no un participante. No hay contraste inferencial por nodo.">descriptivo</span></h3>
        <p class="muted" style="margin:0 0 0.4rem">Cada punto es un canal de la matriz grupal (no un participante). Color = Δ strength. No implica hub diferencial significativo.</p>
        <div class="canvas-wrap"><canvas id="cv-strength" data-h="360"></canvas></div></div>
      <div class="card"><h3>Estadística global de red <span class="pill" title="mean_strength calculado por sujeto en el montaje intersectado; equivale a (n_ch−1) veces su wPLI global medio en ese montaje. Inferencia Mann–Whitney sobre esos valores individuales. Distinto del estimando B (matrices grupales descriptivas) y del A (montaje nativo).">estimando C</span></h3>
        <p class="muted" style="margin:0 0 0.4rem"><code>mean_strength</code> por sujeto en montaje intersectado (= (n_ch−1) × wPLI global medio de ese sujeto en ese montaje). Inferencial. Distinto de B (matrices grupales, descriptivo) y de A (montaje nativo). <code>q_global</code>: FDR-BH entre las 7 bandas.</p>
        <div class="scroll"><table id="tbl-netg"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
  </section>
</main>
<div class="tip" id="tip"></div>
<script src="/static/viewer_common.js"></script>
<script>
function el(id) { return document.getElementById(id); }
var S = null;
var cond = 'EC';
var tab = 'overview';
var edgeLookup = null;

function isCtrl(g) {
  var s = String(g || '').toLowerCase();
  return s === 'control' || s === 'ctrl' || s === 'hc';
}
function isMS(g) {
  var s = String(g || '').toLowerCase();
  return s === 'ms' || s === 'em' || s === 'patient';
}

function setStatus(msg, cls) {
  var node = el('status');
  node.textContent = msg;
  node.className = 'status' + (cls ? ' ' + cls : '');
}

function edgeModeMeta(mode, topn, nEdges) {
  var n = (nEdges != null) ? nEdges : null;
  var suf = (n != null) ? (' — ' + n + ' aristas') : '';
  if (mode === 'fdr') {
    return { title: 'Red FDR' + suf, short: 'FDR q<0.05', cls: 'mode-fdr', explore: false,
      legendFdr: true, legendExplore: false };
  }
  if (mode === 'p05') {
    return { title: 'Red nominal' + suf, short: 'Nominal p<0.05', cls: 'mode-nom', explore: true,
      legendFdr: true, legendExplore: true };
  }
  var tn = topn || 20;
  return { title: 'Red exploratoria Top-' + tn + ' por |efecto|' + suf, short: 'Top-' + tn + ' |efecto|',
    cls: 'mode-explore', explore: true, legendFdr: true, legendExplore: true };
}

function syncEdgeModeUI() {
  var mode = el('edge-mode').value;
  var topn = Number(el('topn').value) || 20;
  var nEdges = (S && S.edges) ? S.edges.length : null;
  var meta = edgeModeMeta(mode, topn, nEdges);
  el('topn-wrap').style.display = mode === 'top_d' ? 'inline-flex' : 'none';
  // Tab keeps default FDR label; interior heading always reflects active filter
  el('tab-net').textContent = 'Red FDR';
  el('net-heading').textContent = meta.title;
  var bm = el('badge-mode');
  bm.textContent = meta.short;
  bm.className = 'badge ' + meta.cls;
  var nmb = el('net-mode-badge');
  nmb.textContent = meta.short + (nEdges != null ? (' · ' + nEdges) : '');
  nmb.className = 'badge ' + meta.cls;
  return meta;
}

function refresh() {
  setStatus('Cargando…');
  syncEdgeModeUI();
  var mode = el('edge-mode').value;
  var topn = Number(el('topn').value) || 20;
  var band = el('band').value || 'ALPHA';
  var url = '/api/state?cond=' + encodeURIComponent(cond) +
            '&band=' + encodeURIComponent(band) +
            '&edge_mode=' + encodeURIComponent(mode) +
            '&topn=' + topn;
  fetch(url).then(function (res) { return res.json(); }).then(function (data) {
    S = data;
    if (!S.ok) { setStatus(S.error || 'Error', 'err'); return; }
    var sel = el('band');
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
    edgeLookup = NMV.edgeLookupFactory(S.edge_stats || []);
    var kpi = S.kpi || {};
    var nCh = S.n_channels != null ? S.n_channels : S.n;
    var nFdr = S.n_sig_band || 0;
    var nNom = S.n_p_uncorr || 0;
    el('badge-n').textContent =
      'MS ' + (kpi.n_ms != null ? kpi.n_ms : '—') +
      ' / Ctrl ' + (kpi.n_ctrl != null ? kpi.n_ctrl : '—') +
      ' · ' + nCh + ' ch wPLI' +
      (S.n_spectral_channels ? (' · ' + S.n_spectral_channels + ' ch potencia') : '') +
      ' · ' + nFdr + ' FDR · ' + nNom + ' nominales (p<0.05)';
    var ban = el('fdr-banner');
    var meta = syncEdgeModeUI();
    if ((S.n_sig_band || 0) === 0) {
      ban.className = 'banner show';
      ban.textContent = 'Sin aristas FDR (q<0.05) en esta banda. Vista exploratoria: usa Top-N |d|. Confirmatorio solo con FDR.';
    } else if (meta.explore) {
      ban.className = 'banner show explore';
      ban.textContent = 'Modo exploratorio activo (' + meta.short + '). Hay ' +
        (S.n_sig_band || 0) + ' aristas FDR en esta banda — cámbialas con «Solo q<0.05» para interpretación confirmatoria.';
    } else {
      ban.className = 'banner';
    }
    renderAll();
    meta = syncEdgeModeUI(); // refresh heading with edge count after render
    setStatus(cond + ' · ' + S.band + ' · ' + nCh + ' ch · ' +
              nFdr + ' FDR · ' + nNom + ' nominales (p<0.05) · ' + meta.short, 'ok');
  }).catch(function (e) { setStatus(String(e), 'err'); });
}

function renderAll() {
  if (!S) return;
  renderFlow();
  renderKPIs();
  renderFdrFamilies();
  renderBars();
  renderInterp();
  renderInclusión();
  renderStrip();
  renderCohenBars();
  var lim = S.shared_lim_ctrl_ms || 0;
  var fdrSet = NMV.buildFdrSet(S.edge_stats || []);
  NMV.drawHeatmap(el('cv-hm-ctrl'), S.matrix_ctrl, {
    n: S.n, channels: S.channels, lim: lim, diverging: false
  });
  NMV.drawHeatmap(el('cv-hm-ms'), S.matrix_ms, {
    n: S.n, channels: S.channels, lim: lim, diverging: false
  });
  NMV.drawHeatmap(el('cv-hm-d'), S.matrix_diff, {
    n: S.n, channels: S.channels, diverging: true, fdrSet: fdrSet, triangle: 'lower'
  });
  drawGraphPanel();
  renderEdges();
  var volcanoX = (el('volcano-x') && el('volcano-x').value) || 'rrb';
  NMV.drawVolcano(el('cv-volcano'), S.edge_stats || [], {
    yKey: 'p', xKey: volcanoX,
    xLabel: volcanoX === 'rrb' ? 'r_rb (MS − Control)' : 'Cohen d pooled (MS − Control)'
  });
  drawTopoPanel();
  renderSpec();
  NMV.drawStrength(el('cv-strength'), S.net_ctrl || [], S.net_ms || [], 'Control', 'MS');
  renderNetGlobal();
  syncEdgeModeUI();
}

function renderFdrFamilies() {
  var bs = (S.band_stats || []).find(function (r) { return r.band === S.band; });
  var nPairs = bs && bs.n_pairs != null ? bs.n_pairs : (S.n ? (S.n * (S.n - 1) / 2) : '—');
  var nSpec = S.n_spectral_channels || '—';
  var fe = el('fdr-fam-edges');
  var fs = el('fdr-fam-spec');
  if (fe) fe.textContent = nPairs + ' aristas dentro de cada banda';
  if (fs) fs.textContent = nSpec + ' canales dentro de cada banda';
}

function renderFlow() {
  var k = S.kpi || {};
  var qcMs = k.n_excluded_qc_ms != null ? k.n_excluded_qc_ms : '—';
  var qcCt = k.n_excluded_qc_ctrl != null ? k.n_excluded_qc_ctrl : '—';
  var steps = [
    { l: 'Diseño', v: (k.n_ms_design || 44) + ' / ' + (k.n_ctrl_design || 40),
      s: 'MS / Ctrl (protocolo)', cls: '' },
    { l: 'Candidatos T1', v: (k.n_ms_candidates != null ? k.n_ms_candidates : '—') + ' / ' +
        (k.n_ctrl_candidates != null ? k.n_ctrl_candidates : '—'),
      s: 'Sin T1: ' + (k.n_ms_no_t1 || 0) + ' MS, ' + (k.n_ctrl_no_t1 || 0) + ' Ctrl', cls: 'warn' },
    { l: 'Excluidos QC', v: qcMs + ' / ' + qcCt,
      s: 'MS / Ctrl · datos: ' + (k.n_excluded_data != null ? k.n_excluded_data : 0), cls: 'warn' },
    { l: 'Incluidos', v: (k.n_ms || 0) + ' / ' + (k.n_ctrl || 0),
      s: 'MS / Ctrl en análisis', cls: 'ok' }
  ];
  var html = '';
  steps.forEach(function (st, i) {
    if (i) html += '<div class="flow-arrow">→</div>';
    html += '<div class="flow-step ' + st.cls + '"><div class="step-l">' + st.l +
      '</div><div class="step-v">' + st.v + '</div><div class="step-s">' + st.s + '</div></div>';
  });
  el('cohort-flow').innerHTML = html;
}

function renderKPIs() {
  var k = S.kpi || {};
  var bestLabel = '— (sin FDR)';
  if (k.best_band && (k.n_total_sig || 0) > 0) {
    bestLabel = k.best_band + ' · ' + (k.best_band_n_sig || '?') + ' FDR';
  }
  var items = [
    ['Resultados arista×banda FDR', k.n_total_sig],
    ['Mayor nº aristas FDR', bestLabel],
    ['Incl. con advertencia', (k.n_warn_ms != null ? k.n_warn_ms : '—') + ' / ' +
      (k.n_warn_ctrl != null ? k.n_warn_ctrl : '—')],
    ['Canales wPLI', S.n_channels != null ? S.n_channels : S.n]
  ];
  el('kpis').innerHTML = items.map(function (pair) {
    return '<div class="kpi"><div class="v">' + (pair[1] != null ? pair[1] : '—') +
           '</div><div class="l">' + pair[0] + '</div></div>';
  }).join('');
}

function renderInterp() {
  var k = S.kpi || {};
  var bs = (S.band_stats || []).find(function (r) { return r.band === S.band; });
  var gb = (S.global_band || []).find(function (r) { return r.band === S.band; });
  var t = 'Condición <b>' + cond + '</b>, banda <b>' + S.band + '</b>: ';
  t += (S.n_sig_band || 0) + ' aristas FDR (q&lt;0.05) y ' + (S.n_p_uncorr || 0) +
       ' con p&lt;0.05 sin corregir (Mann–Whitney). ';
  t += '«Resultados arista×banda FDR» = <b>suma de descubrimientos</b> con BH <b>dentro de cada banda</b> (no FDR conjunta 7×n_pares; una arista puede contar en varias bandas). ';
  if (bs) {
    t += '<br><b>Estimando B</b> (matriz grupal, ' + (bs.n_channels || S.n_channels) +
         ' ch intersectados): Δ = ' + Number(bs.diff_mean).toFixed(4) +
         ' (Ctrl=' + Number(bs.ctrl_mean).toFixed(3) +
         ', MS=' + Number(bs.ms_mean).toFixed(3) + '). ';
  }
  if (gb) {
    t += '<b>Estimando A</b> (wPLI global, montaje nativo por sujeto): Cohen d = ' +
         Number(gb.effect_d_pooled).toFixed(3) +
         ' (q=' + Number(gb.q_value).toFixed(4) +
         ', Δ=' + Number(gb.diff).toFixed(4) + '). ';
  }
  t += '<br>Cohorte: ' + k.n_ms + '/' + k.n_ms_design + ' MS y ' +
       k.n_ctrl + '/' + k.n_ctrl_design + ' controles (T1). EC y EO no se promedian.';
  if ((S.n_sig_band || 0) === 0) {
    t += ' <b>Sin FDR en esta banda</b> — Top-N |d| no es confirmatorio.';
  }
  if ((k.n_warn_ms || 0) + (k.n_warn_ctrl || 0) > 0) {
    t += ' Incluidos con advertencia QC: <b>' + (k.n_warn_ms || 0) + ' MS / ' +
         (k.n_warn_ctrl || 0) + ' Ctrl</b> (ver columna warning_reason). Sensibilidad recomendada: repetir excluyendo include_with_warning.';
  }
  el('interp').innerHTML = t;
}

function renderBars() {
  var rows = NMV.sortBandsPhysio(S.band_stats || []);
  var labs = rows.map(function (r) { return r.band; });
  NMV.barChart(el('cv-nsig'), labs, rows.map(function (r) { return Number(r.n_sig) || 0; }),
    function () { return '#1e3a5f'; }, { height: 230 });
  NMV.barChart(el('cv-diffbar'), labs, rows.map(function (r) { return Number(r.diff_mean) || 0; }),
    function (v) { return v >= 0 ? '#b91c1c' : '#1d4ed8'; }, { height: 230, signed: true });
}

function renderInclusión() {
  NMV.fillTable(el('tbl-incl'),
    ['subject_id', 'group', 'included', 'qc_decision', 'warning_reason', 'n_bands_ok', 'n_epochs_valid', 'excluded_reason'],
    S.inclusion || [], {
      excluded_reason: function (v) {
        if (v == null || v === '' || v === 'missing') return '—';
        return v;
      },
      warning_reason: function (v) {
        if (v == null || v === '' || v === 'missing') return '—';
        return v;
      },
      included: function (v) {
        return (v === true || v === 'true') ? 'sí' : 'no';
      }
    });
}

function renderStrip() {
  var cv = el('cv-strip');
  var sized = NMV.resizeCanvas(cv, { height: 340, fallbackW: 700 });
  var ctx = sized.ctx, w = sized.w, h = sized.h;
  var rows = S.subject_means || [];
  var ctrl = rows.filter(function (r) { return isCtrl(r.group); })
                 .map(function (r) { return Number(r.mean_wpli); })
                 .filter(function (v) { return !isNaN(v); });
  var ms = rows.filter(function (r) { return isMS(r.group); })
               .map(function (r) { return Number(r.mean_wpli); })
               .filter(function (v) { return !isNaN(v); });
  if (!ctrl.length && !ms.length) {
    ctx.fillStyle = '#94a3b8';
    ctx.fillText('Sin subject_band_means', 20, 40);
    return;
  }
  var vals = ctrl.concat(ms);
  var rawMin = Math.min.apply(null, vals), rawMax = Math.max.apply(null, vals);
  // Nice axis from 0 (or below min) with 0.1 steps when plausible
  var ymin = Math.min(0, Math.floor(rawMin * 10) / 10);
  var ymax = Math.ceil(rawMax * 10) / 10;
  if (ymax - ymin < 0.2) ymax = ymin + 0.2;
  var pad = { l: 54, r: 30, t: 32, b: 64 };
  function yscale(v) {
    return pad.t + (1 - (v - ymin) / Math.max(ymax - ymin, 1e-9)) * (h - pad.t - pad.b);
  }
  var xCtrl = w * 0.32, xMS = w * 0.68;

  // Title + Y ticks
  ctx.fillStyle = '#334155';
  ctx.font = '12px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('Mean wPLI por sujeto', w / 2, 16);
  ctx.strokeStyle = '#e2e8f0';
  ctx.fillStyle = '#94a3b8';
  ctx.font = '10px sans-serif';
  ctx.textAlign = 'right';
  var step = 0.1;
  if (ymax - ymin > 0.8) step = 0.2;
  if (ymax - ymin > 1.6) step = 0.25;
  for (var t = ymin; t <= ymax + 1e-9; t += step) {
    var yy = yscale(t);
    ctx.beginPath();
    ctx.moveTo(pad.l, yy);
    ctx.lineTo(w - pad.r, yy);
    ctx.stroke();
    ctx.fillText(t.toFixed(1), pad.l - 6, yy + 3);
  }
  ctx.strokeStyle = '#cbd5e1';
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t);
  ctx.lineTo(pad.l, h - pad.b);
  ctx.lineTo(w - pad.r, h - pad.b);
  ctx.stroke();
  ctx.fillStyle = '#64748b';
  ctx.font = '12px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('Control (n=' + ctrl.length + ')', xCtrl, h - 28);
  ctx.fillText('MS (n=' + ms.length + ')', xMS, h - 28);

  var groupSummary = (S.global_band || []).find(function (r) { return r.band === S.band; });
  function drawGroup(arr, x0, color, prefix) {
    arr.forEach(function (v, i) {
      var j = ((i * 17) % 11 - 5) * 3.2;
      ctx.beginPath();
      ctx.arc(x0 + j, yscale(v), 4, 0, Math.PI * 2);
      ctx.fillStyle = color;
      ctx.globalAlpha = 0.65;
      ctx.fill();
      ctx.globalAlpha = 1;
    });
    var st = {
      m: Number(groupSummary[prefix + '_mean']),
      sem: Number(groupSummary[prefix + '_sem']),
      med: Number(groupSummary[prefix + '_median']),
      q1: Number(groupSummary[prefix + '_q1']),
      q3: Number(groupSummary[prefix + '_q3'])
    };
    if (!isNaN(st.q1) && !isNaN(st.q3)) {
      ctx.fillStyle = color;
      ctx.globalAlpha = 0.12;
      ctx.fillRect(x0 - 22, yscale(st.q3), 44, Math.max(2, yscale(st.q1) - yscale(st.q3)));
      ctx.globalAlpha = 1;
      ctx.strokeStyle = color;
      ctx.globalAlpha = 0.35;
      ctx.strokeRect(x0 - 22, yscale(st.q3), 44, Math.max(2, yscale(st.q1) - yscale(st.q3)));
      ctx.globalAlpha = 1;
    }
    ctx.strokeStyle = color;
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    ctx.moveTo(x0 - 28, yscale(st.m));
    ctx.lineTo(x0 + 28, yscale(st.m));
    ctx.stroke();
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(x0, yscale(st.m - st.sem));
    ctx.lineTo(x0, yscale(st.m + st.sem));
    ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(x0 - 6, yscale(st.m - st.sem));
    ctx.lineTo(x0 + 6, yscale(st.m - st.sem));
    ctx.moveTo(x0 - 6, yscale(st.m + st.sem));
    ctx.lineTo(x0 + 6, yscale(st.m + st.sem));
    ctx.stroke();
    if (!isNaN(st.med)) {
      ctx.setLineDash([3, 2]);
      ctx.beginPath();
      ctx.moveTo(x0 - 16, yscale(st.med));
      ctx.lineTo(x0 + 16, yscale(st.med));
      ctx.stroke();
      ctx.setLineDash([]);
    }
    ctx.fillStyle = color;
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText('μ=' + st.m.toFixed(3) + ' ± ' + st.sem.toFixed(3) + ' SEM', x0 + 34, yscale(st.m) + 4);
    ctx.fillStyle = '#64748b';
    ctx.font = '10px sans-serif';
    ctx.fillText('med=' + st.med.toFixed(3) + '  IQR=[' + st.q1.toFixed(3) + ', ' + st.q3.toFixed(3) + ']',
      x0 + 34, yscale(st.m) + 18);
  }
  drawGroup(ctrl, xCtrl, '#1d4ed8', 'ctrl');
  drawGroup(ms, xMS, '#b91c1c', 'ms');
  var note = el('mw-effect-note');
  if (note && groupSummary) {
    note.textContent = 'Mann–Whitney (banda ' + S.band + '): r_rb = ' +
      Number(groupSummary.effect_rrb).toFixed(3) + ' · P(MS>Ctrl) = ' +
      Number(groupSummary.probability_superiority).toFixed(3) +
      ' · Cohen d = ' + Number(groupSummary.effect_d_pooled).toFixed(3) +
      ' (complementario; el test primario es no paramétrico). q_global en el panel inferior.';
  } else if (note) {
    note.textContent = '';
  }
}

function renderCohenBars() {
  var gb = NMV.sortBandsPhysio(S.global_band || []);
  var note = el('cohen-q-note');
  var cis = gb.map(function (r) {
    return [Number(r.effect_d_pooled_ci_low), Number(r.effect_d_pooled_ci_high)];
  });
  NMV.barChart(el('cv-bandmeans'),
    gb.map(function (r) { return r.band; }),
    gb.map(function (r) { return Number(r.effect_d_pooled) || 0; }),
    function (v, i) {
      var sig = gb[i].is_fdr === true || String(gb[i].is_fdr) === 'true';
      if (sig) return v >= 0 ? '#b91c1c' : '#1d4ed8';
      return v >= 0 ? 'rgba(185,28,28,0.35)' : 'rgba(29,78,216,0.35)';
    },
    { height: 280, signed: true, cis: cis, nTicks: 7, niceLimit: true,
      yLabel: 'Cohen d pooled (MS − Control)' });
  var cur = gb.find(function (r) { return r.band === S.band; });
  note.textContent = 'Banda ' + S.band + ': Cohen d pooled = ' +
    Number(cur.effect_d_pooled).toFixed(3) + ' (IC95% bootstrap [' +
    Number(cur.effect_d_pooled_ci_low).toFixed(3) + ', ' +
    Number(cur.effect_d_pooled_ci_high).toFixed(3) + '])' +
    ' · q_global = ' + Number(cur.q_value).toFixed(4) +
    ' · Δ = ' + Number(cur.diff).toFixed(4) +
    ' · r_rb = ' + Number(cur.effect_rrb).toFixed(3) +
    ' · P(MS>Ctrl) = ' + Number(cur.probability_superiority).toFixed(3) +
    '. Barras intensas = is_fdr de producción (familia: 7 bandas).';
}

function drawGraphPanel() {
  var edges = S.edges || [];
  var empty = el('net-empty');
  var mode = el('edge-mode').value;
  if (!edges.length && mode === 'fdr') {
    empty.style.display = 'block';
  } else {
    empty.style.display = 'none';
  }
  var leg = el('net-legend-explore');
  if (leg) leg.style.display = mode === 'fdr' ? 'none' : '';
  NMV.drawGraph(el('cv-graph'), S, edges, { mode: mode, legend: mode !== 'fdr' });
}

function renderEdges() {
  var fmt = {
    ctrl_mean: function (v) { return Number(v).toFixed(4); },
    ms_mean: function (v) { return Number(v).toFixed(4); },
    diff: function (v) { return Number(v).toFixed(4); },
    p_value: function (v) { return Number(v).toFixed(4); },
    p_mannwhitney: function (v) { return Number(v).toFixed(4); },
    q_value: function (v) { return Number(v).toFixed(4); },
    effect_d_pooled: function (v) { return Number(v).toFixed(3); },
    effect_rrb: function (v) {
      if (v == null || v === '' || isNaN(Number(v))) return '—';
      return Number(v).toFixed(3);
    }
  };
  var cols = ['ch_a', 'ch_b', 'ctrl_mean', 'ms_mean', 'diff', 'p_mannwhitney', 'q_value', 'effect_rrb', 'effect_d_pooled'];
  var hasRrb = (S.edges || []).some(function (e) { return e.effect_rrb != null && e.effect_rrb !== ''; });
  if (!hasRrb) cols = cols.filter(function (c) { return c !== 'effect_rrb'; });
  NMV.fillTable(el('tbl-edges'), cols, S.edges || [], fmt);
}

function drawTopoPanel() {
  NMV.drawTopo(el('cv-topo'), S.spectral || [], S.positions || {});
  var note = el('spec-ch-note');
  var nW = S.n_channels != null ? S.n_channels : S.n;
  var nS = S.n_spectral_channels || 0;
  var fam = nS > 0 ? ('FDR-BH sobre ' + nS + ' canales dentro de la banda seleccionada') : 'FDR-BH por banda entre canales';
  if (nS > nW) {
    note.textContent = 'Absoluta μV² (sin log) · ' + fam +
      '. wPLI usa ' + nW + ' ch (montaje reducido); potencia usa ' + nS + ' — topografías no 1:1 con la matriz wPLI.';
  } else if (nS > 0) {
    note.textContent = 'Absoluta μV² (sin log) · ' + fam + ' · ' + nW + ' canales wPLI.';
  } else {
    note.textContent = 'Absoluta μV² (sin log). ' + fam + '.';
  }
}

function renderSpec() {
  var fmt = {
    diff: function (v) { return Number(v).toFixed(4); },
    p_value: function (v) { return Number(v).toFixed(4); },
    q_value: function (v) { return Number(v).toFixed(4); },
    effect_d_pooled: function (v) { return Number(v).toFixed(3); },
    ms_mean: function (v) { return Number(v).toFixed(3); },
    ctrl_mean: function (v) { return Number(v).toFixed(3); }
  };
  var rows = (S.spectral || []).slice().sort(function (a, b) {
    return Math.abs(b.diff) - Math.abs(a.diff);
  });
  NMV.fillTable(el('tbl-spec'),
    ['channel', 'ms_mean', 'ctrl_mean', 'diff', 'p_value', 'q_value', 'effect_d_pooled', 'n_ms', 'n_ctrl'],
    rows.slice(0, 80), fmt);
  var sn = el('spec-n-note');
  if (sn && rows.length) {
    var nms = rows[0].n_ms, nct = rows[0].n_ctrl;
    var same = rows.every(function (r) { return r.n_ms === nms && r.n_ctrl === nct; });
    sn.textContent = same
      ? ('n_MS = ' + nms + ' · n_Ctrl = ' + nct + ' (constante en todos los canales de esta banda)')
      : 'n_MS / n_Ctrl varían por canal (ver columnas).';
  }
}

function renderNetGlobal() {
  var fmt = {
    ms_mean: function (v) { return Number(v).toFixed(4); },
    ctrl_mean: function (v) { return Number(v).toFixed(4); },
    diff: function (v) { return Number(v).toFixed(4); },
    p_value: function (v) { return Number(v).toFixed(4); },
    q_value: function (v) { return Number(v).toFixed(4); },
    effect_d_pooled: function (v) { return Number(v).toFixed(3); }
  };
  NMV.fillTable(el('tbl-netg'),
    ['band', 'metric', 'ms_mean', 'ctrl_mean', 'diff', 'p_value', 'q_value', 'effect_d_pooled', 'n_ms', 'n_ctrl'],
    S.net_global || [], fmt);
}

['cv-hm-ctrl', 'cv-hm-ms', 'cv-hm-d'].forEach(function (id) {
  NMV.bindHeatmap(el(id), function (a, b) {
    return edgeLookup ? edgeLookup(a, b) : null;
  });
});

document.querySelectorAll('#cond-seg button').forEach(function (b) {
  b.addEventListener('click', function () { cond = b.dataset.cond; refresh(); });
});
el('band').addEventListener('change', refresh);
el('edge-mode').addEventListener('change', function () { syncEdgeModeUI(); refresh(); });
el('topn').addEventListener('change', refresh);
if (el('volcano-x')) {
  el('volcano-x').addEventListener('change', function () {
    if (!S) return;
    var volcanoX = el('volcano-x').value;
    NMV.drawVolcano(el('cv-volcano'), S.edge_stats || [], {
      yKey: 'p', xKey: volcanoX,
      xLabel: volcanoX === 'rrb' ? 'r_rb (MS − Control)' : 'Cohen d pooled (MS − Control)'
    });
  });
}
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

el('btn-png').addEventListener('click', function () {
  setStatus('Exportando PNG…');
  fetch('/api/export_png?cond=' + cond + '&band=' + el('band').value)
    .then(function (r) { return r.json(); })
    .then(function (j) {
      if (!j.ok) throw new Error(j.error || 'fail');
      setStatus('PNG: ' + j.file, 'ok');
    })
    .catch(function (e) { setStatus(String(e.message || e), 'err'); });
});
el('btn-csv').addEventListener('click', function () {
  setStatus('Exportando CSV…');
  var mode = el('edge-mode').value, topn = el('topn').value;
  var url = '/api/export_csv?cond=' + cond + '&band=' + el('band').value +
            '&edge_mode=' + mode + '&topn=' + topn;
  fetch(url)
    .then(function (r) { return r.json(); })
    .then(function (j) {
      if (!j.ok) throw new Error(j.error || 'fail');
      setStatus('CSV: ' + j.file, 'ok');
    })
    .catch(function (e) { setStatus(String(e.message || e), 'err'); });
});

window.addEventListener('resize', function () {
  clearTimeout(window._rt);
  window._rt = setTimeout(renderAll, 150);
});

(function init() {
  syncEdgeModeUI();
  fetch('/api/meta').then(function (r) { return r.json(); }).then(function (meta) {
    if (!meta.ok) {
      el('compat-error').style.display = 'block';
      el('compat-message').textContent = meta.error ||
        'Regenera el análisis transversal con la versión de producción correspondiente.';
      document.querySelector('.toolbar').style.display = 'none';
      document.querySelector('.tabs').style.display = 'none';
      document.querySelector('main').style.display = 'none';
      return;
    }
    cond = meta.default_cond || 'EC';
    refresh();
  });
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

function handle_request(sock, viewer::TransViewer)
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
            return _send(sock, 200, read(COMMON_JS, String);
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
            mode = get(qs, "edge_mode", "fdr")
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
            mode = get(qs, "edge_mode", "fdr")
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
    launch_transversal_viewer(; results_root=nothing, port=8781, cond="EC", open=true)

Carga `results/transversal/{eyesclosed|eyesopen}/` y sirve el visor interactivo caso-control (CLI acepta `EC|EO` como alias).
"""
function launch_transversal_viewer(; results_root=nothing, port::Int=PORT,
                                     cond::String="EC", open::Bool=true)
    viewer = load_viewer(; results_root=results_root, default_cond=uppercase(cond))
    println("NeuroMIND transversal viewer")
    println("  results: $(viewer.results_root)/transversal/")
    for (c, st) in viewer.stores
        println("  $c: $(length(st.bands)) bandas · MS=$(get(st.summary, "n_ms", "?")) · Ctrl=$(get(st.summary, "n_ctrl", "?"))")
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
    launch_transversal_viewer(; results_root=root, cond=cond)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
