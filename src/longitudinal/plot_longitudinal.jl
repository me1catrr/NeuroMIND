#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor interactivo longitudinal (T1 → T2)
# ═══════════════════════════════════════════════════════════════
#
#  Explora results/longitudinal/{EC|EO}/ tras run_longitudinal_analysis.jl
#  · Overview KPIs · spaghetti mean±SEM · heatmaps T1/T2/Δ (escala compartida)
#  · red FDR / top-|dz| · volcano · potencia Δ · métricas de red
#
#  Uso:
#    julia --project=. src/longitudinal/plot_longitudinal.jl
#    julia --project=. src/longitudinal/plot_longitudinal.jl EO
#    # → http://127.0.0.1:8780/
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/longitudinal/plot_longitudinal.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      25-07-2026
# ───────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "viz", "GroupVizCommon.jl"))
using .GroupVizCommon
using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics, TOML

const PROJ = normpath(joinpath(@__DIR__, "..", ".."))
const HOST = "127.0.0.1"
const PORT = 8780
const COMMON_JS = joinpath(@__DIR__, "..", "viz", "group_viewer_common.js")

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
    default_cond::String
end

# ── Loaders ────────────────────────────────────────────────────

function load_condition(dir::String, cond::String)::CondStore
    bands = String[]
    channels = Dict{String,Vector{String}}()
    mat_t1 = Dict{String,Matrix{Float64}}()
    mat_t2 = Dict{String,Matrix{Float64}}()
    mat_diff = Dict{String,Matrix{Float64}}()
    stats = Dict{String,DataFrame}()
    net_t1 = Dict{String,DataFrame}()
    net_t2 = Dict{String,DataFrame}()
    net_delta = Dict{String,DataFrame}()

    for b in BAND_ORDER
        r1 = read_mat_csv(joinpath(dir, "longitudinal_connectivity_t1_$(b).csv"))
        r2 = read_mat_csv(joinpath(dir, "longitudinal_connectivity_t2_$(b).csv"))
        rd = read_mat_csv(joinpath(dir, "longitudinal_difference_$(b).csv"))
        (r1 === nothing || r2 === nothing) && continue
        push!(bands, b)
        channels[b] = r1[1]
        mat_t1[b] = r1[2]
        mat_t2[b] = r2[2]
        mat_diff[b] = rd === nothing ? (r2[2] .- r1[2]) : rd[2]
        sdf = safe_csv(joinpath(dir, "longitudinal_statistics_$(b).csv"))
        if nrow(sdf) > 0
            for c in (:ch_a, :ch_b)
                hasproperty(sdf, c) && (sdf[!, c] = String.(sdf[!, c]))
            end
        end
        stats[b] = sdf
        net_dir = joinpath(dir, "tables", "network")
        net_t1[b] = safe_csv(joinpath(net_dir, "network_metrics_t1_$(b).csv"))
        net_t2[b] = safe_csv(joinpath(net_dir, "network_metrics_t2_$(b).csv"))
        net_delta[b] = safe_csv(joinpath(net_dir, "network_metrics_delta_$(b).csv"))
    end

    CondStore(
        cond, dir, bands, channels, mat_t1, mat_t2, mat_diff, stats,
        safe_csv(joinpath(dir, "band_statistics_longitudinal.csv")),
        safe_csv(joinpath(dir, "subject_band_means.csv")),
        safe_csv(joinpath(dir, "paired_subjects.csv")),
        parse_summary_json(joinpath(dir, "longitudinal_summary.json")),
        safe_csv(joinpath(dir, "tables", "spectral", "band_power_delta_statistics.csv")),
        safe_csv(joinpath(dir, "tables", "network", "network_global_statistics.csv")),
        net_t1, net_t2, net_delta,
    )
end

function load_viewer(; results_root::Union{Nothing,String}=nothing,
                       default_cond::String="EC")::LongViewer
    root = resolve_results_root(PROJ, results_root)
    stores = Dict{String,CondStore}()
    for cond in ("EC", "EO")
        d = joinpath(root, "longitudinal", cond)
        isdir(d) || continue
        st = load_condition(d, cond)
        isempty(st.bands) && continue
        stores[cond] = st
    end
    isempty(stores) && error("No hay datos en $root/longitudinal/{EC|EO}/. Ejecuta primero run_longitudinal_analysis.jl")
    dc = haskey(stores, default_cond) ? default_cond : first(keys(stores))
    return LongViewer(root, stores, dc)
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

function _honest_best_from_summary(summ::Dict{String,Any}, band_stats::DataFrame)::String
    n_total = try
        Int(round(Float64(get(summ, "n_total_sig", 0))))
    catch
        0
    end
    n_total == 0 && return ""
    # Prefer honest pick from band_stats if available
    if nrow(band_stats) > 0 && hasproperty(band_stats, :n_sig) && hasproperty(band_stats, :band)
        n_sigs = Int.(round.(Float64.(band_stats.n_sig)))
        if sum(n_sigs) > 0
            diffs = hasproperty(band_stats, :diff_mean) ?
                    abs.(Float64.(coalesce.(band_stats.diff_mean, 0.0))) :
                    zeros(Float64, nrow(band_stats))
            scores = n_sigs .* 1000.0 .+ diffs
            return string(band_stats.band[argmax(scores)])
        end
    end
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
            edges_out = filter(r -> Float64(r.q_value) < 0.05, edf)
            edges_out = sort(edges_out, :effect_d; by=abs, rev=true)
        else
            edges_out = sort(edf, :effect_d; by=abs, rev=true)
            edges_out = first(edges_out, min(topn, nrow(edges_out)))
        end
    end
    edge_cols = Symbol[:ch_a, :ch_b, :t1_mean, :t2_mean, :diff, :p_value, :q_value,
                       :effect_d, :n]
    edge_cols_full = Symbol[:ch_a, :ch_b, :t1_mean, :t2_mean, :diff, :p_value, :q_value,
                            :effect_d, :n]
    # include p_parametric in filtered edges if present
    if hasproperty(edf, :p_parametric)
        edge_cols = Symbol[:ch_a, :ch_b, :t1_mean, :t2_mean, :diff, :p_value, :q_value,
                           :p_parametric, :effect_d, :n]
    end

    sm = st.subject_means
    sm_band = nrow(sm) > 0 && hasproperty(sm, :band) ?
              filter(r -> string(r.band) == band, sm) : DataFrame()
    sm_cols = Symbol[:subject_id, :timepoint, :band, :mean_wpli]

    sp = st.spectral
    sp_band = nrow(sp) > 0 && hasproperty(sp, :band) ?
              filter(r -> string(r.band) == band, sp) : DataFrame()
    sp_cols = Symbol[:channel, :band, :t1_mean, :t2_mean, :diff, :p_value, :effect_d, :q_value, :n]

    nt1 = get(st.net_t1, band, DataFrame())
    nt2 = get(st.net_t2, band, DataFrame())
    nd  = get(st.net_delta, band, DataFrame())
    for df in (nt1, nt2, nd)
        nrow(df) > 0 && hasproperty(df, :channel) && (df.channel = String.(df.channel))
    end

    summ = st.summary
    band_stats_sorted = _sort_band_stats(st.band_stats)
    best_band = _honest_best_from_summary(summ, band_stats_sorted)
    n_total_sig = try
        Int(round(Float64(get(summ, "n_total_sig", 0))))
    catch
        0
    end

    kpi = """{
      "n_paired": $(get(summ, "n_paired", 0)),
      "n_paired_design": $(get(summ, "n_paired_design", 30)),
      "n_excluded": $(get(summ, "n_excluded", get(summ, "n_loss", 0))),
      "n_excluded_qc": $(get(summ, "n_excluded_qc", get(summ, "n_qc_loss", 0))),
      "n_excluded_data": $(get(summ, "n_excluded_data", get(summ, "n_data_loss", 0))),
      "n_total_sig": $n_total_sig,
      "best_band": "$(json_escape(best_band))",
      "test": "$(json_escape(string(get(summ, "test", "wilcoxon"))))",
      "wpli_method": "$(json_escape(string(get(summ, "wpli_method", ""))))",
      "timestamp": "$(json_escape(string(get(summ, "timestamp", ""))))"
    }"""

    bs_cols = Symbol[:band, :n_channels, :n_pairs, :n_sig, :pct_sig, :t1_mean, :t2_mean,
                     :diff_mean, :mean_p, :mean_d, :n_subjects]
    paired_use = filter(c -> hasproperty(st.paired, c),
        [:subject_id, :session_t1, :session_t2, :n_bands_ok, :included,
         :excluded_reason, :qc_t1, :qc_t2])
    netg_cols = Symbol[:band, :metric, :t1_mean, :t2_mean, :diff, :p_value, :effect_d, :n, :q_value]

    n_sig_band = count_fdr(edf)
    n_p_uncorr = count_uncorr(edf)
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
  "shared_lim_t12": $(round(shared_lim_t12; digits=6)),
  "n_sig_band": $n_sig_band,
  "n_p_uncorr": $n_p_uncorr,
  "edge_mode": "$(json_escape(edge_mode))",
  "edges": $(df_rows_json(edges_out, edge_use)),
  "edge_stats": $(df_rows_json(edf, estat_use)),
  "kpi": $kpi,
  "band_stats": $(df_rows_json(band_stats_sorted, bs_use)),
  "subject_means": $(df_rows_json(sm_band, sm_cols)),
  "subject_means_all": $(df_rows_json(sm, sm_cols)),
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
        out_df = filter(r -> Float64(r.q_value) < 0.05, edf)
        out_df = sort(out_df, :effect_d; by=abs, rev=true)
    else
        out_df = first(sort(edf, :effect_d; by=abs, rev=true), min(topn, nrow(edf)))
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
  --bg: #f4f6f8; --card: #fff; --ink: #0f172a; --muted: #64748b;
  --line: #e2e8f0; --accent: #0f766e; --ms: #b91c1c; --t1: #1d4ed8; --t2: #c2410c;
  --ok: #15803d; --warn: #a16207;
}
* { box-sizing: border-box; }
body { margin:0; font-family: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
  background: var(--bg); color: var(--ink); }
header { background: linear-gradient(120deg, #0f172a 0%, #134e4a 60%, #0f766e 100%);
  color:#fff; padding: 1rem 1.25rem 0.85rem; }
header h1 { margin:0; font-size:1.25rem; font-weight:600; letter-spacing:-0.02em; }
header p { margin:0.25rem 0 0; opacity:0.85; font-size:0.85rem; }
.toolbar { display:flex; flex-wrap:wrap; gap:0.6rem; align-items:center;
  padding:0.75rem 1.25rem; background:var(--card); border-bottom:1px solid var(--line);
  position:sticky; top:0; z-index:20; }
.toolbar label { font-size:0.75rem; color:var(--muted); margin-right:0.2rem; }
.toolbar select, .toolbar button, .seg button {
  font: inherit; font-size:0.85rem; padding:0.35rem 0.65rem; border-radius:6px;
  border:1px solid var(--line); background:#fff; cursor:pointer; }
.seg { display:inline-flex; border:1px solid var(--line); border-radius:8px; overflow:hidden; }
.seg button { border:0; border-radius:0; background:#f8fafc; }
.seg button.active { background:var(--accent); color:#fff; }
.badge { font-size:0.75rem; padding:0.25rem 0.55rem; border-radius:999px;
  background:#ecfdf5; color:#065f46; border:1px solid #a7f3d0; }
.tabs { display:flex; gap:0.25rem; padding:0.5rem 1.25rem 0; overflow-x:auto; }
.tabs button { border:1px solid transparent; background:transparent; padding:0.5rem 0.85rem;
  font:inherit; font-size:0.85rem; color:var(--muted); cursor:pointer; border-radius:8px 8px 0 0; }
.tabs button.active { background:var(--card); color:var(--ink); border-color:var(--line);
  border-bottom-color:var(--card); font-weight:600; }
main { padding:0 1.25rem 2rem; }
.panel { display:none; background:var(--card); border:1px solid var(--line);
  border-radius:0 10px 10px 10px; padding:1rem; min-height:420px; }
.panel.active { display:block; }
.grid { display:grid; gap:1rem; }
.grid-3 { grid-template-columns: repeat(3, 1fr); }
.grid-2 { grid-template-columns: 1fr 1fr; }
@media (max-width: 1100px) { .grid-3, .grid-2 { grid-template-columns: 1fr; } }
.card { border:1px solid var(--line); border-radius:10px; padding:0.75rem; background:#fafbfc; }
.card h3 { margin:0 0 0.5rem; font-size:0.9rem; }
.kpi-row { display:flex; flex-wrap:wrap; gap:0.75rem; margin-bottom:1rem; }
.kpi { min-width:110px; padding:0.65rem 0.85rem; border-radius:10px; background:#f8fafc;
  border:1px solid var(--line); }
.kpi .v { font-size:1.35rem; font-weight:650; letter-spacing:-0.02em; }
.kpi .l { font-size:0.7rem; color:var(--muted); text-transform:uppercase; letter-spacing:0.04em; }
.canvas-wrap { width:100%; }
.canvas-wrap canvas { display:block; background:#fff; border-radius:8px; border:1px solid var(--line); cursor:crosshair; }
table { width:100%; border-collapse:collapse; font-size:0.78rem; }
th, td { padding:0.35rem 0.45rem; border-bottom:1px solid var(--line); text-align:left; }
th { color:var(--muted); font-weight:600; position:sticky; top:0; background:#fff; }
.scroll { max-height:320px; overflow:auto; }
.muted { color:var(--muted); font-size:0.85rem; }
.tip { position:fixed; pointer-events:none; background:#0f172a; color:#fff; font-size:0.75rem;
  padding:0.35rem 0.55rem; border-radius:6px; z-index:50; display:none; max-width:280px; white-space:pre-line; }
.status { font-size:0.8rem; color:var(--muted); margin-left:auto; }
.status.ok { color:var(--ok); } .status.err { color:var(--ms); }
.note { font-size:0.85rem; line-height:1.45; color:#334155; background:#f0fdfa;
  border:1px solid #99f6e4; border-radius:8px; padding:0.75rem; margin-top:0.75rem; }
.banner { display:none; margin:0.75rem 1.25rem 0; padding:0.7rem 0.9rem; border-radius:8px;
  background:#fffbeb; border:1px solid #fcd34d; color:#92400e; font-size:0.85rem; line-height:1.4; }
.banner.show { display:block; }
.empty-cta { text-align:center; padding:2rem 1rem; color:var(--muted); }
.empty-cta button { margin-top:0.75rem; font:inherit; padding:0.45rem 0.9rem; border-radius:6px;
  border:1px solid var(--line); background:var(--accent); color:#fff; cursor:pointer; }
</style>
<script src="/static/group_viewer_common.js"></script>
</head>
<body>
<header>
  <h1>NeuroMIND — Evaluación longitudinal T1 → T2</h1>
  <p>Solo EM con pares completos · EC y EO en paralelo · Wilcoxon + FDR-BH</p>
</header>
<div class="toolbar">
  <div class="seg" id="cond-seg">
    <button type="button" data-cond="EC" class="active">EC</button>
    <button type="button" data-cond="EO">EO</button>
  </div>
  <label>Banda <select id="band"></select></label>
  <label>Edges
    <select id="edge-mode">
      <option value="top_dz" selected>Top-N |dz|</option>
      <option value="fdr">Solo q&lt;0.05</option>
    </select>
  </label>
  <label>N <input id="topn" type="number" min="5" max="100" value="20" style="width:4rem;padding:0.3rem;border:1px solid var(--line);border-radius:6px"/></label>
  <span class="badge" id="badge-n">—</span>
  <button type="button" id="btn-png">Export PNG Δ</button>
  <button type="button" id="btn-csv">Export CSV edges</button>
  <span class="status" id="status">Cargando…</span>
</div>
<div id="fdr-banner" class="banner"></div>
<div class="tabs" id="tabs">
  <button type="button" data-tab="overview" class="active">Overview</button>
  <button type="button" data-tab="global">Cambio global</button>
  <button type="button" data-tab="conn">Conectividad</button>
  <button type="button" data-tab="net">Red significativa</button>
  <button type="button" data-tab="volcano">Volcano</button>
  <button type="button" data-tab="spec">Espectros</button>
  <button type="button" data-tab="hubs">Red / hubs</button>
</div>
<main>
  <section class="panel active" id="panel-overview">
    <div class="kpi-row" id="kpis"></div>
    <div class="grid grid-2">
      <div class="card"><h3>Edges FDR (q&lt;0.05) por banda</h3>
        <div class="canvas-wrap"><canvas id="cv-nsig" data-h="220"></canvas></div></div>
      <div class="card"><h3>Δ mean wPLI (T2−T1) por banda</h3>
        <div class="canvas-wrap"><canvas id="cv-diffbar" data-h="220"></canvas></div></div>
    </div>
    <div class="note" id="interp"></div>
    <div class="card" style="margin-top:1rem"><h3>Inclusión de pares</h3>
      <div class="scroll"><table id="tbl-paired"><thead></thead><tbody></tbody></table></div>
    </div>
  </section>
  <section class="panel" id="panel-global">
    <div class="card"><h3>Mean wPLI por sujeto (T1 → T2) — media ± SEM</h3>
      <div class="canvas-wrap"><canvas id="cv-spaghetti" data-h="320"></canvas></div></div>
    <div class="card" style="margin-top:1rem"><h3>Δ mean wPLI todas las bandas</h3>
      <div class="canvas-wrap"><canvas id="cv-bandmeans" data-h="260"></canvas></div></div>
  </section>
  <section class="panel" id="panel-conn">
    <div class="grid grid-3">
      <div class="card"><h3>T1</h3><div class="canvas-wrap"><canvas id="cv-hm-t1"></canvas></div></div>
      <div class="card"><h3>T2</h3><div class="canvas-wrap"><canvas id="cv-hm-t2"></canvas></div></div>
      <div class="card"><h3>Δ (T2−T1) · contorno FDR</h3><div class="canvas-wrap"><canvas id="cv-hm-d"></canvas></div></div>
    </div>
    <p class="muted" id="hm-tip-hint">Pasa el cursor sobre una celda para ver p, q y Cohen dz (edge_stats completo).</p>
  </section>
  <section class="panel" id="panel-net">
    <div id="net-empty" class="empty-cta" style="display:none">
      <p>Sin aristas FDR (q&lt;0.05) en esta banda.</p>
      <p class="muted">Modo confirmatorio vacío — cambia a Top-N |dz| para exploración de efecto.</p>
      <button type="button" id="btn-to-topn">Usar Top-N |dz|</button>
    </div>
    <div id="net-content" class="grid grid-2">
      <div class="card"><h3>Grafo (edges filtrados)</h3>
        <div class="canvas-wrap"><canvas id="cv-graph"></canvas></div></div>
      <div class="card"><h3>Tabla de edges</h3>
        <div class="scroll"><table id="tbl-edges"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
  </section>
  <section class="panel" id="panel-volcano">
    <div class="card"><h3>Volcano — Cohen dz vs −log10(p) · rojo = q&lt;0.05</h3>
      <div class="canvas-wrap"><canvas id="cv-volcano" data-h="320"></canvas></div>
      <p class="muted">Usa edge_stats completo (todas las aristas de la banda).</p>
    </div>
  </section>
  <section class="panel" id="panel-spec">
    <div class="grid grid-2">
      <div class="card"><h3>Topo Δ band power</h3>
        <div class="canvas-wrap"><canvas id="cv-topo"></canvas></div></div>
      <div class="card"><h3>Canal × banda</h3>
        <div class="scroll"><table id="tbl-spec"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
  </section>
  <section class="panel" id="panel-hubs">
    <div class="grid grid-2">
      <div class="card"><h3>Strength nodal T1 vs T2</h3>
        <div class="canvas-wrap"><canvas id="cv-strength" data-h="360"></canvas></div></div>
      <div class="card"><h3>Estadística global de red</h3>
        <div class="scroll"><table id="tbl-netg"><thead></thead><tbody></tbody></table></div>
      </div>
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

  function updateBanner() {
    var ban = gid('fdr-banner');
    var nSig = (S && S.n_sig_band) || 0;
    var nUnc = (S && S.n_p_uncorr) || 0;
    if (nSig === 0) {
      ban.className = 'banner show';
      ban.innerHTML = 'Sin aristas FDR (q&lt;0.05) en <b>' + ((S && S.band) || '—') +
        '</b>. Hay ' + nUnc + ' con p&lt;0.05 sin corregir. ' +
        'El modo por defecto <b>Top-N |dz|</b> es exploratorio; FDR es confirmatorio.';
    } else {
      ban.className = 'banner';
      ban.textContent = '';
    }
  }

  async function refresh() {
    setStatus('Cargando…');
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
    gid('badge-n').textContent = 'Pares ' + (kpi.n_paired != null ? kpi.n_paired : '—') +
      ' / diseño ' + (kpi.n_paired_design != null ? kpi.n_paired_design : 30) +
      ' · ' + nCh + ' ch';
    updateBanner();
    renderAll();
    setStatus(cond + ' · ' + S.band + ' · ' + nCh + ' ch · FDR ' + (S.n_sig_band || 0) +
      ' · p&lt;0.05 ' + (S.n_p_uncorr || 0), 'ok');
  }

  function renderAll() {
    if (!S) return;
    renderKPIs();
    renderBars();
    renderInterp();
    renderPaired();
    renderSpaghetti();
    renderBandMeans();
    renderHeatmaps();
    renderNet();
    renderVolcano();
    renderTopo();
    renderSpec();
    renderStrength();
    renderNetGlobal();
  }

  function renderKPIs() {
    var k = S.kpi || {};
    var best = k.best_band;
    if (!best) best = '— (sin FDR)';
    var items = [
      ['Pares', k.n_paired],
      ['Diseño', k.n_paired_design],
      ['Excl. QC', k.n_excluded_qc],
      ['Excl. datos', k.n_excluded_data],
      ['Edges FDR total', k.n_total_sig],
      ['Mejor banda', best]
    ];
    gid('kpis').innerHTML = items.map(function (it) {
      return '<div class="kpi"><div class="v">' + (it[1] != null ? it[1] : '—') +
        '</div><div class="l">' + it[0] + '</div></div>';
    }).join('');
  }

  function renderInterp() {
    var k = S.kpi || {};
    var rows = NMV.sortBandsPhysio(S.band_stats || []);
    var bs = null;
    for (var i = 0; i < rows.length; i++) {
      if (rows[i].band === S.band) { bs = rows[i]; break; }
    }
    var mode = gid('edge-mode').value;
    var t = 'Condición <b>' + cond + '</b>, banda <b>' + S.band + '</b>: ';
    t += (S.n_sig_band || 0) + ' aristas FDR (q&lt;0.05) y ' + (S.n_p_uncorr || 0) +
      ' con p&lt;0.05 sin corregir. ';
    if (bs) {
      t += 'Δ mean wPLI = ' + Number(bs.diff_mean).toFixed(4) +
        ' (T1=' + Number(bs.t1_mean).toFixed(3) +
        ', T2=' + Number(bs.t2_mean).toFixed(3) +
        ', |d̄|=' + Number(bs.mean_d).toFixed(3) + '). ';
    }
    t += 'Cohorte: ' + (k.n_paired || '?') + '/' + (k.n_paired_design || 30) +
      ' pares EM · ' + (S.n_channels || S.n || '?') + ' canales. ';
    if (mode === 'fdr') {
      t += '<b>Capa confirmatoria</b> (solo q&lt;0.05). ';
    } else {
      t += '<b>Capa exploratoria</b> (Top-N |dz|): útil para hipótesis, no sustituye FDR. ';
    }
    if ((S.n_sig_band || 0) === 0) {
      t += 'Sin FDR en esta banda: interpreta Top-N con cautela.';
    }
    gid('interp').innerHTML = t;
  }

  function renderBars() {
    var rows = NMV.sortBandsPhysio(S.band_stats || []);
    var labs = rows.map(function (r) { return r.band; });
    NMV.barChart(gid('cv-nsig'), labs, rows.map(function (r) { return Number(r.n_sig) || 0; }),
      function () { return '#0f766e'; }, { height: 220 });
    NMV.barChart(gid('cv-diffbar'), labs, rows.map(function (r) { return Number(r.diff_mean) || 0; }),
      function (v) { return v >= 0 ? '#c2410c' : '#1d4ed8'; }, { height: 220, signed: true });
  }

  function renderPaired() {
    NMV.fillTable(gid('tbl-paired'),
      ['subject_id', 'included', 'qc_t1', 'qc_t2', 'n_bands_ok', 'excluded_reason'],
      S.paired || [], {});
  }

  function renderSpaghetti() {
    var cv = gid('cv-spaghetti');
    var rsz = NMV.resizeCanvas(cv, { height: 320, fallbackW: 700 });
    var ctx = rsz.ctx, w = rsz.w, h = rsz.h;
    var rows = S.subject_means || [];
    var by = {};
    rows.forEach(function (r) {
      var sid = String(r.subject_id);
      by[sid] = by[sid] || {};
      by[sid][String(r.timepoint)] = Number(r.mean_wpli);
    });
    var pairs = Object.keys(by).map(function (sid) {
      return { sid: sid, T1: by[sid].T1, T2: by[sid].T2 };
    }).filter(function (p) { return p.T1 != null && p.T2 != null && !isNaN(p.T1) && !isNaN(p.T2); });
    if (!pairs.length) {
      ctx.fillStyle = '#94a3b8';
      ctx.fillText('Sin subject_band_means', 20, 40);
      return;
    }
    var vals = [];
    pairs.forEach(function (p) { vals.push(p.T1, p.T2); });
    var ymin = Math.min.apply(null, vals), ymax = Math.max.apply(null, vals);
    var pad = { l: 50, r: 30, t: 24, b: 44 };
    var yscale = function (v) {
      return pad.t + (1 - (v - ymin) / Math.max(ymax - ymin, 1e-9)) * (h - pad.t - pad.b);
    };
    var x1 = w * 0.3, x2 = w * 0.7;
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
    var ms1 = NMV.meanSem(pairs.map(function (p) { return p.T1; }));
    var ms2 = NMV.meanSem(pairs.map(function (p) { return p.T2; }));
    ctx.strokeStyle = '#0f766e';
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    ctx.moveTo(x1, yscale(ms1.m));
    ctx.lineTo(x2, yscale(ms2.m));
    ctx.stroke();
    // SEM whiskers
    [[x1, ms1], [x2, ms2]].forEach(function (pair) {
      var x = pair[0], ms = pair[1];
      if (!ms.sem) return;
      ctx.strokeStyle = '#0f766e';
      ctx.lineWidth = 1.5;
      ctx.beginPath();
      ctx.moveTo(x, yscale(ms.m - ms.sem));
      ctx.lineTo(x, yscale(ms.m + ms.sem));
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(x - 6, yscale(ms.m - ms.sem));
      ctx.lineTo(x + 6, yscale(ms.m - ms.sem));
      ctx.moveTo(x - 6, yscale(ms.m + ms.sem));
      ctx.lineTo(x + 6, yscale(ms.m + ms.sem));
      ctx.stroke();
    });
    ctx.fillStyle = '#0f766e';
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText('Media ± SEM (N=' + pairs.length + ')', pad.l, 14);
  }

  function renderBandMeans() {
    var rows = NMV.sortBandsPhysio(S.band_stats || []);
    NMV.barChart(gid('cv-bandmeans'),
      rows.map(function (r) { return r.band; }),
      rows.map(function (r) { return Number(r.t2_mean) - Number(r.t1_mean); }),
      function (v) { return v >= 0 ? '#c2410c' : '#1d4ed8'; },
      { height: 260, signed: true });
  }

  function renderHeatmaps() {
    var lim = S.shared_lim_t12 || null;
    var fdrSet = NMV.buildFdrSet(S.edge_stats || []);
    var base = { n: S.n, channels: S.channels, lim: lim };
    NMV.drawHeatmap(gid('cv-hm-t1'), S.matrix_t1, Object.assign({}, base, { diverging: false }));
    NMV.drawHeatmap(gid('cv-hm-t2'), S.matrix_t2, Object.assign({}, base, { diverging: false }));
    NMV.drawHeatmap(gid('cv-hm-d'), S.matrix_diff, {
      n: S.n, channels: S.channels, diverging: true, fdrSet: fdrSet
    });
  }

  function renderNet() {
    var mode = gid('edge-mode').value;
    var edges = S.edges || [];
    var empty = mode === 'fdr' && edges.length === 0;
    gid('net-empty').style.display = empty ? 'block' : 'none';
    gid('net-content').style.display = empty ? 'none' : '';
    if (empty) return;
    NMV.drawGraph(gid('cv-graph'), S, edges);
    var fmt = {
      t1_mean: function (v) { return Number(v).toFixed(4); },
      t2_mean: function (v) { return Number(v).toFixed(4); },
      diff: function (v) { return Number(v).toFixed(4); },
      p_value: function (v) { return Number(v).toFixed(4); },
      q_value: function (v) { return Number(v).toFixed(4); },
      effect_d: function (v) { return Number(v).toFixed(3); }
    };
    NMV.fillTable(gid('tbl-edges'),
      ['ch_a', 'ch_b', 't1_mean', 't2_mean', 'diff', 'p_value', 'q_value', 'effect_d'],
      edges, fmt);
  }

  function renderVolcano() {
    NMV.drawVolcano(gid('cv-volcano'), S.edge_stats || []);
  }

  function renderTopo() {
    NMV.drawTopo(gid('cv-topo'), S.spectral || [], S.positions || {});
  }

  function renderSpec() {
    var fmt = {
      diff: function (v) { return Number(v).toFixed(4); },
      p_value: function (v) { return Number(v).toFixed(4); },
      q_value: function (v) { return Number(v).toFixed(4); },
      effect_d: function (v) { return Number(v).toFixed(3); }
    };
    var rows = (S.spectral || []).slice().sort(function (a, b) {
      return Math.abs(Number(b.diff) || 0) - Math.abs(Number(a.diff) || 0);
    });
    NMV.fillTable(gid('tbl-spec'),
      ['channel', 'diff', 'p_value', 'q_value', 'effect_d', 'n'],
      rows.slice(0, 80), fmt);
  }

  function renderStrength() {
    NMV.drawStrength(gid('cv-strength'), S.net_t1 || [], S.net_t2 || [], 'T1', 'T2');
  }

  function renderNetGlobal() {
    var fmt = {
      t1_mean: function (v) { return Number(v).toFixed(4); },
      t2_mean: function (v) { return Number(v).toFixed(4); },
      diff: function (v) { return Number(v).toFixed(4); },
      p_value: function (v) { return Number(v).toFixed(4); },
      q_value: function (v) { return Number(v).toFixed(4); },
      effect_d: function (v) { return Number(v).toFixed(3); }
    };
    NMV.fillTable(gid('tbl-netg'),
      ['band', 'metric', 't1_mean', 't2_mean', 'diff', 'p_value', 'q_value', 'effect_d', 'n'],
      S.net_global || [], fmt);
  }

  // bind heatmaps once
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
  gid('edge-mode').addEventListener('change', refresh);
  gid('topn').addEventListener('change', refresh);
  gid('btn-to-topn').addEventListener('click', function () {
    gid('edge-mode').value = 'top_dz';
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
        elseif path == "/static/group_viewer_common.js"
            isfile(COMMON_JS) || return _send(sock, 404, "JS not found")
            return _send_bytes(sock, 200, read(COMMON_JS);
                               content_type="application/javascript; charset=utf-8")
        elseif path == "/api/meta"
            body = "{\"ok\":true,\"default_cond\":\"$(viewer.default_cond)\"," *
                   "\"conditions\":[$(join(["\"$c\"" for c in sort(collect(keys(viewer.stores)))], ","))]}"
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

Carga `results/longitudinal/{EC|EO}/` y sirve el visor interactivo.
"""
function launch_longitudinal_viewer(; results_root=nothing, port::Int=PORT,
                                      cond::String="EC", open::Bool=true)
    viewer = load_viewer(; results_root=results_root, default_cond=uppercase(cond))
    println("NeuroMIND longitudinal viewer")
    println("  results: $(viewer.results_root)/longitudinal/")
    for (c, st) in viewer.stores
        println("  $c: $(length(st.bands)) bandas · pares=$(get(st.summary, "n_paired", "?"))")
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
