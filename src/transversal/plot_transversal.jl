#!/usr/bin/env julia
# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Visor interactivo transversal (MS vs Control)
# ═══════════════════════════════════════════════════════════════
#
#  Explora results/transversal/{EC|EO}/ tras run_transversal_analysis.jl
#  · Overview KPIs · strip MS vs Control · heatmaps Ctrl/MS/Δ
#  · red FDR / top-|d| · volcano · potencia Δ · métricas de red
#
#  Uso:
#    julia --project=. src/transversal/plot_transversal.jl
#    julia --project=. src/transversal/plot_transversal.jl EO
#    # → http://127.0.0.1:8781/
#
# ───────────────────────────────────────────────────────────────
#  Fichero     src/transversal/plot_transversal.jl
#  Autor       Rafael Castro Triguero <me1catrr@uco.es>
#  Creado      25-07-2026
# ───────────────────────────────────────────────────────────────

include(joinpath(@__DIR__, "..", "viz", "GroupVizCommon.jl"))
using .GroupVizCommon
using CSV, DataFrames, CairoMakie, Sockets, Dates, Statistics, TOML

const PROJ = normpath(joinpath(@__DIR__, "..", ".."))
const HOST = "127.0.0.1"
const PORT = 8781
const COMMON_JS = joinpath(@__DIR__, "..", "viz", "group_viewer_common.js")

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
    default_cond::String
end

# ── IO ─────────────────────────────────────────────────────────

function load_condition(dir::String, cond::String)::CondStore
    bands = String[]
    channels = Dict{String,Vector{String}}()
    mat_ctrl = Dict{String,Matrix{Float64}}()
    mat_ms = Dict{String,Matrix{Float64}}()
    mat_diff = Dict{String,Matrix{Float64}}()
    stats = Dict{String,DataFrame}()
    net_ctrl = Dict{String,DataFrame}()
    net_ms = Dict{String,DataFrame}()
    net_diff = Dict{String,DataFrame}()

    for b in BAND_ORDER
        rc = read_mat_csv(joinpath(dir, "group_connectivity_control_$(b).csv"))
        rm = read_mat_csv(joinpath(dir, "group_connectivity_ms_$(b).csv"))
        rd = read_mat_csv(joinpath(dir, "group_difference_$(b).csv"))
        (rc === nothing || rm === nothing) && continue
        push!(bands, b)
        channels[b] = rc[1]
        mat_ctrl[b] = rc[2]
        mat_ms[b] = rm[2]
        mat_diff[b] = rd === nothing ? (rm[2] .- rc[2]) : rd[2]
        sdf = safe_csv(joinpath(dir, "group_statistics_$(b).csv"))
        if nrow(sdf) > 0
            for c in (:ch_a, :ch_b)
                hasproperty(sdf, c) && (sdf[!, c] = String.(sdf[!, c]))
            end
        end
        stats[b] = sdf
        net_dir = joinpath(dir, "tables", "network")
        net_ctrl[b] = safe_csv(joinpath(net_dir, "network_metrics_control_$(b).csv"))
        net_ms[b] = safe_csv(joinpath(net_dir, "network_metrics_ms_$(b).csv"))
        net_diff[b] = safe_csv(joinpath(net_dir, "network_metrics_diff_$(b).csv"))
    end

    CondStore(
        cond, dir, bands, channels, mat_ctrl, mat_ms, mat_diff, stats,
        safe_csv(joinpath(dir, "band_statistics.csv")),
        safe_csv(joinpath(dir, "global_mean_wpli_statistics.csv")),
        safe_csv(joinpath(dir, "subject_band_means.csv")),
        safe_csv(joinpath(dir, "subject_inclusion.csv")),
        parse_summary_json(joinpath(dir, "transversal_summary.json")),
        safe_csv(joinpath(dir, "tables", "spectral", "band_power_group_statistics.csv")),
        safe_csv(joinpath(dir, "tables", "network", "network_global_statistics.csv")),
        net_ctrl, net_ms, net_diff,
    )
end

function load_viewer(; results_root::Union{Nothing,String}=nothing,
                       default_cond::String="EC")::TransViewer
    root = resolve_results_root(PROJ, results_root)
    stores = Dict{String,CondStore}()
    for cond in ("EC", "EO")
        d = joinpath(root, "transversal", cond)
        isdir(d) || continue
        st = load_condition(d, cond)
        isempty(st.bands) && continue
        stores[cond] = st
    end
    isempty(stores) && error("No hay datos en $root/transversal/{EC|EO}/. Ejecuta primero run_transversal_analysis.jl")
    dc = haskey(stores, default_cond) ? default_cond : first(keys(stores))
    return TransViewer(root, stores, dc)
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
                    edge_mode::String="top_d", topn::Int=20)::String
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
    edge_cols = Symbol[:ch_a, :ch_b, :ctrl_mean, :ms_mean, :diff, :p_value,
                       :p_mannwhitney, :p_welch, :q_value, :effect_d, :n_ms, :n_ctrl]

    sm = st.subject_means
    sm_band = nrow(sm) > 0 && hasproperty(sm, :band) ?
              filter(r -> string(r.band) == band, sm) : DataFrame()
    sm_cols = Symbol[:subject_id, :group, :band, :cond, :mean_wpli]

    sp = st.spectral
    sp_band = nrow(sp) > 0 && hasproperty(sp, :band) ?
              filter(r -> string(r.band) == band, sp) : DataFrame()
    sp_cols = Symbol[:channel, :band, :ms_mean, :ctrl_mean, :diff, :p_value, :effect_d, :q_value, :n_ms, :n_ctrl]

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

    kpi = """{
      "n_ms": $(get(summ, "n_ms", 0)),
      "n_ctrl": $(get(summ, "n_ctrl", 0)),
      "n_ms_design": $(get(summ, "n_ms_design", 44)),
      "n_ctrl_design": $(get(summ, "n_ctrl_design", 40)),
      "n_included": $(get(summ, "n_included", 0)),
      "n_excluded": $(get(summ, "n_excluded", 0)),
      "n_excluded_qc": $(get(summ, "n_excluded_qc", 0)),
      "n_excluded_data": $(get(summ, "n_excluded_data", 0)),
      "n_total_sig": $n_total_sig,
      "best_band": "$(json_escape(best_band))",
      "test": "$(json_escape(string(get(summ, "test", "mannwhitney"))))",
      "design": "$(json_escape(string(get(summ, "design", "case_control_T1"))))",
      "wpli_method": "$(json_escape(string(get(summ, "wpli_method", ""))))",
      "timestamp": "$(json_escape(string(get(summ, "timestamp", ""))))"
    }"""

    bs_cols = Symbol[:band, :n_channels, :n_pairs, :n_sig, :pct_sig, :ctrl_mean, :ms_mean,
                     :diff_mean, :mean_p, :mean_d, :n_ms, :n_ctrl]
    gb_cols = Symbol[:band, :ms_mean, :ctrl_mean, :diff, :p_mannwhitney, :p_welch,
                     :effect_d, :n_ms, :n_ctrl, :q_value]
    incl_use = filter(c -> hasproperty(st.inclusion, c),
        [:subject_id, :session_id, :group, :n_bands_ok, :included,
         :excluded_reason, :qc_decision, :n_epochs_valid])
    netg_cols = Symbol[:band, :metric, :ms_mean, :ctrl_mean, :diff, :p_value,
                       :effect_d, :n_ms, :n_ctrl, :q_value]

    n_sig_band = count_fdr(edf)
    n_p_uncorr = count_uncorr(edf)
    Mc = st.mat_ctrl[band]
    Mm = st.mat_ms[band]
    shared_lim = max(maximum(abs, Mc), maximum(abs, Mm), 1e-12)

    bs_sorted = _sort_band_df(st.band_stats)
    gb_sorted = _sort_band_df(st.global_band)

    sp_use = filter(c -> nrow(sp_band) == 0 || hasproperty(sp_band, c), sp_cols)
    ng_use = filter(c -> nrow(st.net_global) == 0 || hasproperty(st.net_global, c), netg_cols)
    nc_use = filter(c -> nrow(nc) == 0 || hasproperty(nc, c), [:channel, :strength, :degree, :norm_strength])
    nm_use = filter(c -> nrow(nm) == 0 || hasproperty(nm, c), [:channel, :strength, :degree, :norm_strength])
    nd_use = filter(c -> nrow(nd) == 0 || hasproperty(nd, c),
                    [:channel, :delta_strength, :delta_degree, :delta_norm_strength])
    gb_use = filter(c -> nrow(gb_sorted) == 0 || hasproperty(gb_sorted, c), gb_cols)
    bs_use = filter(c -> nrow(bs_sorted) == 0 || hasproperty(bs_sorted, c), bs_cols)

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
  "shared_lim_ctrl_ms": $(round(shared_lim; digits=6)),
  "n_sig_band": $n_sig_band,
  "n_p_uncorr": $n_p_uncorr,
  "edges": $(df_rows_json(edges_out, edge_cols)),
  "edge_stats": $(df_rows_json(edf, edge_cols)),
  "kpi": $kpi,
  "band_stats": $(df_rows_json(bs_sorted, bs_use)),
  "global_band": $(df_rows_json(gb_sorted, gb_use)),
  "subject_means": $(df_rows_json(sm_band, sm_cols)),
  "subject_means_all": $(df_rows_json(sm, sm_cols)),
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
    # JS uses string concat (no template literals) to avoid Julia `$` interpolation.
    return """<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>NeuroMIND — Transversal MS vs Control</title>
<style>
:root {
  --bg: #f4f6f8; --card: #fff; --ink: #0f172a; --muted: #64748b;
  --line: #e2e8f0; --accent: #1e3a5f; --ms: #b91c1c; --ctrl: #1d4ed8;
  --ok: #15803d; --warn: #a16207;
}
* { box-sizing: border-box; }
body { margin:0; font-family: "IBM Plex Sans", "Segoe UI", system-ui, sans-serif;
  background: var(--bg); color: var(--ink); }
header { background: linear-gradient(120deg, #0f172a 0%, #1e3a5f 55%, #334155 100%);
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
  background:#eff6ff; color:#1e3a5f; border:1px solid #bfdbfe; }
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
.canvas-wrap canvas { display:block; background:#fff; border-radius:8px;
  border:1px solid var(--line); cursor:crosshair; }
table { width:100%; border-collapse:collapse; font-size:0.78rem; }
th, td { padding:0.35rem 0.45rem; border-bottom:1px solid var(--line); text-align:left; }
th { color:var(--muted); font-weight:600; position:sticky; top:0; background:#fff; }
.scroll { max-height:320px; overflow:auto; }
.muted { color:var(--muted); font-size:0.85rem; }
.tip { position:fixed; pointer-events:none; background:#0f172a; color:#fff; font-size:0.75rem;
  padding:0.35rem 0.55rem; border-radius:6px; z-index:50; display:none; max-width:280px;
  white-space:pre-line; }
.status { font-size:0.8rem; color:var(--muted); margin-left:auto; }
.status.ok { color:var(--ok); } .status.err { color:var(--ms); }
.note { font-size:0.85rem; line-height:1.45; color:#334155; background:#f0f4f8;
  border:1px solid #cbd5e1; border-radius:8px; padding:0.75rem; margin-top:0.75rem; }
.banner { display:none; margin:0.75rem 1.25rem 0; padding:0.65rem 0.85rem; border-radius:8px;
  background:#fffbeb; border:1px solid #fcd34d; color:#92400e; font-size:0.85rem; }
.banner.show { display:block; }
.empty-cta { padding:1.5rem; text-align:center; color:var(--muted); font-size:0.9rem; }
</style>
</head>
<body>
<header>
  <h1>NeuroMIND — Evaluación transversal MS vs Control</h1>
  <p>Caso-control T1 · EC y EO en paralelo · Mann–Whitney + FDR-BH</p>
</header>
<div class="toolbar">
  <div class="seg" id="cond-seg">
    <button type="button" data-cond="EC" class="active">EC</button>
    <button type="button" data-cond="EO">EO</button>
  </div>
  <label>Banda <select id="band"></select></label>
  <label>Edges
    <select id="edge-mode">
      <option value="top_d">Top-N |d|</option>
      <option value="fdr">Solo q&lt;0.05</option>
    </select>
  </label>
  <label>N <input id="topn" type="number" min="5" max="100" value="20" style="width:4rem;padding:0.3rem;border:1px solid var(--line);border-radius:6px"/></label>
  <span class="badge" id="badge-n">—</span>
  <button type="button" id="btn-png">Export PNG Δ</button>
  <button type="button" id="btn-csv">Export CSV edges</button>
  <span class="status" id="status">Cargando…</span>
</div>
<div class="banner" id="fdr-banner">Sin aristas FDR (q&lt;0.05) en esta banda. Vista exploratoria: usa Top-N |d|. Confirmatorio solo con FDR.</div>
<div class="tabs" id="tabs">
  <button type="button" data-tab="overview" class="active">Overview</button>
  <button type="button" data-tab="global">Comparación global</button>
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
      <div class="card"><h3>Δ mean wPLI (MS−Ctrl) por banda</h3>
        <div class="canvas-wrap"><canvas id="cv-diffbar" data-h="220"></canvas></div></div>
    </div>
    <div class="note" id="interp"></div>
    <div class="card" style="margin-top:1rem"><h3>Inclusión de sujetos</h3>
      <div class="scroll"><table id="tbl-incl"><thead></thead><tbody></tbody></table></div>
    </div>
  </section>
  <section class="panel" id="panel-global">
    <div class="card"><h3>Mean wPLI por sujeto — MS vs Control (banda seleccionada)</h3>
      <div class="canvas-wrap"><canvas id="cv-strip" data-h="320"></canvas></div></div>
    <div class="card" style="margin-top:1rem"><h3>Cohen d (mean wPLI global) por banda</h3>
      <div class="canvas-wrap"><canvas id="cv-bandmeans" data-h="260"></canvas></div>
      <p class="muted" id="cohen-q-note"></p></div>
  </section>
  <section class="panel" id="panel-conn">
    <div class="grid grid-3">
      <div class="card"><h3>Control</h3>
        <div class="canvas-wrap"><canvas id="cv-hm-ctrl"></canvas></div></div>
      <div class="card"><h3>MS</h3>
        <div class="canvas-wrap"><canvas id="cv-hm-ms"></canvas></div></div>
      <div class="card"><h3>Δ (MS−Control) · FDR outline</h3>
        <div class="canvas-wrap"><canvas id="cv-hm-d"></canvas></div></div>
    </div>
    <p class="muted" id="hm-tip-hint">Pasa el cursor sobre una celda para ver p, q y Cohen d. Escala Ctrl/MS compartida.</p>
  </section>
  <section class="panel" id="panel-net">
    <div class="grid grid-2">
      <div class="card"><h3>Grafo (edges filtrados)</h3>
        <div class="canvas-wrap"><canvas id="cv-graph"></canvas></div>
        <div class="empty-cta" id="net-empty" style="display:none">Sin edges FDR. Cambia a <b>Top-N |d|</b> para exploración.</div>
      </div>
      <div class="card"><h3>Tabla de edges</h3>
        <div class="scroll"><table id="tbl-edges"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
  </section>
  <section class="panel" id="panel-volcano">
    <div class="card"><h3>Volcano: Cohen d vs −log10(p) (todos los edges)</h3>
      <div class="canvas-wrap"><canvas id="cv-volcano" data-h="300"></canvas></div>
      <p class="muted">Rojo = q&lt;0.05 FDR. Línea discontinua = p=0.05.</p>
    </div>
  </section>
  <section class="panel" id="panel-spec">
    <div class="grid grid-2">
      <div class="card"><h3>Topo Δ band power</h3>
        <div class="canvas-wrap"><canvas id="cv-topo"></canvas></div>
        <p class="muted" id="spec-ch-note"></p></div>
      <div class="card"><h3>Canal × banda</h3>
        <div class="scroll"><table id="tbl-spec"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
  </section>
  <section class="panel" id="panel-hubs">
    <div class="grid grid-2">
      <div class="card"><h3>Strength nodal Control vs MS</h3>
        <div class="canvas-wrap"><canvas id="cv-strength" data-h="360"></canvas></div></div>
      <div class="card"><h3>Estadística global de red</h3>
        <div class="scroll"><table id="tbl-netg"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
  </section>
</main>
<div class="tip" id="tip"></div>
<script src="/static/group_viewer_common.js"></script>
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

function refresh() {
  setStatus('Cargando…');
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
    el('badge-n').textContent =
      'MS ' + (kpi.n_ms != null ? kpi.n_ms : '—') +
      ' / Ctrl ' + (kpi.n_ctrl != null ? kpi.n_ctrl : '—') +
      ' · ' + nCh + ' ch wPLI';
    var ban = el('fdr-banner');
    if ((S.n_sig_band || 0) === 0) ban.classList.add('show');
    else ban.classList.remove('show');
    renderAll();
    setStatus(cond + ' · ' + S.band + ' · ' + nCh + ' ch · ' +
              (S.n_sig_band || 0) + ' FDR · ' + (S.n_p_uncorr || 0) + ' p<0.05', 'ok');
  }).catch(function (e) { setStatus(String(e), 'err'); });
}

function renderAll() {
  if (!S) return;
  renderKPIs();
  renderBars();
  renderInterp();
  renderInclusion();
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
    n: S.n, channels: S.channels, diverging: true, fdrSet: fdrSet
  });
  drawGraphPanel();
  renderEdges();
  NMV.drawVolcano(el('cv-volcano'), S.edge_stats || []);
  drawTopoPanel();
  renderSpec();
  NMV.drawStrength(el('cv-strength'), S.net_ctrl || [], S.net_ms || [], 'Control', 'MS');
  renderNetGlobal();
}

function renderKPIs() {
  var k = S.kpi || {};
  var best = k.best_band;
  if (!best || (k.n_total_sig || 0) === 0) best = '— (sin FDR)';
  var items = [
    ['MS incluidos', k.n_ms], ['Ctrl incluidos', k.n_ctrl],
    ['Diseño MS', k.n_ms_design], ['Diseño Ctrl', k.n_ctrl_design],
    ['Excl. QC', k.n_excluded_qc], ['Excl. datos', k.n_excluded_data],
    ['Edges FDR total', k.n_total_sig], ['Mejor banda', best],
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
  if (bs) {
    t += 'Δ mean wPLI = ' + Number(bs.diff_mean).toFixed(4) +
         ' (Ctrl=' + Number(bs.ctrl_mean).toFixed(3) +
         ', MS=' + Number(bs.ms_mean).toFixed(3) +
         ', |d̄|=' + Number(bs.mean_d).toFixed(3) + '). ';
  }
  if (gb) {
    t += 'Cohen d global = ' + Number(gb.effect_d).toFixed(3) +
         ' (q=' + Number(gb.q_value).toFixed(4) + '). ';
  }
  t += 'Cohorte: ' + k.n_ms + '/' + k.n_ms_design + ' MS y ' +
       k.n_ctrl + '/' + k.n_ctrl_design + ' controles (T1). EC y EO no se promedian. ';
  if ((S.n_sig_band || 0) === 0) {
    t += '<b>Modo exploratorio</b>: sin FDR en esta banda — Top-N |d| no es confirmatorio. ';
  } else {
    t += '<b>Modo confirmatorio</b> disponible vía filtro Solo q&lt;0.05. ';
  }
  if ((k.n_ctrl || 0) < 5) {
    t += ' <b>N control bajo</b>: FDR edge-wise poco fiable.';
  }
  el('interp').innerHTML = t;
}

function renderBars() {
  var rows = NMV.sortBandsPhysio(S.band_stats || []);
  var labs = rows.map(function (r) { return r.band; });
  NMV.barChart(el('cv-nsig'), labs, rows.map(function (r) { return Number(r.n_sig) || 0; }),
    function () { return '#1e3a5f'; }, { height: 220 });
  NMV.barChart(el('cv-diffbar'), labs, rows.map(function (r) { return Number(r.diff_mean) || 0; }),
    function (v) { return v >= 0 ? '#b91c1c' : '#1d4ed8'; }, { height: 220, signed: true });
}

function renderInclusion() {
  NMV.fillTable(el('tbl-incl'),
    ['subject_id', 'group', 'included', 'qc_decision', 'n_bands_ok', 'n_epochs_valid', 'excluded_reason'],
    S.inclusion || [], {});
}

function renderStrip() {
  var cv = el('cv-strip');
  var sized = NMV.resizeCanvas(cv, { height: 320, fallbackW: 700 });
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
  var ymin = Math.min.apply(null, vals), ymax = Math.max.apply(null, vals);
  var pad = { l: 50, r: 30, t: 24, b: 56 };
  function yscale(v) {
    return pad.t + (1 - (v - ymin) / Math.max(ymax - ymin, 1e-9)) * (h - pad.t - pad.b);
  }
  var xCtrl = w * 0.32, xMS = w * 0.68;
  ctx.strokeStyle = '#e2e8f0';
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t);
  ctx.lineTo(pad.l, h - pad.b);
  ctx.lineTo(w - pad.r, h - pad.b);
  ctx.stroke();
  ctx.fillStyle = '#64748b';
  ctx.font = '12px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('Control (n=' + ctrl.length + ')', xCtrl, h - 16);
  ctx.fillText('MS (n=' + ms.length + ')', xMS, h - 16);

  function drawGroup(arr, x0, color) {
    arr.forEach(function (v, i) {
      var j = ((i * 17) % 11 - 5) * 3.2;
      ctx.beginPath();
      ctx.arc(x0 + j, yscale(v), 4, 0, Math.PI * 2);
      ctx.fillStyle = color;
      ctx.globalAlpha = 0.65;
      ctx.fill();
      ctx.globalAlpha = 1;
    });
    var st = NMV.meanSem(arr);
    if (isNaN(st.m)) return;
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
    ctx.fillStyle = color;
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'left';
    ctx.fillText('μ=' + st.m.toFixed(3) + ' ± ' + st.sem.toFixed(3), x0 + 32, yscale(st.m) + 4);
  }
  drawGroup(ctrl, xCtrl, '#1d4ed8');
  drawGroup(ms, xMS, '#b91c1c');
}

function renderCohenBars() {
  var gb = NMV.sortBandsPhysio(S.global_band || []);
  var note = el('cohen-q-note');
  if (gb.length) {
    NMV.barChart(el('cv-bandmeans'),
      gb.map(function (r) { return r.band; }),
      gb.map(function (r) { return Number(r.effect_d) || 0; }),
      function (v) { return v >= 0 ? '#b91c1c' : '#1d4ed8'; },
      { height: 260, signed: true });
    var cur = gb.find(function (r) { return r.band === S.band; });
    if (cur) {
      note.textContent = 'Banda ' + S.band + ': Cohen d = ' + Number(cur.effect_d).toFixed(3) +
        ' · q = ' + Number(cur.q_value).toFixed(4) +
        ' · Δ = ' + Number(cur.diff).toFixed(4) +
        ' (MS ' + Number(cur.ms_mean).toFixed(3) + ' − Ctrl ' + Number(cur.ctrl_mean).toFixed(3) + ')';
    } else {
      note.textContent = 'q FDR del contraste global mean wPLI (MS vs Control) por banda.';
    }
  } else {
    var rows = NMV.sortBandsPhysio(S.band_stats || []);
    NMV.barChart(el('cv-bandmeans'),
      rows.map(function (r) { return r.band; }),
      rows.map(function (r) { return Number(r.diff_mean) || 0; }),
      function (v) { return v >= 0 ? '#b91c1c' : '#1d4ed8'; },
      { height: 260, signed: true });
    note.textContent = 'Sin global_mean_wpli_statistics; mostrando Δ mean de band_statistics.';
  }
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
  NMV.drawGraph(el('cv-graph'), S, edges);
}

function renderEdges() {
  var fmt = {
    ctrl_mean: function (v) { return Number(v).toFixed(4); },
    ms_mean: function (v) { return Number(v).toFixed(4); },
    diff: function (v) { return Number(v).toFixed(4); },
    p_value: function (v) { return Number(v).toFixed(4); },
    p_mannwhitney: function (v) { return Number(v).toFixed(4); },
    q_value: function (v) { return Number(v).toFixed(4); },
    effect_d: function (v) { return Number(v).toFixed(3); }
  };
  NMV.fillTable(el('tbl-edges'),
    ['ch_a', 'ch_b', 'ctrl_mean', 'ms_mean', 'diff', 'p_mannwhitney', 'q_value', 'effect_d'],
    S.edges || [], fmt);
}

function drawTopoPanel() {
  NMV.drawTopo(el('cv-topo'), S.spectral || [], S.positions || {});
  var note = el('spec-ch-note');
  var nW = S.n_channels != null ? S.n_channels : S.n;
  var nS = S.n_spectral_channels || 0;
  if (nS > nW) {
    note.textContent = 'Nota: potencia espectral usa ' + nS +
      ' canales; wPLI/conectividad usa ' + nW +
      ' (montaje reducido). Las topografías no son 1:1 con la matriz wPLI.';
  } else if (nS > 0) {
    note.textContent = nS + ' canales espectrales · ' + nW + ' canales wPLI.';
  } else {
    note.textContent = '';
  }
}

function renderSpec() {
  var fmt = {
    diff: function (v) { return Number(v).toFixed(4); },
    p_value: function (v) { return Number(v).toFixed(4); },
    q_value: function (v) { return Number(v).toFixed(4); },
    effect_d: function (v) { return Number(v).toFixed(3); },
    ms_mean: function (v) { return Number(v).toFixed(3); },
    ctrl_mean: function (v) { return Number(v).toFixed(3); }
  };
  var rows = (S.spectral || []).slice().sort(function (a, b) {
    return Math.abs(b.diff) - Math.abs(a.diff);
  });
  NMV.fillTable(el('tbl-spec'),
    ['channel', 'ms_mean', 'ctrl_mean', 'diff', 'p_value', 'q_value', 'effect_d'],
    rows.slice(0, 80), fmt);
}

function renderNetGlobal() {
  var fmt = {
    ms_mean: function (v) { return Number(v).toFixed(4); },
    ctrl_mean: function (v) { return Number(v).toFixed(4); },
    diff: function (v) { return Number(v).toFixed(4); },
    p_value: function (v) { return Number(v).toFixed(4); },
    q_value: function (v) { return Number(v).toFixed(4); },
    effect_d: function (v) { return Number(v).toFixed(3); }
  };
  NMV.fillTable(el('tbl-netg'),
    ['band', 'metric', 'ms_mean', 'ctrl_mean', 'diff', 'p_value', 'q_value', 'effect_d', 'n_ms', 'n_ctrl'],
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
el('edge-mode').addEventListener('change', refresh);
el('topn').addEventListener('change', refresh);
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
  fetch('/api/meta').then(function (r) { return r.json(); }).then(function (meta) {
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
        elseif path == "/static/group_viewer_common.js"
            isfile(COMMON_JS) || return _send(sock, 404, "JS not found")
            return _send(sock, 200, read(COMMON_JS, String);
                         content_type="application/javascript; charset=utf-8")
        elseif path == "/api/meta"
            body = "{\"ok\":true,\"default_cond\":\"$(viewer.default_cond)\"," *
                   "\"conditions\":[$(join(["\"$c\"" for c in sort(collect(keys(viewer.stores)))], ","))]}"
            return _send(sock, 200, body; content_type="application/json")
        elseif path == "/api/state"
            cond = get(qs, "cond", viewer.default_cond)
            band = get(qs, "band", "ALPHA")
            mode = get(qs, "edge_mode", "top_d")
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
            mode = get(qs, "edge_mode", "top_d")
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

Carga `results/transversal/{EC|EO}/` y sirve el visor interactivo caso-control.
"""
function launch_transversal_viewer(; results_root=nothing, port::Int=PORT,
                                     cond::String="EC", open::Bool=true)
    viewer = load_viewer(; results_root=results_root, default_cond=uppercase(cond))
    println("NeuroMIND transversal viewer")
    println("  results: $(viewer.results_root)/transversal/")
    for (c, st) in viewer.stores
        println("  $c: $(length(st.bands)) bandas · MS=$(get(st.summary, "n_ms", "?")) · Ctrl=$(get(st.summary, "n_ctrl", "?"))")
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
