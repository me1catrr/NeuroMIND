# NeuroMIND/src/report/HTMLReport.jl
# Generación de informes HTML estáticos con figuras embebidas (base64).
# Cada informe es un archivo HTML autocontenido que se abre en cualquier navegador.
# Diseño similar a BrainVision Analyzer Report, adaptado al pipeline del estudio BRAIN.

using Base64, CairoMakie

# ─── Informe por sujeto ───────────────────────────────────────

"""
    generate_report(subj, cfg; condition="EC", verbose=true)

Genera un informe HTML completo por sujeto con todas las etapas del pipeline.
Guarda en: results/{subject_id}/{session_id}/reports/report_{condition}.html
"""
function generate_report(subj::Subject, cfg::PipelineConfig;
                         condition::String="EC", verbose::Bool=true)
    for (sess_id, sess) in subj.sessions
        out_dir = subject_results_dir(cfg, subj.id, sess_id)
        mkpath(joinpath(out_dir, "reports"))
        path = joinpath(out_dir, "reports", "report_$(condition).html")
        html = _build_subject_html(subj, sess, condition, cfg, out_dir)
        write(path, html)
        verbose && println("  ✓ Informe: $path")
    end
end

"""
    generate_group_report(subjects, cfg; session_id="T1", condition="EC", band="ALPHA")

Informe HTML de comparación grupal MS vs Control.
"""
function generate_group_report(
    subjects::Vector{Subject},
    cfg::PipelineConfig;
    session_id::String = "T1",
    condition::String  = "EC",
    band::String       = "ALPHA",
    verbose::Bool      = true
)
    out_dir  = joinpath(results_dir(cfg), "transversal", "reports")
    mkpath(out_dir)
    path     = joinpath(out_dir, "group_$(session_id)_$(condition)_$(band).html")
    html     = _build_group_html(subjects, session_id, condition, band, cfg)
    write(path, html)
    verbose && println("  ✓ Informe grupal: $path")
end

"""
    generate_longitudinal_report(subj, cfg; condition="EC", band="ALPHA")

Informe HTML de evolución longitudinal T1→T2 para un paciente MS.
"""
function generate_longitudinal_report(
    subj::Subject,
    cfg::PipelineConfig;
    condition::String = "EC",
    band::String      = "ALPHA",
    verbose::Bool     = true
)
    subj.group == "MS" || (@warn "Solo pacientes MS tienen seguimiento longitudinal"; return)
    out_dir = joinpath(results_dir(cfg), subj.id, "longitudinal")
    mkpath(out_dir)
    path = joinpath(out_dir, "longitudinal_$(condition)_$(band).html")
    html = _build_longitudinal_html(subj, condition, band, cfg)
    write(path, html)
    verbose && println("  ✓ Informe longitudinal: $path")
end

# ─── Constructores HTML ───────────────────────────────────────

function _build_subject_html(subj, sess, condition, cfg, out_dir)
    band_names = sort(collect(keys(cfg.bands)))

    sections = String[]

    # 1. Cabecera del sujeto
    push!(sections, _section("Información del sujeto", _subject_info_table(subj, sess)))

    # 2. QC de canales
    qc_path = joinpath(out_dir, "tables", "qc_channels_$(condition).csv")
    if isfile(qc_path)
        df = CSV.read(qc_path, DataFrame)
        push!(sections, _section("Quality Control — Canales",
              _html_table(df) * _qc_bad_channels_note(df)))
    end

    # 3. Figuras de espectro (si existen)
    fig_dir = joinpath(out_dir, "figures")
    spec_figs = filter(f -> startswith(f, "psd_") && endswith(f, ".png"), readdir(fig_dir, join=false))
    if !isempty(spec_figs)
        imgs = join([_img64(joinpath(fig_dir, f), basename(f)) for f in spec_figs], "\n")
        push!(sections, _section("Análisis Espectral — PSD", imgs))
    end

    # 4. Matrices wPLI por banda
    wpli_imgs = String[]
    for band in band_names
        fpath = joinpath(fig_dir, "wpli_$(band)_$(condition).png")
        isfile(fpath) && push!(wpli_imgs, _img64(fpath, "wPLI $(band)"))
    end
    if !isempty(wpli_imgs)
        push!(sections, _section("Conectividad wPLI", _gallery(wpli_imgs)))
    end

    # 5. Tabla de edges significativos
    for band in band_names
        tpath = joinpath(out_dir, "tables", "wpli_edges_$(band)_$(condition).csv")
        if isfile(tpath)
            df = CSV.read(tpath, DataFrame)
            sig = df[df.wpli .> 0.0, :]   # mostrar solo los significativos si ya filtrados
            isempty(sig) || push!(sections, _section(
                "Edges significativos — $(band)", _html_table(first(sig, 20))))
        end
    end

    # 6. Graph metrics
    gm_path = joinpath(out_dir, "tables", "graph_metrics_$(condition).csv")
    if isfile(gm_path)
        df = CSV.read(gm_path, DataFrame)
        push!(sections, _section("Métricas de Red (Theory of Graphs)", _html_table(df)))
    end

    _wrap_html(
        "NeuroMIND — Sujeto $(subj.id) / Sesión $(sess.id) / $(condition)",
        join(sections, "\n"),
        cfg
    )
end

function _build_group_html(subjects, session_id, condition, band, cfg)
    ms   = filter(s -> s.group == "MS",      subjects)
    ctrl = filter(s -> s.group == "Control",  subjects)

    sections = String[]
    push!(sections, _section("Diseño del estudio", """
    <table class="info-table">
    <tr><td><b>Sesión</b></td><td>$(session_id)</td></tr>
    <tr><td><b>Condición</b></td><td>$(condition)</td></tr>
    <tr><td><b>Banda</b></td><td>$(band)</td></tr>
    <tr><td><b>N (MS)</b></td><td>$(length(ms))</td></tr>
    <tr><td><b>N (Control)</b></td><td>$(length(ctrl))</td></tr>
    <tr><td><b>Test estadístico</b></td><td>Mann-Whitney U + FDR Benjamini-Hochberg</td></tr>
    </table>
    """))

    # Figura de comparación grupal
    fig_dir = joinpath(results_dir(cfg), "transversal", "figures")
    comp_fig = joinpath(fig_dir, "group_comparison_$(session_id)_$(condition)_$(band).png")
    if isfile(comp_fig)
        push!(sections, _section("Matrices wPLI: MS vs Control", _img64(comp_fig, "Comparación grupal")))
    end

    # Tabla estadística
    stat_path = joinpath(results_dir(cfg), "transversal", "tables",
                         "stats_$(session_id)_$(condition)_$(band).csv")
    if isfile(stat_path)
        df  = CSV.read(stat_path, DataFrame)
        sig = df[df.p_fdr .<= get(cfg.statistics, "fdr_q", 0.05), :]
        push!(sections, _section(
            "Edges significativos (FDR q ≤ $(get(cfg.statistics,"fdr_q",0.05)))",
            isempty(sig) ? "<p>Ningún edge significativo tras corrección FDR.</p>" :
                           _html_table(sig)
        ))
    end

    # Graph metrics comparación
    gm_fig = joinpath(fig_dir, "graph_metrics_$(session_id)_$(condition)_$(band).png")
    if isfile(gm_fig)
        push!(sections, _section("Graph Metrics: MS vs Control", _img64(gm_fig, "Graph metrics")))
    end

    _wrap_html(
        "NeuroMIND — Comparación Grupal $(session_id) / $(condition) / $(band)",
        join(sections, "\n"), cfg
    )
end

function _build_longitudinal_html(subj, condition, band, cfg)
    visits = sort(collect(keys(subj.sessions)))
    sections = String[]

    push!(sections, _section("Paciente y seguimiento", """
    <table class="info-table">
    <tr><td><b>ID</b></td><td>$(subj.id)</td></tr>
    <tr><td><b>Grupo</b></td><td>$(subj.group)</td></tr>
    <tr><td><b>Visitas</b></td><td>$(join(visits, " → "))</td></tr>
    <tr><td><b>Condición</b></td><td>$(condition)</td></tr>
    <tr><td><b>Banda</b></td><td>$(band)</td></tr>
    <tr><td><b>EDSS T1</b></td><td>$(subj.clinical.EDSS)</td></tr>
    </table>
    """))

    # Figuras longitudinales
    long_dir = joinpath(results_dir(cfg), subj.id, "longitudinal")
    evo_fig  = joinpath(long_dir, "evolution_$(condition)_$(band).png")
    if isfile(evo_fig)
        push!(sections, _section("Evolución de la conectividad wPLI", _img64(evo_fig, "Evolución")))
    end

    _wrap_html(
        "NeuroMIND — Longitudinal $(subj.id) / $(condition) / $(band)",
        join(sections, "\n"), cfg
    )
end

# ─── Primitivas HTML ─────────────────────────────────────────

function _wrap_html(title::String, body::String, cfg::PipelineConfig)::String
    """<!DOCTYPE html>
<html lang="es">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>$(title)</title>
<style>
  :root { --accent:#1a6ea0; --bg:#f4f7fb; --card:#fff; --text:#1a1a2e; }
  * { box-sizing:border-box; margin:0; padding:0; }
  body { font-family:'Segoe UI',Arial,sans-serif; background:var(--bg);
         color:var(--text); font-size:14px; }
  header { background:var(--accent); color:#fff; padding:16px 32px;
           display:flex; align-items:center; gap:16px; }
  header h1 { font-size:18px; font-weight:600; }
  header small { font-size:11px; opacity:0.8; }
  nav { background:#0d4f78; padding:8px 32px; display:flex; gap:20px; flex-wrap:wrap; }
  nav a { color:#b8d8f0; text-decoration:none; font-size:13px; }
  nav a:hover { color:#fff; }
  main { max-width:1400px; margin:0 auto; padding:24px 32px; }
  .section { background:var(--card); border-radius:8px; margin-bottom:20px;
             box-shadow:0 1px 4px rgba(0,0,0,.08); overflow:hidden; }
  .section-header { background:var(--accent); color:#fff; padding:10px 18px;
                    font-size:13px; font-weight:600; letter-spacing:.3px; }
  .section-body { padding:16px 18px; }
  table { border-collapse:collapse; width:100%; font-size:12px; }
  th { background:#e8f1f8; padding:7px 10px; text-align:left;
       border-bottom:2px solid #b0cde3; font-weight:600; white-space:nowrap; }
  td { padding:6px 10px; border-bottom:1px solid #e8ecf0; }
  tr:hover td { background:#f0f6fb; }
  .gallery { display:flex; flex-wrap:wrap; gap:12px; padding:8px 0; }
  .gallery figure { margin:0; text-align:center; }
  .gallery img { max-width:380px; border-radius:6px; border:1px solid #dde4ec; }
  .gallery figcaption { font-size:11px; color:#666; margin-top:4px; }
  .info-table td:first-child { font-weight:600; color:#555; width:200px; }
  .badge-bad { background:#fee; color:#c00; border-radius:4px; padding:1px 6px;
               font-size:11px; font-weight:600; }
  .badge-ok  { background:#e8f5e9; color:#2e7d32; border-radius:4px; padding:1px 6px;
               font-size:11px; }
  footer { text-align:center; padding:16px; font-size:11px; color:#888; }
  img.fullwidth { max-width:100%; border-radius:6px; }
</style>
</head>
<body>
<header>
  <div>
    <h1>🧠 NeuroMIND</h1>
    <small>Conectividad funcional EEG · Esclerosis Múltiple · wPLI</small>
  </div>
  <div style="margin-left:auto; text-align:right">
    <div style="font-size:13px; font-weight:600">$(title)</div>
    <small>Generado: $(Dates.format(now(), "dd/mm/yyyy HH:MM"))</small>
  </div>
</header>
<nav>
  <a href="#">Inicio</a>
  <a href="#qc">QC</a>
  <a href="#spectral">Espectral</a>
  <a href="#wpli">wPLI</a>
  <a href="#stats">Estadística</a>
  <a href="#graph">Grafos</a>
</nav>
<main>
$(body)
</main>
<footer>NeuroMIND v0.2 · Rafael Castro Triguero · $(Dates.format(now(), "yyyy"))</footer>
</body>
</html>"""
end

function _section(title::String, body::String)::String
    """<div class="section">
  <div class="section-header">$(title)</div>
  <div class="section-body">$(body)</div>
</div>"""
end

function _html_table(df::DataFrame; max_rows::Int=100)::String
    rows_html = String[]
    for row in eachrow(first(df, max_rows))
        cells = join(["<td>$(v)</td>" for v in values(row)], "")
        push!(rows_html, "<tr>$(cells)</tr>")
    end
    hdrs = join(["<th>$(n)</th>" for n in names(df)], "")
    """<div style="overflow-x:auto"><table>
<thead><tr>$(hdrs)</tr></thead>
<tbody>$(join(rows_html,"\n"))</tbody>
</table></div>"""
end

function _img64(path::String, caption::String="")::String
    isfile(path) || return "<p style='color:#c00'>Figura no encontrada: $(basename(path))</p>"
    b64 = base64encode(read(path))
    ext = lowercase(splitext(path)[2][2:end])
    mime = ext == "png" ? "image/png" : "image/svg+xml"
    """<figure>
  <img src="data:$(mime);base64,$(b64)" alt="$(caption)" class="fullwidth">
  <figcaption>$(caption)</figcaption>
</figure>"""
end

function _gallery(imgs::Vector{String})::String
    "<div class=\"gallery\">$(join(imgs, "\n"))</div>"
end

function _subject_info_table(subj, sess)::String
    """<table class="info-table">
<tr><td>ID sujeto</td><td>$(subj.id)</td></tr>
<tr><td>Grupo</td><td>$(subj.group)</td></tr>
<tr><td>Edad</td><td>$(subj.age)</td></tr>
<tr><td>Sexo</td><td>$(subj.sex)</td></tr>
<tr><td>Sesión</td><td>$(sess.id) (Visita $(sess.visit_number))</td></tr>
<tr><td>EDSS</td><td>$(subj.clinical.EDSS)</td></tr>
<tr><td>Duración enfermedad</td><td>$(subj.clinical.disease_duration_y) años</td></tr>
<tr><td>Medicación</td><td>$(subj.clinical.medication)</td></tr>
<tr><td>Fatiga (MFIS)</td><td>$(subj.clinical.fatigue_score)</td></tr>
<tr><td>Cognición (SDMT)</td><td>$(subj.clinical.cognition_score)</td></tr>
</table>"""
end

function _qc_bad_channels_note(df::DataFrame)::String
    hasproperty(df, :is_bad) || return ""
    bad = df[df.is_bad .== true, :channel]
    isempty(bad) && return "<p style='color:#2e7d32; margin-top:8px'>✓ Sin canales problemáticos detectados.</p>"
    "<p style='margin-top:8px'><span class='badge-bad'>Canales a revisar: $(join(bad, ", "))</span></p>"
end
