module PlutoTemplateBase

using Dates
using Markdown
using CSV
using DataFrames
using Serialization
using Statistics
using DSP
using Plots

# Rutas del proyecto centralizadas en src/modules/paths.jl
include(joinpath(@__DIR__, "..", "src", "modules", "paths.jl"))

const FS_DEFAULT = 500.0
const NOTEBOOK_AUTHOR = "RAFAEL CASTRO TRIGUERO\nALEJANDRO GALVAO\nJAVIER ESPUNY"

function escape_html(s::AbstractString)
    # Escapa caracteres básicos para evitar que el HTML se rompa.
    return replace(replace(replace(s, "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")
end

function published_date_str(d::Date = today())
    month_names_es = [
        "enero", "febrero", "marzo", "abril", "mayo", "junio",
        "julio", "agosto", "septiembre", "octubre", "noviembre", "diciembre",
    ]
    return "$(Dates.day(d)) de $(month_names_es[Dates.month(d)]) de $(Dates.year(d))"
end

function notebook_intro(title::AbstractString; author::AbstractString = NOTEBOOK_AUTHOR, published::Date = today())
    published_str = published_date_str(published)
    author_lines = split(author, '\n')
    author_html = join(escape_html.(author_lines), "<br/>")

    html = """
<h1>$(escape_html(title))</h1>
<table style="border-collapse: collapse; margin: 0 auto; width: 100%; max-width: 760px;">
  <thead>
    <tr>
      <th style="border-bottom: 2px solid #ddd; padding: 6px 12px; text-align: center;">AUTORES</th>
      <th style="border-bottom: 2px solid #ddd; padding: 6px 12px; text-align: center;">PUBLICADO</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 6px 12px; text-align: center; font-weight: 600; line-height: 1.25;">
        $author_html
      </td>
      <td style="padding: 6px 12px; text-align: center;">
        $(escape_html(published_str))
      </td>
    </tr>
  </tbody>
</table>
"""
    return Markdown.HTML(html)
end

export FS_DEFAULT
export NOTEBOOK_AUTHOR, published_date_str, notebook_intro
export project_root, data_root, results_root, stage_dir, bids_root, raw_dir, electrodes_dir

end # module PlutoTemplateBase
