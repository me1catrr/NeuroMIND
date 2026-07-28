# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — ViewerSupport
#  Helpers compartidos por los dos visores de cohorte interactivos
#  (plot_transversal.jl :8781 · plot_longitudinal.jl :8780): IO de
#  CSV/JSON, constantes de banda/electrodos y el único export PNG
#  bajo demanda que hacen (Δ heatmap). No lo usa nada de la capa de
#  generación de figuras (Transversal.jl/Longitudinal.jl tienen su
#  propia copia — ver nota de diseño del plan §6).
# ═══════════════════════════════════════════════════════════════

module ViewerSupport

using CSV, DataFrames, CairoMakie, TOML

export BAND_ORDER, CH_POS, DIVERGING_CMAP
export parse_summary_json, read_mat_csv, safe_csv, resolve_results_root
export json_escape, mat_flat, df_rows_json, pos_json
export count_fdr, count_uncorr, save_heatmap

const BAND_ORDER = ["DELTA", "THETA", "ALPHA", "BETA_LOW", "BETA_MID", "BETA_HIGH", "GAMMA"]

# Positivo (MS−Ctrl / T2−T1) = rojo; negativo = azul.
const DIVERGING_CMAP = Reverse(:RdBu)

const CH_POS = Dict{String,Tuple{Float64,Float64}}(
    "FZ"=>(0.00,0.72),"F3"=>(-0.35,0.55),"F4"=>(0.35,0.55),
    "F7"=>(-0.68,0.42),"F8"=>(0.68,0.42),"FT9"=>(-0.85,0.18),"FT10"=>(0.85,0.18),
    "FC5"=>(-0.50,0.28),"FC1"=>(-0.18,0.28),"FC2"=>(0.18,0.28),"FC6"=>(0.50,0.28),
    "C3"=>(-0.40,0.00),"CZ"=>(0.00,0.00),"C4"=>(0.40,0.00),
    "T7"=>(-0.80,0.00),"T8"=>(0.80,0.00),"TP9"=>(-0.85,-0.22),"TP10"=>(0.85,-0.22),
    "CP5"=>(-0.50,-0.28),"CP1"=>(-0.18,-0.28),"CP2"=>(0.18,-0.28),"CP6"=>(0.50,-0.28),
    "P3"=>(-0.35,-0.55),"PZ"=>(0.00,-0.55),"P4"=>(0.35,-0.55),
    "P7"=>(-0.68,-0.48),"P8"=>(0.68,-0.48),
    "O1"=>(-0.28,-0.82),"OZ"=>(0.00,-0.88),"O2"=>(0.28,-0.82),
    "FP1"=>(-0.25,0.88),"FP2"=>(0.25,0.88),
)

# ── IO ─────────────────────────────────────────────────────────

function parse_summary_json(path::String)::Dict{String,Any}
    out = Dict{String,Any}()
    isfile(path) || return out
    txt = read(path, String)
    for m in eachmatch(r"\"([^\"]+)\"\s*:\s*\"([^\"]*)\"", txt)
        out[String(m.captures[1])] = String(m.captures[2])
    end
    for m in eachmatch(r"\"([^\"]+)\"\s*:\s*(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)", txt)
        k = String(m.captures[1])
        haskey(out, k) && continue
        v = tryparse(Float64, m.captures[2])
        v === nothing || (out[k] = v)
    end
    for m in eachmatch(r"\"([^\"]+)\"\s*:\s*(true|false)", txt)
        k = String(m.captures[1])
        haskey(out, k) && continue
        out[k] = m.captures[2] == "true"
    end
    return out
end

function read_mat_csv(path::String)
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame)
    isempty(df) && return nothing
    ch = String.(df[!, 1])
    n = length(ch)
    M = zeros(Float64, n, n)
    for (j, name) in enumerate(ch)
        col = Symbol(name)
        hasproperty(df, col) || return nothing
        M[:, j] = Float64.(df[!, col])
    end
    return (ch, M)
end

function safe_csv(path::String)::DataFrame
    isfile(path) || return DataFrame()
    try
        return CSV.read(path, DataFrame)
    catch
        return DataFrame()
    end
end

function resolve_results_root(proj::String, explicit::Union{Nothing,String}=nothing)::String
    explicit !== nothing && return abspath(explicit)
    cfg_p = joinpath(proj, "config", "pipeline.toml")
    if isfile(cfg_p)
        cfg = TOML.parsefile(cfg_p)
        r = get(get(cfg, "paths", Dict()), "results", "results")
        return isabspath(r) ? r : joinpath(proj, r)
    end
    return joinpath(proj, "results")
end

# ── JSON helpers ───────────────────────────────────────────────

json_escape(s::AbstractString) =
    replace(replace(replace(String(s), "\\" => "\\\\"), "\"" => "\\\""), "\n" => "\\n")

function mat_flat(M::Matrix{Float64})::String
    n = size(M, 1)
    parts = String[]
    sizehint!(parts, n * n)
    for j in 1:n, i in 1:n
        x = M[i, j]
        push!(parts, isfinite(x) ? repr(x) : "null")
    end
    return "[" * join(parts, ",") * "]"
end

function df_rows_json(df::DataFrame, cols::Vector{Symbol}; limit::Int=0)::String
    nrow(df) == 0 && return "[]"
    n = limit > 0 ? min(limit, nrow(df)) : nrow(df)
    rows = String[]
    for i in 1:n
        fields = String[]
        for c in cols
            hasproperty(df, c) || continue
            v = df[i, c]
            if v isa AbstractString
                push!(fields, "\"$c\":\"$(json_escape(string(v)))\"")
            elseif v isa Bool
                push!(fields, "\"$c\":$(v ? "true" : "false")")
            elseif v isa Integer
                push!(fields, "\"$c\":$v")
            elseif v isa Real
                x = Float64(v)
                push!(fields, "\"$c\":$(isfinite(x) ? repr(x) : "null")")
            else
                push!(fields, "\"$c\":\"$(json_escape(string(v)))\"")
            end
        end
        push!(rows, "{" * join(fields, ",") * "}")
    end
    return "[" * join(rows, ",") * "]"
end

function pos_json()::String
    parts = ["\"$k\":[$(v[1]),$(v[2])]" for (k, v) in CH_POS]
    return "{" * join(parts, ",") * "}"
end

# ── Band / KPI helpers ─────────────────────────────────────────

function count_fdr(edf::DataFrame; qcol=:q_value)::Int
    nrow(edf) == 0 && return 0
    hasproperty(edf, qcol) || return 0
    return count(r -> Float64(r[qcol]) < 0.05, eachrow(edf))
end

function count_uncorr(edf::DataFrame; pcol=:p_value)::Int
    nrow(edf) == 0 && return 0
    hasproperty(edf, pcol) || return 0
    return count(r -> Float64(r[pcol]) < 0.05, eachrow(edf))
end

# ── Export PNG bajo demanda (botón "exportar" del visor) ────────

function save_heatmap(path, W, ch, title; cmap=:viridis, diverging=false,
                      colorrange=nothing, colorbar_label="")
    n = length(ch)
    fig = Figure(size=(720, 640), fontsize=11)
    ax = Axis(fig[1, 1]; title=title, xlabel="Canal", ylabel="Canal",
              xticks=(1:n, ch), yticks=(1:n, ch),
              xticklabelrotation=π/3, xticklabelsize=7, yticklabelsize=7)
    if diverging
        lim = colorrange === nothing ? maximum(abs, W) : maximum(abs, colorrange)
        lim < 1e-12 && (lim = 1.0)
        hm = heatmap!(ax, W; colormap=DIVERGING_CMAP, colorrange=(-lim, lim))
        Colorbar(fig[1, 2], hm; label=isempty(colorbar_label) ? "Δ" : colorbar_label)
    else
        cr = colorrange
        hm = cr === nothing ? heatmap!(ax, W; colormap=cmap) :
             heatmap!(ax, W; colormap=cmap, colorrange=cr)
        Colorbar(fig[1, 2], hm; label=isempty(colorbar_label) ? "wPLI" : colorbar_label)
    end
    mkpath(dirname(path))
    save(path, fig)
    return path
end

end # module
