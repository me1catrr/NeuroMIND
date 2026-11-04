module DefaultConfig

# --- Rutas base del proyecto ---
# @__DIR__ = carpeta "config"
const ROOT_DIR     = normpath(joinpath(@__DIR__, ".."))
const DATA_DIR     = joinpath(ROOT_DIR, "data")
const RESULTS_DIR  = joinpath(ROOT_DIR, "results")
const FIGURES_DIR  = joinpath(RESULTS_DIR, "figures")
const LOGS_DIR     = joinpath(RESULTS_DIR, "logs")
const TABLES_DIR   = joinpath(RESULTS_DIR, "tables")

# --- Definición de Structs de configuración ---

struct FilterConfig
    fs::Float64
    notch_freq::Float64
    bandreject_center::Union{Nothing,Float64}
    bandreject_width::Union{Nothing,Float64}
    hp_cutoff::Float64
    lp_cutoff::Float64
end

struct ICAConfig
    method::Symbol
    max_steps::Int
    tol::Float64
    components_to_zero::Vector{Int}
end

struct SegmentationConfig
    length_s::Float64
    overlap_s::Float64
end

struct ArtifactConfig
    amp_min::Float64
    amp_max::Float64
end

struct FFTConfig
    pad_to::Int
    window::Symbol
    full_spectrum::Bool
end

struct PipelineConfig
    data_dir::String
    output_dir::String
    filter::FilterConfig
    ica::ICAConfig
    segmentation::SegmentationConfig
    artifact::ArtifactConfig
    fft::FFTConfig
end

# --- Configuración por defecto del pipeline ---

const DEFAULT_CONFIG = PipelineConfig(
    DATA_DIR,          # data_dir
    RESULTS_DIR,       # output_dir

    FilterConfig(
        500.0,         # fs
        50.0,          # notch
        100.0,         # bandreject_center
        1.0,           # bandreject_width
        0.5,           # hp
        150.0,         # lp
    ),

    ICAConfig(
        :probabilistic,
        512,
        1e-7,
        [0, 1, 3, 6, 12, 19, 24, 27],
    ),

    SegmentationConfig(1.0, 0.0),
    ArtifactConfig(-70.0, 70.0),
    FFTConfig(512, :hamming, true),
)

# Exportar lo que interesa usar fuera
export ROOT_DIR, DATA_DIR, RESULTS_DIR, FIGURES_DIR, LOGS_DIR, TABLES_DIR
export FilterConfig, ICAConfig, SegmentationConfig, ArtifactConfig, FFTConfig
export PipelineConfig, DEFAULT_CONFIG

end # module
