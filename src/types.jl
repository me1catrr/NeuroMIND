# NeuroMIND/src/types.jl
# Entidades científicas centrales.
# Diseño orientado a objetos: cada struct representa una entidad del estudio.

# ─────────────────────────────────────────────────────────────
# Configuración
# ─────────────────────────────────────────────────────────────

struct PipelineConfig
    project::Dict{String,Any}
    study::Dict{String,Any}
    paths::Dict{String,Any}
    recording::Dict{String,Any}
    filtering::Dict{String,Any}
    segmentation::Dict{String,Any}
    baseline::Dict{String,Any}
    artifact_rejection::Dict{String,Any}
    ica::Dict{String,Any}
    spectral::Dict{String,Any}
    bands::Dict{String,Tuple{Float64,Float64}}
    connectivity::Dict{String,Any}
    surrogates::Dict{String,Any}
    graph::Dict{String,Any}
    clinical::Dict{String,Any}
    longitudinal::Dict{String,Any}
    statistics::Dict{String,Any}
    export_cfg::Dict{String,Any}
    qc::Dict{String,Any}
    montage::Dict{String,Any}
    root::String
end

# ─────────────────────────────────────────────────────────────
# Datos clínicos del sujeto
# ─────────────────────────────────────────────────────────────

struct ClinicalData
    EDSS::Union{Float64,Missing}
    disease_duration_y::Union{Float64,Missing}
    medication::Union{String,Missing}
    fatigue_score::Union{Float64,Missing}
    cognition_score::Union{Float64,Missing}
    lesion_load::Union{Float64,Missing}    # mm³ de lesiones en RM
end

ClinicalData() = ClinicalData(missing, missing, missing, missing, missing, missing)

# ─────────────────────────────────────────────────────────────
# Metadatos de grabación
# ─────────────────────────────────────────────────────────────

struct RecordingMeta
    subject_id::String
    session_id::String
    condition::String        # "EO" | "EC"
    run::Int
    fs::Float64
    n_channels::Int
    channel_names::Vector{String}
    channel_positions::Union{Nothing,Dict{String,Tuple{Float64,Float64}}}
    bids_path::String
end

# ─────────────────────────────────────────────────────────────
# Señal EEG en diferentes etapas
# ─────────────────────────────────────────────────────────────

struct EEGRecording
    meta::RecordingMeta
    data::Matrix{Float64}    # μV, (channels × samples)
    times::Vector{Float64}   # segundos
end

n_channels(r::EEGRecording) = size(r.data, 1)
n_samples(r::EEGRecording)  = size(r.data, 2)
duration(r::EEGRecording)   = n_samples(r) / r.meta.fs

struct EpochSet
    meta::RecordingMeta
    data::Array{Float64,3}   # (channels × samples × epochs)
    epoch_length_s::Float64
    n_valid::Int
    rejected_idx::Vector{Int}
    rejection_reasons::Vector{String}
end

n_epochs(e::EpochSet)        = size(e.data, 3)
n_samples_epoch(e::EpochSet) = size(e.data, 2)
rejection_rate(e::EpochSet)  = length(e.rejected_idx) / (e.n_valid + length(e.rejected_idx))

# Constructor con rejection_reasons opcional (evita ruptura de llamadas existentes)
EpochSet(meta, data, epoch_length_s, n_valid, rejected_idx) =
    EpochSet(meta, data, epoch_length_s, n_valid, rejected_idx, String[])

# ─────────────────────────────────────────────────────────────
# ICA
# ─────────────────────────────────────────────────────────────

struct ICAResult
    meta::RecordingMeta
    mixing_matrix::Matrix{Float64}      # (channels × components)
    unmixing_matrix::Matrix{Float64}    # (components × channels)
    activations::Matrix{Float64}        # (components × samples)
    rejected_components::Vector{Int}
    variance_explained::Vector{Float64}
    diagnostics::Dict{String,Any}       # convergencia FastICA + métricas PCA
end

# Compat: llamadas antiguas sin diagnostics
ICAResult(meta, A, W, S, rej, var_exp) =
    ICAResult(meta, A, W, S, rej, var_exp, Dict{String,Any}())

# ─────────────────────────────────────────────────────────────
# Spectral
# ─────────────────────────────────────────────────────────────

struct SpectralResult
    meta::RecordingMeta
    psd::Matrix{Float64}               # μV², (channels × freqs)
    freqs::Vector{Float64}             # Hz
    band_power::Dict{String,Vector{Float64}}   # banda → potencia por canal
    n_epochs_used::Int
    params::Dict{String,Any}
end

# ─────────────────────────────────────────────────────────────
# Conectividad
# ─────────────────────────────────────────────────────────────

struct ConnectivityMatrix
    meta::RecordingMeta
    method::String                     # "wpli"
    matrices::Dict{String,Matrix{Float64}}     # banda → W[n_ch, n_ch]
    channel_names::Vector{String}
    space::String                      # "sensor" | "CSD"
    n_epochs_used::Int
    params::Dict{String,Any}
end

n_channels(c::ConnectivityMatrix) = length(c.channel_names)

struct SurrogateResult
    connectivity::ConnectivityMatrix
    band::String
    observed::Matrix{Float64}          # wPLI observado
    null_distribution::Array{Float64,3}  # (n_ch × n_ch × n_surrogates)
    p_values::Matrix{Float64}
    sig_mask::BitMatrix
    fdr_threshold::Float64
    n_surrogates::Int
end

# ─────────────────────────────────────────────────────────────
# Graph metrics
# ─────────────────────────────────────────────────────────────

struct GraphMetrics
    band::String
    threshold::Float64
    density::Float64
    strength::Vector{Float64}          # por nodo
    clustering::Vector{Float64}        # por nodo
    path_length::Float64               # global
    efficiency::Float64                # global
    modularity::Float64                # global
    channel_names::Vector{String}
end

# ─────────────────────────────────────────────────────────────
# Resultado estadístico
# ─────────────────────────────────────────────────────────────

struct StatResult
    test_name::String           # "mann_whitney", "wilcoxon", etc.
    statistic::Float64
    p_value::Float64
    p_adjusted::Float64         # tras corrección FDR/Bonferroni
    effect_size::Float64        # Cohen's d, r de Spearman, etc.
    significant::Bool
    group_a_mean::Float64
    group_b_mean::Float64
    n_a::Int
    n_b::Int
end

# ─────────────────────────────────────────────────────────────
# Sujeto y sesión
# ─────────────────────────────────────────────────────────────

mutable struct Session
    id::String
    visit_number::Int
    recordings::Dict{String,EEGRecording}
    epochs::Dict{String,EpochSet}
    ica::Dict{String,ICAResult}
    spectra::Dict{String,SpectralResult}
    connectivity::Dict{String,ConnectivityMatrix}
    surrogates::Dict{String,SurrogateResult}
    graph_metrics::Dict{String,GraphMetrics}    # clave: "condicion_banda"
end

Session(id::String, visit::Int) = Session(
    id, visit, Dict(), Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
)

mutable struct Subject
    id::String
    group::String             # "MS" | "Control"
    age::Union{Int,Missing}
    sex::Union{String,Missing}
    clinical::ClinicalData
    sessions::Dict{String,Session}
end

Subject(id::String, group::String) = Subject(
    id, group, missing, missing, ClinicalData(), Dict()
)

# ─────────────────────────────────────────────────────────────
# Análisis de grupo y longitudinal
# ─────────────────────────────────────────────────────────────

struct GroupAnalysis
    session_id::String
    condition::String
    band::String
    group_ms::String
    group_ctrl::String
    mean_connectivity_ms::Matrix{Float64}
    mean_connectivity_ctrl::Matrix{Float64}
    stat_results::Matrix{StatResult}    # un StatResult por par de canales
    channel_names::Vector{String}
    n_ms::Int
    n_ctrl::Int
end

struct LongitudinalAnalysis
    subject_id::String
    condition::String
    band::String
    visits::Vector{String}
    connectivity_over_time::Vector{Matrix{Float64}}
    graph_metrics_over_time::Vector{GraphMetrics}
    channel_names::Vector{String}
end
