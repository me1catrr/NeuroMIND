# NeuroMIND/src/ica/ICAInspection.jl
# Carga de labels ICA y aplicación de rechazo de componentes.

"""
    load_ica_labels(cfg, subject_id, session_id, condition) -> Vector{Int}

Carga los índices de componentes ICA rechazados desde un CSV de inspección manual.
Formato esperado: CSV con columna `component` (int) y `label` (artifact/brain).
"""
function load_ica_labels(
    cfg::PipelineConfig,
    subject_id::String,
    session_id::String,
    condition::String
)::Vector{Int}
    path = joinpath(
        results_dir(cfg), subject_id, session_id,
        "tables", "ICA_labels_$(condition).csv"
    )
    isfile(path) || return Int[]

    df = CSV.read(path, DataFrame)
    hasproperty(df, :label) || return Int[]

    rejected = df[df.label .== "artifact", :component]
    return convert(Vector{Int}, rejected)
end

"""
    apply_ica_rejection(rec::EEGRecording, ica::ICAResult, rejected::Vector{Int}) -> EEGRecording

Reconstruye la señal EEG eliminando los componentes marcados como artefacto.
"""
function apply_ica_rejection(
    rec::EEGRecording,
    ica::ICAResult,
    rejected::Vector{Int}
)::EEGRecording

    isempty(rejected) && return rec

    n_comp = size(ica.activations, 1)
    keep   = setdiff(1:n_comp, rejected)

    # Reconstrucción: data_clean = A[:, keep] * activations[keep, :]
    A_keep = ica.mixing_matrix[:, keep]
    S_keep = ica.activations[keep, :]

    clean_data = A_keep * S_keep

    # Alinear dimensiones si hay diferencia por la proyección PCA
    n_ch, n_samp = size(rec.data)
    if size(clean_data, 1) != n_ch || size(clean_data, 2) != n_samp
        # Fallback: sustracción de componentes de artefacto
        artifact = ica.mixing_matrix[:, rejected] * ica.activations[rejected, :]
        rows = min(n_ch, size(artifact, 1))
        cols = min(n_samp, size(artifact, 2))
        clean_data = copy(rec.data)
        clean_data[1:rows, 1:cols] .-= artifact[1:rows, 1:cols]
    end

    updated_ica = ICAResult(
        ica.meta, ica.mixing_matrix, ica.unmixing_matrix,
        ica.activations, rejected, ica.variance_explained
    )

    return EEGRecording(rec.meta, clean_data, rec.times)
end
