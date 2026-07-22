# NeuroMIND/src/ica/ICAInspection.jl
# Carga de labels ICA y aplicación de rechazo de componentes.

"""
    load_ica_labels(cfg, subject_id, session_id, condition) -> Vector{Int}

Carga los índices de componentes ICA rechazados desde un CSV de inspección manual.
Formato esperado: CSV con columnas `component` (int) y `label` (artifact/brain).

Busca en múltiples rutas candidatas (pipeline nuevo y legacy):
  1. results/subjects/sub-{id}/ses-{sess}/{cond}/ica_labels.csv
  2. results/subjects/sub-{id}/ses-{sess}/{cond}/ICA_labels.csv
  3. results/subjects/sub-{id}/ses-{sess}/{cond}/tables/ICA_labels_{cond}.csv
  4. results/{id}/{sess}/tables/ICA_labels_{cond}.csv   (legacy)
"""
function load_ica_labels(
    cfg::PipelineConfig,
    subject_id::String,
    session_id::String,
    condition::String
)::Vector{Int}
    res = results_dir(cfg)
    # Árbol ÚNICO: BIDS. task = eyesclosed/eyesopen (no la condición EC/EO).
    task = condition == "EC" ? "eyesclosed" : (condition == "EO" ? "eyesopen" : condition)
    subj_base = joinpath(res, "subjects", "sub-$(subject_id)", "ses-$(session_id)", task)

    candidates = [
        joinpath(subj_base, "ica_labels.csv"),
        joinpath(subj_base, "ICA_labels.csv"),
    ]

    idx = findfirst(isfile, candidates)
    idx === nothing && return Int[]

    df = CSV.read(candidates[idx], DataFrame)
    hasproperty(df, :label)     || return Int[]
    hasproperty(df, :component) || return Int[]

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
