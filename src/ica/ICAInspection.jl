# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Inspección y rechazo ICA
# ═══════════════════════════════════════════════════════════════
#
#  Carga etiquetas de inspección MANUAL (CSV) y reconstruye la
#  señal EEG eliminando los componentes marcados como artefacto.
#  Si no hay CSV manual, SingleSubjectPipeline.jl decide si generar
#  uno automático (write_ica_labels_auto, cfg.ica["auto_reject"]).
#
#  CSV manual esperado
#  ────────────────────
#    columnas: component (Int), label ("artifact" | "brain"|…)
#    rutas candidatas (árbol BIDS):
#      results/subjects/sub-{id}/ses-{sess}/{task}/ica_labels.csv
#      results/subjects/sub-{id}/ses-{sess}/{task}/ICA_labels.csv
#    Tiene SIEMPRE precedencia sobre el CSV automático.
#
#  API pública (exportada por NeuroMIND)
#  ─────────────────────────────────────
#    load_ica_labels(cfg, subject_id, session_id, condition) → Vector{Int}
#    apply_ica_rejection(rec, ica, rejected) → EEGRecording
#    write_ica_labels_auto(export_dir, eval_df) → Vector{Int}
#
# ───────────────────────────────────────────────────────────────
#  Fichero    src/ica/ICAInspection.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 22-07-2026
# ───────────────────────────────────────────────────────────────

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
function _ica_manual_label_candidates(
    cfg::PipelineConfig,
    subject_id::String,
    session_id::String,
    condition::String
)::Vector{String}
    res = results_dir(cfg)
    # Árbol ÚNICO: BIDS. task = eyesclosed/eyesopen (no la condición EC/EO).
    task = condition == "EC" ? "eyesclosed" : (condition == "EO" ? "eyesopen" : condition)
    subj_base = joinpath(res, "subjects", "sub-$(subject_id)", "ses-$(session_id)", task)
    return [
        joinpath(subj_base, "ica_labels.csv"),
        joinpath(subj_base, "ICA_labels.csv"),
    ]
end

"""
    has_manual_ica_labels(cfg, subject_id, session_id, condition) -> Bool

`true` si existe un CSV de inspección MANUAL (`ica_labels.csv` /
`ICA_labels.csv`). Se usa para decidir precedencia manual > auto.
"""
function has_manual_ica_labels(
    cfg::PipelineConfig,
    subject_id::String,
    session_id::String,
    condition::String
)::Bool
    candidates = _ica_manual_label_candidates(cfg, subject_id, session_id, condition)
    return any(isfile, candidates)
end

function load_ica_labels(
    cfg::PipelineConfig,
    subject_id::String,
    session_id::String,
    condition::String
)::Vector{Int}
    candidates = _ica_manual_label_candidates(cfg, subject_id, session_id, condition)

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
        ica.activations, rejected, ica.variance_explained,
        ica.diagnostics,
    )

    return EEGRecording(rec.meta, clean_data, rec.times)
end

"""
    write_ica_labels_auto(export_dir, eval_df) -> Vector{Int}

Genera y persiste `ica_labels_auto.csv` a partir de la clasificación
automática (`evaluate_ica_components`), y devuelve los índices con
`artifact_type != "brain"`.

No toca `ica_labels.csv` (reservado para inspección MANUAL, que
siempre tiene precedencia — ver `load_ica_labels`). Se sobrescribe
en cada ejecución para reflejar el `artifact_threshold` vigente.
"""
function write_ica_labels_auto(export_dir::String, eval_df::DataFrame)::Vector{Int}
    hasproperty(eval_df, :artifact_type) || return Int[]

    types    = String.(eval_df.artifact_type)
    rejected = findall(!=("brain"), types)

    out = DataFrame(
        component     = eval_df.component,
        label         = [t == "brain" ? "brain" : "artifact" for t in types],
        artifact_type = types,
        artifact_score = eval_df.artifact_score,
    )
    mkpath(export_dir)
    CSV.write(joinpath(export_dir, "ica_labels_auto.csv"), out)

    return rejected
end
