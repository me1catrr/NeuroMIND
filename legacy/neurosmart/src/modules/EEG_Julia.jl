module EEG_Julia

"""
Módulo principal del proyecto EEG_Julia.

- Centraliza utilidades reutilizables (rutas, IO binario, logging).
- Expone funciones de alto nivel `run_*` para ejecutar cada etapa de forma individual
  (útil desde REPL y notebooks Pluto.jl).
"""

# Utilidades comunes (deben ir primero)
include(joinpath(@__DIR__, "paths.jl"))
include(joinpath(@__DIR__, "utils.jl"))

# Etapas del pipeline
include(joinpath(@__DIR__, "..", "Setup", "IO.jl"))
include(joinpath(@__DIR__, "..", "Preprocessing", "filtering.jl"))
include(joinpath(@__DIR__, "..", "ICA", "ICA.jl"))
include(joinpath(@__DIR__, "..", "ICA", "ICA_cleaning.jl"))
include(joinpath(@__DIR__, "..", "Processing", "segmentation.jl"))
include(joinpath(@__DIR__, "..", "Processing", "baseline.jl"))
include(joinpath(@__DIR__, "..", "Processing", "artifact_rejection.jl"))
include(joinpath(@__DIR__, "..", "Processing", "baseline_2st.jl"))
include(joinpath(@__DIR__, "..", "Spectral", "FFT.jl"))

include(joinpath(@__DIR__, "..", "Connectivity", "CSD.jl"))
include(joinpath(@__DIR__, "..", "Connectivity", "wPLI.jl"))

# Exports: solo lo "alto nivel" + rutas/IO comunes útiles desde Pluto
export project_root, data_root, results_root, stage_dir, bids_root, raw_dir, electrodes_dir
export ensure_dir, clean_dir_by_ext, run_id, logmsg, open_log, save_bin, load_bin, pick_first_existing

export run_io, run_filtering, run_ica, run_ica_cleaning, run_segmentation
export run_baseline, run_artifact_rejection, run_baseline_2st, run_fft
export run_csd, run_wpli

end # module EEG_Julia

