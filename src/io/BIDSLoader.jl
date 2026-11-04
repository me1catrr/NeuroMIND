# NeuroMIND/src/io/BIDSLoader.jl
# Utilidades del dataset BIDS-EEG (data/bids/). La carga de señal vive en
# BrainVisionLoader.jl (load_eeg_brainvision) y la orquesta
# load_single_subject (SingleSubjectPipeline.jl); las rutas del árbol
# sub-*/ses-*/eeg/ están en Config.jl (bids_vhdr_path, bids_electrodes_path).

# ─── Helpers privados ─────────────────────────────────────────

function _load_electrode_positions(tsv_path::String)
    df  = CSV.read(tsv_path, DataFrame; delim='\t')
    pos = Dict{String,Tuple{Float64,Float64}}()
    for r in eachrow(df)
        type_col = hasproperty(r, :type) ? string(r.type) : "EEG"
        if type_col == "EEG"
            name = uppercase(strip(string(r.name)))
            pos[name] = (Float64(r.x), Float64(r.y))
        end
    end
    return pos
end
