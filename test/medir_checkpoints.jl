# ───────────────────────────────────────────────────────────────
#  Fichero    test/medir_checkpoints.jl
#  Modificado 14-11-2026
# ───────────────────────────────────────────────────────────────
#
# Medición de tamaño y coste de los checkpoints (plan de 2f, apartado 3.6;
# entrega 2f-1). No forma parte de runtests.jl: se ejecuta a mano.
#
#   julia --project=. test/medir_checkpoints.jl [segundos] [canales] [N] [salida.toml]
#
# Genera aquí una grabación sintética representativa (por defecto 31 canales a
# 500 Hz y 300 s, con surrogates activados y N = 20) y mide, por etapa:
#
# - bytes en disco del `.jls` y del `.toml`;
# - el `.bin` que guarda Lab hoy para la misma etapa (`save_result` de lo que
#   guarda `runner/stage.jl`), como línea base;
# - tiempo de cálculo de la etapa;
# - tiempo de `save_checkpoint` solo (escritura segura, sin el inventario de
#   productos del orquestador) y del paso completo del orquestador;
# - tiempo de `load_checkpoint` con verificación e invariantes, y de
#   `verify_checkpoint` sin deserializar;
#
# cinco repeticiones de cada lectura y escritura (mediana y máximo). No lee
# ninguna grabación real.

using NeuroMIND, Random, Statistics, Serialization, Dates, TOML

const CK = NeuroMIND
seconds  = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 300
nchan    = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 31
nsurr    = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20
salida   = length(ARGS) >= 4 ? ARGS[4] : ""
const FS = 500
const REPS = 5

const CHANS32 = ["Fp1", "Fp2", "F7", "F3", "Fz", "F4", "F8", "FC5", "FC1", "FC2", "FC6", "T7", "C3", "Cz",
                 "C4", "T8", "TP9", "CP5", "CP1", "CP2", "CP6", "TP10", "P7", "P3", "Pz", "P4", "P8",
                 "PO9", "O1", "Oz", "O2", "PO10"]

"Carpeta BIDS sintética: `nchan` canales, `seconds` s a 500 Hz; uno anómalo (T8)."
function bids_representativa(root, seconds, nchan)
    chans = CHANS32[1:nchan]
    dir  = joinpath(root, "sub-S01", "ses-T1", "eeg"); mkpath(dir)
    stem = "sub-S01_ses-T1_task-eyesclosed_run-01_eeg"
    n = seconds * FS
    rng = MersenneTwister(20261114)
    data = Matrix{Float32}(undef, length(chans), n)
    for t in 1:n
        x = (t - 1) / FS
        common = 6sin(2π * 50x) + 2sin(2π * 100x) + 14sin(2π * 0.08x)
        for (c, name) in enumerate(chans)
            amp = name[1] in ('O', 'P') ? 16 : 7
            v = amp * sin(2π * 10x + 0.7c) * (0.7 + 0.3sin(2π * 0.21x + c)) + 3sin(2π * 5.5x + 1.3c) +
                common + 6 * (2rand(rng) - 1)
            name == "T8" && (v = 14v + 60 * (2rand(rng) - 1))
            data[c, t] = v
        end
    end
    write(joinpath(dir, stem * ".eeg"), data)
    write(joinpath(dir, stem * ".vhdr"), join([
        "Brain Vision Data Exchange Header File Version 1.0",
        "; Grabación sintética para medir checkpoints. No procede de ningún sujeto.", "",
        "[Common Infos]", "Codepage=UTF-8", "DataFile=$(stem).eeg", "MarkerFile=$(stem).vmrk",
        "DataFormat=BINARY", "DataOrientation=MULTIPLEXED", "NumberOfChannels=$(length(chans))",
        "SamplingInterval=$(1e6 / FS)", "", "[Binary Infos]", "BinaryFormat=IEEE_FLOAT_32", "",
        "[Channel Infos]", ["Ch$(i)=$(nm),,1,µV" for (i, nm) in enumerate(chans)]..., ""], "\n"))
    write(joinpath(dir, stem * ".vmrk"), join([
        "Brain Vision Data Exchange Marker File, Version 1.0", "", "[Common Infos]", "Codepage=UTF-8",
        "DataFile=$(stem).eeg", "", "[Marker Infos]", "Mk1=New Segment,,1,1,0,20261114000000000000", ""], "\n"))
    open(joinpath(dir, "sub-S01_ses-T1_electrodes.tsv"), "w") do io
        println(io, "name\tx\ty\tz")
        for (i, nm) in enumerate(chans)
            a = 2π * (i - 1) / length(chans)
            println(io, nm, "\t", round(0.08cos(a), digits=4), "\t", round(0.08sin(a), digits=4), "\t0.05")
        end
    end
    return root
end

"Protocolo del repositorio con el sujeto sintético, surrogates activados y resultados en `root`."
function config_representativa(root, bids, nsurr)
    lines = split(read(joinpath(@__DIR__, "..", "config", "pipeline.toml"), String), '\n')
    function set!(section, key, value)
        inside = false
        for (i, l) in enumerate(lines)
            startswith(strip(l), "[") && (inside = startswith(strip(l), "[$(section)]"))
            inside && occursin(Regex("^\\s*$(key)\\s*="), l) && (lines[i] = "$(key) = $(value)"; return)
        end
        error("El protocolo no tiene [$(section)] $(key)")
    end
    set!("subject", "subject_id", "\"S01\""); set!("subject", "session_id", "\"T1\"")
    set!("subject", "task", "\"eyesclosed\""); set!("paths", "bids_root", repr(bids))
    set!("surrogates", "enabled", "true"); set!("surrogates", "n_surrogates", string(nsurr))
    path = joinpath(root, "config", "pipeline.toml"); mkpath(dirname(path))
    write(path, join(lines, '\n'))
    return path
end

"Bytes del `.bin` que guarda Lab hoy (envoltorio de `save_result`) para lo que guarda en esa etapa."
function bytes_lab(obj)
    obj === nothing && return 0
    io = IOBuffer()
    serialize(io, Dict("version" => "0.1.0", "stage" => "x", "subject" => "S01", "session" => "T1",
                       "saved_at" => string(now()), "data" => obj))
    return position(io)
end

tmp  = mktempdir()
bids = bids_representativa(joinpath(tmp, "bids"), seconds, nchan)
cfg  = config_representativa(joinpath(tmp, "run"), bids, nsurr)
println("Sintética: $(nchan) canales · $(seconds) s · $(FS) Hz · surrogates N = $(nsurr)")

o = Dict{Symbol,Any}(); t_stage = Dict{Symbol,Float64}(); t_ck = Dict{Symbol,Float64}()
ctx = redirect_stdout(devnull) do; pipeline_context(cfg); end
ck  = CK._checkpoint_run(ctx, nothing)
function step!(s, f; kw...)
    t_stage[s] = @elapsed (o[s] = redirect_stdout(f, devnull))
    t_ck[s]    = @elapsed CK._checkpoint_after!(ck, s, o[s]; kw...)
end
step!(:load,   () -> stage_load(ctx))
step!(:qc,     () -> stage_qc(ctx, o[:load]))
step!(:filter, () -> stage_filter(ctx, o[:load]))
step!(:ica,    () -> stage_ica(ctx, o[:filter]); filter=o[:filter])
step!(:montage, () -> stage_montage(ctx, o[:ica], o[:qc]); qc=o[:qc])
step!(:segmentation, () -> stage_segmentation(ctx, o[:montage], o[:ica]))
step!(:spectral, () -> stage_spectral(ctx, o[:segmentation]))
step!(:connectivity, () -> stage_connectivity(ctx, o[:segmentation], o[:montage]); segmentation=o[:segmentation])
step!(:surrogates, () -> stage_surrogates(ctx, o[:connectivity]))
step!(:finalize, () -> stage_finalize(ctx, o[:load], o[:qc], o[:filter], o[:montage], o[:segmentation],
                                      o[:spectral], o[:connectivity]))
close(ctx.log_io)

store = CK.results_dir(ctx.cfg)
dir   = joinpath(ctx.export_dir, "checkpoints")
lab_obj = Dict(:load => o[:load].rec, :qc => nothing, :filter => o[:filter].rec, :ica => o[:ica].rec,
               :montage => o[:montage].rec, :segmentation => o[:segmentation].epochs, :spectral => o[:spectral],
               :connectivity => nothing, :surrogates => nothing, :finalize => nothing)

# Escritura sola: la misma salida, guardada otra vez en carpetas nuevas del almacén
function save_only(s, rep)
    prov = TOML.parsefile(joinpath(dir, "$(s).toml"))
    ex = joinpath(store, "medida", "$(s)-$(rep)")
    c2 = PipelineContext(ctx.cfg, ctx.config_path, ctx.subject_id, ctx.session_id, ctx.task, ctx.run,
                         ctx.condition, ex, ctx.dpi, devnull, ctx.t0)
    refs = Dict{Any,Any}()
    for r in prov["references"]
        if r["kind"] == "ica_decomposition"
            refs[:ica_decomposition] = (path = joinpath(store, r["path"]), identity = r["identity"])
        else
            # La nula tiene que estar dentro de la carpeta de la etapa: se copia
            p = joinpath(ex, r["path"]); mkpath(dirname(p))
            cp(joinpath(ctx.export_dir, r["path"]), p; force=true)
            refs[r["band"]] = p
        end
    end
    out = s == :finalize ? nothing : o[s]
    kw = s == :connectivity ? (epochs_from_segmentation = true,) : NamedTuple()
    return @elapsed save_checkpoint(c2, s, out; key=prov["key"], store,
                                    consumed=Dict(c => dir for c in CK.CHECKPOINT_SPEC[s].consumes),
                                    references=refs, facts=prov["facts"], kw...)
end

stats(v) = (mediana = median(v), maximo = maximum(v))
filas = Dict{String,Any}[]
for s in CHECKPOINT_STAGES
    w = [save_only(s, r) for r in 1:REPS]
    l = [@elapsed(load_checkpoint(dir, s; store)) for _ in 1:REPS]
    v = [@elapsed(verify_checkpoint(dir, s; store)) for _ in 1:REPS]
    push!(filas, Dict{String,Any}(
        "etapa" => String(s), "jls_bytes" => filesize(joinpath(dir, "$(s).jls")),
        "toml_bytes" => filesize(joinpath(dir, "$(s).toml")), "lab_bin_bytes" => bytes_lab(lab_obj[s]),
        "calculo_s" => t_stage[s], "orquestador_s" => t_ck[s],
        "escritura_s" => stats(w).mediana, "escritura_max_s" => stats(w).maximo,
        "lectura_s" => stats(l).mediana, "lectura_max_s" => stats(l).maximo,
        "verificacion_s" => stats(v).mediana, "verificacion_max_s" => stats(v).maximo))
end
null_bytes = sum(filesize(joinpath(ctx.export_dir, "cache", "surrogate", f))
                 for f in readdir(joinpath(ctx.export_dir, "cache", "surrogate")) if endswith(f, ".jls"); init=0)
ica_cache_bytes = filesize(joinpath(ctx.export_dir, "cache", "ica_result.jls"))

# Arranque de Julia por proceso: en caliente (ya precompilado), cinco veces
arranque = [@elapsed(run(pipeline(`$(Base.julia_cmd()) --project=$(dirname(@__DIR__)) -e "using NeuroMIND"`;
                                  stdout=devnull, stderr=devnull))) for _ in 1:REPS]

f(x) = round(x; digits=3)
mb(x) = round(x / 2^20; digits=2)
println()
println("| Etapa | .jls (MB) | .toml (kB) | Lab hoy (MB) | Cálculo (s) | Escritura (s) | Orquestador (s) | Lectura verificada (s) | Verificación sin deserializar (s) |")
println("|---|---|---|---|---|---|---|---|---|")
for r in filas
    println("| $(r["etapa"]) | $(mb(r["jls_bytes"])) | $(round(r["toml_bytes"]/1024; digits=1)) | $(mb(r["lab_bin_bytes"])) | " *
            "$(f(r["calculo_s"])) | $(f(r["escritura_s"])) (máx. $(f(r["escritura_max_s"]))) | $(f(r["orquestador_s"])) | " *
            "$(f(r["lectura_s"])) (máx. $(f(r["lectura_max_s"]))) | $(f(r["verificacion_s"])) (máx. $(f(r["verificacion_max_s"]))) |")
end
tot_ck  = sum(r["jls_bytes"] + r["toml_bytes"] for r in filas)
tot_lab = sum(r["lab_bin_bytes"] for r in filas)
println()
println("Total checkpoints: $(mb(tot_ck)) MB · Lab hoy (.bin de las mismas etapas): $(mb(tot_lab)) MB · cociente $(round(tot_ck / tot_lab; digits=3))")
println("Archivos referenciados, que ya escribe el motor: caché de ICA $(mb(ica_cache_bytes)) MB · distribuciones nulas $(mb(null_bytes)) MB")
println("Arranque de Julia con `using NeuroMIND`: mediana $(f(median(arranque))) s · máx. $(f(maximum(arranque))) s")
println("Máquina: $(Sys.MACHINE) · Julia $(VERSION) · hilos BLAS $(NeuroMIND.LinearAlgebra.BLAS.get_num_threads()) · hilos Julia $(Threads.nthreads())")

if !isempty(salida)
    open(salida, "w") do io
        TOML.print(io, Dict("parametros" => Dict("segundos" => seconds, "canales" => nchan, "fs" => FS,
                                                 "n_surrogates" => nsurr, "repeticiones" => REPS),
                            "etapas" => filas,
                            "totales" => Dict("checkpoints_bytes" => tot_ck, "lab_bin_bytes" => tot_lab,
                                              "ica_cache_bytes" => ica_cache_bytes, "null_bytes" => null_bytes),
                            "arranque_s" => arranque,
                            "entorno" => Dict("machine" => String(Sys.MACHINE), "julia" => string(VERSION),
                                              "blas_threads" => NeuroMIND.LinearAlgebra.BLAS.get_num_threads(),
                                              "threads" => Threads.nthreads())); sorted=true)
    end
end
