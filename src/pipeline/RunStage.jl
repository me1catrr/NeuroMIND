# ───────────────────────────────────────────────────────────────
#  Fichero    src/pipeline/RunStage.jl
#  Modificado 14-11-2026
# ───────────────────────────────────────────────────────────────
#
# Ejecución por etapas (PR 3, entrega 2f-2; plan de 2f, apartado 3.8).
#
#   stage_graph()                         qué consume y qué lee cada etapa
#   stage_key(stage, config_path; …)      identidad de una etapa sin ejecutarla
#   run_stage(stage, config_path; …)      carga y verifica sus entradas, ejecuta
#                                         la etapa y guarda su checkpoint
#
# La identidad tiene una sola fuente, `_stage_identity`, que usan `stage_key`,
# `run_stage` y el orquestador con checkpoints (E-2). Nada de lo que hay aquí
# cambia una función de etapa: `run_stage` llama a las mismas `stage_*` que el
# pipeline completo, con entradas cargadas de checkpoints verificados.

"""
    StageError(stage, code, message)

`run_stage` o `stage_key` no pueden ejecutar o identificar la etapa con lo que
han recibido. No es un error de los datos de un checkpoint (eso es
`CheckpointError`), sino de la petición. Códigos:

- `:CHECKPOINT_EXISTE`: la carpeta ya tiene el checkpoint de la etapa; no se
  sobrescribe (apartado 3.2);
- `:ENTRADA_FALTA`: falta una entrada que la etapa necesita;
- `:ICA_LABELS_OTRA_DESCOMPOSICION`: la decisión manual se tomó sobre otra
  descomposición; no se aplica nada (decisiones 065 y 084);
- `:ICA_LABELS_SIN_IDENTIDAD`: no se puede calcular la identidad de la
  descomposición, así que la decisión no se puede comprobar;
- `:ICA_LABELS_NO_ENTREGADAS`: hay un archivo de etiquetas en la carpeta que
  no se ha entregado como decisión verificada;
- `:ICA_CACHE_NO_VERIFICADA`: la carpeta ya tiene una caché de ICA que no se
  ha entregado ni verificado;
- `:SEGMENTACION_ROTULO`: el rótulo de entrada de la segmentación no
  corresponde a la lista de rechazados de la ICA consumida.
"""
struct StageError <: Exception
    stage::Symbol
    code::Symbol
    message::String
end

Base.showerror(io::IO, e::StageError) = print(io, "Etapa $(e.stage): $(e.code) · $(e.message)")

"""
Entradas que `run_stage` carga para cada etapa: lo que consume y, en el
montaje, el QC, cuya lista de canales marcados entra en su identidad por su
contenido (apartado 3.4) sin que el montaje consuma el checkpoint del QC.
"""
const STAGE_INPUTS = Dict{Symbol,Vector{Symbol}}(
    s => (s == :montage ? [:ica, :qc] : copy(CHECKPOINT_SPEC[s].consumes)) for s in CHECKPOINT_STAGES)

"""
    stage_graph() -> Vector{NamedTuple}

Las diez etapas en orden de ejecución, con lo que consume cada una
(`consumes`), las secciones de la configuración que lee (`sections`, incluida
`output` donde cambia sus figuras), sus entradas propias de identidad
(`inputs`) y los checkpoints que `run_stage` necesita recibir (`needs`). Es la
única fuente para Lab: `schemas/etapas.json` se contrasta con ella.
"""
stage_graph() = [(stage = s,
                  consumes = copy(CHECKPOINT_SPEC[s].consumes),
                  sections = copy(CHECKPOINT_SPEC[s].sections),
                  inputs = copy(get(CHECKPOINT_INPUTS, s, Symbol[])),
                  needs = copy(STAGE_INPUTS[s])) for s in CHECKPOINT_STAGES]

# ── Identidad: una sola fuente ───────────────────────────────

"""
Identidad de `stage` y las entradas que la forman. La usan `stage_key`,
`run_stage` y el orquestador, así que las tres calculan lo mismo con lo mismo.

- `consumed_keys`: etapa consumida ⇒ su identidad;
- `filter`: la salida del filtrado (ICA: identidad de la descomposición);
- `qc`: la salida del QC (montaje: canales marcados, por su contenido);
- `labels_file`: el archivo de la decisión manual de ICA, o `nothing`;
- `toml_sha256`: cierre, la huella de los bytes del TOML si ya se calculó al
  empezar (el orquestador la calcula una vez); si no, se calcula aquí.

Devuelve `(key, inputs, decomposition, manual)`.
"""
function _stage_identity(stage::Symbol, cfg::PipelineConfig, ids, config_path::String,
                         consumed_keys::AbstractDict; filter=nothing, qc=nothing, labels_file=nothing,
                         toml_sha256=nothing)
    decomposition = nothing; manual = nothing
    inputs = if stage == :load
        checkpoint_load_inputs(cfg, ids.subject_id, ids.session_id, ids.task, ids.run)
    elseif stage == :ica
        decomposition = try _ica_cache_identity(filter.rec, cfg) catch; nothing end
        manual = labels_file === nothing ? nothing : _read_ica_labels_csv(labels_file)
        (decomposition = decomposition === nothing ? "sin identidad" : decomposition.key,
         manual_rejected = manual === nothing ? "sin decisión manual" : manual)
    elseif stage == :montage
        (qc_bad_channels = qc.bad_channels,)
    elseif stage == :finalize
        (toml_sha256 = toml_sha256 === nothing ? _file_sha256(config_path) : String(toml_sha256),)
    else
        NamedTuple()
    end
    key = checkpoint_key(stage, cfg, consumed_keys; inputs...)
    return (key = key, inputs = inputs, decomposition = decomposition, manual = manual)
end

# ── Entradas ─────────────────────────────────────────────────

"Normaliza `inputs` (etapa ⇒ carpeta `checkpoints/`) y valida cada carpeta dentro del almacén."
function _input_dirs(stage::Symbol, inputs, store::String)
    dirs = Dict{Symbol,String}()
    for (k, v) in inputs
        s = Symbol(k)
        s in CHECKPOINT_STAGES || throw(ArgumentError("inputs: etapa desconocida $(k)"))
        dirs[s] = _check_inside(stage, store, String(v); what="entrada $(s)")
    end
    missing_ = setdiff(STAGE_INPUTS[stage], keys(dirs))
    isempty(missing_) ||
        throw(StageError(stage, :ENTRADA_FALTA, "faltan las entradas $(join(sort(collect(missing_)), ", "))"))
    return dirs
end

"""
Carga y verifica las entradas de `stage`. Las carpetas de `dirs` que no son
entradas directas sirven para encontrar lo que consumen las entradas (E-1): se
pasan como `consumed_dirs`, y todo se sigue comprobando por identidad y huellas.

Las entradas de `deserialize` se cargan enteras con `load_checkpoint`; las demás
solo se verifican, sin deserializar. Devuelve `(outputs, keys)`.
"""
function _load_inputs(stage::Symbol, dirs::Dict{Symbol,String}, store::String; deserialize)
    outputs = Dict{Symbol,Any}(); keys_ = Dict{Symbol,String}()
    for s in STAGE_INPUTS[stage]
        cdirs = Dict(c => dirs[c] for c in CHECKPOINT_SPEC[s].consumes if haskey(dirs, c))
        if s in deserialize
            outputs[s] = load_checkpoint(dirs[s], s; store, consumed_dirs=cdirs)
        else
            ok, reason, msg = verify_checkpoint(dirs[s], s; store, consumed_dirs=cdirs)
            ok || throw(CheckpointError(reason, s, msg))
        end
        keys_[s] = TOML.parsefile(_prov(dirs[s], s))["key"]
    end
    return outputs, keys_
end

"""
P-10: la identidad de la descomposición que registra una ICA consumida se
recalcula desde la señal filtrada que esa ICA consumió, con la configuración de
quien la consume. Si no coincide, la ICA se calculó con otra configuración de
ICA o en otra máquina (apartado 3.3: un cambio de máquina invalida la ICA y lo
posterior): `OTRA_CONFIGURACION`.
"""
function _check_consumed_decomposition(stage::Symbol, cfg::PipelineConfig, ica_dir::String,
                                       dirs::Dict{Symbol,String}, store::String)
    prov = TOML.parsefile(_prov(ica_dir, :ica))
    registered = String(get(get(prov, "facts", Dict{String,Any}()), "decomposition_identity", ""))
    isempty(registered) && return nothing
    entry = only(c for c in prov["consumed"] if c["stage"] == "filter")
    fdir = haskey(dirs, :filter) ? dirs[:filter] :
           dirname(_inside_store(stage, store, store, entry["path"]; what="filtrado de la ICA consumida"))
    fprov = TOML.parsefile(_prov(fdir, :filter))
    fprov["key"] == entry["key"] ||
        _alterado(stage, "el filtrado entregado no es el que consumió la ICA")
    cdirs = Dict(c => dirs[c] for c in CHECKPOINT_SPEC[:filter].consumes if haskey(dirs, c))
    filt = load_checkpoint(fdir, :filter; store, consumed_dirs=cdirs)
    now_ = try _ica_cache_identity(filt.rec, cfg).key catch; "" end
    now_ == registered ||
        throw(CheckpointError(:OTRA_CONFIGURACION, stage,
            "la ICA consumida se calculó con otra configuración de ICA o en otra máquina " *
            "(identidad de la descomposición $(first(registered, 12))…, aquí $(first(now_, 12))…)"))
    return nothing
end

"""
Lo que `run_stage` y `stage_key` comparten antes de ejecutar: configuración,
sujeto, entradas verificadas y decisión manual comprobada. Sin efectos sobre el
disco.
"""
function _prepare_stage(stage::Symbol, config_path::String, inputs, ica_labels, store::String)
    stage in CHECKPOINT_STAGES || throw(ArgumentError("Etapa sin checkpoint: $(stage)"))
    store = String(store)
    cp_ = _check_inside(stage, store, abspath(config_path); what="configuración")
    ica_labels === nothing || stage == :ica || throw(ArgumentError("ica_labels solo se entrega a :ica"))
    cfg = load_ss_config(cp_)
    ids = _context_ids(cfg, cp_)
    dirs = _input_dirs(stage, inputs, store)
    deser = stage == :ica ? (:filter,) : stage == :montage ? (:qc,) : ()
    outputs, keys_ = _load_inputs(stage, dirs, store; deserialize = deser)
    stage in (:montage, :segmentation) && _check_consumed_decomposition(stage, cfg, dirs[:ica], dirs, store)

    labels_file = nothing
    if ica_labels !== nothing
        path, identity = ica_labels
        labels_file = _check_inside(stage, store, abspath(String(path)); what="decisión de ICA")
        isfile(labels_file) || throw(ArgumentError("ica_labels: no existe $(path)"))
        decomposition = try _ica_cache_identity(outputs[:filter].rec, cfg) catch; nothing end
        decomposition === nothing &&
            throw(StageError(:ica, :ICA_LABELS_SIN_IDENTIDAD,
                "no se puede calcular la identidad de la descomposición: la decisión no se puede comprobar"))
        String(identity) == decomposition.key ||
            throw(StageError(:ica, :ICA_LABELS_OTRA_DESCOMPOSICION,
                "la decisión se tomó sobre la descomposición $(first(String(identity), 12))… " *
                "y esta etapa calcula $(first(decomposition.key, 12))…; no se aplica nada"))
    end
    consumed_keys = Dict(c => keys_[c] for c in CHECKPOINT_SPEC[stage].consumes)
    id = _stage_identity(stage, cfg, ids, cp_, consumed_keys;
                         filter = get(outputs, :filter, nothing), qc = get(outputs, :qc, nothing), labels_file)
    return (; cfg, ids, config_path = cp_, dirs, keys = keys_, id, labels_file, store)
end

"""
    stage_key(stage, config_path; inputs=Dict(), ica_labels=nothing, store) -> String

La identidad que `run_stage` registrará para `stage` con esta configuración y
estas entradas, sin ejecutar la etapa ni escribir nada (E-2). Verifica las
entradas como `run_stage` y calcula la identidad con la misma función. No
depende de dónde estén las carpetas: solo de la configuración efectiva, de las
identidades consumidas y de las entradas propias de la etapa.

Sirve para reutilizar: quien ya tiene un checkpoint de la etapa lo carga con
`load_checkpoint(dir, stage; store, expected_key = stage_key(…))`, y si no es el
cálculo pedido obtiene `OTRA_CONFIGURACION`.
"""
function stage_key(stage::Symbol, config_path::AbstractString; inputs=Dict{Symbol,String}(),
                   ica_labels=nothing, store::AbstractString)
    return _prepare_stage(stage, String(config_path), inputs, ica_labels, String(store)).id.key
end

# ── Caché de ICA entregada ───────────────────────────────────

"""
Decide si se usa la caché de ICA de otra ejecución (`src`, la carpeta
`checkpoints/` de una ICA anterior, dentro del almacén). Solo se usa si su
checkpoint registra la huella de bytes de `ica_result.jls`, la identidad es la
de esta señal y los bytes coinciden con la huella. No deserializa nada.
Devuelve `(use, reason, bytes, key)`; si no se usa, la ICA se recalcula y
`reason` dice por qué.
"""
function _check_ica_cache(src::String, store::String, decomposition)
    dir = _check_inside(:ica, store, abspath(src); what="caché de ICA entregada")
    no(reason) = (use = false, reason = reason, bytes = UInt8[], key = "")
    decomposition === nothing && return no("no se puede calcular la identidad de esta descomposición")
    prov = try
        _read_provenance(dir, :ica, store)
    catch e
        e isa CheckpointError || rethrow()
        return no("su checkpoint no se puede usar: $(e.reason) · $(e.message)")
    end
    r = _ref(prov, "ica_decomposition")
    r === nothing && return no("su checkpoint no registra la huella de bytes de la caché")
    get(r, "identity", "") == decomposition.key || return no("es de otra señal o de otra configuración de ICA")
    path = _inside_store(:ica, store, dirname(dir), r["path"]; what="caché de ICA entregada")
    bytes = try
        _read_checked(:ica, path, r["file_sha256"]; what="la caché de ICA entregada")
    catch e
        e isa CheckpointError || rethrow()
        return no(e.message)
    end
    return (use = true, reason = "", bytes = bytes, key = decomposition.key)
end

# ── run_stage ────────────────────────────────────────────────

"""
    run_stage(stage, config_path; inputs=Dict(), ica_labels=nothing, ica_cache=nothing,
              store, commit=nothing) -> (key, dir, ica_cache)

Ejecuta una etapa y guarda su checkpoint (plan de 2f, apartado 3.8). Pensada
para un proceso nuevo por etapa.

- `config_path`: el TOML de la etapa. Su carpeta de salida (`export_dir`) es la
  de la etapa.
- `inputs`: etapa ⇒ carpeta `checkpoints/` de cada entrada (`stage_graph()`,
  campo `needs`). Puede incluir otras etapas: sirven para encontrar lo que
  consumen las entradas si se han trasladado (E-1).
- `ica_labels = (ruta, identidad)`: la decisión del investigador y la identidad
  de la descomposición sobre la que se tomó. Si no es la de esta ejecución, o no
  se puede comprobar, la etapa falla sin aplicar nada.
- `ica_cache`: la carpeta `checkpoints/` de una ICA anterior. Su caché solo se
  usa si su identidad es la de esta señal y la huella de sus bytes es la
  registrada; si no, la ICA se recalcula y se dice.
- `store`: el almacén autorizado. Todas las rutas recibidas tienen que estar
  dentro.

Antes de abrir el registro comprueba que la carpeta no tiene ya el checkpoint
de la etapa. El registro se cierra también ante una excepción, que se vuelve a
lanzar sin cambiarla, y entonces no queda ningún checkpoint (apartado 3.9).
"""
function run_stage(stage::Symbol, config_path::AbstractString; inputs=Dict{Symbol,String}(),
                   ica_labels=nothing, ica_cache=nothing, store::AbstractString, commit=nothing)
    ica_cache === nothing || stage == :ica || throw(ArgumentError("ica_cache solo se entrega a :ica"))
    p = _prepare_stage(stage, String(config_path), inputs, ica_labels, String(store))
    ckdir = joinpath(p.ids.export_dir, "checkpoints")
    _check_inside(stage, p.store, ckdir; what="carpeta de la etapa")
    # Antes de cualquier escritura: no se sobrescribe un checkpoint ni se
    # aplica nada que no se haya entregado y verificado
    (ispath(_jls(ckdir, stage)) || ispath(_prov(ckdir, stage))) &&
        throw(StageError(stage, :CHECKPOINT_EXISTE, "ya hay un checkpoint de $(stage) en $(ckdir): no se sobrescribe"))
    cache = nothing
    if stage == :ica
        present = filter(isfile, _ica_manual_label_candidates(p.cfg, p.ids.subject_id, p.ids.session_id, p.ids.condition))
        isempty(present) || throw(StageError(:ica, :ICA_LABELS_NO_ENTREGADAS,
            "la carpeta ya tiene $(join(basename.(present), ", ")) y no se ha entregado como decisión verificada"))
        stale = [f for f in _ICA_CACHE_FILES if ispath(joinpath(p.ids.export_dir, split(f, '/')...))]
        isempty(stale) || throw(StageError(:ica, :ICA_CACHE_NO_VERIFICADA,
            "la carpeta ya tiene $(join(stale, ", ")), que no se ha entregado ni verificado"))
        ica_cache === nothing || (cache = _check_ica_cache(String(ica_cache), p.store, p.id.decomposition))
    end

    ctx = pipeline_context(p.config_path)
    failed = true
    try
        run = _CheckpointRun(ctx, p.store, ckdir, Dict{Symbol,String}(p.keys), _file_sha256(p.config_path),
                             Dict{String,Tuple{Int,Float64}}(), Dict{String,String}(),
                             Dict{Symbol,String}(p.dirs), commit === nothing ? nothing : String(commit),
                             Dict{String,Any}())
        run.inventory = _inventory(run)
        if p.labels_file !== nothing
            cp(p.labels_file, joinpath(ctx.export_dir, "ica_labels.csv"))
        end
        if cache !== nothing
            if cache.use
                write(joinpath(ctx.export_dir, "cache", "ica_result.jls"), cache.bytes)
                write(joinpath(ctx.export_dir, "cache", "ica_config.hash"), cache.key)
                run.notes["ica_cache"] = "entregada y usada"
                _log(ctx.log_io, "Caché de ICA entregada: verificada y usada")
            else
                run.notes["ica_cache"] = "entregada y no usada: $(cache.reason)"
                _log(ctx.log_io, "Caché de ICA entregada: no se usa ($(cache.reason)); se recalcula")
                println("  ⚠ Caché de ICA entregada: no se usa ($(cache.reason)); se recalcula")
            end
        end
        L = Dict{Symbol,Any}()
        for s in STAGE_INPUTS[stage]
            cdirs = Dict(c => p.dirs[c] for c in CHECKPOINT_SPEC[s].consumes if haskey(p.dirs, c))
            L[s] = load_checkpoint(p.dirs[s], s; store=p.store, consumed_dirs=cdirs)
        end
        out = _call_stage(stage, ctx, L)
        kw = stage == :ica ? (filter = L[:filter],) :
             stage == :montage ? (qc = L[:qc],) :
             stage == :segmentation ? (ica = L[:ica],) :
             stage == :connectivity ? (segmentation = L[:segmentation],) : NamedTuple()
        key = _checkpoint_after!(run, stage, out; kw...)
        key == p.id.key || error("run_stage($(stage)): la identidad registrada no es la calculada antes de ejecutar")
        failed = false
        return (key = key, dir = ckdir, ica_cache = cache === nothing ? nothing : (used = cache.use, reason = cache.reason))
    catch e
        _log_failure(ctx.log_io, stage, e)
        rethrow()
    finally
        failed || _log(ctx.log_io, "\n✓ Etapa $(stage) completada")
        isopen(ctx.log_io) && close(ctx.log_io)
    end
end

function _call_stage(stage::Symbol, ctx::PipelineContext, L::Dict{Symbol,Any})
    stage == :load         && return stage_load(ctx)
    stage == :qc           && return stage_qc(ctx, L[:load])
    stage == :filter       && return stage_filter(ctx, L[:load])
    stage == :ica          && return stage_ica(ctx, L[:filter])
    stage == :montage      && return stage_montage(ctx, L[:ica], L[:qc])
    stage == :segmentation && return stage_segmentation(ctx, L[:montage], L[:ica])
    stage == :spectral     && return stage_spectral(ctx, L[:segmentation])
    stage == :connectivity && return stage_connectivity(ctx, L[:segmentation], L[:montage])
    stage == :surrogates   && return stage_surrogates(ctx, L[:connectivity])
    stage == :finalize     && return stage_finalize(ctx, L[:load], L[:qc], L[:filter], L[:montage],
                                                    L[:segmentation], L[:spectral], L[:connectivity])
    error("Etapa sin checkpoint: $(stage)")
end

"Anota en el registro una etapa que falla: la etapa y el mensaje (apartado 3.9). No lanza."
function _log_failure(io::IO, stage, e)
    isopen(io) || return nothing
    try
        _log(io, "\n✗ La etapa $(stage) ha fallado: $(sprint(showerror, e))")
        _log(io, "  No se ha guardado ningún checkpoint de esta etapa.")
    catch
    end
    return nothing
end
