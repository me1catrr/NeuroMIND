# ───────────────────────────────────────────────────────────────
#  Fichero    src/pipeline/Checkpoints.jl
#  Modificado 14-11-2026
# ───────────────────────────────────────────────────────────────
#
# Checkpoints comunes de las etapas del pipeline individual (PR 3, entrega 2f-1).
#
# Un checkpoint es la salida completa de una etapa, con su tipo del PR 1, y se
# puede verificar antes de reutilizarla. Dos archivos por etapa, en
# `<export_dir>/checkpoints/`:
#
#   <etapa>.jls   la salida, en el envoltorio de `save_result` (formato, etapa,
#                 sujeto, sesión, datos) y sin hora
#   <etapa>.toml  la procedencia, sin fecha: identidad, compatibilidad,
#                 implementación, configuración efectiva, checkpoints
#                 consumidos, archivos referenciados, hechos y productos
#
# Tres conceptos separados (plan de 2f, apartado 3.3):
#
#   identidad (`key`)        ¿es el mismo cálculo?   → OTRA_CONFIGURACION
#   compatibilidad           ¿se puede leer aquí?    → INCOMPATIBLE
#   procedencia del entorno  ¿dónde y con qué?       → se registra, no invalida
#
# Ningún `.jls` se deserializa antes de comprobar la huella SHA-256 de sus
# bytes, tampoco los archivos externos que referencia: la descomposición de ICA
# (`cache/ica_result.jls`) y las distribuciones nulas de surrogates
# (`cache/surrogate/null_distribution_<banda>.jls`, decisión 130, F6). Las
# huellas verifican la integridad respecto a lo registrado; no hacen seguro ni
# auténtico un archivo que alguien haya podido reescribir junto con su
# procedencia (apartado 3.5).
#
# Este archivo no cambia ninguna función de etapa, ni `save_result` ni
# `load_result`. El orquestador solo escribe checkpoints si se le pide.

using SHA: sha256, SHA256_CTX, update!, digest!

"Versión del formato de los checkpoints. Un checkpoint de otra versión es `INCOMPATIBLE`."
const CHECKPOINT_FORMAT = 1

"Versión del envoltorio del `.jls`."
const CHECKPOINT_WRAPPER = 1

"Etapas con checkpoint, en orden de ejecución."
const CHECKPOINT_STAGES = (:load, :qc, :filter, :ica, :montage, :segmentation,
                           :spectral, :connectivity, :surrogates, :finalize)

"""
Lo que consume cada etapa y las secciones de la configuración que lee
(plan de 2f, apartado 3.4). Las secciones son las que leen la etapa y sus
auxiliares en el código; `output` es `[output]` (`figure_dpi`), que cambia las
figuras de las etapas que las escriben. Además:

- carga: huellas de `.vhdr`, `.vmrk`, `.eeg` y `electrodes.tsv`; sujeto,
  sesión, tarea y `run`;
- ICA: identidad de la descomposición (PR 2) y la decisión manual consumida,
  por los componentes que rechaza (decisión 129), o su ausencia;
- montaje: la lista de canales marcados por el QC, por su contenido, no la
  identidad del QC;
- cierre: la huella de los bytes completos del TOML enviado (decisión 130, F12).
"""
const CHECKPOINT_SPEC = Dict{Symbol,NamedTuple{(:consumes, :sections),Tuple{Vector{Symbol},Vector{String}}}}(
    :load         => (consumes = Symbol[],             sections = String[]),
    :qc           => (consumes = [:load],              sections = ["qc", "bands"]),
    :filter       => (consumes = [:load],              sections = ["filtering", "output"]),
    :ica          => (consumes = [:filter],            sections = ["ica"]),
    :montage      => (consumes = [:ica],               sections = ["montage"]),
    :segmentation => (consumes = [:montage, :ica],     sections = ["segmentation", "baseline", "artifact_rejection"]),
    :spectral     => (consumes = [:segmentation],      sections = ["spectral", "bands", "output"]),
    :connectivity => (consumes = [:segmentation, :montage], sections = ["connectivity", "bands", "graph", "output"]),
    :surrogates   => (consumes = [:connectivity],      sections = ["surrogates", "connectivity", "bands", "output"]),
    :finalize     => (consumes = [:load, :qc, :filter, :montage, :segmentation, :spectral, :connectivity],
                      sections = ["segmentation", "artifact_rejection", "output"]),
)

"Entradas que, además de configuración y consumidos, forman la identidad de cada etapa."
const CHECKPOINT_INPUTS = Dict{Symbol,Vector{Symbol}}(
    :load     => [:input_files, :subject, :session, :task, :run],
    :ica      => [:decomposition, :manual_rejected],
    :montage  => [:qc_bad_channels],
    :finalize => [:toml_sha256],
)

# ── Errores ──────────────────────────────────────────────────

"""
    CheckpointError(reason, stage, message)

Rechazo de un checkpoint. `reason` es `:INCOMPATIBLE` (otra versión de formato,
de envoltorio o de Julia), `:OTRA_CONFIGURACION` (la identidad no es la
pedida) o `:ALTERADO` (falta, ruta fuera del almacén, enlace simbólico, huella
que no coincide, archivo ilegible, consumido distinto del declarado o
invariante que falla). Solo `:ALTERADO` es un error de los datos.
"""
struct CheckpointError <: Exception
    reason::Symbol
    stage::Symbol
    message::String
end

Base.showerror(io::IO, e::CheckpointError) =
    print(io, "Checkpoint $(e.stage): $(e.reason) · $(e.message)")

_alterado(stage, msg)    = throw(CheckpointError(:ALTERADO, stage, msg))
_incompatible(stage, msg) = throw(CheckpointError(:INCOMPATIBLE, stage, msg))

# ── Huella canónica ──────────────────────────────────────────
#
# La misma regla que la clave de la caché de ICA (`_ica_put`): no pasa por
# texto; enteros como Int64, reales como sus bits IEEE-754, cadenas como
# longitud y bytes UTF-8, todo en little-endian. Cada valor lleva una etiqueta
# de tipo, así que dos estructuras distintas no se confunden.
#
# Fuera de la huella de contenido quedan solo:
# - `RecordingMeta.bids_path`: dónde estaba la entrada, no lo calculado (la
#   identidad de la carga ya lleva la huella de los archivos);
# - las claves de reloj de los diccionarios (`_CLOCK_KEYS`), como
#   `optimize_s` en los diagnósticos de ICA: dos cálculos iguales las
#   escriben distintas.

const _CLOCK_KEYS = ("optimize_s", "duration_s", "timestamp", "saved_at")

mutable struct _Hasher
    ctx::SHA256_CTX
end
_Hasher() = _Hasher(SHA256_CTX())
_hex(h::_Hasher) = bytes2hex(digest!(h.ctx))

_raw!(h::_Hasher, b::AbstractVector{UInt8}) = (update!(h.ctx, b); h)
_tag!(h::_Hasher, t::String) = _raw!(h, codeunits(t))
_i64!(h::_Hasher, x::Integer) = _raw!(h, reinterpret(UInt8, [htol(Int64(x))]))
_str!(h::_Hasher, s::AbstractString) = (b = codeunits(String(s)); _i64!(h, length(b)); _raw!(h, b))

"Bytes de un arreglo de tipo elemental, en little-endian, sin copiarlo si la máquina ya lo es."
function _bits!(h::_Hasher, A::Array{T}) where {T<:Union{Float64,Float32,Int64,Int32,UInt8,Bool}}
    isempty(A) && return h
    if ENDIAN_BOM == 0x04030201
        GC.@preserve A update!(h.ctx, unsafe_wrap(Array, Ptr{UInt8}(pointer(A)), sizeof(A)))
    else
        _raw!(h, reinterpret(UInt8, vec(htol.(A))))
    end
    return h
end

_canon!(h::_Hasher, ::Nothing) = _tag!(h, "0")
_canon!(h::_Hasher, ::Missing) = _tag!(h, "M")
_canon!(h::_Hasher, x::Bool) = (_tag!(h, "B"); _raw!(h, UInt8[x]))
_canon!(h::_Hasher, x::Integer) = (_tag!(h, "I"); _i64!(h, x))
_canon!(h::_Hasher, x::AbstractFloat) = (_tag!(h, "F"); _raw!(h, reinterpret(UInt8, [htol(reinterpret(UInt64, Float64(x)))])))
_canon!(h::_Hasher, s::AbstractString) = (_tag!(h, "S"); _str!(h, s))
_canon!(h::_Hasher, s::Symbol) = (_tag!(h, "Y"); _str!(h, String(s)))

function _canon!(h::_Hasher, A::AbstractArray)
    _tag!(h, "A"); _str!(h, string(eltype(A))); _i64!(h, ndims(A))
    foreach(n -> _i64!(h, n), size(A))
    if A isa Array && eltype(A) <: Union{Float64,Float32,Int64,Int32,UInt8,Bool}
        _bits!(h, A)
    elseif A isa BitArray
        _bits!(h, Array{Bool}(A))
    else
        foreach(x -> _canon!(h, x), A)
    end
    return h
end

function _canon!(h::_Hasher, t::Tuple)
    _tag!(h, "T"); _i64!(h, length(t)); foreach(x -> _canon!(h, x), t); h
end

function _canon!(h::_Hasher, t::NamedTuple)
    _tag!(h, "N"); _i64!(h, length(t))
    for (k, v) in pairs(t)
        _str!(h, String(k)); _canon!(h, v)
    end
    h
end

function _canon!(h::_Hasher, d::AbstractDict)
    ks = sort!([string(k) for k in keys(d) if !(string(k) in _CLOCK_KEYS)])
    byname = Dict(string(k) => v for (k, v) in d)
    _tag!(h, "D"); _i64!(h, length(ks))
    for k in ks
        _str!(h, k); _canon!(h, byname[k])
    end
    h
end

function _canon!(h::_Hasher, df::DataFrame)
    _tag!(h, "DF"); _i64!(h, ncol(df)); _i64!(h, nrow(df))
    for n in names(df)
        _str!(h, n); _canon!(h, collect(df[!, n]))
    end
    h
end

# Tipos del motor: por campos, en su orden, con el nombre del tipo.
function _canon!(h::_Hasher, x::T) where {T}
    isstructtype(T) || error("Huella canónica: tipo no previsto $(T)")
    _tag!(h, "R"); _str!(h, string(nameof(T))); _i64!(h, fieldcount(T))
    for f in fieldnames(T)
        T === RecordingMeta && f === :bids_path && continue
        _str!(h, String(f)); _canon!(h, getfield(x, f))
    end
    h
end

"""
    checkpoint_content_sha256(x) -> String

Huella canónica del contenido de una salida de etapa: lo calculado y lo
aplicado. Dos salidas iguales campo a campo (con los NaN en su sitio) dan la
misma huella, aunque procedan de procesos o carpetas distintos.
"""
checkpoint_content_sha256(x) = _hex(_canon!(_Hasher(), x))

_file_sha256(path::AbstractString) = bytes2hex(open(sha256, path))

# ── Implementación ───────────────────────────────────────────

const _IMPLEMENTATION = Ref{Any}(nothing)

"""
Huella del código del motor (`src/`, todos los archivos, por ruta y contenido)
y de `Manifest.toml`. Entra en la identidad: cualquier cambio de código o de
dependencias da otro cálculo. Se calcula una vez por proceso.
"""
function checkpoint_implementation()
    _IMPLEMENTATION[] === nothing || return _IMPLEMENTATION[]
    src  = dirname(@__DIR__)
    root = dirname(src)
    files = sort!([relpath(joinpath(d, f), src) for (d, _, fs) in walkdir(src) for f in fs])
    h = _Hasher()
    _tag!(h, "neuromind-src")
    for f in files
        _str!(h, join(splitpath(f), '/'))
        _raw!(h, read(joinpath(src, f)))
    end
    manifest = joinpath(root, "Manifest.toml")
    version = try
        String(TOML.parsefile(joinpath(root, "Project.toml"))["version"])
    catch
        "desconocida"
    end
    _IMPLEMENTATION[] = (
        src_sha256      = _hex(h),
        manifest_sha256 = isfile(manifest) ? _file_sha256(manifest) : "ausente",
        neuromind_version = version,
    )
    return _IMPLEMENTATION[]
end

# ── Configuración efectiva ───────────────────────────────────

"""
Sección de la configuración efectiva tal como la recibe la etapa: la de
`PipelineConfig`, con los valores por defecto de `load_ss_config` aplicados.
Las claves que un perfil ignora siguen dentro: con ellas se puede recalcular
de más, nunca reutilizar de menos.
"""
function checkpoint_section(cfg::PipelineConfig, name::String)
    name == "output" && return cfg.export_cfg
    name == "bands"  && return Dict{String,Any}(k => [v[1], v[2]] for (k, v) in cfg.bands)
    return getfield(cfg, Symbol(name))
end

# Para la procedencia: TOML no admite `nothing`, símbolos ni tuplas.
_tomlable(x::AbstractDict) = Dict{String,Any}(string(k) => _tomlable(v) for (k, v) in x if v !== nothing)
_tomlable(x::Tuple) = Any[_tomlable(v) for v in x]
_tomlable(x::AbstractVector) = Any[_tomlable(v) for v in x]
_tomlable(x::Symbol) = String(x)
_tomlable(x) = x

# ── Identidad ────────────────────────────────────────────────

"""
    checkpoint_key(stage, cfg_or_ctx, consumed=Dict(); inputs...) -> String

Identidad del cálculo de `stage` (plan de 2f, apartados 3.3 y 3.4): etapa,
configuración efectiva de las secciones que lee, identidades de los
checkpoints que consume (`consumed`, etapa ⇒ clave), entradas propias de la
etapa (`CHECKPOINT_INPUTS`) e implementación (`src/` y `Manifest.toml`). La
máquina y la versión de Julia no entran aquí: la máquina llega a ICA y a lo
posterior a través de la identidad de la descomposición (PR 2), y Julia está
en la compatibilidad.

Falla si falta un consumido o una entrada de la etapa, o si sobra alguno:
una entrada que no está en la tabla no puede influir.
"""
checkpoint_key(stage::Symbol, ctx::PipelineContext, consumed::AbstractDict=Dict{Symbol,String}(); inputs...) =
    checkpoint_key(stage, ctx.cfg, consumed; inputs...)

function checkpoint_key(stage::Symbol, cfg::PipelineConfig, consumed::AbstractDict=Dict{Symbol,String}();
                        implementation = checkpoint_implementation(), inputs...)
    spec = get(CHECKPOINT_SPEC, stage, nothing)
    spec === nothing && error("Etapa sin checkpoint: $(stage)")
    want = Set(spec.consumes)
    Set(Symbol.(keys(consumed))) == want ||
        error("checkpoint_key($(stage)): consume $(sort(collect(want))), recibido $(sort(Symbol.(collect(keys(consumed)))))")
    expected_inputs = Set(get(CHECKPOINT_INPUTS, stage, Symbol[]))
    Set(keys(inputs)) == expected_inputs ||
        error("checkpoint_key($(stage)): entradas $(sort(collect(expected_inputs))), recibidas $(sort(collect(keys(inputs))))")

    h = _Hasher()
    _tag!(h, "neuromind-checkpoint-key"); _i64!(h, CHECKPOINT_FORMAT); _str!(h, String(stage))
    _str!(h, implementation.src_sha256); _str!(h, implementation.manifest_sha256)
    for s in spec.sections
        _str!(h, s); _canon!(h, checkpoint_section(cfg, s))
    end
    for c in sort!(collect(spec.consumes))
        _str!(h, String(c)); _str!(h, String(consumed[c]))
    end
    for k in sort!(collect(expected_inputs))
        _str!(h, String(k)); _canon!(h, inputs[k])
    end
    return _hex(h)
end

"""
    checkpoint_load_inputs(cfg, subject, session, task, run) -> NamedTuple

Entradas de la identidad de la carga: huellas de bytes de `.vhdr`, de los
`.eeg` y `.vmrk` que declara y de `electrodes.tsv` («ausente» si no existe),
sujeto, sesión, tarea y `run`. La ruta no entra: mover los datos no cambia el
cálculo; cambiar un byte, sí.
"""
function checkpoint_load_inputs(cfg::PipelineConfig, sid::String, sess::String, task::String, run::Int)
    root = bids_root_dir(cfg)
    vhdr = bids_vhdr_path(root, sid, sess, task, run)
    files = Dict{String,String}("vhdr" => isfile(vhdr) ? _file_sha256(vhdr) : "ausente")
    if isfile(vhdr)
        for line in eachline(vhdr)
            m = match(r"^\s*(DataFile|MarkerFile)\s*=\s*(.+?)\s*$", line)
            m === nothing && continue
            p = joinpath(dirname(vhdr), m.captures[2])
            files[m.captures[1] == "DataFile" ? "eeg" : "vmrk"] = isfile(p) ? _file_sha256(p) : "ausente"
        end
    end
    elec = bids_electrodes_path(root, sid, sess)
    files["electrodes"] = isfile(elec) ? _file_sha256(elec) : "ausente"
    return (input_files = files, subject = sid, session = sess, task = task, run = run)
end

# ── Rutas dentro del almacén ─────────────────────────────────

"""
Resuelve `rel` (relativa a `base`) y comprueba que queda dentro de `store`,
sin `..` que salga y sin enlaces simbólicos en ningún punto del recorrido
desde `store`. Una ruta absoluta en la procedencia se rechaza.
"""
function _inside_store(stage::Symbol, store::String, base::String, rel::AbstractString; what="ruta",
                       where="del almacén autorizado")
    isabspath(rel) && _alterado(stage, "$(what) absoluta en la procedencia: $(rel)")
    return _check_inside(stage, store, normpath(joinpath(base, rel)); what, where)
end

function _check_inside(stage::Symbol, store::String, path::AbstractString; what="ruta",
                       where="del almacén autorizado")
    st = normpath(abspath(store))
    st = endswith(st, "/") ? st[1:end-1] : st
    p  = normpath(abspath(path))
    (p == st || startswith(p, st * "/")) || _alterado(stage, "$(what) fuera $(where): $(path)")
    islink(st) && _alterado(stage, "el almacén es un enlace simbólico: $(store)")
    cur = st
    for part in splitpath(relpath(p, st))
        part == "." && continue
        cur = joinpath(cur, part)
        islink(cur) && _alterado(stage, "$(what) atraviesa un enlace simbólico: $(relpath(cur, st))")
    end
    return p
end

"Ruta de `path` relativa a `base`, con `/`. Puede empezar por `..` si `path` está fuera de `base`."
_rel(base::String, path::String) = join(splitpath(relpath(normpath(abspath(path)), normpath(abspath(base)))), '/')

# ── Lectura verificada ───────────────────────────────────────

"Contador de deserializaciones de checkpoints y archivos referenciados (lo usan las pruebas)."
const CHECKPOINT_DESERIALIZATIONS = Ref(0)

"Deserializa bytes ya comprobados: nunca lee el archivo otra vez."
function _deserialize_checked(bytes::Vector{UInt8})
    CHECKPOINT_DESERIALIZATIONS[] += 1
    return Serialization.deserialize(IOBuffer(bytes))
end

"Lee un archivo y comprueba la huella de sus bytes antes de devolverlos."
function _read_checked(stage::Symbol, path::String, expected::AbstractString; what="archivo")
    isfile(path) || _alterado(stage, "falta $(what): $(path)")
    bytes = try
        read(path)
    catch e
        _alterado(stage, "no se puede leer $(what): $(sprint(showerror, e))")
    end
    bytes2hex(sha256(bytes)) == expected ||
        _alterado(stage, "huella de bytes que no coincide en $(what) ($(basename(path)))")
    return bytes
end

_jls(dir, stage)  = joinpath(dir, "$(stage).jls")
_prov(dir, stage) = joinpath(dir, "$(stage).toml")

"""
Pasos 1 a 3 del apartado 3.5 para un checkpoint: rutas, procedencia y
compatibilidad. Devuelve la procedencia leída. No deserializa nada.
"""
function _read_provenance(dir::String, stage::Symbol, store::String)
    _check_inside(stage, store, dir; what="carpeta del checkpoint")
    jls = _check_inside(stage, store, _jls(dir, stage); what="checkpoint")
    tml = _check_inside(stage, store, _prov(dir, stage); what="procedencia")
    isfile(tml) || _alterado(stage, isfile(jls) ? "falta la procedencia: el checkpoint no está completo" :
                                                  "no hay checkpoint en $(dir)")
    isfile(jls) || _alterado(stage, "falta el archivo del checkpoint")
    prov = try
        TOML.parsefile(tml)
    catch e
        _alterado(stage, "procedencia ilegible: $(sprint(showerror, e))")
    end
    get(prov, "stage", "") == String(stage) || _alterado(stage, "la procedencia es de otra etapa: $(get(prov, "stage", "?"))")
    comp = get(prov, "compatibility", Dict{String,Any}())
    get(prov, "format", -1) == CHECKPOINT_FORMAT && get(comp, "checkpoint_format", -1) == CHECKPOINT_FORMAT ||
        _incompatible(stage, "formato $(get(prov, "format", "?")); este motor lee el $(CHECKPOINT_FORMAT)")
    get(comp, "wrapper", -1) == CHECKPOINT_WRAPPER ||
        _incompatible(stage, "envoltorio $(get(comp, "wrapper", "?")); este motor lee el $(CHECKPOINT_WRAPPER)")
    get(comp, "julia", "") == string(VERSION) ||
        _incompatible(stage, "escrito con Julia $(get(comp, "julia", "?")); este proceso es $(VERSION)")
    for k in ("key", "file_sha256", "content_sha256")
        get(prov, k, "") isa String && !isempty(prov[k]) || _alterado(stage, "procedencia sin $(k)")
    end
    return prov
end

"""
Paso 5 sobre los checkpoints consumidos: dentro del almacén, con la
procedencia y la clave declaradas y con los bytes registrados. No los
deserializa. Devuelve etapa ⇒ (carpeta, procedencia).

`consumed_dirs` (etapa ⇒ carpeta `checkpoints/`) da la ubicación efectiva de
un consumido que se ha trasladado desde que se registró, como la carpeta
temporal de una tanda de Lab renombrada a su destino (E-1). Solo cambia
dónde se busca: la carpeta se valida dentro del almacén y el consumido tiene
que ser el declarado, por identidad, huella de contenido y huella de bytes.
"""
function _verify_consumed(stage::Symbol, prov, store::String; consumed_dirs=nothing)
    out = Dict{Symbol,Tuple{String,Dict{String,Any}}}()
    entries = get(prov, "consumed", Any[])
    declared = Set(Symbol(get(c, "stage", "")) for c in entries)
    declared == Set(CHECKPOINT_SPEC[stage].consumes) ||
        _alterado(stage, "consume $(sort(collect(declared))), y la etapa consume $(sort(CHECKPOINT_SPEC[stage].consumes))")
    given = consumed_dirs === nothing ? Dict{Symbol,String}() :
            Dict{Symbol,String}(Symbol(k) => String(v) for (k, v) in consumed_dirs)
    extra = setdiff(keys(given), declared)
    isempty(extra) ||
        throw(ArgumentError("consumed_dirs: $(stage) no consume $(join(sort(collect(extra)), ", "))"))
    for c in entries
        cst  = Symbol(c["stage"])
        jls  = if haskey(given, cst)
            _check_inside(stage, store, _jls(given[cst], cst); what="consumido $(cst) entregado")
        else
            p = _inside_store(stage, store, store, c["path"]; what="consumido $(cst)")
            (ispath(p) || ispath(_prov(dirname(p), cst))) ||
                _alterado(stage, "el consumido $(cst) no está en la ruta registrada ($(c["path"])). " *
                                 "Si se ha trasladado, su ubicación se da con `consumed_dirs`")
            p
        end
        cdir = dirname(jls)
        basename(jls) == "$(cst).jls" || _alterado(stage, "el consumido $(cst) no es un checkpoint de esa etapa")
        cprov = try
            _read_provenance(cdir, cst, store)
        catch e
            e isa CheckpointError || rethrow()
            _alterado(stage, "consumido $(cst): $(e.reason) · $(e.message)")
        end
        cprov["key"] == c["key"] || _alterado(stage, "el consumido $(cst) no es el declarado: otra identidad")
        cprov["content_sha256"] == c["content_sha256"] || _alterado(stage, "el consumido $(cst) no es el declarado: otro contenido")
        cprov["file_sha256"] == c["file_sha256"] || _alterado(stage, "el consumido $(cst) no es el declarado: otra huella")
        _file_sha256(jls) == c["file_sha256"] || _alterado(stage, "la huella de bytes del consumido $(cst) no coincide")
        out[cst] = (cdir, cprov)
    end
    return out
end

"Carga del envoltorio de un checkpoint cuya huella ya se ha comprobado."
function _payload(stage::Symbol, bytes::Vector{UInt8})
    w = try
        _deserialize_checked(bytes)
    catch e
        _alterado(stage, "el checkpoint no se puede deserializar: $(sprint(showerror, e))")
    end
    (w isa Dict && get(w, "wrapper", -1) == CHECKPOINT_WRAPPER && get(w, "stage", "") == String(stage) &&
     haskey(w, "data")) || _alterado(stage, "el envoltorio no es el de un checkpoint de $(stage)")
    return w["data"]
end

"Datos de un consumido, ya verificado por `_verify_consumed`, comprobando sus bytes otra vez."
_consumed_payload(stage::Symbol, cdir::String, cprov) =
    _payload(Symbol(cprov["stage"]), _read_checked(stage, _jls(cdir, Symbol(cprov["stage"])), cprov["file_sha256"];
                                                    what="consumido $(cprov["stage"])"))

"""
    load_checkpoint(dir, stage; store, expected_key=nothing, invariants=true, consumed_dirs=nothing) -> salida

Carga el checkpoint de `stage` de la carpeta `dir` (la carpeta `checkpoints/`)
con la verificación del apartado 3.5, en su orden:

1. rutas: todas dentro de `store`, sin rutas absolutas, sin `..` que salga y
   sin enlaces simbólicos;
2. procedencia (texto: leerla no la hace fiable);
3. compatibilidad: formato, envoltorio y versión de Julia → `INCOMPATIBLE`;
4. identidad frente a `expected_key`, si se da → `OTRA_CONFIGURACION`;
5. huellas de bytes del `.jls`, de los archivos referenciados y de los
   checkpoints consumidos → `ALTERADO`;
6. solo entonces deserializa, reconstruye la salida completa y comprueba su
   huella de contenido y, con `invariants=true`, las invariantes de la etapa.

`consumed_dirs` (etapa consumida ⇒ su carpeta `checkpoints/`) da la ubicación
efectiva de los consumidos que se han trasladado desde que se registraron
(E-1). No omite ninguna comprobación: cada consumido entregado se valida dentro
del almacén y tiene que ser el declarado en identidad y huellas, o el
resultado es `ALTERADO`. Una etapa que no se consume da `ArgumentError`.

Devuelve exactamente la salida que la etapa devolvió en memoria. Lanza
`CheckpointError` con el motivo.
"""
function load_checkpoint(dir::String, stage::Symbol; store::String, expected_key=nothing, invariants::Bool=true,
                         consumed_dirs=nothing)
    prov = _read_provenance(dir, stage, store)
    expected_key === nothing || prov["key"] == expected_key ||
        throw(CheckpointError(:OTRA_CONFIGURACION, stage, "la identidad guardada no es la pedida"))
    stage_root = dirname(normpath(abspath(dir)))
    refs = _verify_references(stage, prov, store, stage_root)
    consumed = _verify_consumed(stage, prov, store; consumed_dirs)
    bytes = _read_checked(stage, _jls(dir, stage), prov["file_sha256"]; what="el checkpoint")
    data = _payload(stage, bytes)
    out = _reconstruct(stage, data, prov, refs, consumed)
    checkpoint_content_sha256(out) == prov["content_sha256"] ||
        _alterado(stage, "la huella canónica del contenido no coincide")
    invariants && _check_invariants(stage, out, data, prov, consumed, store, stage_root)
    return out
end

"""
    verify_checkpoint(dir, stage; store, expected_key=nothing, deserialize=false, consumed_dirs=nothing)
        -> (ok, reason, message)

La verificación de `load_checkpoint` para auditar una carpeta. Por defecto
recorre los cinco primeros pasos sin deserializar nada; con `deserialize=true`
hace también el sexto y las invariantes. `consumed_dirs`, como en
`load_checkpoint`. No lanza por un checkpoint rechazado: devuelve el motivo.
"""
function verify_checkpoint(dir::String, stage::Symbol; store::String, expected_key=nothing, deserialize::Bool=false,
                           consumed_dirs=nothing)
    try
        if deserialize
            load_checkpoint(dir, stage; store, expected_key, invariants=true, consumed_dirs)
        else
            prov = _read_provenance(dir, stage, store)
            expected_key === nothing || prov["key"] == expected_key ||
                throw(CheckpointError(:OTRA_CONFIGURACION, stage, "la identidad guardada no es la pedida"))
            _verify_references(stage, prov, store, dirname(normpath(abspath(dir))))
            _verify_consumed(stage, prov, store; consumed_dirs)
            _file_sha256(_jls(dir, stage)) == prov["file_sha256"] ||
                _alterado(stage, "huella de bytes que no coincide en el checkpoint ($(stage).jls)")
            # Los productos de la etapa, relativos a su carpeta (decisión 133, P-1)
            stage_root = dirname(normpath(abspath(dir)))
            for (rel, h) in get(get(prov, "products", Dict{String,Any}()), "files", Dict{String,Any}())
                p = _inside_store(stage, store, stage_root, rel; what="producto")
                isfile(p) || _alterado(stage, "falta el producto $(rel)")
                _file_sha256(p) == h || _alterado(stage, "la huella del producto $(rel) no coincide")
            end
        end
        return (true, nothing, "")
    catch e
        e isa CheckpointError || rethrow()
        return (false, e.reason, e.message)
    end
end

# ── Referencias: descomposición de ICA y distribuciones nulas ──

"""
Paso 5 sobre los archivos referenciados: rutas dentro del almacén y huella de
bytes. Devuelve los bytes comprobados, sin deserializarlos. Las dos rutas,
`ica_result.jls` y las distribuciones nulas, son relativas a la carpeta de la
etapa (decisión 133, P-1), la que contiene `checkpoints/`, y se validan dentro
del almacén. Las nulas, además, no pueden salir de esa carpeta (decisión 130).
"""
function _verify_references(stage::Symbol, prov, store::String, stage_root::String)
    out = Dict{String,Vector{UInt8}}()
    for r in get(prov, "references", Any[])
        kind = get(r, "kind", "")
        p = if kind == "ica_decomposition"
            _inside_store(stage, store, stage_root, r["path"]; what="referencia")
        elseif kind == "null_distribution"
            # Más estricto que «dentro del almacén»: la distribución nula no
            # puede salir de la carpeta de su etapa, ni por `..`, ni por una
            # ruta absoluta, ni por un enlace simbólico.
            _inside_store(stage, stage_root, stage_root, r["path"]; what="distribución nula",
                          where="de la carpeta de la etapa")
            _inside_store(stage, store, stage_root, r["path"]; what="referencia")
        else
            _alterado(stage, "referencia de tipo desconocido: $(kind)")
        end
        out[r["path"]] = _read_checked(stage, p, r["file_sha256"]; what="la referencia $(r["path"])")
    end
    return out
end

function _ref(prov, kind; band=nothing)
    for r in get(prov, "references", Any[])
        get(r, "kind", "") == kind && (band === nothing || get(r, "band", "") == band) && return r
    end
    return nothing
end

# ── Contenido por etapa ──────────────────────────────────────

# Lo que se guarda en el `.jls` (la salida, sin lo que no se duplica) …
_checkpoint_data(::Val{S}, out) where {S} = out

function _checkpoint_data(::Val{:ica}, out::ICAOutput; referenced::Bool)
    Dict{String,Any}(
        "rec" => out.rec, "status" => out.status, "rejection_source" => out.rejection_source,
        "classification_ok" => out.classification_ok,
        "rejected" => out.result === nothing ? Int[] : copy(out.result.rejected_components),
        "has_decomposition" => out.result !== nothing,
        # Sin caché de la que referenciarla, la descomposición va dentro: sigue
        # habiendo una sola copia.
        "decomposition" => (out.result === nothing || referenced) ? nothing : out.result,
    )
end

function _checkpoint_data(::Val{:connectivity}, out::ConnectivityOutput; epochs_from_segmentation::Bool)
    Dict{String,Any}("conn" => out.conn, "n_edges" => out.n_edges,
                     "epochs_from" => epochs_from_segmentation ? "segmentation" : "inline",
                     "epochs" => epochs_from_segmentation ? nothing : out.epochs)
end

function _checkpoint_data(::Val{:surrogates}, out::SurrogatesOutput; referenced_bands::Set{String})
    Dict{String,Any}(
        "failed_bands" => copy(out.failed_bands),
        "results" => [Dict{String,Any}(
            "band" => sr.band, "p_values" => sr.p_values, "sig_mask" => sr.sig_mask,
            "fdr_threshold" => sr.fdr_threshold, "n_surrogates" => sr.n_surrogates,
            "null_distribution" => sr.band in referenced_bands ? nothing : sr.null_distribution,
        ) for sr in out.results],
    )
end

_checkpoint_data(::Val{:finalize}, ::Nothing) = Dict{String,Any}()

# … y cómo se reconstruye la salida completa.
function _reconstruct(stage::Symbol, data, prov, refs, consumed)
    T = Dict(:load => LoadOutput, :qc => QCOutput, :filter => FilterOutput, :montage => MontageOutput,
             :segmentation => SegmentationOutput, :spectral => SpectralResult)
    if haskey(T, stage)
        data isa T[stage] || _alterado(stage, "el checkpoint no contiene un $(T[stage])")
        return data
    elseif stage == :finalize
        data == Dict{String,Any}() || _alterado(stage, "el cierre no tiene contenido propio")
        return nothing
    elseif stage == :ica
        return _reconstruct_ica(data, prov, refs)
    elseif stage == :connectivity
        return _reconstruct_connectivity(data, consumed)
    elseif stage == :surrogates
        return _reconstruct_surrogates(data, prov, refs, consumed)
    end
    error("Etapa sin checkpoint: $(stage)")
end

function _reconstruct_ica(data, prov, refs)
    data isa Dict || _alterado(:ica, "contenido no reconocido")
    rec = data["rec"]; rejected = Vector{Int}(data["rejected"])
    result = if !data["has_decomposition"]
        nothing
    else
        base = if data["decomposition"] !== nothing
            data["decomposition"]
        else
            r = _ref(prov, "ica_decomposition")
            r === nothing && _alterado(:ica, "falta la referencia a la descomposición")
            obj = try
                _deserialize_checked(refs[r["path"]])
            catch e
                _alterado(:ica, "la descomposición referenciada no se puede deserializar: $(sprint(showerror, e))")
            end
            obj isa ICAResult || _alterado(:ica, "la referencia no contiene una descomposición")
            obj
        end
        n_ch, n_samp = size(rec.data)
        n_comp = size(base.activations, 1)
        (size(base.mixing_matrix) == (n_ch, n_comp) && size(base.activations, 2) == n_samp &&
         base.meta.channel_names == rec.meta.channel_names) ||
            _alterado(:ica, "la descomposición no corresponde a la señal del checkpoint")
        ICAResult(base.meta, base.mixing_matrix, base.unmixing_matrix, base.activations,
                  rejected, base.variance_explained, base.diagnostics)
    end
    return ICAOutput(result, rec, data["status"], data["rejection_source"], data["classification_ok"])
end

function _reconstruct_connectivity(data, consumed)
    data isa Dict || _alterado(:connectivity, "contenido no reconocido")
    epochs = if data["epochs_from"] == "segmentation"
        cdir, cprov = consumed[:segmentation]
        seg = _consumed_payload(:connectivity, cdir, cprov)
        seg isa SegmentationOutput || _alterado(:connectivity, "la segmentación consumida no contiene épocas")
        seg.epochs
    else
        data["epochs"]
    end
    return ConnectivityOutput(data["conn"], epochs, data["n_edges"])
end

function _reconstruct_surrogates(data, prov, refs, consumed)
    data isa Dict || _alterado(:surrogates, "contenido no reconocido")
    results = SurrogateResult[]
    conn = nothing
    if !isempty(data["results"])
        cdir, cprov = consumed[:connectivity]
        cdata = _consumed_payload(:surrogates, cdir, cprov)
        cdata isa Dict && haskey(cdata, "conn") || _alterado(:surrogates, "la conectividad consumida no es reconocible")
        conn = cdata["conn"]
    end
    for r in data["results"]
        band = r["band"]
        haskey(conn.matrices, band) || _alterado(:surrogates, "la banda $(band) no está en la conectividad consumida")
        null = if r["null_distribution"] !== nothing
            r["null_distribution"]
        else
            ref = _ref(prov, "null_distribution"; band)
            ref === nothing && _alterado(:surrogates, "falta la referencia a la distribución nula de $(band)")
            payload = try
                _deserialize_checked(refs[ref["path"]])
            catch e
                _alterado(:surrogates, "la distribución nula de $(band) no se puede deserializar: $(sprint(showerror, e))")
            end
            (payload isa Dict && get(payload, "band", "") == band &&
             get(payload, "channels", nothing) == conn.channel_names) ||
                _alterado(:surrogates, "la distribución nula de $(band) no corresponde a la conectividad")
            payload["null_distribution"]
        end
        n = length(conn.channel_names)
        (null isa Array{Float64,3} && size(null) == (n, n, r["n_surrogates"])) ||
            _alterado(:surrogates, "la distribución nula de $(band) no tiene canales × canales × N")
        push!(results, SurrogateResult(conn, band, conn.matrices[band], null, r["p_values"], r["sig_mask"],
                                       r["fdr_threshold"], r["n_surrogates"]))
    end
    return SurrogatesOutput(results, Vector{String}(data["failed_bands"]))
end

# ── Invariantes ──────────────────────────────────────────────

function _check_invariants(stage::Symbol, out, data, prov, consumed, store, stage_root)
    if stage == :ica
        cdir, cprov = consumed[:filter]
        filt = _consumed_payload(:ica, cdir, cprov)
        filt isa FilterOutput || _alterado(:ica, "el filtrado consumido no es reconocible")
        expected = if out.result === nothing || isempty(out.result.rejected_components)
            filt.rec.data
        else
            apply_ica_rejection(filt.rec, out.result, out.result.rejected_components).data
        end
        isequal(expected, out.rec.data) ||
            _alterado(:ica, "retirar los rechazados de la descomposición no da la señal guardada")
        if out.result !== nothing
            facts = get(prov, "facts", Dict{String,Any}())
            r = _ref(prov, "ica_decomposition")
            r === nothing || get(r, "identity", "") == get(facts, "decomposition_identity", "") ||
                _alterado(:ica, "la identidad de la descomposición referenciada no es la registrada")
        end
    elseif stage == :montage
        cdir, cprov = consumed[:ica]
        ica = _consumed_payload(:montage, cdir, cprov)
        ica isa Dict && haskey(ica, "rec") || _alterado(:montage, "la ICA consumida no es reconocible")
        src = ica["rec"]
        idx = [findfirst(==(c), src.meta.channel_names) for c in out.rec.meta.channel_names]
        (all(!isnothing, idx) && isequal(src.data[Int.(idx), :], out.rec.data)) ||
            _alterado(:montage, "la señal del montaje no son filas de la señal de ICA consumida")
    elseif stage == :segmentation
        out.n_valid + out.n_rejected == out.n_total ||
            _alterado(:segmentation, "los recuentos de épocas no cuadran")
        size(out.epochs.data, 3) == out.n_valid ||
            _alterado(:segmentation, "el número de épocas no es el de válidas")
    elseif stage == :finalize
        facts = get(prov, "facts", Dict{String,Any}())
        files = get(get(prov, "products", Dict{String,Any}()), "files", Dict{String,Any}())
        # `config_snapshot.toml` está en la carpeta de la etapa, y sus productos
        # son relativos a ella (decisión 133, P-1)
        for k in (k for k in keys(files) if k == "config_snapshot.toml")
            files[k] == get(facts, "toml_sha256", "") ||
                _alterado(:finalize, "la huella de config_snapshot.toml no es la del TOML de la identidad")
            p = _inside_store(:finalize, store, stage_root, k; what="producto")
            isfile(p) && _file_sha256(p) == files[k] ||
                _alterado(:finalize, "config_snapshot.toml no es el registrado")
        end
    end
    return nothing
end

# ── Escritura segura ─────────────────────────────────────────

struct _SimulatedKill <: Exception
    step::Symbol
end

"""
Pruebas: simula que el proceso muere en un punto de la escritura
(`:after_jls`, `:after_provenance_tmp`, `:after_jls_rename`). Sin uso fuera de
las pruebas.
"""
const _CHECKPOINT_FAULT = Ref{Symbol}(:none)
_fault(step::Symbol) = _CHECKPOINT_FAULT[] === step && throw(_SimulatedKill(step))

function _fsync(io::IOStream)
    flush(io)
    ccall(:fsync, Cint, (RawFD,), fd(io)) == 0 || error("fsync falló: $(Libc.strerror())")
end

"""
    save_checkpoint(ctx, stage, output; key, store, consumed=Dict(), references=Dict(),
                    facts=Dict(), products=Dict(), absent=String[], commit=nothing) -> carpeta

Escribe el checkpoint de `stage` en `<export_dir>/checkpoints/` (plan de 2f,
apartado 3.2):

1. en una carpeta temporal `checkpoints/.tmp-<etapa>-<aleatorio>/`, el `.jls`,
   cerrado y volcado a disco;
2. la huella de los bytes de ese archivo definitivo;
3. la procedencia, con esa huella;
4. se renombra el `.jls`;
5. se renombra la procedencia, la última: es la marca de checkpoint completo.

No sobrescribe nunca: si ya hay un checkpoint de la etapa, falla. Una escritura
que falla no deja ningún checkpoint cargable.

- `consumed`: etapa consumida ⇒ su carpeta de checkpoints.
- `references`: para ICA, `:ica_decomposition => (path, identity)`; para
  surrogates, `banda => path` de cada distribución nula. Rutas absolutas;
  se registran relativas a la carpeta de la etapa (`export_dir`, la que
  contiene `checkpoints/`; decisión 133, P-1).
- `products`: archivos que escribió la etapa (ruta relativa a la carpeta de la
  etapa ⇒ huella; puede empezar por `..` si el archivo está fuera de ella y
  dentro del almacén, como las tablas globales del cierre); `absent`: los que
  debía escribir y no están.

Los consumidos se registran relativos a la raíz del almacén.
"""
function save_checkpoint(ctx::PipelineContext, stage::Symbol, output; key::String, store::String,
                         consumed::AbstractDict = Dict{Symbol,String}(),
                         references::AbstractDict = Dict{Any,Any}(),
                         facts::AbstractDict = Dict{String,Any}(),
                         products::AbstractDict = Dict{String,String}(),
                         absent::AbstractVector = String[],
                         epochs_from_segmentation::Union{Nothing,Bool} = nothing,
                         commit = nothing)
    spec = CHECKPOINT_SPEC[stage]
    Set(Symbol.(keys(consumed))) == Set(spec.consumes) ||
        error("save_checkpoint($(stage)): consume $(spec.consumes), recibido $(collect(keys(consumed)))")
    dir = joinpath(ctx.export_dir, "checkpoints")
    stage_root = ctx.export_dir
    _check_inside(stage, store, dir; what="carpeta del checkpoint")
    mkpath(dir)
    (ispath(_jls(dir, stage)) || ispath(_prov(dir, stage))) &&
        error("Ya existe un checkpoint de $(stage) en $(dir): no se sobrescribe")

    # Contenido del .jls, sin lo que no se duplica
    refs_toml = Dict{String,Any}[]
    data = if stage == :ica
        r = get(references, :ica_decomposition, nothing)
        if r !== nothing
            p = _check_inside(stage, store, r.path; what="referencia")
            push!(refs_toml, Dict{String,Any}("kind" => "ica_decomposition", "path" => _rel(stage_root, p),
                                              "identity" => r.identity, "file_sha256" => _file_sha256(p)))
        end
        _checkpoint_data(Val(:ica), output; referenced = r !== nothing)
    elseif stage == :connectivity
        from_seg = epochs_from_segmentation === nothing ?
                   _same_epochs(output.epochs, consumed[:segmentation], store) : epochs_from_segmentation
        _checkpoint_data(Val(:connectivity), output; epochs_from_segmentation = from_seg)
    elseif stage == :surrogates
        bands = Set{String}()
        for sr in output.results
            p = get(references, sr.band, nothing)
            p === nothing && continue
            _check_inside(stage, stage_root, p; what="distribución nula", where="de la carpeta de la etapa")
            p = _check_inside(stage, store, p; what="referencia")
            push!(refs_toml, Dict{String,Any}("kind" => "null_distribution", "band" => sr.band,
                                              "path" => _rel(stage_root, p),
                                              "dims" => collect(size(sr.null_distribution)),
                                              "file_sha256" => _file_sha256(p)))
            push!(bands, sr.band)
        end
        _checkpoint_data(Val(:surrogates), output; referenced_bands = bands)
    else
        _checkpoint_data(Val(stage), output)
    end

    consumed_toml = Dict{String,Any}[]
    for c in sort!(collect(spec.consumes))
        cprov = TOML.parsefile(_prov(consumed[c], c))
        push!(consumed_toml, Dict{String,Any}(
            "stage" => String(c), "path" => _rel(store, _jls(consumed[c], c)), "key" => cprov["key"],
            "file_sha256" => cprov["file_sha256"], "content_sha256" => cprov["content_sha256"]))
    end

    impl = checkpoint_implementation()
    wrapper = Dict{String,Any}("wrapper" => CHECKPOINT_WRAPPER, "format" => CHECKPOINT_FORMAT,
                               "stage" => String(stage), "subject" => ctx.subject_id,
                               "session" => ctx.session_id, "task" => ctx.task, "data" => data)

    tmp = mktempdir(dir; prefix=".tmp-$(stage)-", cleanup=false)
    try
        jls_tmp = joinpath(tmp, "$(stage).jls")
        open(jls_tmp, "w") do io
            Serialization.serialize(io, wrapper)
            _fsync(io)
        end
        _fault(:after_jls)
        file_sha = _file_sha256(jls_tmp)
        prov = Dict{String,Any}(
            "format" => CHECKPOINT_FORMAT, "stage" => String(stage), "key" => key,
            "file_sha256" => file_sha, "content_sha256" => checkpoint_content_sha256(output),
            "compatibility" => Dict{String,Any}("checkpoint_format" => CHECKPOINT_FORMAT,
                                                "wrapper" => CHECKPOINT_WRAPPER, "julia" => string(VERSION)),
            "implementation" => Dict{String,Any}("src_sha256" => impl.src_sha256,
                                                 "manifest_sha256" => impl.manifest_sha256,
                                                 "neuromind_version" => impl.neuromind_version),
            "environment" => Dict{String,Any}("machine" => String(Sys.MACHINE),
                                              (commit === nothing ? () : ("commit" => String(commit),))...),
            "config" => Dict{String,Any}(s => _tomlable(checkpoint_section(ctx.cfg, s)) for s in spec.sections),
            "consumed" => consumed_toml,
            "references" => refs_toml,
            "facts" => _tomlable(facts),
            "products" => Dict{String,Any}("files" => Dict{String,Any}(products),
                                           "absent" => sort!(String.(collect(absent)))),
        )
        prov_tmp = joinpath(tmp, "$(stage).toml")
        open(prov_tmp, "w") do io
            TOML.print(io, prov; sorted=true)
            _fsync(io)
        end
        _fault(:after_provenance_tmp)
        mv(jls_tmp, _jls(dir, stage))
        _fault(:after_jls_rename)
        mv(prov_tmp, _prov(dir, stage))
        rm(tmp; force=true, recursive=true)
    catch e
        e isa _SimulatedKill || rm(tmp; force=true, recursive=true)
        rethrow()
    end
    return dir
end

"¿Son las épocas de la conectividad las de la segmentación consumida? Por contenido."
function _same_epochs(epochs::EpochSet, sdir::String, store::String)
    cprov = TOML.parsefile(_prov(sdir, :segmentation))
    seg = _payload(:segmentation, _read_checked(:connectivity, _jls(sdir, :segmentation), cprov["file_sha256"];
                                                 what="segmentación consumida"))
    seg isa SegmentationOutput || return false
    return checkpoint_content_sha256(seg.epochs) == checkpoint_content_sha256(epochs)
end

# ── Hechos por etapa (apartado 3.7) ──────────────────────────
#
# Lo que la etapa ya calcula y lo que se le aplicó, sin reglas de validez: la
# validez científica, su agregación y su presentación son de Lab.

function checkpoint_facts(stage::Symbol, out, cfg::PipelineConfig; extra...)
    f = Dict{String,Any}()
    if stage == :load
        f["n_channels"] = out.rec.meta.n_channels; f["n_samples"] = n_samples(out.rec)
        f["fs"] = out.rec.meta.fs; f["channels"] = out.rec.meta.channel_names
        f["electrodes_ok"] = Bool(get(out.electrodes, :ok, false))
    elseif stage == :qc
        f["bad_channels"] = out.bad_channels; f["amplitude_warning"] = out.amplitude_warning
        f["raw_sigma_mean"] = out.raw_sigma_mean
    elseif stage == :filter
        f["preview_written"] = out.preview_written
    elseif stage == :ica
        f["status"] = String(out.status); f["rejection_source"] = out.rejection_source
        f["classification_ok"] = out.classification_ok
        f["rejected"] = out.result === nothing ? Int[] : out.result.rejected_components
        if out.result !== nothing
            d = out.result.diagnostics
            for k in ("converged", "n_iter", "max_iter", "final_error", "n_attempts", "max_attempts",
                      "restart_needed", "seed", "reference_valid")
                haskey(d, k) && (f[k] = d[k])
            end
        end
    elseif stage == :montage
        f["excluded"] = out.excluded; f["fp2_excluded"] = out.fp2_excluded
        f["bad_channels_non_fp2"] = out.bad_channels_non_fp2; f["n_channels"] = out.n_channels
        f["channels"] = out.rec.meta.channel_names
    elseif stage == :segmentation
        f["n_total"] = out.n_total; f["n_valid"] = out.n_valid; f["n_rejected"] = out.n_rejected
        f["epoch_length_s"] = out.epochs.epoch_length_s
        f["min_epochs"] = Int(get(cfg.segmentation, "min_epochs", 10))
    elseif stage == :spectral
        f["n_epochs_used"] = out.n_epochs_used
        f["bands"] = sort!(collect(keys(out.band_power)))
        f["n_nonfinite_band_power"] = sum(count(!isfinite, v) for v in values(out.band_power); init=0)
    elseif stage == :connectivity
        c = out.conn
        f["method"] = c.method; f["space"] = c.space; f["channels"] = c.channel_names
        f["n_edges"] = out.n_edges; f["n_epochs_used"] = c.n_epochs_used
        f["epoch_length_s"] = out.epochs.epoch_length_s
        f["bands_computed"] = sort!(collect(keys(c.matrices)))
        f["bands_omitted"] = sort!([b for b in keys(cfg.bands) if !haskey(c.matrices, b)])
        f["wpli_method"] = String(get(cfg.connectivity, "wpli_method", "hilbert"))
        f["use_dwpli"] = Bool(get(cfg.connectivity, "use_dwpli", false))
    elseif stage == :surrogates
        s = cfg.surrogates
        f["enabled"] = Bool(get(s, "enabled", false))
        f["bands_computed"] = [sr.band for sr in out.results]; f["failed_bands"] = out.failed_bands
        f["n_surrogates"] = Int(get(s, "n_surrogates", 0)); f["seed"] = Int(get(s, "seed", 42))
        f["alpha"] = Float64(get(s, "alpha", 0.05)); f["fdr_method"] = String(get(s, "fdr_method", "bh"))
        # Lo que hace el motor, como hecho: no es una elección metodológica validada
        # (decisiones 120 y 128). El `method` del TOML es solo un metadato.
        f["method_applied"] = "circular_shift"
        f["correction_family"] = "por banda, triángulo superior"
        f["null_dims"] = Dict(sr.band => collect(size(sr.null_distribution)) for sr in out.results)
    end
    for (k, v) in extra
        f[String(k)] = v
    end
    return f
end

# ── Orquestador: checkpoints opcionales ──────────────────────

"""
Estado de una ejecución con checkpoints, del orquestador o de `run_stage`:
carpeta, almacén, claves de las etapas ya guardadas o recibidas, carpetas de
lo consumido (por defecto, la propia), inventario de archivos para saber qué
escribió cada etapa, commit del motor si quien llama lo conoce y notas que
pasan a los hechos.
"""
mutable struct _CheckpointRun
    ctx::PipelineContext
    store::String
    dir::String
    keys::Dict{Symbol,String}
    toml_sha256::String
    inventory::Dict{String,Tuple{Int,Float64}}
    hashes::Dict{String,String}
    consumed_dirs::Dict{Symbol,String}
    commit::Union{Nothing,String}
    notes::Dict{String,Any}
end

# El registro es de toda la ejecución, no de una etapa; la decisión manual de
# ICA es una entrada, no un producto.
const _EXCLUDED_FROM_PRODUCTS = ("pipeline_log.txt", "ica_labels.csv", "ICA_labels.csv")

function _inventory(run::_CheckpointRun)
    inv = Dict{String,Tuple{Int,Float64}}()
    exp = run.ctx.export_dir
    for (d, _, fs) in walkdir(exp), f in fs
        p = joinpath(d, f)
        r = relpath(p, exp)
        (startswith(r, "checkpoints") || r in _EXCLUDED_FROM_PRODUCTS) && continue
        st = stat(p); inv[p] = (st.size, st.mtime)
    end
    rd = results_dir(run.ctx.cfg)
    for p in (joinpath(rd, "subjects_index.csv"), joinpath(rd, "qc", "qc_decision_table.csv"))
        isfile(p) && (st = stat(p); inv[p] = (st.size, st.mtime))
    end
    return inv
end

"Archivos nuevos o reescritos desde el inventario anterior, con su huella."
function _products!(run::_CheckpointRun)
    now_inv = _inventory(run)
    prods = Dict{String,String}()
    for (p, sm) in now_inv
        get(run.inventory, p, nothing) == sm && continue
        prods[_rel(run.ctx.export_dir, p)] = _file_sha256(p)
    end
    run.inventory = now_inv
    return prods
end

"Checkpoints que ya hay en `dir` (`.jls` y `.toml`)."
_existing_checkpoints(dir::String) =
    isdir(dir) ? [f for f in readdir(dir) if endswith(f, ".jls") || endswith(f, ".toml")] : String[]

"""
P-5: el orquestador con checkpoints no usa una carpeta que ya los tiene. Se
comprueba antes de abrir el registro, sin escribir nada.
"""
function _refuse_existing_checkpoints(config_path::String)
    cfg = load_ss_config(config_path)
    dir = joinpath(_context_ids(cfg, config_path).export_dir, "checkpoints")
    prev = _existing_checkpoints(dir)
    isempty(prev) || error("La carpeta ya tiene checkpoints ($(join(prev, ", "))): no se sobrescriben. " *
                           "Use otra carpeta de resultados.")
    return nothing
end

function _checkpoint_run(ctx::PipelineContext, store)
    st = store === nothing ? results_dir(ctx.cfg) : String(store)
    dir = joinpath(ctx.export_dir, "checkpoints")
    prev = _existing_checkpoints(dir)
    isempty(prev) || error("La carpeta ya tiene checkpoints ($(join(prev, ", "))): no se sobrescriben. " *
                           "Use otra carpeta de resultados.")
    run = _CheckpointRun(ctx, st, dir, Dict{Symbol,String}(), _file_sha256(ctx.config_path),
                         Dict{String,Tuple{Int,Float64}}(), Dict{String,String}(),
                         Dict{Symbol,String}(), nothing, Dict{String,Any}())
    run.inventory = _inventory(run)
    return run
end

"""
Guarda el checkpoint de la etapa que acaba de terminar el orquestador. No
cambia nada de lo que la etapa escribió: solo lee la carpeta y escribe en
`checkpoints/`.
"""
function _checkpoint_after!(run::_CheckpointRun, stage::Symbol, out; consumed_outputs...)
    ctx = run.ctx; cfg = ctx.cfg
    spec = CHECKPOINT_SPEC[stage]
    consumed_keys = Dict(c => run.keys[c] for c in spec.consumes)
    consumed_dirs = Dict(c => get(run.consumed_dirs, c, run.dir) for c in spec.consumes)
    products = _products!(run)
    absent = String[]
    references = Dict{Any,Any}()
    extra = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in run.notes)
    extra_kw = NamedTuple()

    # La identidad, con la misma función que `stage_key` y `run_stage` (E-2).
    # La decisión manual de ICA es el archivo que la etapa ha leído.
    ids = (subject_id = ctx.subject_id, session_id = ctx.session_id, task = ctx.task, run = ctx.run)
    labels = stage == :ica ?
        (i = findfirst(isfile, _ica_manual_label_candidates(cfg, ctx.subject_id, ctx.session_id, ctx.condition));
         i === nothing ? nothing : _ica_manual_label_candidates(cfg, ctx.subject_id, ctx.session_id, ctx.condition)[i]) :
        nothing
    id = _stage_identity(stage, cfg, ids, ctx.config_path, consumed_keys;
                         filter = get(consumed_outputs, :filter, nothing), qc = get(consumed_outputs, :qc, nothing),
                         labels_file = labels, toml_sha256 = run.toml_sha256)

    if stage == :ica
        ident = id.decomposition
        extra[:decomposition_identity] = ident === nothing ? "" : ident.key
        extra[:manual_decision] = id.manual === nothing ? "none" :
            (out.rejection_source == "manual" ? "applied" : "not_applicable")
        cache = joinpath(ctx.export_dir, "cache", "ica_result.jls")
        hashf = joinpath(ctx.export_dir, "cache", "ica_config.hash")
        if out.result !== nothing && ident !== nothing && isfile(cache) && isfile(hashf) &&
           _ica_stored_key(hashf) == ident.key
            references[:ica_decomposition] = (path = cache, identity = ident.key)
        end
    elseif stage == :segmentation && haskey(consumed_outputs, :ica)
        # P-9: el motor rotula la entrada a partir de la lista de rechazados de
        # la ICA consumida; el rótulo escrito tiene que ser el que corresponde
        ica = consumed_outputs[:ica]
        want = (ica.result !== nothing && !isempty(ica.result.rejected_components)) ? "ICA-limpiada" : "filtrada"
        summary = joinpath(ctx.export_dir, "json", "segmentation_summary.json")
        m = isfile(summary) ? match(r"\"signal_input\":\s*\"([^\"]*)\"", read(summary, String)) : nothing
        got = m === nothing ? "" : String(m.captures[1])
        got == want || throw(StageError(:segmentation, :SEGMENTACION_ROTULO,
            "el motor rotula la entrada como «$(got)» y a la lista de rechazados le corresponde «$(want)»"))
        extra[:signal_input] = got
    elseif stage == :connectivity
        # Sin CSD, las épocas de la conectividad son las de la segmentación: se
        # referencian a ese checkpoint en lugar de guardarlas otra vez.
        seg = consumed_outputs[:segmentation]
        extra_kw = (epochs_from_segmentation = out.epochs === seg.epochs ||
                    checkpoint_content_sha256(out.epochs) == checkpoint_content_sha256(seg.epochs),)
    elseif stage == :surrogates
        # Solo se referencia una distribución nula que esta ejecución acaba de
        # escribir: una de una ejecución anterior no describe este cálculo.
        for sr in out.results
            p = joinpath(ctx.export_dir, "cache", "surrogate", "null_distribution_$(sr.band).jls")
            isfile(p) && haskey(products, _rel(ctx.export_dir, p)) && (references[sr.band] = p)
        end
    elseif stage == :finalize
        snap = joinpath(ctx.export_dir, "config_snapshot.toml")
        rel = _rel(ctx.export_dir, snap)
        if !(haskey(products, rel) && products[rel] == run.toml_sha256)
            delete!(products, rel); push!(absent, rel)
        end
        extra[:toml_sha256] = run.toml_sha256
        extra[:config_snapshot] = isempty(absent) ? "presente" : "ausente"
    end
    key = id.key
    facts = checkpoint_facts(stage, out, cfg; extra...)
    save_checkpoint(ctx, stage, out; key, store=run.store, consumed=consumed_dirs, references,
                    facts, products, absent, commit=run.commit, extra_kw...)
    run.keys[stage] = key
    run.inventory = _inventory(run)
    return key
end
