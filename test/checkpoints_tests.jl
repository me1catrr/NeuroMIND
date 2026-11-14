# ───────────────────────────────────────────────────────────────
#  Fichero    test/checkpoints_tests.jl
#  Modificado 14-11-2026
# ───────────────────────────────────────────────────────────────
#
# Conjuntos de la entrega 2f-1 (plan de 2f, apartado 7.2): checkpoints comunes
# de las etapas. Datos sintéticos: la grabación de `synthetic_bids`. Lo incluye
# runtests.jl después de los conjuntos existentes, que no cambian.
#
# `run_stage`, la caché de ICA entregada y la decisión manual verificada contra
# su descomposición son de 2f-2 y no se prueban aquí.

using SHA: sha256

const CK = NeuroMIND

# ── Ayudantes ────────────────────────────────────────────────

"Igualdad campo a campo, con los NaN en su sitio y los tipos exactos."
deepeq(a, b) = _deq(a, b)
_deq(a::AbstractArray, b::AbstractArray) =
    typeof(a) == typeof(b) && size(a) == size(b) &&
    (eltype(a) <: Union{Number,Bool,Missing} ? isequal(a, b) : all(_deq(x, y) for (x, y) in zip(a, b)))
_deq(a::AbstractDict, b::AbstractDict) =
    typeof(a) == typeof(b) && Set(keys(a)) == Set(keys(b)) && all(_deq(a[k], b[k]) for k in keys(a))
_deq(a::DataFrame, b::DataFrame) =
    names(a) == names(b) && all(typeof(a[!, c]) == typeof(b[!, c]) && isequal(a[!, c], b[!, c]) for c in names(a))
_deq(a::NamedTuple, b::NamedTuple) = keys(a) == keys(b) && all(_deq(a[k], b[k]) for k in keys(a))
_deq(a::Tuple, b::Tuple) = length(a) == length(b) && all(_deq(x, y) for (x, y) in zip(a, b))
function _deq(a, b)
    typeof(a) == typeof(b) || return false
    (a isa Number || a isa AbstractString || a isa Symbol || a === nothing) && return isequal(a, b)
    isstructtype(typeof(a)) || return isequal(a, b)
    return all(_deq(getfield(a, f), getfield(b, f)) for f in fieldnames(typeof(a)))
end

"Carpeta de checkpoints de la grabación sintética bajo `root/results`."
ckdir(root) = joinpath(root, "results", "subjects", "sub-S01", "ses-T1", "eyesclosed", "checkpoints")
ckstore(root) = joinpath(root, "results")

"""
Ejecuta las etapas a mano, como el orquestador, guardando el checkpoint de cada
una, y devuelve las salidas en memoria. `until` corta la cadena; `store` es el
almacén (por defecto, la carpeta de resultados).
"""
function chain_with_checkpoints(cfg; until::Symbol=:finalize, before_ica=nothing, before_finalize_ck=nothing,
                                store=nothing)
    redirect_stdout(devnull) do
        ctx = pipeline_context(cfg)
        ck  = CK._checkpoint_run(ctx, store)
        o = Dict{Symbol,Any}(:ctx => ctx)
        try
            o[:load] = stage_load(ctx);                 CK._checkpoint_after!(ck, :load, o[:load])
            if until != :ica
                o[:qc] = stage_qc(ctx, o[:load]);       CK._checkpoint_after!(ck, :qc, o[:qc])
            end
            o[:filter] = stage_filter(ctx, o[:load]);   CK._checkpoint_after!(ck, :filter, o[:filter])
            before_ica === nothing || before_ica(ctx)
            o[:ica] = stage_ica(ctx, o[:filter]);       CK._checkpoint_after!(ck, :ica, o[:ica]; filter=o[:filter])
            until == :ica && return o
            o[:montage] = stage_montage(ctx, o[:ica], o[:qc]); CK._checkpoint_after!(ck, :montage, o[:montage]; qc=o[:qc])
            o[:segmentation] = stage_segmentation(ctx, o[:montage], o[:ica])
            CK._checkpoint_after!(ck, :segmentation, o[:segmentation])
            o[:spectral] = stage_spectral(ctx, o[:segmentation]); CK._checkpoint_after!(ck, :spectral, o[:spectral])
            o[:connectivity] = stage_connectivity(ctx, o[:segmentation], o[:montage])
            CK._checkpoint_after!(ck, :connectivity, o[:connectivity]; segmentation=o[:segmentation])
            o[:surrogates] = stage_surrogates(ctx, o[:connectivity]); CK._checkpoint_after!(ck, :surrogates, o[:surrogates])
            stage_finalize(ctx, o[:load], o[:qc], o[:filter], o[:montage], o[:segmentation], o[:spectral], o[:connectivity])
            before_finalize_ck === nothing || before_finalize_ck(ctx)
            o[:finalize] = nothing;                     CK._checkpoint_after!(ck, :finalize, nothing)
        finally
            close(ctx.log_io)
        end
        o
    end
end

"Lanza `f` y devuelve el motivo del `CheckpointError`, o `:ninguno`."
function reason_of(f)
    try
        f(); return :ninguno
    catch e
        e isa CheckpointError || rethrow()
        return e.reason
    end
end

"Copia un resultado entero a otra raíz, para alterarlo sin tocar el original."
function copy_results(src_root, dst_root)
    mkpath(dst_root); cp(joinpath(src_root, "results"), joinpath(dst_root, "results")); dst_root
end

"Cambia un byte de un archivo (el del medio)."
function flip_byte!(path)
    b = read(path); i = cld(length(b), 2); b[i] = b[i] ⊻ 0xff
    chmod(path, 0o644); write(path, b)
end

"Reescribe la procedencia de un checkpoint aplicando `f!` al diccionario."
function edit_provenance!(f!, dir, stage)
    p = joinpath(dir, "$(stage).toml")
    d = CK.TOML.parsefile(p); f!(d)
    open(io -> CK.TOML.print(io, d; sorted=true), p, "w")
end

const SURR20 = (("surrogates", "enabled") => "true", ("surrogates", "n_surrogates") => "20")

# Ejecuciones compartidas por los conjuntos
const CK_TMP  = mktempdir()
const CK_BIDS = synthetic_bids(joinpath(CK_TMP, "bids"))
# A: orquestador con checkpoints; D: el mismo sin ellos; C: cadena a mano con checkpoints
const CK_A = joinpath(CK_TMP, "A"); const CK_D = joinpath(CK_TMP, "D"); const CK_C = joinpath(CK_TMP, "C")
redirect_stdout(devnull) do
    run_single_subject_pipeline(synthetic_config(CK_A, CK_BIDS, SURR20...); checkpoints=true)
    run_single_subject_pipeline(synthetic_config(CK_D, CK_BIDS, SURR20...))
end
const CK_MEM = chain_with_checkpoints(synthetic_config(CK_C, CK_BIDS, SURR20...))
# E: casos límite. ICA que falla, ninguna época válida (potencias no numéricas, D-23)
# y surrogates desactivados, que es lo del protocolo
const CK_E = joinpath(CK_TMP, "E")
const CK_EDGE = chain_with_checkpoints(synthetic_config(CK_E, CK_BIDS,
    ("ica", "profile") => "\"default\"", ("ica", "n_components") => "-3",
    ("artifact_rejection", "max_amplitude_uv") => "0.001"))

# ── Checkpoints_RoundTrip ────────────────────────────────────

@testset "Checkpoints_RoundTrip" begin
    @testset "Los diez tipos: lo leído es lo guardado" begin
        for s in CHECKPOINT_STAGES
            out = load_checkpoint(ckdir(CK_C), s; store=ckstore(CK_C))
            @test deepeq(out, CK_MEM[s])
            @test checkpoint_content_sha256(out) == CK.TOML.parsefile(joinpath(ckdir(CK_C), "$(s).toml"))["content_sha256"]
        end
    end

    @testset "El orquestador escribe lo mismo que hay en memoria" begin
        for s in CHECKPOINT_STAGES
            out = load_checkpoint(ckdir(CK_A), s; store=ckstore(CK_A))
            @test checkpoint_content_sha256(out) == checkpoint_content_sha256(CK_MEM[s])
        end
    end

    @testset "La opción no cambia nada de lo que escribe el pipeline" begin
        a = ckstore(CK_A); d = ckstore(CK_D)
        ta = [f for f in tree(a) if !occursin("/checkpoints/", f)]
        @test ta == tree(d)
        @test length(tree(a)) - length(ta) == 2 * length(CHECKPOINT_STAGES)
        # Tablas, figuras y distribuciones nulas, byte a byte (la caché de ICA
        # lleva el tiempo de optimización; las tablas globales, la fecha)
        compared = [f for f in ta if (endswith(f, ".csv") || endswith(f, ".png") || occursin("cache/surrogate/", f)) &&
                                     !occursin("subjects_index", f) && !occursin("qc_decision_table", f)]
        @test length(compared) > 40
        @test isempty([f for f in compared if read(joinpath(a, f)) != read(joinpath(d, f))])
    end

    @testset "Surrogates: la nula referenciada y la conectividad reconstruida" begin
        out = load_checkpoint(ckdir(CK_C), :surrogates; store=ckstore(CK_C))
        @test length(out.results) == 7
        for (sr, m) in zip(out.results, CK_MEM[:surrogates].results)
            @test deepeq(sr.observed, m.observed)
            @test deepeq(sr.connectivity, m.connectivity)
            @test deepeq(sr.null_distribution, m.null_distribution)
            @test sr.observed === sr.connectivity.matrices[sr.band]
        end
        prov = CK.TOML.parsefile(joinpath(ckdir(CK_C), "surrogates.toml"))
        refs = [r for r in prov["references"] if r["kind"] == "null_distribution"]
        @test length(refs) == 7
        @test all(startswith(r["path"], "cache/surrogate/null_distribution_") for r in refs)
        @test all(r["dims"] == [15, 15, 20] for r in refs)
        # Dentro del .jls no hay ninguna distribución nula
        w = deserialize(joinpath(ckdir(CK_C), "surrogates.jls"))
        @test all(r["null_distribution"] === nothing for r in w["data"]["results"])
    end

    @testset "Conectividad: las épocas son las de la segmentación, referenciadas" begin
        w = deserialize(joinpath(ckdir(CK_C), "connectivity.jls"))
        @test w["data"]["epochs_from"] == "segmentation"
        @test w["data"]["epochs"] === nothing
        out = load_checkpoint(ckdir(CK_C), :connectivity; store=ckstore(CK_C))
        @test deepeq(out.epochs, CK_MEM[:segmentation].epochs)
    end

    @testset "Casos límite: ICA fallida, cero épocas, potencias no numéricas, surrogates desactivados" begin
        @test CK_EDGE[:ica].status == :failed
        @test CK_EDGE[:segmentation].n_valid == 0
        @test any(v -> any(!isfinite, v), values(CK_EDGE[:spectral].band_power))
        @test isempty(CK_EDGE[:surrogates].results)
        for s in CHECKPOINT_STAGES
            out = load_checkpoint(ckdir(CK_E), s; store=ckstore(CK_E))
            @test deepeq(out, CK_EDGE[s])
        end
        prov = CK.TOML.parsefile(joinpath(ckdir(CK_E), "ica.toml"))
        @test isempty(prov["references"])
        @test prov["facts"]["status"] == "failed"
        @test CK.TOML.parsefile(joinpath(ckdir(CK_E), "surrogates.toml"))["facts"]["enabled"] == false
    end
end

# ── Checkpoints_Write ────────────────────────────────────────

@testset "Checkpoints_Write" begin
    function load_only(name)
        root = joinpath(CK_TMP, "w-$(name)")
        ctx = redirect_stdout(devnull) do
            pipeline_context(synthetic_config(root, CK_BIDS))
        end
        loaded = redirect_stdout(devnull) do
            stage_load(ctx)
        end
        close(ctx.log_io)
        key = checkpoint_key(:load, ctx; checkpoint_load_inputs(ctx.cfg, "S01", "T1", "eyesclosed", 1)...)
        return ctx, loaded, key, ckstore(root)
    end

    for step in (:after_jls, :after_provenance_tmp, :after_jls_rename)
        @testset "Interrupción $(step)" begin
            ctx, loaded, key, store = load_only(step)
            dir = joinpath(ctx.export_dir, "checkpoints")
            CK._CHECKPOINT_FAULT[] = step
            try
                @test_throws CK._SimulatedKill save_checkpoint(ctx, :load, loaded; key, store)
            finally
                CK._CHECKPOINT_FAULT[] = :none
            end
            # Nunca queda un checkpoint cargable
            @test !isfile(joinpath(dir, "load.toml"))
            @test reason_of(() -> load_checkpoint(dir, :load; store)) == :ALTERADO
            @test verify_checkpoint(dir, :load; store)[2] == :ALTERADO
            @test isfile(joinpath(dir, "load.jls")) == (step == :after_jls_rename)
            if step == :after_jls_rename
                # Un .jls sin procedencia no se carga, y no se sobrescribe
                @test_throws ErrorException save_checkpoint(ctx, :load, loaded; key, store)
            else
                # Los temporales que quedan se ignoran
                @test any(startswith(".tmp-load-"), readdir(dir))
                save_checkpoint(ctx, :load, loaded; key, store)
                @test deepeq(load_checkpoint(dir, :load; store), loaded)
            end
        end
    end

    @testset "La huella registrada es la del archivo publicado; no se sobrescribe" begin
        ctx, loaded, key, store = load_only("huella")
        dir = save_checkpoint(ctx, :load, loaded; key, store)
        prov = CK.TOML.parsefile(joinpath(dir, "load.toml"))
        before = read(joinpath(dir, "load.jls"))
        @test prov["file_sha256"] == bytes2hex(sha256(before))
        @test_throws ErrorException save_checkpoint(ctx, :load, loaded; key, store)
        @test read(joinpath(dir, "load.jls")) == before
        @test !any(startswith(".tmp-"), readdir(dir))
    end

    @testset "El orquestador no escribe sobre checkpoints anteriores" begin
        @test_throws ErrorException redirect_stdout(devnull) do
            run_single_subject_pipeline(joinpath(CK_A, "config", "pipeline.toml"); checkpoints=true)
        end
    end
end

# ── Checkpoints_Untrusted ────────────────────────────────────

@testset "Checkpoints_Untrusted" begin
    # Cada caso altera una copia del resultado A y espera ALTERADO sin
    # deserializar nada.
    function untrusted(name, stage, alter!)
        root = copy_results(CK_A, joinpath(CK_TMP, "u-$(name)"))
        dir = ckdir(root); store = ckstore(root)
        alter!(dir, store)
        CK.CHECKPOINT_DESERIALIZATIONS[] = 0
        m = try
            load_checkpoint(dir, stage; store); nothing
        catch e
            e isa CheckpointError || rethrow()
            e
        end
        n = CK.CHECKPOINT_DESERIALIZATIONS[]
        v = verify_checkpoint(dir, stage; store)
        return (m === nothing ? :ninguno : m.reason), (m === nothing ? "" : m.message), n, v
    end
    eyes(dir) = dirname(dir)

    cases = [
        ("byte del checkpoint", :segmentation, (d, s) -> flip_byte!(joinpath(d, "segmentation.jls"))),
        ("byte de ica_result.jls", :ica, (d, s) -> flip_byte!(joinpath(eyes(d), "cache", "ica_result.jls"))),
        ("referencia absoluta", :ica, (d, s) -> edit_provenance!(p -> (p["references"][1]["path"] =
            joinpath(eyes(d), "cache", "ica_result.jls")), d, :ica)),
        # La referencia es relativa a la carpeta de la etapa (decisión 133): la
        # ruta sale del almacén y apunta a una copia con los mismos bytes, para
        # que el rechazo sea por la ruta y no porque falte el archivo
        ("referencia que sale del almacén", :ica, function (d, s)
            out = joinpath(dirname(s), "fuera", "ica_result.jls"); mkpath(dirname(out))
            cp(joinpath(eyes(d), "cache", "ica_result.jls"), out)
            edit_provenance!(p -> (p["references"][1]["path"] = CK._rel(eyes(d), out)), d, :ica)
        end),
        ("consumido que sale del almacén", :montage, (d, s) -> edit_provenance!(p -> (p["consumed"][1]["path"] =
            "../ica.jls"), d, :montage)),
        ("enlace simbólico en el recorrido", :ica, function (d, s)
            c = joinpath(eyes(d), "cache"); mv(c, joinpath(eyes(d), "cache_real")); symlink(joinpath(eyes(d), "cache_real"), c)
        end),
        ("checkpoint que es un enlace simbólico", :qc, function (d, s)
            p = joinpath(d, "qc.jls"); mv(p, joinpath(eyes(d), "qc_real.jls")); symlink(joinpath(eyes(d), "qc_real.jls"), p)
        end),
        ("distribución nula que falta", :surrogates, (d, s) ->
            rm(joinpath(eyes(d), "cache", "surrogate", "null_distribution_ALPHA.jls"))),
        ("byte de una distribución nula", :surrogates, (d, s) ->
            flip_byte!(joinpath(eyes(d), "cache", "surrogate", "null_distribution_BETA_LOW.jls"))),
        ("byte de un consumido", :spectral, (d, s) -> flip_byte!(joinpath(d, "segmentation.jls"))),
        # La distribución nula no puede salir de la carpeta de su etapa, aunque
        # se quede dentro del almacén y sus bytes sean los registrados
        ("nula fuera de la carpeta de la etapa, con ../", :surrogates, function (d, s)
            src = joinpath(eyes(d), "cache", "surrogate", "null_distribution_ALPHA.jls")
            out = joinpath(dirname(eyes(d)), "otra-etapa", "null_distribution_ALPHA.jls")
            mkpath(dirname(out)); cp(src, out)
            edit_provenance!(p -> (r = only(x for x in p["references"] if x["band"] == "ALPHA");
                                   r["path"] = "../otra-etapa/null_distribution_ALPHA.jls"), d, :surrogates)
        end),
        ("nula con ruta absoluta", :surrogates, (d, s) -> edit_provenance!(p -> (
            r = only(x for x in p["references"] if x["band"] == "ALPHA");
            r["path"] = joinpath(eyes(d), "cache", "surrogate", "null_distribution_ALPHA.jls")), d, :surrogates)),
        ("nula a través de un enlace simbólico", :surrogates, function (d, s)
            c = joinpath(eyes(d), "cache", "surrogate"); mv(c, joinpath(eyes(d), "cache", "surrogate_real"))
            symlink(joinpath(eyes(d), "cache", "surrogate_real"), c)
        end),
    ]
    # Los rechazos por la ruta tienen que serlo por la ruta, no porque falte
    # el archivo
    motive = Dict(
        "referencia absoluta" => "absoluta",
        "referencia que sale del almacén" => "fuera del almacén autorizado",
        "consumido que sale del almacén" => "fuera del almacén autorizado",
        "enlace simbólico en el recorrido" => "enlace simbólico",
        "checkpoint que es un enlace simbólico" => "enlace simbólico",
        "nula fuera de la carpeta de la etapa, con ../" => "fuera de la carpeta de la etapa",
        "nula con ruta absoluta" => "absoluta",
        "nula a través de un enlace simbólico" => "enlace simbólico",
    )
    for (name, stage, alter!) in cases
        @testset "$(name)" begin
            r, msg, n, v = untrusted(replace(name, " " => "-"), stage, alter!)
            @test r == :ALTERADO
            @test n == 0
            @test v[1] == false && v[2] == :ALTERADO
            haskey(motive, name) && @test occursin(motive[name], msg)
        end
    end

    @testset "Un almacén distinto no da acceso a los archivos de otro" begin
        other = mktempdir()
        @test reason_of(() -> load_checkpoint(ckdir(CK_A), :load; store=other)) == :ALTERADO
    end
end

# ── Checkpoints_Load ─────────────────────────────────────────

@testset "Checkpoints_Load" begin
    function altered(name, stage, f!)
        root = copy_results(CK_A, joinpath(CK_TMP, "l-$(name)"))
        f!(ckdir(root))
        return reason_of(() -> load_checkpoint(ckdir(root), stage; store=ckstore(root)))
    end

    @testset "INCOMPATIBLE: otra versión de formato, de envoltorio o de Julia" begin
        @test altered("julia", :qc, d -> edit_provenance!(p -> (p["compatibility"]["julia"] = "1.10.0"), d, :qc)) == :INCOMPATIBLE
        @test altered("formato", :qc, d -> edit_provenance!(p -> (p["format"] = 2), d, :qc)) == :INCOMPATIBLE
        @test altered("envoltorio", :qc, d -> edit_provenance!(p -> (p["compatibility"]["wrapper"] = 9), d, :qc)) == :INCOMPATIBLE
    end

    @testset "OTRA_CONFIGURACION: la identidad no es la pedida" begin
        store = ckstore(CK_A); dir = ckdir(CK_A)
        key = CK.TOML.parsefile(joinpath(dir, "qc.toml"))["key"]
        @test deepeq(load_checkpoint(dir, :qc; store, expected_key=key), CK_MEM[:qc])
        other = load_ss_config(synthetic_config(joinpath(CK_TMP, "l-otra"), CK_BIDS, ("qc", "bad_channel_zscore_threshold") => "4.0"))
        load_key = CK.TOML.parsefile(joinpath(dir, "load.toml"))["key"]
        @test reason_of(() -> load_checkpoint(dir, :qc; store,
            expected_key=checkpoint_key(:qc, other, Dict(:load => load_key)))) == :OTRA_CONFIGURACION
        @test verify_checkpoint(dir, :qc; store, expected_key="otra")[2] == :OTRA_CONFIGURACION
    end

    @testset "ALTERADO: contenido, procedencia, consumido e invariante" begin
        @test altered("contenido", :spectral, d -> edit_provenance!(p -> (p["content_sha256"] = "0"^64), d, :spectral)) == :ALTERADO
        @test altered("sin-procedencia", :qc, d -> rm(joinpath(d, "qc.toml"))) == :ALTERADO
        @test altered("consumido-otro", :montage, d -> edit_provenance!(p -> (p["consumed"][1]["key"] = "0"^64), d, :montage)) == :ALTERADO
        @test altered("consumido-falta", :montage, d -> rm(joinpath(d, "ica.toml"))) == :ALTERADO
        @test altered("otra-etapa", :qc, d -> edit_provenance!(p -> (p["stage"] = "load"), d, :qc)) == :ALTERADO

        # Invariante de ICA: un checkpoint bien escrito cuya señal no es la que
        # resulta de retirar los rechazados. Las huellas coinciden; la
        # invariante, no.
        root = joinpath(CK_TMP, "l-invariante")
        o = chain_with_checkpoints(synthetic_config(root, CK_BIDS); until=:ica)
        dir = ckdir(root); store = ckstore(root)
        ica = o[:ica]
        bad = ICAOutput(ica.result, EEGRecording(ica.rec.meta, ica.rec.data .+ 1e-9, ica.rec.times),
                        ica.status, ica.rejection_source, ica.classification_ok)
        mv(joinpath(dir, "ica.toml"), joinpath(root, "ica.toml.orig")); mv(joinpath(dir, "ica.jls"), joinpath(root, "ica.jls.orig"))
        save_checkpoint(o[:ctx], :ica, bad; key="k", store, consumed=Dict(:filter => dir),
                        references=Dict(:ica_decomposition => (path=joinpath(dirname(dir), "cache", "ica_result.jls"),
                                                                identity=CK._ica_cache_identity(o[:filter].rec, o[:ctx].cfg).key)),
                        facts=Dict("decomposition_identity" => CK._ica_cache_identity(o[:filter].rec, o[:ctx].cfg).key))
        @test reason_of(() -> load_checkpoint(dir, :ica; store)) == :ALTERADO
        @test reason_of(() -> load_checkpoint(dir, :ica; store, invariants=false)) == :ninguno
    end

    @testset "Un checkpoint ALTERADO se conserva y nunca se selecciona; el cálculo nuevo es otro, en otro sitio" begin
        r1 = joinpath(CK_TMP, "l-alterado"); c1 = synthetic_config(r1, CK_BIDS)
        redirect_stdout(devnull) do; run_single_subject_pipeline(c1; checkpoints=true); end
        bad = joinpath(ckdir(r1), "segmentation.jls")
        flip_byte!(bad)
        evidence = read(bad); evidence_prov = read(joinpath(ckdir(r1), "segmentation.toml"))
        @test reason_of(() -> load_checkpoint(ckdir(r1), :segmentation; store=ckstore(r1))) == :ALTERADO
        # Lo que lo consume tampoco se puede cargar: no se reutiliza por la cadena
        @test reason_of(() -> load_checkpoint(ckdir(r1), :spectral; store=ckstore(r1))) == :ALTERADO
        # No se puede escribir encima: la ejecución nueva no usa esa carpeta
        @test_throws ErrorException redirect_stdout(devnull) do
            run_single_subject_pipeline(c1; checkpoints=true)
        end
        # El cálculo nuevo, en otra carpeta: el mismo cálculo (misma identidad),
        # otra ubicación y otro archivo, cargable
        r2 = joinpath(CK_TMP, "l-recalculo")
        redirect_stdout(devnull) do; run_single_subject_pipeline(synthetic_config(r2, CK_BIDS); checkpoints=true); end
        p1 = CK.TOML.parsefile(joinpath(ckdir(r1), "segmentation.toml"))
        p2 = CK.TOML.parsefile(joinpath(ckdir(r2), "segmentation.toml"))
        @test p1["key"] == p2["key"]
        @test ckdir(r1) != ckdir(r2)
        @test bytes2hex(sha256(read(joinpath(ckdir(r2), "segmentation.jls")))) == p2["file_sha256"]
        @test reason_of(() -> load_checkpoint(ckdir(r2), :segmentation; store=ckstore(r2), expected_key=p2["key"])) == :ninguno
        # El alterado sigue ahí, sin tocar, y sigue rechazándose con la misma identidad pedida
        @test read(bad) == evidence
        @test read(joinpath(ckdir(r1), "segmentation.toml")) == evidence_prov
        @test reason_of(() -> load_checkpoint(ckdir(r1), :segmentation; store=ckstore(r1), expected_key=p2["key"])) == :ALTERADO
    end
end

# ── Checkpoints_Dependencies ─────────────────────────────────

@testset "Checkpoints_Dependencies" begin
    # Quién consume a quién (apartado 3.4), escrito aquí aparte del código:
    # el montaje no consume la identidad del QC, solo su lista.
    DEPENDENTS = Dict(
        :load => [:qc, :filter, :finalize], :qc => [:finalize], :filter => [:ica, :finalize],
        :ica => [:montage, :segmentation], :montage => [:segmentation, :connectivity, :finalize],
        :segmentation => [:spectral, :connectivity, :finalize], :spectral => [:finalize],
        :connectivity => [:surrogates, :finalize], :surrogates => Symbol[], :finalize => Symbol[])
    closure(direct) = (seen = Set{Symbol}(direct); todo = collect(direct);
                       while !isempty(todo); for d in DEPENDENTS[pop!(todo)]; d in seen || (push!(seen, d); push!(todo, d)); end; end; seen)

    files = Dict("vhdr" => "a", "eeg" => "b", "vmrk" => "c", "electrodes" => "d")
    function chain_keys(cfg; files=files, decomposition="descomposición", manual="sin decisión manual",
                        qc_bad=["T8"], toml="toml")
        k = Dict{Symbol,String}()
        k[:load] = checkpoint_key(:load, cfg; input_files=files, subject="S01", session="T1", task="eyesclosed", run=1)
        k[:qc] = checkpoint_key(:qc, cfg, Dict(:load => k[:load]))
        k[:filter] = checkpoint_key(:filter, cfg, Dict(:load => k[:load]))
        k[:ica] = checkpoint_key(:ica, cfg, Dict(:filter => k[:filter]); decomposition, manual_rejected=manual)
        k[:montage] = checkpoint_key(:montage, cfg, Dict(:ica => k[:ica]); qc_bad_channels=qc_bad)
        k[:segmentation] = checkpoint_key(:segmentation, cfg, Dict(:montage => k[:montage], :ica => k[:ica]))
        k[:spectral] = checkpoint_key(:spectral, cfg, Dict(:segmentation => k[:segmentation]))
        k[:connectivity] = checkpoint_key(:connectivity, cfg, Dict(:segmentation => k[:segmentation], :montage => k[:montage]))
        k[:surrogates] = checkpoint_key(:surrogates, cfg, Dict(:connectivity => k[:connectivity]))
        k[:finalize] = checkpoint_key(:finalize, cfg, Dict(c => k[c] for c in
            (:load, :qc, :filter, :montage, :segmentation, :spectral, :connectivity)); toml_sha256=toml)
        return k
    end
    changed(a, b) = Set(s for s in keys(a) if a[s] != b[s])
    cfgat(name, changes...) = load_ss_config(synthetic_config(joinpath(CK_TMP, "dep-$(name)"), CK_BIDS, changes...))
    base = chain_keys(cfgat("base"))

    # Sección ⇒ etapas que la leen directamente
    perturb = [
        (("qc", "bad_channel_zscore_threshold") => "4.0",        [:qc]),
        (("bands", "ALPHA") => "[8.0, 12.0]",                     [:qc, :spectral, :connectivity, :surrogates]),
        (("filtering", "notch_hz") => "60.0",                     [:filter]),
        (("output", "figure_dpi") => "100",                       [:filter, :spectral, :connectivity, :surrogates, :finalize]),
        (("ica", "artifact_threshold") => "2.0",                  [:ica]),
        (("montage", "exclude_fp2") => "true",                    [:montage]),
        (("segmentation", "min_segments") => "12",                [:segmentation, :finalize]),
        (("baseline", "baseline_end_s") => "0.2",                 [:segmentation]),
        (("artifact_rejection", "max_amplitude_uv") => "60.0",    [:segmentation, :finalize]),
        (("spectral", "window_pct") => "20.0",                    [:spectral]),
        (("connectivity", "use_dwpli") => "true",                 [:connectivity, :surrogates]),
        (("graph", "density") => "0.2",                           [:connectivity]),
        (("surrogates", "n_surrogates") => "30",                  [:surrogates]),
    ]
    for (i, (change, direct)) in enumerate(perturb)
        @testset "$(change.first[1]).$(change.first[2])" begin
            @test changed(base, chain_keys(cfgat("p$(i)", change))) == closure(direct)
        end
    end

    @testset "Entradas propias de cada etapa" begin
        cfg = cfgat("base")
        @test changed(base, chain_keys(cfg; files=merge(files, Dict("eeg" => "otro")))) == closure([:load])
        @test changed(base, chain_keys(cfg; manual=[1])) == closure([:ica])
        # El umbral del QC cambia sin que cambie la lista: el montaje no cambia
        @test !(:montage in changed(base, chain_keys(cfgat("p1", perturb[1][1]))))
        # La lista cambia: el montaje y lo posterior
        @test changed(base, chain_keys(cfg; qc_bad=["T8", "F7"])) == closure([:montage])
    end

    @testset "La máquina cambia ICA y lo posterior; no carga, QC ni filtrado" begin
        cfg = cfgat("base")
        rec = mock_recording(8, 3000)
        here  = CK._ica_cache_identity(rec, cfg).key
        other = CK._ica_cache_identity(rec, cfg; implementation=merge(CK._ica_implementation(), (machine="otra-máquina",))).key
        @test here != other
        a = chain_keys(cfg; decomposition=here); b = chain_keys(cfg; decomposition=other)
        @test changed(a, b) == closure([:ica])
        @test isempty(intersect(changed(a, b), Set([:load, :qc, :filter])))
    end

    @testset "Decisión manual: cuenta lo que rechaza, no la fecha (decisión 129)" begin
        # Mismas etiquetas, otra fecha y otra justificación: la misma identidad
        function key_with(labels, name)
            root = joinpath(CK_TMP, "manual-$(name)")
            o = chain_with_checkpoints(synthetic_config(root, CK_BIDS); until=:ica,
                                       before_ica = ctx -> write(joinpath(ctx.export_dir, "ica_labels.csv"), labels))
            return CK.TOML.parsefile(joinpath(ckdir(root), "ica.toml"))["key"]
        end
        k1 = key_with("component,label,fecha,justificacion\n1,artifact,2026-11-13T10:00,parpadeo\n2,brain,2026-11-13T10:00,\n", "a")
        k2 = key_with("component,label,fecha,justificacion\n1,artifact,2026-11-14T18:30,ojos\n2,brain,2026-11-14T18:30,otra\n", "b")
        k3 = key_with("component,label,fecha,justificacion\n1,artifact,2026-11-13T10:00,parpadeo\n2,artifact,2026-11-13T10:00,\n", "c")
        @test k1 == k2
        @test k1 != k3
    end

    @testset "Cierre: los bytes del TOML (decisión 130, F12)" begin
        p1 = synthetic_config(joinpath(CK_TMP, "f12-a"), CK_BIDS)
        p2 = synthetic_config(joinpath(CK_TMP, "f12-b"), CK_BIDS)
        write(p2, read(p2, String) * "\n# Un comentario: mismo protocolo, otro texto\n")
        c1 = load_ss_config(p1); c2 = load_ss_config(p2)
        sections = ["qc", "bands", "filtering", "output", "ica", "montage", "segmentation", "baseline",
                    "artifact_rejection", "spectral", "connectivity", "graph", "surrogates"]
        @test all(deepeq(CK.checkpoint_section(c1, s), CK.checkpoint_section(c2, s)) for s in sections)
        h1 = bytes2hex(open(sha256, p1)); h2 = bytes2hex(open(sha256, p2))
        @test h1 != h2
        a = chain_keys(c1; toml=h1); b = chain_keys(c2; toml=h2)
        @test changed(a, b) == Set([:finalize])
        @test chain_keys(c1; toml=h1) == a
    end

    @testset "Cierre: los bytes del TOML solo entran en su identidad" begin
        # La separación está en la tabla de entradas y la función la hace
        # cumplir: ninguna etapa de 01 a 09 acepta la huella del TOML
        @test [s for (s, ins) in CK.CHECKPOINT_INPUTS if :toml_sha256 in ins] == [:finalize]
        cfg = cfgat("base"); k = chain_keys(cfg)
        for s in CHECKPOINT_STAGES
            s == :finalize && continue
            consumed = Dict(c => k[c] for c in CHECKPOINT_SPEC[s].consumes)
            own = Dict(:load => (input_files=files, subject="S01", session="T1", task="eyesclosed", run=1),
                       :ica => (decomposition="d", manual_rejected="sin decisión manual"),
                       :montage => (qc_bad_channels=["T8"],))
            @test_throws ErrorException checkpoint_key(s, cfg, consumed; get(own, s, NamedTuple())..., toml_sha256="h")
        end
        # Y el cierre no se calcula sin ella
        @test_throws ErrorException checkpoint_key(:finalize, cfg,
            Dict(c => k[c] for c in CHECKPOINT_SPEC[:finalize].consumes))
    end

    @testset "Cierre: la copia del TOML es la de la identidad, o su ausencia queda registrada" begin
        prov = CK.TOML.parsefile(joinpath(ckdir(CK_A), "finalize.toml"))
        h = bytes2hex(open(sha256, joinpath(CK_A, "config", "pipeline.toml")))
        snap = "config_snapshot.toml"     # relativo a la carpeta de la etapa (decisión 133)
        @test prov["facts"]["toml_sha256"] == h
        @test prov["products"]["files"][snap] == h
        @test isempty(prov["products"]["absent"])
        # Las tablas globales están fuera de la carpeta de la etapa y dentro del almacén
        @test "../../../../subjects_index.csv" in keys(prov["products"]["files"])
        @test "../../../../qc/qc_decision_table.csv" in keys(prov["products"]["files"])

        # La copia falla: se registra su ausencia y no se da por escrita
        root = joinpath(CK_TMP, "f12-ausente")
        chain_with_checkpoints(synthetic_config(root, CK_BIDS);
            before_finalize_ck = ctx -> rm(joinpath(ctx.export_dir, "config_snapshot.toml")))
        p = CK.TOML.parsefile(joinpath(ckdir(root), "finalize.toml"))
        @test p["products"]["absent"] == [snap]
        @test !haskey(p["products"]["files"], snap)
        @test p["facts"]["config_snapshot"] == "ausente"
        @test load_checkpoint(ckdir(root), :finalize; store=ckstore(root)) === nothing

        # Una copia que no es la del TOML de la identidad: ALTERADO
        r2 = copy_results(CK_A, joinpath(CK_TMP, "f12-otra"))
        sp = joinpath(dirname(ckdir(r2)), "config_snapshot.toml")
        chmod(sp, 0o644); write(sp, read(sp, String) * "\n# cambiado\n")
        @test reason_of(() -> load_checkpoint(ckdir(r2), :finalize; store=ckstore(r2))) == :ALTERADO
    end
end

# ── Checkpoints_ICA ──────────────────────────────────────────

@testset "Checkpoints_ICA" begin
    function ica_with(labels, name)
        root = joinpath(CK_TMP, "ica-$(name)")
        o = chain_with_checkpoints(synthetic_config(root, CK_BIDS); until=:ica,
            before_ica = labels === nothing ? nothing : ctx -> write(joinpath(ctx.export_dir, "ica_labels.csv"), labels))
        return o, load_checkpoint(ckdir(root), :ica; store=ckstore(root)), ckdir(root)
    end
    for (name, labels, source, rejected) in [
            ("propuesta automática", nothing, "auto", nothing),
            ("decisión que retira uno", "component,label\n3,artifact\n", "manual", [3]),
            ("decisión que no retira ninguno", "component,label\n" * join(["$(i),brain" for i in 1:16], "\n") * "\n", "manual", Int[])]
        @testset "$(name)" begin
            o, out, dir = ica_with(labels, replace(name, " " => "-"))
            @test out.rejection_source == source
            @test deepeq(out, o[:ica])
            rejected === nothing || @test out.result.rejected_components == rejected
            prov = CK.TOML.parsefile(joinpath(dir, "ica.toml"))
            @test prov["facts"]["rejected"] == out.result.rejected_components
            # La descomposición se lee de la caché, comprobada; no está en el checkpoint
            @test length(prov["references"]) == 1
            @test prov["references"][1]["identity"] == prov["facts"]["decomposition_identity"]
            @test deserialize(joinpath(dir, "ica.jls"))["data"]["decomposition"] === nothing
            # El .jls lleva la señal tras el rechazo y poco más: no las activaciones
            @test filesize(joinpath(dir, "ica.jls")) < sizeof(out.rec.data) + sizeof(out.result.activations) ÷ 2
        end
    end
end

# ── Checkpoints_Rename (decisión 133, P-1) ───────────────────

@testset "Checkpoints_Rename" begin
    # Como Lab: la etapa se calcula en una carpeta temporal dentro del almacén
    # y se renombra al terminar. Las referencias y los productos son relativos
    # a la carpeta de la etapa, así que sobreviven al renombrado.
    store = joinpath(CK_TMP, "almacen-p1")
    tmp   = joinpath(store, ".en-curso", "m1", "etapa")
    final = joinpath(store, "etapa")
    chain_with_checkpoints(synthetic_config(tmp, CK_BIDS, SURR20...); store)
    mv(tmp, final)
    dir = ckdir(final)

    prov(s) = CK.TOML.parsefile(joinpath(dir, "$(s).toml"))
    @testset "Rutas registradas relativas a la carpeta de la etapa" begin
        @test only(prov(:ica)["references"])["path"] == "cache/ica_result.jls"
        @test all(startswith(r["path"], "cache/surrogate/") for r in prov(:surrogates)["references"])
        fk = collect(keys(prov(:filter)["products"]["files"]))
        @test !isempty(fk) && all(k -> !startswith(k, "..") && !startswith(k, "subjects/"), fk)
        @test haskey(prov(:finalize)["products"]["files"], "config_snapshot.toml")
        # Los consumidos, relativos a la raíz del almacén: apuntan a la carpeta temporal
        @test all(startswith(c["path"], ".en-curso/m1/etapa/") for c in prov(:qc)["consumed"])
    end

    @testset "Una etapa sin consumidos carga y verifica tras el renombrado" begin
        @test verify_checkpoint(dir, :load; store) == (true, nothing, "")
        @test load_checkpoint(dir, :load; store) isa LoadOutput
    end

    @testset "Tras el renombrado, lo único que no se encuentra es lo consumido" begin
        for s in (:ica, :surrogates, :finalize)
            ok, reason, msg = verify_checkpoint(dir, s; store)
            @test !ok && reason == :ALTERADO
            @test occursin("consumido", msg)
        end
    end
end

