# ───────────────────────────────────────────────────────────────
#  Fichero    test/run_stage_tests.jl
#  Modificado 15-11-2026
# ───────────────────────────────────────────────────────────────
#
# Conjuntos de la entrega 2f-2 (plan de 2f, apartado 7.2): ejecución por
# etapas con `run_stage`, identidad con `stage_key` y grafo de etapas, y el
# alcance de las exclusiones de la huella de contenido (P-4). Datos
# sintéticos. Lo incluye runtests.jl después de checkpoints_tests.jl, cuyas
# ejecuciones compartidas (CK_A, CK_BIDS, …) reutiliza.

# ── Ayudantes ────────────────────────────────────────────────

"Carpeta raíz de la etapa `s` en un almacén con una carpeta por etapa, como Lab."
stage_root(store, s; tag="") = joinpath(store, "$(s)$(tag)")

"""
Ejecuta las etapas de `stages` con `run_stage`, cada una en su carpeta, con la
configuración que da `config(root, s)`. Devuelve etapa ⇒ carpeta `checkpoints/`.
"""
function staged_chain(store, config; stages=CHECKPOINT_STAGES, dirs=Dict{Symbol,String}(), tag="")
    dirs = copy(dirs)
    for s in stages
        r = redirect_stdout(devnull) do
            run_stage(s, config(stage_root(store, s; tag), s); inputs=dirs, store)
        end
        dirs[s] = r.dir
    end
    return dirs
end

"Igual, pero cada etapa en un proceso de Julia nuevo. Devuelve también el tiempo de cada proceso."
function staged_chain_processes(store, config)
    dirs = Dict{Symbol,String}(); times = Dict{Symbol,Float64}()
    for s in CHECKPOINT_STAGES
        cfg = config(stage_root(store, s), s)
        code = """
            using NeuroMIND
            inputs = Dict{Symbol,String}($(join(["$(repr(k)) => $(repr(v))" for (k, v) in dirs], ", ")))
            r = redirect_stdout(devnull) do
                run_stage($(repr(s)), $(repr(cfg)); inputs, store=$(repr(store)))
            end
            print(r.dir)
            """
        cmd = `$(Base.julia_cmd()) --project=$(pkgdir(NeuroMIND)) --startup-file=no -e $code`
        t = @elapsed out = read(pipeline(cmd; stderr=devnull), String)
        dirs[s] = String(strip(out)); times[s] = t
    end
    return dirs, times
end

stage_prov(dir, s) = CK.TOML.parsefile(joinpath(dir, "$(s).toml"))

"Lanza `f` y devuelve el `StageError` o `CheckpointError`, o `nothing`."
function stage_error(f)
    try
        f(); return nothing
    catch e
        (e isa StageError || e isa CheckpointError) || rethrow()
        return e
    end
end

surr_config(root, _) = synthetic_config(root, CK_BIDS, SURR20...)

# Ejecución por etapas compartida, en el mismo proceso
const RS_STORE = joinpath(CK_TMP, "rs")
const RS_DIRS  = staged_chain(RS_STORE, surr_config)
# La ICA con el perfil «default», sobre el mismo filtrado: otra descomposición
default_ica_config(root, _) = synthetic_config(root, CK_BIDS, SURR20..., ("ica", "profile") => "\"default\"")
const RS_ICA_DEFAULT = staged_chain(RS_STORE, default_ica_config; stages=(:ica,),
                                    dirs=Dict(:load => RS_DIRS[:load], :filter => RS_DIRS[:filter]), tag="-default")[:ica]

# ── Stage_Graph ──────────────────────────────────────────────

@testset "Stage_Graph" begin
    g = stage_graph()
    @test [n.stage for n in g] == collect(CHECKPOINT_STAGES)
    for n in g
        @test n.consumes == CHECKPOINT_SPEC[n.stage].consumes
        @test n.sections == CHECKPOINT_SPEC[n.stage].sections
        @test issubset(n.consumes, n.needs)
    end
    byname = Dict(n.stage => n for n in g)
    # P-3: [output] cambia las figuras y entra en la identidad de quien las escribe
    @test all("output" in byname[s].sections for s in (:filter, :spectral, :connectivity, :surrogates, :finalize))
    # El montaje lee la lista del QC por su contenido, sin consumir su checkpoint
    @test byname[:montage].needs == [:ica, :qc] && byname[:montage].consumes == [:ica]
    @test byname[:montage].inputs == [:qc_bad_channels]
    @test byname[:finalize].inputs == [:toml_sha256]
end

# ── Stage_Key (E-2) ──────────────────────────────────────────

@testset "Stage_Key" begin
    @testset "stage_key da la identidad que run_stage registra, en las diez etapas" begin
        prev = Dict{Symbol,String}()
        for s in CHECKPOINT_STAGES
            cfg = surr_config(stage_root(RS_STORE, s), s)
            @test stage_key(s, cfg; inputs=prev, store=RS_STORE) == stage_prov(RS_DIRS[s], s)["key"]
            prev[s] = RS_DIRS[s]
        end
    end

    @testset "La identidad no depende de dónde estén las carpetas" begin
        # La misma configuración y los mismos datos que el orquestador (CK_A),
        # en otro almacén y con una carpeta por etapa: mismas identidades
        for s in CHECKPOINT_STAGES
            @test stage_prov(RS_DIRS[s], s)["key"] == stage_prov(ckdir(CK_A), s)["key"]
        end
        # Y stage_key, con las entradas en otro sitio del almacén, tampoco cambia
        moved = joinpath(RS_STORE, "trasladado"); mkpath(moved)
        cp(stage_root(RS_STORE, :filter), joinpath(moved, "filter"))
        dirs = Dict(:filter => ckdir(joinpath(moved, "filter")), :load => RS_DIRS[:load])
        cfg = surr_config(joinpath(CK_TMP, "rs-cfg-ica"), :ica)
        @test stage_key(:ica, cfg; inputs=dirs, store=CK_TMP) == stage_prov(RS_DIRS[:ica], :ica)["key"]
    end

    @testset "stage_key no escribe nada" begin
        root = joinpath(RS_STORE, "sin-efectos"); cfg = surr_config(root, :montage)
        before = tree(root)
        k = stage_key(:montage, cfg; inputs=RS_DIRS, store=RS_STORE)
        @test tree(root) == before
        @test k == stage_prov(RS_DIRS[:montage], :montage)["key"]
    end

    @testset "Reutilizar: cargar con la identidad esperada" begin
        cfg = surr_config(joinpath(RS_STORE, "reuso"), :spectral)
        k = stage_key(:spectral, cfg; inputs=RS_DIRS, store=RS_STORE)
        @test load_checkpoint(RS_DIRS[:spectral], :spectral; store=RS_STORE, expected_key=k) isa SpectralResult
        # Otra configuración del espectro: otra identidad, y no se reutiliza
        cfg2 = synthetic_config(joinpath(RS_STORE, "reuso-2"), CK_BIDS, SURR20..., ("spectral", "nfft") => "2048")
        k2 = stage_key(:spectral, cfg2; inputs=RS_DIRS, store=RS_STORE)
        @test k2 != k
        @test reason_of(() -> load_checkpoint(RS_DIRS[:spectral], :spectral; store=RS_STORE, expected_key=k2)) ==
              :OTRA_CONFIGURACION
    end

    @testset "Entradas que faltan o están fuera del almacén" begin
        cfg = surr_config(joinpath(RS_STORE, "faltan"), :montage)
        e = stage_error(() -> stage_key(:montage, cfg; inputs=Dict(:ica => RS_DIRS[:ica]), store=RS_STORE))
        @test e isa StageError && e.code == :ENTRADA_FALTA
        fuera = joinpath(CK_TMP, "rs-fuera"); cp(stage_root(RS_STORE, :qc), fuera)
        e = stage_error(() -> stage_key(:montage, cfg; inputs=Dict(:ica => RS_DIRS[:ica], :qc => ckdir(fuera)),
                                        store=RS_STORE))
        @test e isa CheckpointError && e.reason == :ALTERADO && occursin("fuera del almacén", e.message)
        e = stage_error(() -> stage_key(:load, surr_config(joinpath(CK_TMP, "rs-cfg-fuera"), :load); store=RS_STORE))
        @test e isa CheckpointError && e.reason == :ALTERADO && occursin("configuración", e.message)
    end
end

# ── Pipeline_Staged_Equivalence ──────────────────────────────

# Lo que el pipeline escribe con reloj o fecha, o con rutas de la ejecución
const RS_CLOCK = r"\"(timestamp|duration_s|optimize_s|saved_at)\"\s*:\s*[^,\n}]*"
skip_compare(f) = occursin("/checkpoints/", f) || endswith(f, "pipeline_log.txt") || endswith(f, ".toml") ||
                  endswith(f, "subjects_index.csv") || endswith(f, "qc_decision_table.csv") ||
                  occursin("cache/ica_result.jls", f) || occursin("cache/ica_config.hash", f) ||
                  endswith(f, "ica_labels.csv")

@testset "Pipeline_Staged_Equivalence" begin
    store = joinpath(CK_TMP, "rs-procesos")
    dirs, times = staged_chain_processes(store, surr_config)
    @info "Un proceso de Julia por etapa: " * join(["$(s) $(round(times[s], digits=1)) s" for s in CHECKPOINT_STAGES], " · ") *
          " · total $(round(sum(values(times)), digits=1)) s"

    @testset "Mismos checkpoints por contenido e identidad que el orquestador en memoria" begin
        for s in CHECKPOINT_STAGES
            a = stage_prov(ckdir(CK_A), s); b = stage_prov(dirs[s], s)
            @test b["content_sha256"] == a["content_sha256"]
            @test b["key"] == a["key"]
            @test b["content_sha256"] == checkpoint_content_sha256(CK_MEM[s])
        end
    end

    @testset "Tablas y figuras idénticas byte a byte; JSON salvo el reloj" begin
        a = ckstore(CK_A)
        staged = Dict{String,String}()       # archivo ⇒ carpeta de etapa que lo escribió
        for s in CHECKPOINT_STAGES, f in tree(joinpath(stage_root(store, s), "results"))
            skip_compare(f) && continue
            staged[f] = joinpath(stage_root(store, s), "results", f)
        end
        orch = [f for f in tree(a) if !skip_compare(f)]
        @test length(orch) > 40
        @test isempty(setdiff(orch, keys(staged)))
        differ = String[]
        for f in orch
            x = read(joinpath(a, f)); y = read(staged[f])
            if endswith(f, ".json")
                x = replace(String(x), RS_CLOCK => ""); y = replace(String(y), RS_CLOCK => "")
            end
            x == y || push!(differ, f)
        end
        @test isempty(differ)
    end
end

# ── Checkpoints_ICA (2f-2): decisión manual y caché entregadas ──

@testset "Checkpoints_ICA_RunStage" begin
    base = Dict(:load => RS_DIRS[:load], :filter => RS_DIRS[:filter])
    ident = stage_prov(RS_DIRS[:ica], :ica)["facts"]["decomposition_identity"]
    labels = joinpath(RS_STORE, "decision", "ica_labels.csv")
    mkpath(dirname(labels)); write(labels, "component,label\n3,artifact\n")
    cfg_for(name) = surr_config(joinpath(RS_STORE, "ica-$(name)"), :ica)
    exp_dir(name) = joinpath(RS_STORE, "ica-$(name)", "results", "subjects", "sub-S01", "ses-T1", "eyesclosed")

    @testset "Decisión verificada: se aplica y entra en la identidad" begin
        cfg = cfg_for("decision")
        k = stage_key(:ica, cfg; inputs=base, ica_labels=(labels, ident), store=RS_STORE)
        r = redirect_stdout(devnull) do
            run_stage(:ica, cfg; inputs=base, ica_labels=(labels, ident), store=RS_STORE)
        end
        @test r.key == k
        @test k != stage_prov(RS_DIRS[:ica], :ica)["key"]
        out = load_checkpoint(r.dir, :ica; store=RS_STORE, consumed_dirs=Dict(:filter => RS_DIRS[:filter]))
        @test out.rejection_source == "manual" && out.result.rejected_components == [3]
        @test stage_prov(r.dir, :ica)["facts"]["manual_decision"] == "applied"
    end

    @testset "Decisión sobre otra descomposición: falla sin aplicar ni escribir nada" begin
        cfg = cfg_for("otra")
        e = stage_error(() -> run_stage(:ica, cfg; inputs=base, ica_labels=(labels, "0"^64), store=RS_STORE))
        @test e isa StageError && e.code == :ICA_LABELS_OTRA_DESCOMPOSICION
        @test !isdir(exp_dir("otra"))
        e = stage_error(() -> stage_key(:ica, cfg; inputs=base, ica_labels=(labels, "0"^64), store=RS_STORE))
        @test e isa StageError && e.code == :ICA_LABELS_OTRA_DESCOMPOSICION
    end

    @testset "Etiquetas en la carpeta que no se han entregado: no se aplican" begin
        cfg = cfg_for("previas"); mkpath(exp_dir("previas"))
        write(joinpath(exp_dir("previas"), "ica_labels.csv"), "component,label\n1,artifact\n")
        e = stage_error(() -> run_stage(:ica, cfg; inputs=base, store=RS_STORE))
        @test e isa StageError && e.code == :ICA_LABELS_NO_ENTREGADAS
        @test !isfile(joinpath(exp_dir("previas"), "pipeline_log.txt"))
    end

    ica_cache_of(name; inputs=base, cache) = redirect_stdout(devnull) do
        run_stage(:ica, cfg_for(name); inputs, ica_cache=cache, store=RS_STORE)
    end

    @testset "Caché de la misma señal, con su huella: se usa" begin
        r = ica_cache_of("cache-ok"; cache=RS_DIRS[:ica])
        @test r.ica_cache.used
        @test stage_prov(r.dir, :ica)["facts"]["ica_cache"] == "entregada y usada"
        # La descomposición, la señal y el rechazo son los calculados, y la
        # identidad no cambia. La etapa dice que la leyó de caché (`status`).
        fd = Dict(:filter => RS_DIRS[:filter])
        a = load_checkpoint(RS_DIRS[:ica], :ica; store=RS_STORE)
        b = load_checkpoint(r.dir, :ica; store=RS_STORE, consumed_dirs=fd)
        @test a.status == :computed && b.status == :cache
        @test all(checkpoint_content_sha256(getfield(a, f)) == checkpoint_content_sha256(getfield(b, f))
                  for f in (:result, :rec, :rejection_source, :classification_ok))
        @test r.key == stage_prov(RS_DIRS[:ica], :ica)["key"]
    end

    @testset "Caché sin huella de bytes, alterada o de otra señal: se recalcula y se dice" begin
        # Sin huella: un checkpoint de ICA sin referencia a la caché
        sin = joinpath(RS_STORE, "cache-sin-huella"); cp(stage_root(RS_STORE, :ica), sin)
        d = ckdir(sin)
        p = stage_prov(d, :ica); p["references"] = Any[]
        chmod(joinpath(d, "ica.toml"), 0o644)
        open(io -> CK.TOML.print(io, p; sorted=true), joinpath(d, "ica.toml"), "w")
        r = ica_cache_of("cache-sin"; cache=d)
        @test !r.ica_cache.used && occursin("huella", r.ica_cache.reason)
        @test startswith(stage_prov(r.dir, :ica)["facts"]["ica_cache"], "entregada y no usada")
        # Un byte cambiado en ica_result.jls
        alt = joinpath(RS_STORE, "cache-alterada"); cp(stage_root(RS_STORE, :ica), alt)
        flip_byte!(joinpath(alt, "results", "subjects", "sub-S01", "ses-T1", "eyesclosed", "cache", "ica_result.jls"))
        r = ica_cache_of("cache-alt"; cache=ckdir(alt))
        @test !r.ica_cache.used && occursin("huella", r.ica_cache.reason)
        # De otra descomposición: la ICA del perfil «default» sobre la misma señal
        r = ica_cache_of("cache-otra"; cache=RS_ICA_DEFAULT)
        @test !r.ica_cache.used && occursin("otra señal o de otra configuración de ICA", r.ica_cache.reason)
    end

    @testset "Caché entregada fuera del almacén: ALTERADO" begin
        fuera = joinpath(CK_TMP, "rs-cache-fuera"); cp(stage_root(RS_STORE, :ica), fuera)
        e = stage_error(() -> ica_cache_of("cache-fuera"; cache=ckdir(fuera)))
        @test e isa CheckpointError && e.reason == :ALTERADO
    end
end

# ── P-9 y P-10 ───────────────────────────────────────────────

@testset "Stage_Consumed_ICA" begin
    @testset "P-10: una ICA calculada con otra configuración de ICA no se consume" begin
        cfg = default_ica_config(joinpath(RS_STORE, "montaje-otra-ica"), :montage)
        e = stage_error(() -> run_stage(:montage, cfg; inputs=RS_DIRS, store=RS_STORE))
        @test e isa CheckpointError && e.reason == :OTRA_CONFIGURACION
        @test occursin("otra configuración de ICA o en otra máquina", e.message)
    end

    @testset "P-9: el rótulo de entrada de la segmentación corresponde a la ICA consumida" begin
        root = joinpath(CK_TMP, "p9")
        ctx = redirect_stdout(devnull) do; pipeline_context(synthetic_config(root, CK_BIDS)); end
        try
            run = CK._checkpoint_run(ctx, nothing)
            ld, qc, flt, ica, mont, seg = redirect_stdout(devnull) do
                ld = stage_load(ctx);            CK._checkpoint_after!(run, :load, ld)
                qc = stage_qc(ctx, ld);          CK._checkpoint_after!(run, :qc, qc)
                flt = stage_filter(ctx, ld);     CK._checkpoint_after!(run, :filter, flt)
                ica = stage_ica(ctx, flt);       CK._checkpoint_after!(run, :ica, ica; filter=flt)
                mont = stage_montage(ctx, ica, qc); CK._checkpoint_after!(run, :montage, mont; qc=qc)
                ld, qc, flt, ica, mont, stage_segmentation(ctx, mont, ica)
            end
            # Una ICA con otra lista de rechazados que la que rotuló el motor
            other = isempty(ica.result.rejected_components) ? [1] : Int[]
            r = ica.result
            fake = ICAOutput(ICAResult(r.meta, r.mixing_matrix, r.unmixing_matrix, r.activations, other,
                                       r.variance_explained, r.diagnostics), ica.rec, ica.status, ica.rejection_source,
                             ica.classification_ok)
            e = stage_error(() -> CK._checkpoint_after!(run, :segmentation, seg; ica=fake))
            @test e isa StageError && e.code == :SEGMENTACION_ROTULO
            @test !isfile(joinpath(ckdir(root), "segmentation.toml"))
            # Con la ICA verdadera, se guarda
            CK._checkpoint_after!(run, :segmentation, seg; ica=ica)
            @test isfile(joinpath(ckdir(root), "segmentation.toml"))
        finally
            close(ctx.log_io)
        end
    end
end

# ── Checkpoints_Finalize_Reuse (decisión 130, F12) ───────────

@testset "Checkpoints_Finalize_Reuse" begin
    root = stage_root(RS_STORE, :finalize)
    cfg = joinpath(root, "config", "pipeline.toml")
    old_key = stage_prov(RS_DIRS[:finalize], :finalize)["key"]
    log = joinpath(dirname(RS_DIRS[:finalize]), "pipeline_log.txt")
    log_before = read(log)

    # Otro texto del TOML (un comentario), misma configuración efectiva
    cfg2 = joinpath(RS_STORE, "finalize-otro", "config", "pipeline.toml")
    mkpath(dirname(cfg2)); write(cfg2, read(cfg, String) * "\n# otro texto\n")
    k2 = stage_key(:finalize, cfg2; inputs=RS_DIRS, store=RS_STORE)
    @test k2 != old_key
    # El cierre existente no es el pedido: no se reutiliza
    @test reason_of(() -> load_checkpoint(RS_DIRS[:finalize], :finalize; store=RS_STORE, expected_key=k2)) ==
          :OTRA_CONFIGURACION
    # En la misma carpeta, run_stage no lo sobrescribe, y no toca su registro
    chmod(cfg, 0o644); write(cfg, read(cfg2))
    e = stage_error(() -> run_stage(:finalize, cfg; inputs=RS_DIRS, store=RS_STORE))
    @test e isa StageError && e.code == :CHECKPOINT_EXISTE
    @test read(log) == log_before
    @test stage_prov(RS_DIRS[:finalize], :finalize)["key"] == old_key
    # En una carpeta nueva se calcula otro cierre, con la identidad pedida
    r = redirect_stdout(devnull) do
        run_stage(:finalize, cfg2; inputs=RS_DIRS, store=RS_STORE)
    end
    @test r.key == k2
    @test stage_prov(r.dir, :finalize)["facts"]["toml_sha256"] == bytes2hex(open(sha256, cfg2))
end

# ── Stage_Log_On_Error (apartado 3.9) y P-5 ──────────────────

@testset "Stage_Log_On_Error" begin
    csd(root) = synthetic_config(root, CK_BIDS, ("connectivity", "use_csd") => "true")
    logof(root) = joinpath(root, "results", "subjects", "sub-S01", "ses-T1", "eyesclosed", "pipeline_log.txt")

    @testset "Pipeline completo: la excepción original, el registro cerrado con el fallo" begin
        root = joinpath(CK_TMP, "log-pipeline")
        e = try
            redirect_stdout(devnull) do; run_single_subject_pipeline(csd(root)); end; nothing
        catch err
            err
        end
        @test e isa DimensionMismatch          # D-21: la de la etapa, sin envolver
        txt = read(logof(root), String)
        @test occursin("La etapa connectivity ha fallado", txt)
        @test endswith(rstrip(txt), "No se ha guardado ningún checkpoint de esta etapa.")
    end

    @testset "Con checkpoints: tampoco queda el de la etapa que falla" begin
        root = joinpath(CK_TMP, "log-pipeline-ck")
        @test_throws DimensionMismatch redirect_stdout(devnull) do
            run_single_subject_pipeline(csd(root); checkpoints=true)
        end
        d = ckdir(root)
        @test isfile(joinpath(d, "segmentation.toml"))
        @test !ispath(joinpath(d, "connectivity.toml")) && !ispath(joinpath(d, "connectivity.jls"))
        @test isempty([f for f in readdir(d) if startswith(f, ".tmp-")])
    end

    @testset "run_stage: la excepción original, el registro cerrado y sin checkpoint" begin
        root = joinpath(RS_STORE, "log-conectividad")
        @test_throws DimensionMismatch redirect_stdout(devnull) do
            run_stage(:connectivity, csd(root); inputs=RS_DIRS, store=RS_STORE)
        end
        txt = read(logof(root), String)
        @test occursin("La etapa connectivity ha fallado", txt)
        d = joinpath(dirname(logof(root)), "checkpoints")
        @test !ispath(joinpath(d, "connectivity.toml")) && !ispath(joinpath(d, "connectivity.jls"))
    end

    @testset "P-5: rechazar una carpeta con checkpoints no borra su registro" begin
        root = joinpath(CK_TMP, "p5"); cp(CK_A, root)
        before = read(logof(root))
        @test_throws ErrorException redirect_stdout(devnull) do
            run_single_subject_pipeline(joinpath(root, "config", "pipeline.toml"); checkpoints=true)
        end
        @test read(logof(root)) == before
    end
end

# ── Checkpoints_Content_Scope (P-4, INS-003) ─────────────────

"Copia de `x` con el campo `f` cambiado por `v`."
with_field(x::T, f::Symbol, v) where {T} = T((g === f ? v : getfield(x, g) for g in fieldnames(T))...)

@testset "Checkpoints_Content_Scope" begin
    H = checkpoint_content_sha256

    @testset "Solo dos datos quedan fuera de la huella, por su posición" begin
        @test CK._EXCLUDED_FIELDS == ((RecordingMeta, :bids_path),)
        @test CK._EXCLUDED_DICT_KEYS == ((ICAResult, :diagnostics, ("optimize_s",)),)
        @test !isdefined(CK, :_CLOCK_KEYS)
    end

    @testset "Las claves de reloj cuentan en cualquier otro diccionario" begin
        # El mismo `_canon!` de diccionarios da la identidad de las secciones de
        # configuración: tampoco ahí se excluye ninguna clave por su nombre
        for k in ("optimize_s", "duration_s", "timestamp", "saved_at", "bids_path")
            @test H(Dict{String,Any}(k => 1.0)) != H(Dict{String,Any}(k => 2.0))
            @test H(Dict{String,Any}("x" => 0)) != H(Dict{String,Any}("x" => 0, k => 0))
            @test H(Dict(Symbol(k) => 1.0)) != H(Dict(Symbol(k) => 2.0))
            @test H(Dict("a" => Dict{String,Any}(k => 1.0))) != H(Dict("a" => Dict{String,Any}(k => 2.0)))
            @test H(NamedTuple{(Symbol(k),)}((1.0,))) != H(NamedTuple{(Symbol(k),)}((2.0,)))
        end
        # Ni en el diccionario de un tipo del motor que no es la ICA
        sp = load_checkpoint(ckdir(CK_A), :spectral; store=ckstore(CK_A))
        for k in ("optimize_s", "duration_s", "timestamp", "saved_at")
            p1 = merge(sp.params, Dict{String,Any}(k => 1.0)); p2 = merge(sp.params, Dict{String,Any}(k => 2.0))
            @test H(with_field(sp, :params, p1)) != H(with_field(sp, :params, p2))
        end
    end

    @testset "ICAResult.diagnostics: solo optimize_s, y solo en su raíz" begin
        out = load_checkpoint(ckdir(CK_A), :ica; store=ckstore(CK_A))
        r = out.result; d = r.diagnostics
        @test haskey(d, "optimize_s") && !isempty(d["attempts"])
        with_diag(d2) = with_field(out, :result, with_field(r, :diagnostics, d2))
        set(k, v) = (d2 = deepcopy(d); d2[k] = v; with_diag(d2))
        @test H(with_diag(deepcopy(d))) == H(out)
        @test H(set("optimize_s", d["optimize_s"] + 100.0)) == H(out)
        @test H(set("n_iter", d["n_iter"] + 1)) != H(out)
        @test H(set("final_error", d["final_error"] * 2)) != H(out)
        for k in ("duration_s", "timestamp", "saved_at")
            @test H(set(k, 1.0)) != H(out)
        end
        # Dentro de un intento, `optimize_s` ya no es la clave excluida
        nested(v) = (d2 = deepcopy(d); d2["attempts"] = [Dict{String,Any}(a) for a in d["attempts"]];
                     d2["attempts"][1]["optimize_s"] = v; with_diag(d2))
        @test H(nested(1.0)) != H(nested(2.0))
        # El mismo diccionario fuera de su sitio entra completo
        @test H(d) != H(merge(d, Dict{String,Any}("optimize_s" => d["optimize_s"] + 100.0)))
    end

    @testset "RecordingMeta: solo bids_path queda fuera" begin
        a = load_checkpoint(ckdir(CK_A), :load; store=ckstore(CK_A))
        m = a.rec.meta
        rec(m2) = with_field(a, :rec, with_field(a.rec, :meta, m2))
        @test H(rec(with_field(m, :bids_path, "/otro/sitio"))) == H(a)
        @test H(rec(with_field(m, :fs, m.fs + 1))) != H(a)
        @test H(rec(with_field(m, :subject_id, m.subject_id * "x"))) != H(a)
        @test H(rec(with_field(m, :channel_names, reverse(m.channel_names)))) != H(a)
    end

    @testset "bids_path: la ubicación del BIDS no cuenta; su contenido, sí" begin
        # El mismo conjunto BIDS en otro sitio, las diez etapas con run_stage
        moved = joinpath(CK_TMP, "bids-trasladado"); cp(CK_BIDS, moved)
        moved_config(root, _) = synthetic_config(root, moved, SURR20...)
        dirs = staged_chain(RS_STORE, moved_config; tag="-bids")
        a = load_checkpoint(RS_DIRS[:load], :load; store=RS_STORE)
        b = load_checkpoint(dirs[:load], :load; store=RS_STORE)
        @test a.rec.meta.bids_path != b.rec.meta.bids_path
        @test occursin("bids-trasladado", b.rec.meta.bids_path)
        for s in CHECKPOINT_STAGES
            @test stage_prov(dirs[s], s)["content_sha256"] == stage_prov(RS_DIRS[s], s)["content_sha256"]
        end
        # La identidad tampoco cambia, salvo la del cierre, que lleva el texto del TOML
        for s in CHECKPOINT_STAGES
            same = stage_prov(dirs[s], s)["key"] == stage_prov(RS_DIRS[s], s)["key"]
            @test same == (s != :finalize)
        end

        # Otro contenido en la misma ubicación: otra huella y otra identidad
        alt = joinpath(CK_TMP, "bids-alterado"); cp(CK_BIDS, alt)
        alt_config(root, _) = synthetic_config(root, alt, SURR20...)
        before = staged_chain(RS_STORE, alt_config; stages=(:load,), tag="-bids-antes")[:load]
        @test stage_prov(before, :load)["content_sha256"] == stage_prov(RS_DIRS[:load], :load)["content_sha256"]
        eeg = only(joinpath(d, f) for (d, _, fs) in walkdir(alt) for f in fs if endswith(f, ".eeg"))
        flip_byte!(eeg)
        after = staged_chain(RS_STORE, alt_config; stages=(:load,), tag="-bids-despues")[:load]
        @test stage_prov(after, :load)["content_sha256"] != stage_prov(before, :load)["content_sha256"]
        @test stage_prov(after, :load)["key"] != stage_prov(before, :load)["key"]
        c = load_checkpoint(after, :load; store=RS_STORE)
        @test c.rec.meta.bids_path == load_checkpoint(before, :load; store=RS_STORE).rec.meta.bids_path
    end
end
