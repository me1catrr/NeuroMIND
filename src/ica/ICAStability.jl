# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Repetibilidad de los componentes del perfil «estable»
# ═══════════════════════════════════════════════════════════════
#
#  Diagnóstico del perfil EXPERIMENTAL `estable` (ICAStable.jl). Es un
#  producto distinto de la descomposición, con su propia identidad
#  (`_ica_stability_identity`): depende de la descomposición de referencia,
#  de K, de las semillas, del criterio de emparejamiento y de este archivo.
#  Cambiar cualquiera de ellos invalida este diagnóstico, no la ICA de
#  referencia, cuya caché no incluye este archivo.
#
#    K = `stable_runs` descomposiciones más, con semillas seed+1 … seed+K,
#    sobre la misma señal blanqueada. Cada una se empareja con la de
#    referencia por asignación global uno a uno (algoritmo húngaro) que
#    maximiza la suma de correlaciones absolutas entre activaciones.
#    Solo las ejecuciones convergidas cuentan en el resumen; las demás se
#    emparejan y se guardan, sin contar. Si la referencia no converge, todas
#    las medidas son un diagnóstico exploratorio.
#
#    Es una aproximación simplificada INSPIRADA en ICASSO (Himberg,
#    Hyvärinen y Esposito, NeuroImage 22:1214-1222, 2004). ICASSO agrupa los
#    componentes de muchas estimaciones y mide la compacidad de cada grupo;
#    aquí no hay agrupamiento, ni índice de calidad, ni centrotipo.
#
#  Repetibilidad no es identificación: un componente puede repetirse en
#  todas las semillas sin ser la fuente, y una fuente bien recuperada puede
#  repetirse poco (test/validacion_ica_estable.jl). La estabilidad no cambia
#  la clasificación automática ni el rechazo.
#
#  API
#  ───
#    ica_stability(rec, cfg, ref)        → (components, runs, similarity, identity, …)
#
# ───────────────────────────────────────────────────────────────
#  Fichero    src/ica/ICAStability.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 12-11-2026
# ───────────────────────────────────────────────────────────────

"Versión del análisis de estabilidad. Forma parte de su identidad, no de la de la descomposición."
const ICA_STABILITY_VERSION = 1

"Umbrales de correlación absoluta que se informan. Exploratorios: ninguno es un criterio de validez."
const ICA_STABLE_THRESHOLDS = (0.90, 0.95, 0.99)

"Criterio de emparejamiento, tal como queda registrado en las tablas."
const ICA_STABLE_MATCHING = "hungarian_max_sum_abs_corr"

"""
    _ica_hungarian(C) -> Vector{Int}

Asignación uno a uno que MAXIMIZA `sum(C[i, p[i]])` en una matriz cuadrada
(algoritmo húngaro, O(n³)). `p[i]` es la columna asignada a la fila `i`.
"""
function _ica_hungarian(C::AbstractMatrix{<:Real})
    n = size(C, 1)
    size(C, 2) == n || throw(ArgumentError("la matriz de similitud debe ser cuadrada"))
    M   = maximum(C) .- Float64.(C)                  # coste a minimizar
    u   = zeros(n + 1); v = zeros(n + 1)
    p   = zeros(Int, n + 1); way = zeros(Int, n + 1)
    for i in 1:n
        p[1] = i
        j0   = 1
        minv = fill(Inf, n + 1)
        used = falses(n + 1)
        while true
            used[j0] = true
            i0 = p[j0]; delta = Inf; j1 = 0
            for j in 2:n+1
                used[j] && continue
                cur = M[i0, j-1] - u[i0+1] - v[j]
                if cur < minv[j]; minv[j] = cur; way[j] = j0; end
                if minv[j] < delta; delta = minv[j]; j1 = j; end
            end
            for j in 1:n+1
                if used[j]
                    u[p[j]+1] += delta; v[j] -= delta
                else
                    minv[j] -= delta
                end
            end
            j0 = j1
            p[j0] == 0 && break
        end
        while true
            j1 = way[j0]; p[j0] = p[j1]; j0 = j1
            j0 == 1 && break
        end
    end
    assignment = zeros(Int, n)
    for j in 2:n+1
        p[j] > 0 && (assignment[p[j]] = j - 1)
    end
    return assignment
end

"Curtosis en exceso de una serie (0 para una gaussiana). Descriptiva: no clasifica nada."
function _ica_excess_kurtosis(x::AbstractVector{<:Real})
    y  = x .- mean(x)
    m2 = mean(abs2, y)
    return m2 == 0 ? 0.0 : mean(y .^ 4) / m2^2 - 3
end


"Huella de este archivo: la implementación del análisis de estabilidad."
_ica_stability_sha() = bytes2hex(sha256(read(joinpath(@__DIR__, "ICAStability.jl"))))

"""
    _ica_stability_identity(rec, cfg, runs) -> NamedTuple

Identidad del análisis de estabilidad: la clave de la descomposición de
referencia, K, las semillas, el criterio de emparejamiento, los umbrales que
se informan, la versión y la huella de este archivo. Su `key` es un SHA-256
con la misma serialización canónica que la caché de ICA.
"""
function _ica_stability_identity(rec::EEGRecording, cfg::PipelineConfig, runs::Int)
    p     = _ica_effective_params(cfg, rec.meta.n_channels)
    seeds = collect(p.seed .+ (1:runs))
    decomposition = _ica_cache_key(rec, cfg)
    impl  = _ica_stability_sha()
    io = IOBuffer()
    _ica_put(io, "neuromind-ica-stability"); _ica_put(io, ICA_STABILITY_VERSION); _ica_put(io, impl)
    _ica_put(io, decomposition); _ica_put(io, runs); foreach(s -> _ica_put(io, s), seeds)
    _ica_put(io, ICA_STABLE_MATCHING); foreach(t -> _ica_put(io, t), ICA_STABLE_THRESHOLDS)
    return (key = bytes2hex(sha256(take!(io))), version = ICA_STABILITY_VERSION, sha256 = impl,
            files = ["ica/ICAStability.jl"], decomposition = decomposition, runs, seeds,
            matching = ICA_STABLE_MATCHING, thresholds = collect(ICA_STABLE_THRESHOLDS))
end

"Componentes que el clasificador automático propone rechazar (los que no son «brain»)."
function _ica_auto_proposal(A, S, fs, names, threshold)
    ev = evaluate_ica_components(compute_ica_features(Float64.(A), Float64.(S), fs, names);
                                 artifact_thresh=threshold)
    return String.(ev.artifact_type) .!= "brain"
end

"""
    ica_stability(rec, cfg, ref; runs) -> NamedTuple

Repetibilidad de cada componente de `ref` frente a `runs` descomposiciones con
las semillas `seed+1 … seed+runs`. Devuelve tres tablas:

- `components`: una fila por componente de referencia.
- `runs`: una fila por ejecución adicional.
- `similarity`: las matrices de similitud completas, en formato largo.

Las medidas se marcan como `"exploratory"` si la referencia no ha convergido.
No modifica `ref` ni la clasificación.
"""
function ica_stability(rec::EEGRecording, cfg::PipelineConfig, ref::ICAResult;
                       runs::Int = Int(get(cfg.ica, "stable_runs", 5)),
                       _run = _ica_stable_run)    # solo para pruebas: provocar el fallo de una ejecución
    p = _ica_effective_params(cfg, rec.meta.n_channels)
    p.profile == "estable" || throw(ArgumentError("ica_stability requiere ica.profile = \"estable\""))
    runs >= 1 || throw(ArgumentError("stable_runs debe ser al menos 1"))
    Z, V, _, _ = _ica_stable_whiten(rec)
    S_ref  = Float64.(ref.activations)
    n      = size(S_ref, 1)
    fs     = rec.meta.fs
    names  = rec.meta.channel_names
    thr    = Float64(get(cfg.ica, "artifact_threshold", 1.5))
    ref_ok = Bool(get(ref.diagnostics, "converged", false))
    measures = ref_ok ? "reference_converged" : "exploratory"
    step   = Float64.(get(ref.diagnostics, "component_step", fill(NaN, n)))
    ref_proposed = _ica_auto_proposal(ref.mixing_matrix, S_ref, fs, names, thr)

    run_rows = DataFrame(run=Int[], seed=Int[], status=String[], converged=Bool[], n_iter=Int[],
                         final_error=Float64[], matching=String[], assignment_score=Float64[],
                         measures=String[], error=String[])
    sim_rows = DataFrame(run=Int[], reference_component=Int[], run_component=Int[],
                         abs_corr=Float64[], paired=Bool[])
    r_match  = fill(NaN, n, runs)          # correlación con la pareja, por ejecución
    proposed = falses(n, runs)             # el clasificador propone la pareja en esa ejecución
    conv     = falses(runs)
    for k in 1:runs
        seed = p.seed + k
        # Una ejecución que falla (p. ej., la decorrelación de ICACore no se puede
        # calcular) no es una ejecución no convergida: se registra como fallida.
        r = try
            _run(Z, V, p.max_iter, p.tol, seed)
        catch e
            push!(run_rows, (k, seed, "failed", false, 0, NaN, ICA_STABLE_MATCHING, NaN, measures,
                             first(sprint(showerror, e), 300)))
            continue
        end
        C = abs.(cor(S_ref', r.S'))         # referencia × ejecución
        a = _ica_hungarian(C)
        prop_k = _ica_auto_proposal(r.A, r.S, fs, names, thr)
        conv[k] = r.converged
        for i in 1:n
            r_match[i, k]  = C[i, a[i]]
            proposed[i, k] = prop_k[a[i]]
            for j in 1:n
                push!(sim_rows, (k, i, j, C[i, j], a[i] == j))
            end
        end
        push!(run_rows, (k, seed, r.converged ? "converged" : "not_converged", r.converged, r.n_iter,
                         r.final_error, ICA_STABLE_MATCHING, sum(C[i, a[i]] for i in 1:n), measures, ""))
    end

    used = findall(conv)
    comp = DataFrame(
        component              = 1:n,
        reference_converged    = fill(ref_ok, n),
        measures               = fill(measures, n),
        runs                   = fill(runs, n),
        runs_converged         = fill(length(used), n),
        evaluable              = fill(!isempty(used), n),
        r_min                  = [isempty(used) ? missing : minimum(r_match[i, used]) for i in 1:n],
        r_median               = [isempty(used) ? missing : median(r_match[i, used]) for i in 1:n],
        n_r_ge_090             = [count(>=(0.90), r_match[i, used]) for i in 1:n],
        n_r_ge_095             = [count(>=(0.95), r_match[i, used]) for i in 1:n],
        n_r_ge_099             = [count(>=(0.99), r_match[i, used]) for i in 1:n],
        # Todas las ejecuciones que han producido una descomposición, convergidas o no
        r_min_all_runs         = [(v = filter(!isnan, r_match[i, :]); isempty(v) ? missing : minimum(v)) for i in 1:n],
        auto_proposed_reference = ref_proposed,
        auto_proposed_runs     = [count(proposed[i, used]) for i in 1:n],
        excess_kurtosis        = [_ica_excess_kurtosis(S_ref[i, :]) for i in 1:n],
        last_step_change       = step,
    )
    return (; components = comp, runs = run_rows, similarity = sim_rows,
              identity = _ica_stability_identity(rec, cfg, runs),
              measures, reference_converged = ref_ok, runs_converged = length(used), n_runs = runs,
              runs_failed = count(==("failed"), run_rows.status),
              thresholds = ICA_STABLE_THRESHOLDS, matching = ICA_STABLE_MATCHING,
              seeds = collect(p.seed .+ (1:runs)))
end
