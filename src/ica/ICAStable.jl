# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — Perfil experimental de ICA «estable»
# ═══════════════════════════════════════════════════════════════
#
#  Perfil EXPERIMENTAL. No es el protocolo vigente: el perfil por defecto
#  sigue siendo "eeg_julia", cuyo código (ICACore.jl) este archivo no toca.
#
#  Separa tres cosas que no se deducen unas de otras:
#
#    convergencia     ¿terminó FastICA según su criterio de parada?
#                     (por ejecución, y cuánto se movería cada componente
#                     en una iteración más)
#    repetibilidad    ¿aparece cada componente en descomposiciones
#                     calculadas desde otras semillas?
#    interpretación   qué origen fisiológico tiene un componente. Este
#                     archivo NO la evalúa: ni la convergencia ni la
#                     repetibilidad la certifican.
#
#  Descomposición de referencia
#  ────────────────────────────
#    FastICA simétrico (tanh) de ICACore.jl, un intento, con la semilla de
#    la configuración y hasta `stable_max_iter` iteraciones. Si no converge,
#    la descomposición se conserva con sus diagnósticos, pero no es válida
#    para análisis definitivo y la etapa no aplica rechazos.
#
#  Repetibilidad
#  ─────────────
#    En ICAStability.jl, con su propia identidad: cambiar K o el
#    emparejamiento invalida ese diagnóstico, no la descomposición.
#
#  La estabilidad no cambia la clasificación automática ni el rechazo.
#
#  API
#  ───
#    run_ica_stable(rec, cfg)            → ICAResult
#
# ───────────────────────────────────────────────────────────────
#  Fichero    src/ica/ICAStable.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 12-11-2026
# ───────────────────────────────────────────────────────────────

"Versión del centrado del perfil `estable`. Forma parte de su implementación."
const ICA_STABLE_CENTERING = 1

"""
    _ica_stable_center(X) -> Matrix{Float64}

Centrado del perfil `estable`: a cada canal se le resta su media, calculada
como `sum(fila) / n` sobre `Float64`. Es la única forma de centrar del perfil;
`eeg_julia` conserva la suya (ICACore.jl).
"""
function _ica_stable_center(X::AbstractMatrix{<:Real})
    Xf = Matrix{Float64}(X)
    n  = size(Xf, 2)
    return Xf .- (sum(Xf; dims=2) ./ n)
end

"""
    _ica_component_step(W, Z) -> Vector{Float64}

Cuánto se movería cada componente en una iteración más de FastICA desde `W`:
`|1 − |diag(W_nueva·W')||`, el mismo término que usa el criterio de parada de
`_fastica_attempt`, pero componente a componente. Describe dónde está la falta
de convergencia; no dice por qué.
"""
function _ica_component_step(W::AbstractMatrix{<:Real}, Z::AbstractMatrix{<:Real})
    n  = size(Z, 2)
    GY = tanh.(W * Z)
    Gp = 1 .- GY .^ 2
    Wn = try
        _sym_decorr((GY * Z') / n - Diagonal(vec(mean(Gp; dims=2))) * W)
    catch
        return fill(NaN, size(W, 1))      # la decorrelación no se puede calcular: sin diagnóstico
    end
    return abs.(1 .- abs.(diag(Wn * W')))
end

"Señal centrada y blanqueada del perfil `estable`."
function _ica_stable_whiten(rec::EEGRecording)
    Xc = _ica_stable_center(rec.data)
    Z, V, λ_k, λ_all = _whiten_pca(Xc, size(Xc, 1))
    return Z, V, λ_k, λ_all
end

"Una ejecución de FastICA del perfil, con su diagnóstico de convergencia por componente."
function _ica_stable_run(Z, V, max_iter::Int, tol::Float64, seed::Int)
    n_comp = size(Z, 1)
    W, n_iter, err, ok = _fastica_attempt(Z, n_comp, max_iter, tol, 1.0, seed)
    S = W * Z
    A = inv(W * V)
    return (; W, S, A, n_iter, final_error = err, converged = ok, seed,
              step = _ica_component_step(W, Z))
end

"""
    run_ica_stable(rec, cfg) -> ICAResult

Descomposición de referencia del perfil `estable`. Mismo formato que `run_ica`.
`diagnostics["reference_valid"]` es `false` si no ha convergido: la
descomposición se conserva para su estudio, pero no es válida para análisis
definitivo.
"""
function run_ica_stable(rec::EEGRecording, cfg::PipelineConfig)::ICAResult
    t0 = time()
    p  = _ica_effective_params(cfg, rec.meta.n_channels)
    p.profile == "estable" || throw(ArgumentError("run_ica_stable requiere ica.profile = \"estable\""))
    Z, V, λ_k, λ_all = _ica_stable_whiten(rec)
    r = _ica_stable_run(Z, V, p.max_iter, p.tol, p.seed)

    n_comp    = size(r.S, 1)
    total_var = sum(var(rec.data; dims=2))
    var_exp   = [sum(r.A[:, i] .^ 2) * var(r.S[i, :]) for i in 1:n_comp] ./ total_var
    λ_pos     = filter(>(0), λ_all)
    λ_k_pos   = filter(>(0), λ_k)
    WW        = r.W * r.W'

    diagnostics = Dict{String,Any}(
        "converged"        => r.converged,
        "n_iter"           => r.n_iter,
        "max_iter"         => p.max_iter,
        "final_error"      => r.final_error,
        "tol"              => p.tol,
        "seed"             => p.seed,
        "seed0"            => p.seed,
        "n_attempts"       => 1,
        "max_attempts"     => 1,
        "restart_needed"   => false,
        "optimize_s"       => round(time() - t0; digits=3),
        "whitening_cond"   => isempty(λ_k_pos) ? NaN : maximum(λ_k_pos) / max(minimum(λ_k_pos), eps(Float64)),
        "pca_rank"         => n_comp,
        "pca_rank_full"    => rec.meta.n_channels,
        "pca_var_pct"      => isempty(λ_pos) ? NaN : 100.0 * sum(λ_k) / sum(λ_pos),
        "orthogonality"    => sqrt(sum(abs2, WW .- Matrix{Float64}(I, n_comp, n_comp))) / n_comp,
        "profile"          => "estable",
        "centering"        => ICA_STABLE_CENTERING,
        "reference_valid"  => r.converged,
        "component_step"   => r.step,
        "attempts"         => Dict{String,Any}[Dict{String,Any}(
            "attempt" => 1, "seed" => p.seed, "n_iter" => r.n_iter,
            "final_error" => r.final_error, "converged" => r.converged)],
    )
    return ICAResult(rec.meta, r.A, r.W * V, r.S, Int[], var_exp, diagnostics)
end

