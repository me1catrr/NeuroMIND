# ═══════════════════════════════════════════════════════════════
#  NeuroMIND — FastICA (núcleo)
# ═══════════════════════════════════════════════════════════════
#
#  FastICA simétrico con blanqueo PCA.
#  Implementación propia en Julia puro (LinearAlgebra + Random);
#  no usa MultivariateStats.
#  Portado de EEG_Julia/Pluto/ICA/ICA.jl.
#
#  Perfiles
#  ────────
#    "eeg_julia"  n_comp = n_ch, max_iter=512, tol=1e-7, seed=1234
#                 A = inv(W_total); max_attempts forzado a 1
#    "default"    parámetros desde config/pipeline.toml [ica]
#                 A = pinv(W_total)  (n_comp ≤ n_ch)
#
#  Diagnósticos (ICAResult.diagnostics)
#  ────────────────────────────────────
#    converged, n_iter, max_iter, final_error, tol, seed,
#    n_attempts, max_attempts, restart_needed, optimize_s,
#    whitening_cond, pca_rank, pca_var_pct, orthogonality,
#    attempts::Vector{Dict}  (detalle por intento)
#
#  API pública (exportada por NeuroMIND)
#  ─────────────────────────────────────
#    run_ica(rec, cfg) → ICAResult
#
#  Helpers internos
#  ────────────────
#    _whiten_pca(X, k)   blanqueo PCA → (Z, V_whit, λ_k)
#    _sym_decorr(W)      decorrelación simétrica de W
#    _fastica_attempt(…) un intento FastICA → (W, n_iter, err, ok)
#
# ───────────────────────────────────────────────────────────────
#  Fichero    src/ica/ICACore.jl
#  Autor      Rafael Castro Triguero <me1catrr@uco.es>
#  Modificado 22-07-2026
# ───────────────────────────────────────────────────────────────

# ─── Blanqueo PCA ─────────────────────────────────────────────
"""
    _whiten_pca(X, k) -> (Z, V_whit, λ_k)

Blanqueo PCA reteniendo los `k` autovalores mayores.
Devuelve también `λ_k` para métricas (condición, varianza).
"""
function _whiten_pca(X::AbstractMatrix{<:Real}, k::Int)
    n_ch, n_samp = size(X)
    CovX = (X * X') / n_samp
    F    = eigen(Symmetric(CovX))
    idx  = sortperm(F.values; rev=true)
    λ_all = F.values[idx]
    λ_k  = λ_all[1:k]
    E_k  = F.vectors[:, idx][:, 1:k]
    V    = Diagonal(1 ./ sqrt.(λ_k)) * E_k'   # (k × n_ch)
    return V * X, V, λ_k, λ_all               # Z, V_whit, λ_k, λ_all
end

# ─── Decorrelación simétrica ──────────────────────────────────
function _sym_decorr(W::AbstractMatrix{<:Real})
    F   = eigen(Symmetric(W * W'))
    return (F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors') * W
end

# ─── Un intento FastICA ───────────────────────────────────────
"""
    _fastica_attempt(Z, n_comp, max_iter, tol, a, seed)
        -> (W, n_iter, final_error, converged)

Ejecuta FastICA simétrico (tanh) desde una semilla dada.
`final_error` = max|1 − |diag(W_new·W')|| en la última iteración.
"""
function _fastica_attempt(
    Z::AbstractMatrix{<:Real},
    n_comp::Int,
    max_iter::Int,
    tol::Float64,
    a::Float64,
    seed::Int,
)
    Random.seed!(seed)
    W = _sym_decorr(randn(n_comp, n_comp))
    n_samp = size(Z, 2)
    final_err = Inf
    n_iter = 0
    converged = false

    for iter in 1:max_iter
        n_iter = iter
        Y     = W * Z
        GY    = tanh.(a .* Y)
        Gp    = a .* (1 .- GY .^ 2)
        D     = Diagonal(vec(mean(Gp, dims=2)))
        W_new = _sym_decorr((GY * Z') / n_samp - D * W)
        M     = W_new * W'
        final_err = maximum(abs.(1 .- abs.(diag(M))))
        W = W_new
        if final_err < tol
            converged = true
            break
        end
    end
    return W, n_iter, final_err, converged
end

# ─── run_ica ──────────────────────────────────────────────────
"""
    run_ica(rec::EEGRecording, cfg::PipelineConfig) -> ICAResult

FastICA simétrico sobre la señal continua filtrada.
Usa sólo LinearAlgebra y Random; no requiere MultivariateStats.

Perfil `"eeg_julia"` (reproducibilidad estricta con EEG_Julia/src/ICA/ICA.jl):
  - n_components = n_channels (todos)
  - max_iter = 512, tol = 1e-7, seed = 1234, a = 1.0
  - max_attempts = 1 (sin reinicios; no altera M05)
  - A = inv(W_total)

Perfil `"default"` (modo reducido, configurable):
  - n_components / max_iter / tol / seed / max_attempts desde config
  - A = pinv(W_total)

Si un intento no alcanza `tol` y `max_attempts > 1`, se reintenta
con `seed + intento − 1` hasta converger o agotar intentos.
"""
function run_ica(rec::EEGRecording, cfg::PipelineConfig)::ICAResult
    ica_cfg = cfg.ica
    profile = String(get(ica_cfg, "profile", "default"))
    n_ch    = rec.meta.n_channels
    t0      = time()

    # Parámetros según perfil
    if profile == "eeg_julia"
        n_comp       = n_ch
        max_iter     = 512
        tol          = 1e-7
        seed0        = 1234
        max_attempts = 1   # reproducibilidad estricta — sin reinicios
    else
        n_comp_raw = get(ica_cfg, "n_components", 30)
        n_comp = if n_comp_raw == 0 || n_comp_raw == "auto"
            n_ch
        else
            min(Int(n_comp_raw), n_ch)
        end
        max_iter     = Int(get(ica_cfg, "max_iter", 500))
        tol          = Float64(get(ica_cfg, "tol",  1e-5))
        seed0        = Int(get(ica_cfg, "seed",  42))
        max_attempts = max(1, Int(get(ica_cfg, "max_attempts", 1)))
    end
    a = 1.0   # parámetro tanh supergaussiano (fijo)

    # 1) Centrado por canal
    Xc = copy(rec.data)
    for i in 1:size(Xc, 1)
        Xc[i, :] .-= mean(Xc[i, :])
    end

    # 2) Blanqueo PCA
    Z, V_whit, λ_k, λ_all = _whiten_pca(Xc, n_comp)
    λ_pos = filter(>(0), λ_all)
    λ_k_pos = filter(>(0), λ_k)
    whitening_cond = isempty(λ_k_pos) ? NaN :
        maximum(λ_k_pos) / max(minimum(λ_k_pos), eps(Float64))
    pca_var_pct = isempty(λ_pos) ? NaN :
        100.0 * sum(λ_k) / sum(λ_pos)

    # 3) Intentos FastICA (reinicio con seed+k si no converge)
    attempts = Dict{String,Any}[]
    W = Matrix{Float64}(undef, 0, 0)
    n_iter = 0
    final_err = Inf
    converged = false
    seed_used = seed0

    for attempt in 1:max_attempts
        seed_k = seed0 + (attempt - 1)
        W_k, n_iter_k, err_k, ok_k = _fastica_attempt(Z, n_comp, max_iter, tol, a, seed_k)
        push!(attempts, Dict{String,Any}(
            "attempt"     => attempt,
            "seed"        => seed_k,
            "n_iter"      => n_iter_k,
            "final_error" => err_k,
            "converged"   => ok_k,
        ))
        W = W_k
        n_iter = n_iter_k
        final_err = err_k
        converged = ok_k
        seed_used = seed_k
        ok_k && break
    end

    # 4) Matrices de salida
    S       = W * Z
    W_total = W * V_whit
    A = (n_comp == n_ch) ? inv(W_total) : pinv(W_total)

    # Ortogonalidad de W en espacio blanqueado: ||W W' − I||_F / n
    WW = W * W'
    I_n = Matrix{Float64}(I, n_comp, n_comp)
    orthogonality = sqrt(sum(abs2, WW .- I_n)) / n_comp

    # Varianza explicada por componente
    total_var = sum(var(rec.data, dims=2))
    var_exp   = [sum(A[:, i] .^ 2) * var(S[i, :]) for i in 1:n_comp] ./ total_var

    optimize_s = round(time() - t0, digits=3)
    n_attempts = length(attempts)

    diagnostics = Dict{String,Any}(
        "converged"       => converged,
        "n_iter"          => n_iter,
        "max_iter"        => max_iter,
        "final_error"     => final_err,
        "tol"             => tol,
        "seed"            => seed_used,
        "seed0"           => seed0,
        "n_attempts"      => n_attempts,
        "max_attempts"    => max_attempts,
        "restart_needed"  => n_attempts > 1,
        "optimize_s"      => optimize_s,
        "whitening_cond"  => whitening_cond,
        "pca_rank"        => n_comp,
        "pca_rank_full"   => n_ch,
        "pca_var_pct"     => pca_var_pct,
        "orthogonality"   => orthogonality,
        "profile"         => profile,
        "attempts"        => attempts,
    )

    return ICAResult(rec.meta, A, W_total, S, Int[], var_exp, diagnostics)
end
