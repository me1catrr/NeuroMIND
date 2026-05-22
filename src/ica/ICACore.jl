# NeuroMIND/src/ica/ICACore.jl
# FastICA simétrico con PCA whitening.
# Implementación propia en Julia puro (LinearAlgebra + Random).
# Portada desde EEG_Julia/Pluto/ICA/ICA.jl.

# ─── Blanqueo PCA ─────────────────────────────────────────────
function _whiten_pca(X::AbstractMatrix{<:Real}, k::Int)
    n_ch, n_samp = size(X)
    CovX = (X * X') / n_samp
    F    = eigen(Symmetric(CovX))
    idx  = sortperm(F.values; rev=true)
    λ_k  = F.values[idx][1:k]
    E_k  = F.vectors[:, idx][:, 1:k]
    V    = Diagonal(1 ./ sqrt.(λ_k)) * E_k'   # (k × n_ch)
    return V * X, V                             # Z, V_whit
end

# ─── Decorrelación simétrica ──────────────────────────────────
function _sym_decorr(W::AbstractMatrix{<:Real})
    F   = eigen(Symmetric(W * W'))
    return (F.vectors * Diagonal(1 ./ sqrt.(F.values)) * F.vectors') * W
end

# ─── run_ica ──────────────────────────────────────────────────
"""
    run_ica(rec::EEGRecording, cfg::PipelineConfig) -> ICAResult

FastICA simétrico sobre la señal continua filtrada.
Usa sólo LinearAlgebra y Random; no requiere MultivariateStats.

Perfil `"eeg_julia"` (reproducibilidad estricta con EEG_Julia/src/ICA/ICA.jl):
  - n_components = n_channels (todos)
  - max_iter = 512, tol = 1e-7, seed = 1234, a = 1.0
  - A = inv(W_total)  (exactamente como EEG_Julia)

Perfil `"default"` (modo reducido, configurable):
  - n_components: de config (0 / "auto" → todos los canales)
  - max_iter, tol, seed: de config
  - A = pinv(W_total)
"""
function run_ica(rec::EEGRecording, cfg::PipelineConfig)::ICAResult
    ica_cfg = cfg.ica
    profile = String(get(ica_cfg, "profile", "default"))
    n_ch    = rec.meta.n_channels

    # Parámetros según perfil
    if profile == "eeg_julia"
        n_comp   = n_ch     # todos los canales, igual que EEG_Julia
        max_iter = 512
        tol      = 1e-7
        seed     = 1234
    else
        n_comp_raw = get(ica_cfg, "n_components", 30)
        n_comp = if n_comp_raw == 0 || n_comp_raw == "auto"
            n_ch
        else
            min(Int(n_comp_raw), n_ch)
        end
        max_iter = Int(get(ica_cfg, "max_iter", 500))
        tol      = Float64(get(ica_cfg, "tol",  1e-5))
        seed     = Int(get(ica_cfg, "seed",  42))
    end
    a = 1.0   # parámetro tanh supergaussiano (fijo)

    # 1) Centrado por canal
    Xc = copy(rec.data)
    for i in 1:size(Xc, 1)
        Xc[i, :] .-= mean(Xc[i, :])
    end

    # 2) Blanqueo PCA
    Z, V_whit = _whiten_pca(Xc, n_comp)   # Z: (n_comp × n_samp)

    # 3) Inicialización de W
    Random.seed!(seed)
    W = _sym_decorr(randn(n_comp, n_comp))

    # 4) Iteraciones FastICA
    n_samp = size(Z, 2)
    for iter in 1:max_iter
        Y    = W * Z                                # (n_comp × n_samp)
        GY   = tanh.(a .* Y)
        Gp   = a .* (1 .- GY .^ 2)
        D    = Diagonal(vec(mean(Gp, dims=2)))
        W_new = _sym_decorr((GY * Z') / n_samp - D * W)
        M    = W_new * W'
        if maximum(abs.(1 .- abs.(diag(M)))) < tol
            W = W_new; break
        end
        W = W_new
    end

    # 5) Matrices de salida
    S       = W * Z                          # activaciones  (n_comp × n_samp)
    W_total = W * V_whit                     # unmixing total (n_comp × n_ch)
    # inv() cuando n_comp == n_ch (perfil eeg_julia, exactamente como EEG_Julia);
    # pinv() en modo reducido (n_comp < n_ch).
    A = (n_comp == n_ch) ? inv(W_total) : pinv(W_total)

    # Varianza explicada: contribución de cada componente a la varianza total de la señal
    # var_i = ||A[:,i]||² × var(S[i,:]) / Σ var(canales originales)
    total_var = sum(var(rec.data, dims=2))
    var_exp   = [sum(A[:, i] .^ 2) * var(S[i, :]) for i in 1:n_comp] ./ total_var

    return ICAResult(rec.meta, A, W_total, S, Int[], var_exp)
end
