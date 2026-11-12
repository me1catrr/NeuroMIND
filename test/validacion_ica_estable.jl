# ═══════════════════════════════════════════════════════════════
#  Validación del perfil experimental de ICA «estable»
# ═══════════════════════════════════════════════════════════════
#
#  NO forma parte de runtests.jl. Muestra, con fuentes sintéticas conocidas,
#  comportamientos que dependen de la muestra y del redondeo y que, por eso,
#  no se pueden fijar como pruebas portables entre entornos:
#
#    · mezclas con dos fuentes gaussianas: la referencia converge o no, y los
#      componentes que corresponden a las gaussianas se repiten o no entre
#      semillas, sin que repetirse signifique recuperar las fuentes;
#    · el subespacio que generan esos componentes, que sí coincide con el de
#      las fuentes gaussianas;
#    · ejecuciones que fallan dentro de ICACore (`_sym_decorr`).
#
#    julia --project=. test/validacion_ica_estable.jl
#
#  Solo señales generadas aquí. No lee ninguna grabación.
# ───────────────────────────────────────────────────────────────

using NeuroMIND, Statistics, Random, LinearAlgebra, Printf, DataFrames

function fuentes(rng, n; fs = 500.0)
    t = (0:n-1) ./ fs
    S = [(randn(rng, n) .* sign.(randn(rng, n))) .* abs.(randn(rng, n)), rand(rng, n) .- 0.5,
         sin.(2π * 10 .* t) .* (1 .+ 0.5 .* sin.(2π * 0.2 .* t)),
         [rand(rng) < 0.002 ? 8.0 : 0.0 for _ in 1:n] .+ 0.1 .* randn(rng, n),
         mod.(3 .* t, 1.0) .- 0.5, sin.(2π * 50 .* t)]
    return permutedims(hcat([(x .- mean(x)) ./ std(x) for x in S]...))
end
function grabacion(X; fs = 500.0)
    meta = NeuroMIND.RecordingMeta("S01", "T1", "EC", 1, fs, size(X, 1), ["X$(i)" for i in 1:size(X, 1)], nothing, "")
    return NeuroMIND.EEGRecording(meta, Matrix{Float64}(X), collect((0:size(X, 2)-1) ./ fs))
end
cfg = load_ss_config(joinpath(@__DIR__, "..", "config", "pipeline.toml")); cfg.ica["profile"] = "estable"

n  = 50_000
S6 = fuentes(MersenneTwister(20261112), n)
println("| Mezcla | Referencia | Ejecuciones convergidas · fallidas | Recuperación de las no gaussianas | Recuperación de las gaussianas | Repetibilidad de las gaussianas (mín. en todas las ejecuciones) | Cosenos principales del subespacio gaussiano |")
println("|---|---|---|---|---|---|---|")
for semilla in 2:9
    rng = MersenneTwister(semilla)
    S = vcat(S6, randn(rng, 2, n)); rec = grabacion(randn(rng, 8, 8) * S)
    ref = try run_ica_stable(rec, cfg) catch e
        println("| $(semilla) | falla: $(typeof(e)) | | | | | |"); continue
    end
    st = ica_stability(rec, cfg, ref)
    C = abs.(cor(S', Float64.(ref.activations)')); a = NeuroMIND._ica_hungarian(C); g = a[7:8]
    Q1 = Matrix(qr(S[7:8, :]').Q)[:, 1:2]; Q2 = Matrix(qr(Float64.(ref.activations[g, :])').Q)[:, 1:2]
    cs = svdvals(Q1' * Q2)
    rep = [ismissing(x) ? "—" : @sprintf("%.3f", x) for x in st.components.r_min_all_runs[g]]
    @printf("| %d | %s | %d · %d | %.4f | %.3f / %.3f | %s | %.4f / %.4f |\n", semilla,
        ref.diagnostics["converged"] ? "converge" : "no converge", st.runs_converged, st.runs_failed,
        minimum(C[i, a[i]] for i in 1:6), C[7, a[7]], C[8, a[8]], join(rep, " / "), cs...)
end
println("\nEntorno: Julia $(VERSION) · $(Sys.MACHINE) · BLAS $(BLAS.get_num_threads()) hilos")
