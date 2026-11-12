# ═══════════════════════════════════════════════════════════════
#  Experimento de validación del perfil de ICA «estable»: ocho mezclas
# ═══════════════════════════════════════════════════════════════
#
#  NO forma parte de runtests.jl: sus resultados dependen de la muestra y
#  del redondeo, y no se pueden fijar como pruebas portables entre entornos.
#  Lo reproducible es el procedimiento, que se especifica aquí completo.
#  Solo señales generadas aquí. No lee ninguna grabación.
#
#    julia --project=. test/validacion_ica_estable.jl [salida.md]
#
#  ── Especificación ──────────────────────────────────────────────
#
#  Datos
#    · n = 50 000 muestras, fs = 500 Hz.
#    · Seis fuentes no gaussianas, generadas con MersenneTwister(20261112) en
#      este orden: (1) supergaussiana, randn·signo(randn)·|randn|; (2)
#      subgaussiana, uniforme en [−0,5, 0,5]; (3) alfa de 10 Hz modulada al
#      0,2 Hz con profundidad 0,5; (4) picos de amplitud 8 con probabilidad
#      0,002 más ruido gaussiano de σ = 0,1; (5) diente de sierra de 3 Hz;
#      (6) seno de 50 Hz. Cada una, centrada y con varianza 1.
#    · Mezcla m ∈ 2…9: con MersenneTwister(m), primero dos fuentes gaussianas
#      (randn, 2 × n) y después una matriz de mezcla 8 × 8 (randn). X = A·S.
#
#  ICA
#    · Perfil `estable` con sus valores por defecto: stable_seed = 1234,
#      stable_max_iter = 2048, stable_tol = 1e-7, stable_runs = 5 (semillas
#      1235…1239), contraste tanh, n_comp = 8.
#
#  Métricas
#    · Recuperación de una fuente: |r de Pearson| entre la fuente verdadera y
#      el componente que le asigna la asignación global uno a uno (húngaro)
#      sobre la matriz |r| fuentes × componentes de la referencia.
#    · Repetibilidad de un componente: |r| con su pareja en cada ejecución
#      adicional (asignación global, como en el perfil). Se dan dos mínimos:
#      sobre las ejecuciones convergidas (el que usa el perfil) y sobre todas
#      las que produjeron una descomposición (exploratorio).
#    · Subespacio gaussiano: cosenos de los ángulos principales entre el
#      subespacio que generan las dos fuentes gaussianas y el que generan
#      sus dos componentes asignados, como series temporales: valores
#      singulares de Q₁ᵀQ₂, con Q₁ y Q₂ bases ortonormales (QR) de cada par.
#      1 significa el mismo plano.
#
#  Lectura
#    · Las cifras son descriptivas: no se compara ninguna con un umbral.
#      Sirven para delimitar conclusiones, no para validar componentes.
#    · Cambian con el entorno numérico (versión de Julia, biblioteca de
#      álgebra lineal, hilos): la salida registra el entorno y el commit.
# ───────────────────────────────────────────────────────────────

using NeuroMIND, Statistics, Random, LinearAlgebra, Printf, DataFrames, Dates

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
fmt(x) = ismissing(x) ? "—" : @sprintf("%.3f", x)

cfg = load_ss_config(joinpath(@__DIR__, "..", "config", "pipeline.toml")); cfg.ica["profile"] = "estable"
n  = 50_000
S6 = fuentes(MersenneTwister(20261112), n)
commit = try strip(read(`git -C $(joinpath(@__DIR__, "..")) rev-parse --short HEAD`, String)) catch; "desconocido" end

lineas = String[
    "# Validación del perfil «estable» · ocho mezclas",
    "",
    "$(Dates.format(now(), "dd-mm-yyyy HH:MM")) · NeuroMIND $(commit) · Julia $(VERSION) · $(Sys.MACHINE) · BLAS $(BLAS.get_config().loaded_libs[1].libname |> basename) con $(BLAS.get_num_threads()) hilos",
    "",
    "| Mezcla | Referencia | Ejecuciones convergidas · fallidas | Recuperación mínima de las no gaussianas | Recuperación de las gaussianas | Repetibilidad de las gaussianas, convergidas | Repetibilidad de las gaussianas, todas | Cosenos principales del subespacio gaussiano |",
    "|---|---|---|---|---|---|---|---|",
]
for m in 2:9
    rng = MersenneTwister(m)
    S = vcat(S6, randn(rng, 2, n)); rec = grabacion(randn(rng, 8, 8) * S)
    ref = try run_ica_stable(rec, cfg) catch e
        push!(lineas, "| $(m) | falla: $(typeof(e)) | | | | | | |"); continue
    end
    st = ica_stability(rec, cfg, ref)
    C = abs.(cor(S', Float64.(ref.activations)')); a = NeuroMIND._ica_hungarian(C); g = a[7:8]
    Q1 = Matrix(qr(S[7:8, :]').Q)[:, 1:2]; Q2 = Matrix(qr(Float64.(ref.activations[g, :])').Q)[:, 1:2]
    cs = svdvals(Q1' * Q2)
    push!(lineas, @sprintf("| %d | %s | %d · %d | %.4f | %.3f / %.3f | %s | %s | %.4f / %.4f |", m,
        ref.diagnostics["converged"] ? "converge" : "no converge", st.runs_converged, st.runs_failed,
        minimum(C[i, a[i]] for i in 1:6), C[7, a[7]], C[8, a[8]],
        join(fmt.(st.components.r_min[g]), " / "), join(fmt.(st.components.r_min_all_runs[g]), " / "), cs...))
end
salida = join(lineas, "\n") * "\n"
print(salida)
isempty(ARGS) || write(ARGS[1], salida)
