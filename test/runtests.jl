# test/runtests.jl
# Wrapper estándar para Pkg.test() (Julia convention: test/runtests.jl).
# La suite real vive en tests/runtests.jl para preservar el flujo existente de
# `julia --project=. tests/runtests.jl`.
#
# Pkg.test() en Julia 1.12 crea un env aislado con LOAD_PATH = ["@", tmpdir],
# omitiendo "@stdlib". Lo restauramos para que `using Test/Statistics/LinearAlgebra`
# funcionen correctamente dentro del runner.
"@stdlib" in Base.LOAD_PATH || push!(Base.LOAD_PATH, "@stdlib")

include(joinpath(@__DIR__, "..", "tests", "runtests.jl"))
