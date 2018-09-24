using Simplex, Base.Test
Pkg.add("JuMP")
Pkg.add("Cbc")
include("test_simplexluup.jl")
include("test_simplexinv.jl")
