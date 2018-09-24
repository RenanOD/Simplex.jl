Pkg.add("JuMP")
Pkg.add("Cbc")

using Simplex, Base.Test, JuMP, Cbc

include("test_simplexluup.jl")
include("test_simplexinv.jl")
