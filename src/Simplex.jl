module Simplex

using Base.Test, JuMP, Cbc, Base.SparseArrays.halfperm!
include("simplexluup.jl")
include("simplexinv.jl")

end
