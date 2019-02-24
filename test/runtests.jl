using Pkg

Pkg.add("JuMP")
Pkg.add("GLPK")

using Main.Simplex, Test, JuMP, GLPK, Random, LinearAlgebra, SparseArrays

  include("test_simplexluup.jl")
  #include("test_simplexinv.jl")
