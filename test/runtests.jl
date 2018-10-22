Pkg.add("JuMP")
Pkg.add("GLPK")

using Simplex, Base.Test, JuMP, GLPKMathProgInterface

  include("test_simplexluup.jl")
  include("test_simplexinv.jl")
