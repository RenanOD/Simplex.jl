using Base.Test

function test_simplex_mps()

  problems = ["basic.mps", "farm.mps","kleemin3.mps","kleemin4.mps","kleemin5.mps","kleemin6.mps","kleemin7.mps","kleemin8.mps","afiro.mps","refine.mps","sc50a.mps","sc50b.mps", "share2b.mps","scagr7.mps","p0282.mps","p0548.mps","p0291.mps","agg2.mps","agg3.mps","nsic1.mps"]
  @testset "Simplex mps" begin
    for i in 1:length(problems)
	    mpsfile = problems[i]
	    c, A, b, status, zchange = mpstostd(mpsfile)
	    m, n = size(A)

	    mod = Model(solver=GLPKSolverLP())
	    internal_model = MathProgBase.LinearQuadraticModel(GLPKSolverLP())
        MathProgBase.loadproblem!(internal_model, mpsfile)
        status = solvelp(internal_model).status
        zorig = solvelp(internal_model).objval

        zsimp = simplex(c, A, b)[2]

        @test zsimp ≈ (zorig - zchange)
      end
    end
end

test_mpstostd()
