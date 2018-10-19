using Base.Test

function test_mpstostd()

	problems = ["basic.mps", "farm.mps","kleemin3.mps","kleemin4.mps","kleemin5.mps","kleemin6.mps","kleemin7.mps","kleemin8.mps","afiro.mps","refine.mps","sc50a.mps","sc50b.mps", "share2b.mps","scagr7.mps","p0282.mps","p0548.mps","p0291.mps","agg2.mps","agg3.mps","nsic1.mps"]
	cd("mps")

	@testset "Many tests" begin
	  for i in 1:length(problems)
	    mpsfile = problems[i]
	    c, A, b, status, zchange = mpstostd(mpsfile)
	    m, n = size(A)

	    mod = Model(solver=ClpSolver())
	    internal_model = MathProgBase.LinearQuadraticModel(ClpSolver())
        MathProgBase.loadproblem!(internal_model, mpsfile)
        status = solvelp(internal_model).status
        zorig = solvelp(internal_model).objval

        mo = Model(solver=ClpSolver())
        @variable(mo, x[1:length(c)] >= 0)
        @objective(mo, Min, dot(c,x))
        @constraint(mo, full(A)*x .== b)
        statusstd = solve(mo)
        zstd = getobjectivevalue(mo)

        @test zstd ≈ (zorig - zchange)
      end
    end
end

test_mpstostd()
cd("../")
