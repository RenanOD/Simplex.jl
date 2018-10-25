# run 'runtests.jl' before benchmarking

using Base.Test, BenchmarkTools

@testset "Testing benchmark problems" begin
  open("benchmark.log", "w") do logfile
    for (m, n) in [(90,10) (450,50) (900,100)]
      for t = 1:1
        # construct random transport problem given m and n
        A, b, c = transport_instance(m, n)
        model = Model(solver = GLPKSolverLP(presolve=false))
        @variable(model, x[1:m*n] >= 0)
        @objective(model, Min, dot(x, c))
        @constraint(model, A * x .== b)

        #benchmarking GLPK
        benchmarkGLPK  = @benchmark status = solve($model)             samples=10 evals=1
        xj = getvalue(x)
        zj = getobjectivevalue(model)

        #benchmarking simplexluup
        x = Array{Float64,1}(m*n)
        benchmarkluup = @benchmark $x .= simplexluup($c, $A, $b, max_iter = 100*$m*$n)[1]             samples=10 evals=1
        #benchmarklu   = @benchmark simplexluup($c, $A, $b, maxups=0)   samples=10 evals=1

        println(logfile, "transport problem: $(m+n) rows, $(m*n) cols, $(nnz(A)) NNZ's luup/GLPK: ")
        println(logfile, "  luup          = $(median(benchmarkluup))")
        println(logfile, "  GLPK          = $(median(benchmarkGLPK))")
        #println(logfile, "  lu no update  = $(median(benchmarklu))")
        println(logfile, "  ", judge(median(benchmarkluup), median(benchmarkGLPK)))
        #println(logfile, "  ", judge(median(benchmarkluup), median(benchmarklu)))

        # test if simplexluup successfuly solved the benchmark problems
        @test dot(c, x) ≈ zj atol = zj*1e-9
      end
    end
  end
end
