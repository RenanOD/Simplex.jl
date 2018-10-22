# run 'runtests.jl' before benchmarking

using Base.Test, BenchmarkTools

@testset "Testing benchmark problems" begin
  open("benchmark.log", "w") do logfile
    for (m, n) in [(90,10) (450,50)]
      for t = 1:2
        # construct random transport problem given m and n
        A, b, c = transport_instance(m, n)
        model = Model(solver = CbcSolver(logLevel=0))
        @variable(model, x[1:m*n] >= 0)
        @objective(model, Min, dot(x, c))
        @constraint(model, A * x .== b)
        status = solve(model)
        xj = getvalue(x)
        zj = getobjectivevalue(model)

        # benchmark
        #benchmarkCbc  = @benchmark status = solve($model)             samples=10 evals=1
        benchmarkluup = @benchmark simplexluup($c, $A, $b)             samples=10 evals=1
        benchmarklu   = @benchmark simplexluup($c, $A, $b, maxups=0)   samples=10 evals=1

        println(logfile, "transport problem: $(m+n) rows, $(m*n) cols, $(nnz(A)) NNZ's luup/Clp: ")
        println(logfile, "  luup          = $(median(benchmarkluup))")
        #println(logfile, "  Cbc          = $(median(benchmarkCbc))")
        println(logfile, "  lu no update  = $(median(benchmarklu))")
        #println(logfile, "  ", judge(median(benchmarkluup), median(benchmarkCbc)))
        println(logfile, "  ", judge(median(benchmarkluup), median(benchmarklu)))

        # test if simplexluup successfuly solved the benchmark problems
        x, z, status = simplexluup(c, A, b, max_iter = 100*m*n)
        @test dot(c, x) ≈ zj atol = 1e-3
      end
    end
  end
end
