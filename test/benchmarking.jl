# run 'runtests.jl' before benchmarking

using Base.Test, BenchmarkTools

  @testset "Testing benchmark problems" begin
    for (m,n) in [(90,10) (450,50)]

      CbcSolver(print_level=0)

      # construct random transport problem given m and n
      A, b, c = transport_instance(m, n)
      model = Model(solver = CbcSolver())
      @variable(model, x[1:m*n] >= 0)
      @objective(model, Min, dot(x, c))
      @constraint(model, A * x .== b)
      status = solve(model)
      xj = getvalue(x)
      zj = getobjectivevalue(model)

      # benchmark
      benchmarkCbc = @benchmark status = solve(model)
      benchmarkluup = @benchmark simplexluup(c, A, b)

      print("transport problem: $(m+n) rows, $(m*n) cols, $(nnz(A)) NNZ's luup/Clp: ")
      println(judge(median(benchmarkluup), median(benchmarkCbc)))

      # test if simplexluup successfuly solved the benchmark problems
      x, z, status = simplexluup(c, A, b, max_iter = 100*m*n)
      @test dot(c, x) ≈ zj atol = 1e-3

    end
  end
