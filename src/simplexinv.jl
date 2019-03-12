export simplexinv

function simplexinv(c, A, b, IB=0, invB=0; max_iter = 20000)
  m, n = size(A)
  iter = 0
  if IB == 0 # construct artificial problem
    artificial = true
    signb = sign.(b)
    A = [A sparse(Diagonal(signb))]
    IB = collect(n+1:n+m) # indexes of basic variables
    IN = collect(1:n)
    co = collect(c)
    c = [zeros(n); ones(m)]
    invB = sparse(Diagonal(signb))
    r = -A[:,IN]'*signb # artificial relative costs
    xB = collect(abs.(b)) # solution in current basis
  else
    artificial = false
    IN = setdiff(1:n, IB)
    r = c[IN] - A[:,IN]' * (invB' * c[IB])
    xB = invB*b
  end
  q = findfirst(r .< 0) # Bland's Rule

  status = :Optimal

  while !(q == nothing || iter >= max_iter) # relative variable changes to directioner >= max_iter)
    iter += 1
    d = invB * A[:,IN[q]] # viable direction

    xq = Inf
    for k in 1:m # find min xB/d s.t. d .> 0
      if d[k] >= 2e-16
        dfrac = xB[k]/d[k]
        if dfrac < xq
          xq = dfrac
          p = k
        end
      end
    end
    if xq == Inf
      status = :Unbounded; break
    end

    xB -= xq * d; xB[p] = xq # update solution
    IB[p], IN[q] = IN[q], IB[p] # update indexes
    #update of inverse of B
    E = one(zeros(m,m))
    dp = d[p]
    d[p] = -1
    E[:, p] = -d / dp
    invB = E*invB #no need to create E
    r = c[IN] - A[:,IN]' * (invB' * c[IB])
    q = findfirst(r .< 0) # Bland's Rule
  end

  if iter >= max_iter
    status = :UserLimit
  end

  x = zeros(n)
  if !artificial
    x[IB] = xB
    z = dot(c, x)
  else
    if dot(xB, c[IB]) > 0
      status = :Infeasible
      I = findall(IB .<= n - m)
      x[I] = xB[I]
      z = dot(co, x)
    elseif maximum(IB) > n # check for artificial variables in basis
      deleteat!(IN, findall(IN .> n))
      Irows = collect(1:m)
      p = findfirst(IB .> n)
      while p != nothing
        q = findfirst(invB[p,:]' * A[Irows,IN] .!= 0)
        if q == nothing
          deleteat!(Irows, p)
          deleteat!(IB, p)
          invB = inv(Matrix(A[Irows,IB]))
        else
          IB[p] = IN[q]
          deleteat!(IN, q)
          d = invB * A[Irows,IN[q]]
          E = one(zeros(m,m))
          dp = d[p]
          d[p] = -1
          E[:, p] = -d / dp
          invB = Q*invB
        end
        p = findfirst(IB .> n)
      end
      x, z, status = simplexinv(co, A[Irows,1:n], b[Irows], IB, invB)
    else
      x, z, status = simplexinv(co, A[:,1:n], b, IB, invB)
    end
  end

  return x, z, status
end
