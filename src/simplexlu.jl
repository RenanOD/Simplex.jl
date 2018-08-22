export simplexlu

function simplexlu(c, A, b, IB=0, L=0, U=0, prow=0, Rs=0; max_iter = 4000) #no need to pass xB
  m, n = size(A)
  A = sparse(A)*1.0 #only needed if IB == 0
  iter = 0
  if IB == 0 # construct artificial problem
    artificial = true
    signb = sign.(b)
    A = [A spdiagm(signb)]
    IB = collect(n+1:n+m) # indexes of basic variables
    IN = collect(1:n)
    co = collect(c)
    c = [zeros(n); ones(m)]
    U = A[:,IB]
    r = -A[:,IN]'*signb # artificial relative costs
    xB = collect(abs.(b)) # solution in current basis
  else
    artificial = false
    IN = setdiff(1:n, IB)
    if L == 0 #only if user gives IB
      F = lufact(A[:,IB])
      xB = F\b
      L, U, prow, pcol, Rs = F[:(:)] # (Rs.*A)[prow,pcol] * x[pcol] = b[prow]
      IB, xB = IB[pcol], xB[pcol]
      r = c[IN] - A[:,IN]'*((speye(m)[:,prow]*(L'\(U'\(c[IB])))).*Rs)
    else
      r = c[IN] - A[:,IN]'*((speye(m)[:,prow]*(L'\(U'\(c[IB])))).*Rs)
      xB = (U\(L\((b.*Rs)[prow])))
    end
  end
  q = findfirst(r .< -1e-12) # Bland's Rule

  status = :Optimal

  while !(q == 0 || iter >= max_iter)
    iter += 1
    @assert all(xB .>= 0)
    d = (L==0)? U\A[:,IN[q]] : U\(L\((A[:,IN[q]].*Rs)[prow]))
    apfrac = xB ./ d # relative variable changes to direction
    indpos = find(d .> 1e-12) # variables that decrease in d direction
    if length(indpos) == 0
      status = :Unbounded
      break
    end
    indxq = indmin(apfrac[indpos])
    xq = apfrac[indpos[indxq]]
    @assert xq >= 0
    @assert xq < Inf
    #update basis
    p = findfirst(apfrac, xq) # Bland's Rule
    xB -= xq * d; xB[p] = xq # update solution
    IB[p], IN[q] = IN[q], IB[p] # update indexes
    #reset LU
    F = lufact(A[:,IB])
    L, U, prow, pcol, Rs = F[:(:)] # (Rs.*A)[prow,pcol] * x[pcol] = b[prow]
    IB, xB = IB[pcol], xB[pcol]

    if L == 0
      r = c[IN] - A[:,IN]' * U'\((c[IB]))
    else
      r = c[IN] - A[:,IN]'*((speye(m)[:,prow]*(L'\(U'\((c[IB]))))).*Rs)
    end
    q = findfirst(r .< -1e-12) # Bland's Rule
  end

  if iter >= max_iter
    status = :UserLimit
  end

  x = zeros(n)
  if !artificial
    x[IB] = xB
    z = dot(c, x)
  else
    if dot(xB, c[IB])/norm(xB) > 1e-12 #use epsilon, not zero
      status = :Infeasible
      I = find(IB .<= n - m)
      x[I] = xB[I]
      z = dot(co, x)
    elseif maximum(IB) > n # check for artificial variables in basis
      deleteat!(IN, find(IN .> n))
      Irows = collect(1:m)
      p = findfirst(IB .> n)

      while p != 0
        q = 1
        while q <= length(IN) #searching for columns to substitute artificials in basis
          di = (L==0)? U\A[Irows,IN[q]] : U\(L\((A[Irows,IN[q]].*Rs)[prow]))
          dipivot = (L==0)? di[findfirst((A[Irows,IB[p]]))] : di[findfirst((A[Irows,IB[p]][prow]))]
          if abs(dipivot) > 1e-12
            break
          end
          q += 1
        end
        if q > length(IN)
          deleteat!(Irows, findfirst(A[Irows,IB[p]]))
          deleteat!(IB, p)
          F = lufact(A[Irows,IB])
          L, U, prow, pcol, Rs = F[:(:)]
          IB = IB[pcol]
        else
          IB[p] = IN[q]
          deleteat!(IN, q)
          F = lufact(A[Irows,IB])
          L, U, prow, pcol, Rs = F[:(:)]
          IB = IB[pcol]
        end
        p = findfirst(IB .> n)
      end
      @assert length(Irows) > 0
      x, z, status = simplexlu(co, A[Irows,1:n], b[Irows], IB, L, U, prow, Rs)
    else
      x, z, status = simplexlu(co, A[:,1:n], b, IB, L, U, prow, Rs)
    end
  end

  return x, z, status
end

