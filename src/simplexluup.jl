export simplexluup

function simplexluup(c, A, b, IB=0, xB=0, L=0, U=0, MP=[], P=[], prow=0, Rs=0; max_iter = 4000) #no need to pass xB
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
    if (MP == [] && L == 0) #only if user gives IB
      F = lufact(A[:,IB])
      xB = F\b
      L, U, prow, pcol, Rs = F[:(:)] # (Rs.*A)[prow,pcol] * x[pcol] = b[prow]
      IB, xB = IB[pcol], xB[pcol]
      r = c[IN] - A[:,IN]'*((speye(m)[:,prow]*(L'\(U'\(c[IB])))).*Rs)
    elseif MP==[] #there is L, no MP
      r = c[IN] - A[:,IN]'*((speye(m)[:,prow]*(L'\(U'\(c[IB])))).*Rs)
    else #there is MP, maybe L
      lambda = U'\((c[IB]))
      for j in 1:length(MP)
        lambda -= lambda[end]*MP[end-j+1]
        lambda[P[end-j+1]] = lambda
      end
      r = (L==0)? c[IN] - A[:,IN]'*lambda : c[IN] - A[:,IN]'*((speye(m)[:,prow]*(L'\lambda)).*Rs)
    end
  end
  q = findfirst(r .< -1e-12) # Bland's Rule

  status = :Optimal

  while !(q == 0 || iter >= max_iter)
    iter += 1
    @assert all(xB .>= 0)
    w = (L==0)? A[:,IN[q]] : L\((A[:,IN[q]].*Rs)[prow])
    for j in 1:length(MP)
      w = w[P[j]]
      w[end] -= dot(MP[j], w)
    end
    d = (U\w)
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

    p = findfirst(apfrac, xq) # Bland's Rule
    xB -= xq * d; xB[p] = xq # update solution
    IB[p], IN[q] = IN[q], IB[p] # update indexes
    if length(MP) >= 25 #reset LU
      F = lufact(A[:,IB])
      L, U, prow, pcol, Rs = F[:(:)] # (Rs.*A)[prow,pcol] * x[pcol] = b[prow]
      IB, xB = IB[pcol], xB[pcol]
      MP, P = [], []
    else #update LU
      U[:,p] = w
      push!(P, reverse(reverse(1:m,p,m),p,m-1))
      U = U[:,P[end]]
      mp = spzeros(m)
      for i in 1:m-p
        mp[p+i-1] = U[p,p+i-1]/U[p+i,p+i-1]
        U[p,i:m] -= mp[p+i-1]*U[p+i,i:m]
      end
      push!(MP, mp)
      U = U[P[end],:]
      IB, xB = IB[P[end]], xB[P[end]]
    end

    lambda = U'\((c[IB]))
    for j in 1:length(MP)
      lambda -= lambda[end]*MP[end-j+1]
      lambda[P[end-j+1]] = lambda
    end
    if L == 0
      r = c[IN] - A[:,IN]'*lambda
    else
      r = c[IN] - A[:,IN]'*((speye(m)[:,prow]*(L'\lambda)).*Rs)
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
    if dot(xB, c[IB])/norm(xB) > 1e-12 #use eps, not zero
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
        Ap = (prow==0)? A[Irows,IB[p]] : A[Irows,IB[p]][prow]
        for pj in P
          Ap = Ap[pj]
        end
        PivotAp = findfirst(Ap .> 0)
        while q <= length(IN) #searching for columns to substitute artificials in basis
          w = (L==0)? A[Irows,IN[q]] : L\((A[Irows,IN[q]].*Rs)[prow])
          for j in 1:length(MP)
            w = w[P[j]]
            w[end] -= dot(MP[j], w)
          end
          if abs((U\w)[PivotAp]) > 1e-12
            break
          end
          q += 1
        end
        if q > length(IN)
          deleteat!(Irows, findfirst(A[Irows,IB[p]]))
          deleteat!(IB, p)
          deleteat!(xB, p)
          F = lufact(A[Irows,IB])
          L, U, prow, pcol, Rs = F[:(:)]
          IB, xB = IB[pcol], xB[pcol]
          MP, P = [], []
        else
          IB[p] = IN[q]
          deleteat!(IN, q)
          #update LU
          U[:,p] = w
          push!(P, reverse(reverse(1:m,p,m),p,m-1))
          U = U[:,P[end]]
          mp = spzeros(m)
          for i in 1:m-p
            mp[p+i-1] = U[p,p+i-1]/U[p+i,p+i-1]
            U[p,:] -= mp[p+i-1]*U[p+i,:]
          end
          push!(MP, mp)
          U = U[P[end],:]
          IB, xB = IB[P[end]], xB[P[end]]
        end
        p = findfirst(IB .> n)
      end
      @assert length(Irows) > 0
      x, z, status = simplexluup(co, A[Irows,1:n], b[Irows], IB, xB, L, U, MP, P, prow, Rs)
    else
      x, z, status = simplexluup(co, A[:,1:n], b, IB, xB, L, U, MP, P, prow, Rs)
    end
  end

  return x, z, status
end

