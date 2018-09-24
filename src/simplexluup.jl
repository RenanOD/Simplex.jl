export simplexluup

#export simplexluup

function simplexluup(c, A, b, 𝔹=0, L=0, U=0, prow=0, Rs=0, xB=0; max_iter = 10000, maxups = 15)
  # preparations
  m, n = size(A)
  iter = 0; ups = 0
  Ufiller = spzeros(0)
  P, MP = Vector{Vector{Int64}}(maxups), Vector{SparseVector{Float64,Int64}}(maxups)
  if 𝔹 == 0 # construct artificial problem
    artificial = true
    signb = sign.(b)
    A = [A spdiagm(signb)] # try to save A memory to return
    𝔹 = collect(n+1:n+m); ℕ = collect(1:n) # artificial indexes
    co = collect(c); c = [zeros(n); ones(m)] # try to save c to save memory
    U = A[:,𝔹]; AN = A[:,ℕ]; cN = @view c[ℕ]
    r = -(signb'*AN)' # artificial relative costs
    xB = abs.(b) # solution in current basis
  else
    artificial = false
    ℕ = setdiff(1:n, 𝔹)
    AN = A[:,ℕ]; cN = @view c[ℕ]
    if L == 0 # if user gives 𝔹
      F = lufact(A[:,𝔹])
      xB = F\b
      L, U, prow, pcol, Rs = F[:(:)] # (Rs.*A)[prow,pcol] * x[pcol] = b[prow]
      𝔹, xB = 𝔹[pcol], xB[pcol]
      r = cN - ((ipermute!(L'\(U'\c[𝔹]),prow).*Rs)'*AN)'
    else
      r = cN - ((ipermute!(L'\(U'\c[𝔹]),prow).*Rs)'*AN)'
    end
  end
  q = findfirst(r .< -1e-12) # Bland's Rule
  status = :Optimal

  # simplex search
  apfrac = Array{Float64, 1}(m)
  w = Array{Float64, 1}(m); d = Array{Float64, 1}(m)
  λ = Array{Float64, 1}(m); Up = Array{Float64, 1}(m)
  while !(q == 0 || iter >= max_iter)
    # finding viable columns to enter basis
    iter += 1
    #@assert all(xB .>= 0)
    (L == 0) ? w .= A[:,ℕ[q]] : w .= L\((A[:,ℕ[q]].*Rs)[prow])
    for j in 1:ups
      permute!(w, P[j])
      w[end] -= dot(MP[j], w)
    end
    d .= U\w
    apfrac .= xB ./ d # relative variable changes to direction
    indpos = find(d .> 1e-12) # variables that decrease in d direction
    if length(indpos) == 0
      status = :Unbounded; break
    end
    indxq = indmin(apfrac[indpos])
    xq = apfrac[indpos[indxq]]
    #@assert xq >= 0
    p = findfirst(apfrac, xq) # Bland's Rule

    # column change
    xB -= xq * d; xB[p] = xq # upfiz minha inscrição no SIICUSP/2018 poisdate solution
    𝔹[p], ℕ[q] = ℕ[q], 𝔹[p] # update indexes
    if ups >= maxups # reset LU
      F = lufact(A[:,𝔹])
      L, U, prow, pcol, Rs = F[:(:)] # (Rs.*A)[prow,pcol]*x = b[prow]
      𝔹, xB = 𝔹[pcol], xB[pcol]
      ups = 0
    else # update LU
      ups += 1
      U[:,p] .= w
      if nnz(Ufiller) < nnz(U)
        nnz(Ufiller) == 0 ? Ufiller = similar(U) : copy!(Ufiller, U)
      end
      P[ups] = reverse(reverse(1:m,p,m),p,m-1)
      MP[ups] = spzeros(m)
      halfperm!(Ufiller, U, P[ups])
      Up .= Ufiller[:,p]
      for i in p:m-1
        (MP[ups])[i] = Up[i]/Ufiller[i,i+1]
        for j in nzrange(Ufiller,i+1)
          Up[Ufiller.rowval[j]] -= (MP[ups])[i]*Ufiller.nzval[j]
        end
        Ufiller[i, p] = 0
      end
      Ufiller[end,p] = Up[end]
      halfperm!(U, Ufiller, P[ups])
      𝔹, xB = 𝔹[P[ups]], xB[P[ups]]
    end

    # check optimality and choose variable to leave basis if necessary
    λ .= U'\(@view c[𝔹])
    for j in 1:ups
      λ .= (-).(λ, λ[end]*MP[ups-j+1])
      λ[P[ups-j+1]] = λ
  end
    AN[:,q] = A[:,ℕ[q]] # do something similar to 𝔹, use a for
    (L==0) ? r .= (-).(cN, (λ'*AN)') : r .= (-).(cN, ((ipermute!(L\λ, prow).*Rs)'*AN)')
    q = findfirst(r .< -1e-12) # Bland's Rule
  end
  if iter >= max_iter
    status = :UserLimit
  end

  # finalization
  x = zeros(n)
  if !artificial
    x[𝔹] = xB
    z = dot(c, x)
  else
    Irows = collect(1:m)
    if dot(xB, c[𝔹])/norm(xB) > 1e-12
      status = (iter >= max_iter) ? :UserLimit : :Infeasible
      I = find(𝔹 .<= n - m)
      x[𝔹[I]] = xB[I]
      z = dot(co, x)
    elseif maximum(𝔹) > n # check for artificial variables in basis
      # remove artificial variables from basis
      deleteat!(ℕ, find(ℕ .> n))
      p = findfirst(𝔹 .> n)
      if ups != 0
        F = lufact(A[Irows,𝔹])
        L, U, prow, pcol, Rs = F[:(:)]
        𝔹, xB = 𝔹[pcol], xB[pcol]
      end
      Ap = Array{Float64, 1}(m)
      while p != 0
        q = 1
        (prow == 0) ? Ap .= A[Irows,𝔹[p]] : Ap .= A[Irows,𝔹[p]][prow]
        PivotAp = findfirst(Ap .> 0)
        while q <= length(ℕ) # searching for columns to substitute artificials ℕ basis
          d .= U\(L\((A[Irows,ℕ[q]].*Rs)[prow]))
          (abs(d[PivotAp]) > 1e-12) ? break : q += 1
        end
        if q > length(ℕ)
          deleteat!(Irows, findfirst(A[Irows,𝔹[p]]))
          deleteat!(𝔹, p); deleteat!(xB, p)
          deleteat!(Ap, p); deleteat!(d, p)
          F = lufact(A[Irows,𝔹])
          L, U, prow, pcol, Rs = F[:(:)]
          𝔹, xB = 𝔹[pcol], xB[pcol]
        else # falta trocar a coluna
          𝔹[p] = ℕ[q]
          F = lufact(A[Irows,𝔹])
          L, U, prow, pcol, Rs = F[:(:)]
          𝔹, xB = 𝔹[pcol], xB[pcol]
        end
        p = findfirst(𝔹 .> n)
      end
      return simplexluup(co, A[Irows,1:n], b[Irows], 𝔹, L, U, prow, Rs, xB) # stop creating matrix to return
    else
      return simplexluup(co, A[:,1:n], b, 𝔹, L, U, prow, Rs, xB)
    end
  end
  return x, z, status
end
