export simplexluup

function simplexluup(c, A, b, 𝔹=0, L=0, U=0, prow=0, Rs=0, xB=0; max_iter = 10000, maxups = 15)
  # preparations
  m, n = size(A)
  iter = 0; ups = 0
  Ufiller = spzeros(0)
  P, MP = Vector{Vector{Int64}}(maxups), Vector{SparseVector{Float64,Int64}}(maxups)
  lnz = Ref{Int64}(); unz = Ref{Int64}(); nz_diag = Ref{Int64}()
  n_row = Ref{Int64}(); n_col = Ref{Int64}()
  Lp = Vector{Int64}(m + 1); Up = Vector{Int64}(m + 1)
  pcol = Vector{Int64}(m)
  tempperm = Vector{Int64}(m)
  if Rs == 0
    Rs = Vector{Float64}(m)
    prow = Vector{Int64}(m)
  end

  if 𝔹 == 0 # construct artificial problem
    artificial = true
    signb = sign.(b*1.)
    AN = copy(A); U = spdiagm(signb)
    Ao = A
    A = [A U] # try to save A memory to return
    𝔹 = collect(n+1:n+m); ℕ = collect(1:n) # artificial indexes
    ca = [zeros(n); ones(m)]; cN = @view ca[ℕ]
    r = -(signb'*AN)' # artificial relative costs
    xB = abs.(b) # solution in current basis
  else
    artificial = false
    ℕ = setdiff(1:n, 𝔹)
    AN = A[:,ℕ]; cN = @view c[ℕ]
    if L == 0 # if user gives 𝔹
      F = lufact(A[:,𝔹]) # (Rs.*A)[prow,pcol] * x[pcol] = b[prow]
      xB = F\b
      L, U, prow, pcol, Rs = F[:(:)]
      copy!(tempperm, pcol); permute!!(𝔹, tempperm)
      copy!(tempperm, pcol); permute!!(xB, tempperm)
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
  λ = Array{Float64, 1}(m); Ucolp = Array{Float64, 1}(m)
  while !(q == 0 || iter >= max_iter)
    # finding viable columns to enter basis
    iter += 1
    (L == 0) ? w .= A[:,ℕ[q]] : w .= L\((A[:,ℕ[q]].*Rs)[prow])
    for j in 1:ups
      copy!(tempperm, P[j])
      permute!!(w, tempperm)
      w[end] -= dot(MP[j], w)
    end
    d .= U\w
    apfrac .= xB ./ d # relative variable changes to direction
    apfracpos = apfrac[find(d .> 1e-12)] # variables that decrease in d direction
    if length(apfracpos) == 0
      status = :Unbounded; break
    end
    xq = minimum(apfracpos)
    p = findfirst(apfrac, xq) # Bland's Rule

    # column change
    xB -= xq*d; xB[p] = xq # update solution
    𝔹[p], ℕ[q] = ℕ[q], 𝔹[p] # update indexes
    if ups >= maxups # reset LU
      F = lufact(A[:,𝔹])
      ccall(("umfpack_dl_get_lunz",:libumfpack), Int64,(Ptr{Int64},Ptr{Int64},
            Ptr{Int64},Ptr{Int64},Ptr{Int64},Ptr{Void}),
            lnz,unz,n_row,n_col,nz_diag,F.numeric)
      Lj = Vector{Int64}(lnz[]); Lx = Vector{Float64}(lnz[])
      Ui = Vector{Int64}(unz[]); Ux = Vector{Float64}(unz[])
      ccall(("umfpack_dl_get_numeric",:libumfpack),Int64, (Ptr{Int64},
             Ptr{Int64},Ptr{Float64},Ptr{Int64},Ptr{Int64},Ptr{Float64},
             Ptr{Int64},Ptr{Int64},Ptr{Void},Ref{Int64},Ptr{Float64},
             Ptr{Void}),Lp,Lj,Lx,Up,Ui,Ux,prow,pcol,C_NULL,0, Rs, F.numeric)
      if L == 0
        L = transpose(SparseMatrixCSC(m, m, increment!(Lp), increment!(Lj), Lx))
      else
        copy!(L, transpose(SparseMatrixCSC(m, m, increment!(Lp), increment!(Lj), Lx)))
      end
      copy!(U, SparseMatrixCSC(m, m, increment!(Up), increment!(Ui), Ux))
      increment!(prow); increment!(pcol)
      copy!(tempperm, pcol); permute!!(𝔹, tempperm)
      copy!(tempperm, pcol); permute!!(xB, tempperm)
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
      Ucolp .= Ufiller[:,p]
      for i in p:m-1
        if Ucolp[i] != 0
          (MP[ups])[i] = Ucolp[i]/Ufiller[i,i+1]
          for j in nzrange(Ufiller,i+1)
            Ucolp[Ufiller.rowval[j]] -= (MP[ups])[i]*Ufiller.nzval[j]
          end
        end
      end
      Ufiller.nzval[nzrange(Ufiller,p)] = 0
      Ufiller[end,p] = Ucolp[end]
      halfperm!(U, Ufiller, P[ups])
      𝔹, xB = 𝔹[P[ups]], xB[P[ups]]
    end

    # check optimality and choose variable to leave basis if necessary
    artificial ? λ .= U'\(@view ca[𝔹]) : λ .= U'\(@view c[𝔹])
    for j in 1:ups
      λ .= (-).(λ, λ[end]*MP[ups-j+1])
      copy!(tempperm, P[ups-j+1])
      ipermute!!(λ, tempperm)
    end
    AN[:,q] = A[:,ℕ[q]] # do something similar to 𝔹, use a for
    (L==0) ? r .= (-).(cN, (λ'*AN)') : r .= (-).(cN, ((ipermute!(L'\λ, prow).*Rs)'*AN)')
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
    if dot(xB, ca[𝔹])/norm(xB) > 1e-12
      status = (iter >= max_iter) ? :UserLimit : :Infeasible
      I = find(𝔹 .<= n - m)
      x[𝔹[I]] = xB[I]
      z = dot(c, x)
    elseif maximum(𝔹) > n # check for artificial variables in basis
      # remove artificial variables from basis
      deleteat!(ℕ, find(ℕ .> n))
      p = findfirst(𝔹 .> n)
      Ap = Array{Float64, 1}(m)
      while p != 0
        q = 1
        (L == 0) ? Ap .= A[Irows,𝔹[p]] : Ap .= A[Irows,𝔹[p]][prow]
        for j in 1:ups
          copy!(tempperm, P[j])
          permute!!(Ap, tempperm)
        end
        PivotAp = findfirst(Ap .> 0)
        while q <= length(ℕ) # searching for columns to substitute artificials ℕ basis
          (L == 0) ? d .= A[:,ℕ[q]] : d .= L\((A[:,ℕ[q]].*Rs)[prow])
          for j in 1:ups
            copy!(tempperm, P[j])
            permute!!(d, tempperm)
            d[end] -= dot(MP[j], d)
          end
          d .= U\d
          (abs(d[PivotAp]) > 1e-12) ? break : q += 1
        end
        if q > length(ℕ)
          deleteat!(Irows, findfirst(A[Irows,𝔹[p]]))
          deleteat!(𝔹, p); deleteat!(xB, p)
          deleteat!(Ap, p); deleteat!(d, p)
        else
          𝔹[p] = ℕ[q]
        end
        F = lufact(A[Irows,𝔹])
        L, U, prow, pcol, Rs = F[:(:)]
        ups = 0
        copy!(tempperm, pcol); permute!!(𝔹, tempperm)
        copy!(tempperm, pcol); permute!!(xB, tempperm)
        p = findfirst(𝔹 .> n)
      end
      return simplexluup(c, Ao[Irows,:], b[Irows], 𝔹, L, U, prow, Rs, xB)
    else
      return simplexluup(c, Ao, b, 𝔹, L, U, prow, Rs, xB) # stop creating matrix to return
    end
  end
  return x, z, status
end
