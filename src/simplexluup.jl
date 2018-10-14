export simplexluup

# myprod! is equivalent to xD = (l'*A[:,N])'
# avoids creating A[:,N]
function myprod!(xD::Vector{Float64}, l::Vector{Float64},
                 A::SparseMatrixCSC{Float64,Int64}, N::Vector{Int64})
  xD .= 0.
  for (colN,colA) in enumerate(N)
    for i in nzrange(A,colA)
      if l[A.rowval[i]] != 0
        xD[colN] += l[A.rowval[i]]*A.nzval[i]
      end
    end
  end
end

# getcX! is equivalent to l .= c[X]
# avoids creating a vector when calling c[X]
function getcX!(l::Vector{Float64}, c::Vector{Float64}, X::Vector{Int64})
  l .= 0.
  for (i,j) in enumerate(X)
    l[i] = c[j]
  end
end

# getAcol! is equivalent to l .= A[:,col]
# avoids creating a vector when calling A[:,col]
function getAcol!(l::Vector{Float64},
                  A::SparseMatrixCSC{Float64,Int64}, col::Int64)
  l .= 0.
  for i in nzrange(A,col)
    l[A.rowval[i]] = A.nzval[i]
  end
end
function getAcol!(l::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64},
                  col::Int64, Irows::Vector{Int64})
  l .= 0.; j = 0
  for i in nzrange(A,col)
    j = findfirst(Irows,A.rowval[i])
    if j != 0
      l[j] = A.nzval[i]
    end
  end
end

# same as l -= z*a
function subdot!(l::Vector{Float64},a,z::Float64,m::Int64)
  for i in 1:m
    if a[i] != 0
      l[i] -= z*a[i]
    end
  end
end

#returns l[perm] without allocating unnecessary memory
#returns l[invperm(perm)] if inv = true
function savepermute!(tempperm, perm, l; inv::Bool = false)
  copy!(tempperm, perm)
  !inv ? permute!!(l, tempperm) : ipermute!!(l, tempperm)
end


function simplexluup(c::Vector{Float64}, A::SparseMatrixCSC{Float64,Int64},
                     b::Vector{Float64}, 𝔹=0, L=0, U=0, prow::Vector{Int64}=Int64[],
                     Rs::Vector{Float64}=Float64[], xB::Vector{Float64}=Float64[];
                     max_iter::Int64 = 10000, maxups::Int64 = 15)

  m, n = size(A) # preparations
  iter = 0; ups = 0; maxed = false
  Ufiller = spzeros(0)
  P, MP = Vector{Vector{Int64}}(maxups), Vector{SparseVector{Float64,Int64}}(maxups)
  cB = Vector{Float64}(m)
  lnz = Ref{Int64}(); unz = Ref{Int64}(); nz_diag = Ref{Int64}()
  n_row = Ref{Int64}(); n_col = Ref{Int64}()
  Lp = Vector{Int64}(m + 1); Up = Vector{Int64}(m + 1)
  pcol = Vector{Int64}(m); tempperm = Vector{Int64}(m)
  if L == 0
    Rs = Vector{Float64}(m)
    prow = Vector{Int64}(m)
  end

  if 𝔹 == 0 # construct artificial problem
    artificial = true
    signb = sign.(b)
    Ao = A
    U = spdiagm(signb); A = [A U]
    Utri = UpperTriangular(U)
    𝔹 = collect(n+1:n+m); ℕ = collect(1:n) # artificial indexes
    ca = [zeros(n); ones(m)]; getcX!(cB,ca,𝔹)
    cN = Vector{Float64}(n); getcX!(cN,ca,ℕ)
    xD = Vector{Float64}(n)
    myprod!(xD,signb,A,ℕ)
    r = -xD # artificial relative costs
    xB = abs.(b) # solution in current basis
  else
    artificial = false
    ℕ = setdiff(1:n, 𝔹)
    xD = Vector{Float64}(n-m)
    cN = Vector{Float64}(n-m); getcX!(cN,c,ℕ); getcX!(cB,c,𝔹)
    if L == 0
      F = lufact(A[:,𝔹]) # (Rs.*A)[prow,pcol] * x[pcol] = b[prow]
      xB = F\b
      L, U, prow, pcol, Rs = F[:(:)]
      Utri = UpperTriangular(U); Ltri = LowerTriangular(L)
      savepermute!(tempperm, pcol, 𝔹)
      permute!!(xB, pcol) # pcol is lost
      At_ldiv_B!(Utri,cB); At_ldiv_B!(Ltri,cB)
      savepermute!(tempperm, prow, cB, inv = true)
      cB .= cB.*Rs
      myprod!(xD,cB,A,ℕ)
      r = (-).(cN, xD)
    else
      Utri = UpperTriangular(U); Ltri = LowerTriangular(L)
      At_ldiv_B!(Utri,cB); At_ldiv_B!(Ltri,cB)
      savepermute!(tempperm, prow, cB, inv = true)
      cB .= cB.*Rs
      myprod!(xD,cB,A,ℕ)
      r = (-).(cN, xD)
    end
  end
  q = findfirst(r .< -1e-12) # Bland's Rule
  status = :Optimal

  # simplex search
  apfrac = Array{Float64,1}(m)
  w = Array{Float64,1}(m); d = Array{Float64,1}(m)
  λ = Array{Float64,1}(m); Ucolp = Array{Float64,1}(m)
  while !(q == 0 || iter >= max_iter)
    # finding viable columns to enter basis
    iter += 1
    if L == 0
      getAcol!(w,A,ℕ[q])
    else
      getAcol!(w,A,ℕ[q])
      w .= w.*Rs
      savepermute!(tempperm, prow, w)
      A_ldiv_B!(Ltri,w)
    end
    for j in 1:ups
      savepermute!(tempperm, P[j], w)
      w[end] -= dot(MP[j], w)
    end
    d .= w
    A_ldiv_B!(Utri,d)
    apfrac .= xB ./ d # relative variable changes to direction
    xq = Inf
    for i in 1:m # find min(apfrac[d .> 0])
      if d[i] > 1e-12
        if apfrac[i] < xq
          xq = apfrac[i]
          p = i
        end
      end
    end
    if xq == Inf
      status = :Unbounded; break
    end

    # column change
    subdot!(xB,d,xq,m)
    xB[p] = xq # update solution
    𝔹[p], ℕ[q] = ℕ[q], 𝔹[p] # update indexes
    if ups >= maxups # reset LU
      maxed = true
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
        Ltri = LowerTriangular(L)
      else
        copy!(L, transpose(SparseMatrixCSC(m, m, increment!(Lp), increment!(Lj), Lx)))
      end
      copy!(U, SparseMatrixCSC(m, m, increment!(Up), increment!(Ui), Ux))
      increment!(prow); increment!(pcol)
      savepermute!(tempperm, pcol, 𝔹)
      permute!!(xB, pcol) # pcol is lost
      ups = 0
    else # update LU
      U[:,p] .= w
      if findlast(w) > p
        ups += 1
        if maxed
          P[ups] .= 1:m; reverse!(reverse!(P[ups],p,m),p,m-1)
        else
          P[ups] = 1:m; reverse!(reverse!(P[ups],p,m),p,m-1)
        end
        maxed? MP[ups] .= spzeros(m) : MP[ups] = spzeros(m)
        if nnz(Ufiller) < nnz(U)
          nnz(Ufiller) == 0 ? Ufiller = similar(U) : copy!(Ufiller, U)
        end
        halfperm!(Ufiller, U, P[ups])
        getAcol!(Ucolp,Ufiller,p)
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
        dropzeros!(Ufiller)
        halfperm!(U, Ufiller, P[ups])
        savepermute!(tempperm, P[ups], 𝔹)
        savepermute!(tempperm, P[ups], xB)
      end
    end

    # check optimality and choose variable to leave basis if necessary
    artificial? getcX!(λ,ca,𝔹) : getcX!(λ,c,𝔹)
    At_ldiv_B!(Utri,λ)
    for j in 1:ups
      subdot!(λ,MP[ups-j+1],λ[end],m)
      savepermute!(tempperm, P[ups-j+1], λ, inv = true)
    end
    if L == 0
      myprod!(xD,λ,A,ℕ)
      artificial? getcX!(cN,ca,ℕ) : getcX!(cN,c,ℕ)
      r .= (-).(cN, xD)
    else
      At_ldiv_B!(Ltri,λ)
      savepermute!(tempperm, prow, λ, inv = true)
      λ .= λ.*Rs
      myprod!(xD,λ,A,ℕ)
      artificial? getcX!(cN,ca,ℕ) : getcX!(cN,c,ℕ)
      r .= (-).(cN, xD)
    end
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
    artificial? getcX!(λ,ca,𝔹) : getcX!(λ,c,𝔹)
    if dot(xB, λ)/norm(xB) > 1e-12
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
        if L ==0
          getAcol!(Ap,A,𝔹[p],Irows)
        else
          getAcol!(Ap,A,𝔹[p],Irows)
          savepermute!(tempperm, prow, Ap)
        end
        for j in 1:ups
          savepermute!(tempperm, P[j], Ap)
        end
        PivotAp = findfirst(Ap)
        while q <= length(ℕ) # searching for columns to substitute artificials ℕ basis
          if L == 0
            getAcol!(d,A,ℕ[q],Irows)
          else
            getAcol!(d,A,ℕ[q],Irows)
            d .= d.*Rs
            savepermute!(tempperm, prow, d)
            A_ldiv_B!(Ltri,d)
          end
          for j in 1:ups
            savepermute!(tempperm, P[j], d)
            d[end] -= dot(MP[j], d)
          end
          A_ldiv_B!(Utri,d)
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
        Utri = UpperTriangular(U); Ltri = LowerTriangular(L)
        ups = 0
        savepermute!(tempperm, pcol, 𝔹)
        permute!!(xB, pcol) # pcol is lost
        p = findfirst(𝔹 .> n)
      end
      return simplexluup(c, A[Irows,1:n], b[Irows], 𝔹, L, U, prow, Rs, xB)
    else
      return simplexluup(c, Ao, b, 𝔹, L, U, prow, Rs, xB)
    end
  end
  return x, z, status
end

