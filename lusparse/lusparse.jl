#lu without pivoting
function lubasic(A)
  Al = copy(A)
  n = size(Al,1)
  for k in 1:n-1
    rows = k+1:n
    Al[rows,k] /= Al[k,k]
    for j in rows
      Al[j,rows] -= Al[j,k]*Al[k,rows]
    end
  end
  U = zeros(n,n)
  L = eye(n)
  for j in 1:n
    U[j,j:n] = Al[j,j:n]
  end
  for j in 2:n
    L[j,1:j-1] = Al[j,1:j-1]
  end
  return L, U
end

#lu with partial pivoting
function lupartial(A)
  Al = copy(A)
  n = size(Al,1)
  p = collect(1:n)
  for k in 1:n-1
    rows = k+1:n
    indpivot = indmax(abs.(Al[k:n,k])) + k - 1
    p[k], p[indpivot] = p[indpivot], p[k]
    Al[k,:], Al[indpivot,:] = Al[indpivot,:], Al[k,:]
    Al[rows,k] /= Al[k,k]
    for j in rows
      Al[j,rows] -= Al[j,k]*Al[k,rows]
    end
  end
  U = spzeros(n,n)
  L = speye(n)
  for j in 1:n
    U[j,j:n] = Al[j,j:n]
  end
  for j in 2:n
    L[j,1:j-1] = Al[j,1:j-1]
  end
  return L, U, p
end

#lu with partial pivoting weighting stability and sparsity

function lusparse(Al; α = 100, β = 1.2) #α is how much bigger the pivot is to be worth the change,
#and β how much sparser the pivot's line is to be worth the change
  n = size(Al,1)
  p = collect(1:n)
  for k in 1:n-1
    rows = k+1:n
    abspivot = abs(Al[p[k],k])
    nzpivot = countnz(Al[p[k],k:n])
    abss = abs.(Al[p[rows],k]) #getting all pivot candidate's absolute value
    nzs = collect(rows)
    for i in rows #counting nonzeros in each pivot candidate's line
      nzs[i-k] = countnz(Al[p[i],k:n])
    end
    for i in rows #pivoting based on stability and sparsity
      if (abss[i-k] > α*abspivot && β*nzpivot > nzs[i-k] && abss[i-k] > 1e-14) || (β*nzs[i-k] < nzpivot && abspivot < α*abss[i-k] && abss[i-k] > 1e-14) || (abspivot < 1e-14 && abss[i-k] > 1e-14)
        p[k], p[i] = p[i], p[k]
        abspivot = abss[i-k]
        nzpivot = nzs[i-k]
      end
    end
    #applying Gaussian reduction
    Al[p[rows],k] /= Al[p[k],k]
    for j in rows
      Al[p[j],rows] -= Al[p[j],k]*Al[p[k],rows]
    end
  end
  U = spzeros(n,n) #L and U construction for tests
  for j in 1:n
    U[j,j:n] = Al[p[j],j:n]
    Al[p[j],j:n] = 0
  end
  dropzeros(Al)
  return Al[p,:] + speye(n), U, p
end
