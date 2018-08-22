include("lusparse.jl")

function lutest()
  density = 0.075
  n = 100 # dimension of the square matrix A
  α = 100 # how much bigger the pivot is to be worth the change
  β = 1.1 # how much sparser the pivot's line is to be worth the change

  # building A with given density, ~90% nonzero random values from (-3:2:3) and <~10% either -1.5*α or 1.5*α
  elements = [zeros(round(n*(1-density))); rand(-3:2:3,round(Int,n*density*0.9)); rand([-1.5*α; 1.5*α],max(round(Int,n*density*0.1),1))]
  A = rand(elements, n, n)
  if det(A) == 0
    return error("A is singular, increase density")
  end

  # checking lufact solution
  F = lufact(A)
  LUfact = F[:L]*F[:U]
  pfact = F[:p]
  if norm(A[pfact,:] - LUfact) < 1
    println("lufact solved!")
  else
    println("ill problem!")
  end

  # checking our solution
  L, U, p = lusparse(A, α = α, β = β)
  ourLU = L*U
  if NaN in Set(ourLU)
    println("lusparse failed! NaN found.")
  elseif norm(A[p,:] - ourLU) < 1
    println("lusparse solved!")
    extrasparsity = 100*round((1-countnz(U)/countnz(F[:U])),4)
    println("U obtained by lusparse is $extrasparsity","% sparser!")
  else
    println("lusparse failed! unknown reason!")
  end

  return L, U, p, A
end

L, U, p, A = lutest()
