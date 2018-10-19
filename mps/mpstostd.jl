#The function 'mpstostd()' receives a MPS file and returns its standard form
#min [c; zeros]'x
#s.a. [A C]x = b
#     x >= 0,
#in the order [c; zeros], [A C], b, status.
#
#The GLPK solver is used to convert the MPS file to a problem with form 
#min c'x
#s.a. l <= Ax <= u
#     xl <= x <= xu,
#which is then explicitly standardized.

using GLPKMathProgInterface, JuMP, MathProgBase

function mpstostd(mpsfile)
    mod = Model(solver=GLPKSolverLP())
	internal_model = MathProgBase.LinearQuadraticModel(GLPKSolverLP())
    MathProgBase.loadproblem!(internal_model, mpsfile)

	c = MathProgBase.getobj(internal_model)
    A = MathProgBase.getconstrmatrix(internal_model)
    m, n = size(A)
    xl = MathProgBase.getvarLB(internal_model)
    xu = MathProgBase.getvarUB(internal_model)
    l = MathProgBase.getconstrLB(internal_model)
    u = MathProgBase.getconstrUB(internal_model)
    zchange = 0
    status = :Successful
#Notation: A = [a1 ... ai ... an].
#For simplicity problems that have free equations such as
#-Inf <= aixi <= Inf, or fixed unbounded values
#such as -Inf <= aixi <= -Inf or Inf <= aixi <= Inf,
#will be considered ill or unbounded problems and be discarded,
#even though the problem could still possibly be solved.
#
#Here we'll transform xl <= x <= xu into x >= 0.
	for i = 1:n
	  if xl[i] == -Inf
	    if xu[i] == -Inf
	      status = :UnboundedConstraint
	      break
	    elseif xu[i] == Inf #free variable -Inf <= x <= Inf
	      A = [A -A[:,i]]
	      push!(c, -c[i])
	    else # -Inf < x <= u
	      u -= xu[i]*A[:,i]
	      l -= xu[i]*A[:,i]
	      A[:,i] = -A[:,i]
	      zchange += c[i]*xu[i]
	      c[i] = -c[i]
	    end
	  elseif xl[i] == Inf
	    status = :UnboundedConstraint
	    break
	  elseif xu[i] == -Inf
	    status = :UnboundedConstraint
	    break
	  elseif xu[i] == Inf
	    u -= xl[i]*A[:,i]
	    l -= xl[i]*A[:,i]
	    zchange += c[i]*xl[i]
	  elseif xu[i] == xl[i]
        status = :FixedVariables #Can't we just remove these variables?
        break
      else
        u -= xl[i]*A[:,i]
	    l -= xl[i]*A[:,i]
	    zchange += c[i]*xl[i]
	    xu[i] -= xl[i]
	    A = [A zeros(m)]
	    slackrow = zeros(n+1)'
        slackrow[i] = 1
        slackrow[end] = 1
        A = [A; slackrow]
        push!(l, xu[i])
        push!(u, xu[i])
        push!(c, 0)
        m, n = size(A)
	  end
	end

#Next we'll transform l <= Ax <= u into Ax = b.
    b = zeros(m);
	new_variables = 0
	for i = 1:m
	  if status != :Successful
	    break
	  end
	  slack = zeros(size(A,1))
	  if l[i] == -Inf
	    if u[i] == Inf
	      status = :FreeEquations #Can't we just remove these lines?
	      break
	    elseif u[i] == -Inf
	      	status = :UnboundedConstraint
	      	break
	    else # -Inf <= aixi <= ui < Inf
	      slack[i] = 1.0
	      b[i] = u[i]
	      A = [A slack]
	      new_variables += 1
	    end
	  elseif l[i] == Inf
	    status = :UnboundedConstraint
	    break
	  elseif u[i] == -Inf
	    status = :UnboundedConstraint
	    break
	  elseif u[i] == Inf # -Inf < li <= aixi <= Inf
	    slack[i] = -1.0
	    b[i] = l[i]
	    A = [A slack]
	    new_variables += 1
	  elseif l[i] == u[i] # bi <= aixi <= bi -> aixi = bi
	      b[i] = l[i]
	  else # -Inf < li <= aixi <= ui < Inf
	    A = [A; A[i,:]']
	    slack = zeros(size(A,1))
	    slack[i] = 1; A = [A slack]
	    slack = zeros(size(A,1))
	    slack[end] = -1; A = [A slack]
	    b[i] = u[i]
	    push!(b, l[i])
	    new_variables += 2
	  end
	end
	c = [c; zeros(new_variables)]

	return c, A, b, status, zchange
end
