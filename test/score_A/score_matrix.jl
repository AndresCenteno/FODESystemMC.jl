# checking, or trying to do, score function over a matrix

using FODESystemMC, Random, Statistics

Random.seed!(1)
nnodes = 4
problem = myrand(randFODESystem(),nnodes)
init_state = 1

# deterministic sensitivities
using FiniteDifferences
uTD = L1Solver(problem,init_state;Nt=500)
f(α,A) = L1Solver(A,problem.u0,α,problem.T,init_state;Nt=500)
duTdαD = grad(central_fdm(2, 1), α->f(α,problem.A), problem.α)[1]
duTdAD = grad(central_fdm(2, 1), A->f(problem.α,A), problem.A)[1]

uTS, duTdαS, duTdAS = @time MCSolver(problem,init_state;nsims=Int(1e5))

@show uTS, uTD
display(duTdαD); display(duTdαS);

display(duTdAD)
display(duTdAS);