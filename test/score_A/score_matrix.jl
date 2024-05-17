# checking, or trying to do, score function over a matrix

using FODESystemMC, Random, Statistics

Random.seed!(0)
nnodes = 30; init_node = 1
problem = myrand(randFODESystem(),nnodes)
DETsol = FD_L1Solver(problem,init_node;Nt=1000)

nsims = Int(1e5)
# testing no save
# MCsol = MCSolver(problem,init_node;nsims=nsims)
# tests_passed = compare(DETsol,MCsol)
# tests_passed

# testing save
MCsol = MCSolver(problem,init_node,SaveSamples();nsims=nsims)
tests_passed = compare(DETsol,MCsol)
tests_passed
@show tests_passed
@show fieldnames(FODESystemMC.forwback)

s = getmeans(MCsol)
compare(DETsol,s)