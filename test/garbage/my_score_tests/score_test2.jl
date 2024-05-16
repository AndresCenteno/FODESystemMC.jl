using FODESystemMC
using HypothesisTests
using Random
using Statistics

Random.seed!(1)

n = 5
problem = myrand(randFODESystem(),n)
init_node = 1

MCsolNoBranching = MCSolver(problem,init_node,SaveSamplesNoBranching();nsims=100000)
DETsol = FD_L1Solver(problem,init_node;Nt=2000)

OneSampleTTest(MCsolNoBranching[1],DETsol[1])
OneSampleHotellingT2Test(MCsolNoBranching[2]',DETsol[2])