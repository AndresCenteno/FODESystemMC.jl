using FODESystemMC
using HypothesisTests
using Random
using Statistics
using NaNStatistics

Random.seed!(0)

n = 10
problem = myrand(randFODESystem(),n)
init_node = 1

# I suddenly get a domain error after writing the -diagA
MCsol = MCSolver(problem,init_node,SaveSamples();nsims=100000)
DETsol = FD_L1Solver(problem,init_node;Nt=1000)

# forward
OneSampleTTest(MCsol[1],DETsol[1])
OneSampleHotellingT2Test(MCsol[2]',DETsol[2])
mean(MCsol[2],dims=2)
DETsol[2]
cov(MCsol[2],dims=2)