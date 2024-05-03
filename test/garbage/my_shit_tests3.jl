using FODESystemMC
using HypothesisTests
using Random
using Statistics

Random.seed!(0)

n = 10
problem = myrand(randFODESystem(),n)
init_node = 1

# I suddenly get a domain error after writing the -diagA
# I wrote too much branching
MCsol = MCSolver(problem,init_node,SaveSamples();nsims=100000)
MCsolNoBranching = MCSolver(problem,init_node,SaveSamplesNoBranching();nsims=100000)

DETsol = FD_L1Solver(problem,init_node;Nt=1000)

# forward
OneSampleTTest(MCsol[1],DETsol[1])
OneSampleHotellingT2Test(MCsol[2]',DETsol[2])

# forward no branching
OneSampleTTest(MCsolNoBranching[1],DETsol[1])
OneSampleHotellingT2Test(MCsolNoBranching[2]',DETsol[2])