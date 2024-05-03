using FODESystemMC
using HypothesisTests
using Random
using Statistics
using NaNStatistics

Random.seed!(0)

n = 5
problem = myrand(randFODESystem(),n)
init_node = 1

MCsol = MCSolver(problem,init_node,SaveSamples();nsims=100000)
DETsol = FD_L1Solver(problem,init_node;Nt=1000)
OneSampleTTest(MCsol[1],DETsol[1])

# getting mfking nans
stack([filter(!isnan,col) for col in eachcol(MCsol[2])])
soj = sojourn(problem.A,problem.α)
mean_rowwise = [mean(filter(!isnan,row)) for row in eachrow(MCsol[2])]

DETsol[2]
MCder = MCsol[2][:,(all(!isnan,MCsol[2],dims=1))[:]]
OneSampleHotellingT2Test(MCder',DETsol[2])