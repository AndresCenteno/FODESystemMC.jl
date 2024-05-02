using FODESystemMC
using HypothesisTests
using Random

Random.seed!(0)

n = 5
problem = myrand(randFODESystem(),n)