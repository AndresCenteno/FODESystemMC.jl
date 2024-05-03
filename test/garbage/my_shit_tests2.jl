# comparing
using Random
using LinearAlgebra
using Statistics
using FODESystemMC
using HypothesisTests

Random.seed!(0)

n = 5
problem = myrand(randFODESystem(),n)
init_state = 1

MCSol = MCSolver(problem,init_state, SaveSamples(); nsims=100000)
DETSol = L1Solver(problem; Nt=2000)

OneSampleTTest(MCSol,DETSol[init_state,end])
