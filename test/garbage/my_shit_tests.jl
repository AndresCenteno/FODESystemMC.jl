using FODESystemMC
using LinearAlgebra
using HypothesisTests
using Random
include("shit_aux.jl")
Random.seed!(0)
# problem = FODESystemMC.myrand(randFODESystem(),10)
n = 5
L = random_Laplacian(n); u0 = rand(n); α = rand(n)*0.8 .+ 0.2
T = 1; init_state = 1
problem = FODESystem(-L,u0,α,T)

Random.seed!(0)
# another huge mistake, I have two methods
MCSol = MCSolver(problem,init_node,SaveSamples(),nsims= 100000)
##
DETSol = L1Solver(problem,Nt=8000) # be careful with stability here
OneSampleTTest(MCSol,DETSol[init_node,end]) # still doesn't work, but in the other one works!!