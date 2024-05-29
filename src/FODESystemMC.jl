module FODESystemMC

using LinearAlgebra, SparseArrays, SharedArrays
using Statistics: mean
using UnPack
using Distributed
using QuadGK

export FODESystem, randFODESystem, MCSolver, L1Solver, FD_L1Solver
export SaveSamples, SaveSamplesNoBranching, NoSave, MCDecomposition, myrand, SharedSamples
export SaveSamplesThreaded
export LaplaceInv, StronglyTyped, Zyg
# API for tests, will probably delete next line when everything works
export sojourn, score
# API for comparisons
export compare, getmeans, forwback

include("types.jl")
include("utils.jl")
include("solvers/deterministic_solvers.jl"), include("solvers/stochastic_solvers.jl")
include("solvers/strongly_typed_solver.jl"); include("solvers/parallel_solvers.jl")
include("solvers/thread_solvers.jl")
include("aux/mittag_leffler.jl"); include("aux/matlab.jl")

# for benchmarking 
include("solvers/one_sim.jl")
export one_sim!, MCDecomposition2
end
