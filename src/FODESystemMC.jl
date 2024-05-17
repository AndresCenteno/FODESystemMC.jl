module FODESystemMC

using LinearAlgebra
using Statistics: mean
using UnPack
using MittagLeffler # to store Malliavin weights efficiently
using FiniteDifferences

export FODESystem, randFODESystem, MCSolver, L1Solver, FD_L1Solver
export SaveSamples, SaveSamplesNoBranching, NoSave, MCDecomposition, myrand
# API for tests, will probably delete next line when everything works
export sojourn, score
# API for comparisons
export compare, getmeans

include("types.jl")
include("utils.jl")
include("solvers/deterministic_solvers.jl"), include("solvers/stochastic_solvers.jl")


end
