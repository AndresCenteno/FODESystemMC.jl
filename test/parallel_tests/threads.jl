# this one will attempt to do the tests in parallel using threads for ulam
using Pkg
Pkg.activate(".")
using FODESystemMC
using Random
using DelimitedFiles

ntests = 3
nnodes = 10
test_passed = zeros(Bool,5,ntests)

Threads.@threads for test = 1:ntests
    @show test
    Random.seed!(test)
    init_node = rand(1:nnodes)
    problem = myrand(randFODESystem(),nnodes)
    DETsol = FD_L1Solver(problem,init_node; Nt = 1000)
    nsims = rand(Int(1e4):Int(1e5))
    MCsol = MCSolver(problem,init_node,SharedSamples();nsims = nsims)
    test_passed[:,test] = compare(DETsol,MCsol)
    @show test_passed
end
writedlm("test/parallel_tests/test_passed.csv",test_passed)
