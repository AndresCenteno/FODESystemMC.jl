using SharedArrays
using Distributed

addprocs(5)

@everywhere begin

using FODESystemMC
using Random
using DelimitedFiles

end

ntests = 30
max_nnodes = 15
test_passed = SharedArray{Bool}(5,ntests)

@distributed for test = 1:ntests
    @show test
    Random.seed!(test)
    nnodes = rand(3:max_nnodes)
    init_node = rand(1:nnodes)
    problem = myrand(randFODESystem(),nnodes;sparse=rand(0:1),density = 0.05)
    DETsol = FD_L1Solver(problem,init_node; Nt = 1000)
    nsims = rand(Int(1e4):Int(1e5))
    MCsol = MCSolver(problem,init_node,SharedSamples();nsims = nsims)
    test_passed[:,test] = compare(DETsol,MCsol)
    @show test_passed
end

writedlm("test/parallel_test")