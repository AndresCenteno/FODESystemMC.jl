# this one will attempt to do the tests in parallel using threads for ulam
@info "Executing threads.jl"
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using FODESystemMC
using Random
using DelimitedFiles
@show Threads.nthreads()
ntests = 100
nnodes = 5
test_passed = zeros(Int,5,ntests)

Threads.@threads for test = 1:ntests
    @show test
    Random.seed!(test)
    try
        init_node = rand(1:nnodes)
        problem = myrand(randFODESystem(),nnodes)
        DETsol = FD_L1Solver(problem,init_node; Nt = 1000)
        nsims = Int(1e5)
        MCsol = MCSolver(problem,init_node,SharedSamples();nsims = nsims)
        test_passed[:,test] = compare(DETsol,MCsol) # why do I get a 0 at uT
    catch
        test_passed[:,test] = -ones(5)
    end
end

@show test_passed
writedlm("test/parallel_tests/test_passed3_executed_in_REPL.csv",test_passed)