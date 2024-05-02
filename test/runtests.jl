using FODESystemMC
using Test
using HypothesisTests
using Random

@testset "Forward" begin
    # Write your tests here.
    Random.seed!(0)
    for test = 1:10
        display("Forward test #$(test)")
        nnodes = 10; nsims = Int(1e5); init_node = rand(1:nnodes)
        problem = myrand(randFODESystem(),nnodes)
        MCsol = MCSolver(problem,init_node,SaveSamples(); nsims=nsims)
        DETsol = L1Solver(problem; Nt = 4000)
        tst = OneSampleTTest(MCsol,DETsol[init_node,end])
        @test pvalue(tst) > 0.05 # just as they do in the docs!
    end
end