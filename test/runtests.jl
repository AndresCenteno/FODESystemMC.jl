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

@testset "ForwardBackward" begin
    Random.seed!(0)
    display("Testing null-hypothesis H₀: E(X) = μ₀ for forward and backward")
    for test=1:5
        display("Forward-backward test #$(test)")
        nnodes = rand(2:10); nsims = Int(1e6); init_node = rand(1:nnodes)
        problem = myrand(randFODESystem(),nnodes)
        MCsol = MCSolver(problem,init_node,SaveSamplesNoBranching(); nsims=nsims)
        DETsol = FD_L1Solver(problem,init_node; Nt = 2000)
        tst_forward = OneSampleTTest(MCsol[1],DETsol[1])
        tst_backward = OneSampleHotellingT2Test(MCsol[2]',DETsol[2])
        boolforw = (pvalue(tst_forward) > 0.05)
        boolback = (pvalue(tst_backward) > 0.05)
        @test boolforw
        @test boolback
        display("Forward returned $(boolforw), backward returned $(boolback)")
    end
end

# ForwardBackward |   10      1     11  34m47.8s
# I had one error and then it breaks so for test=1:5 it should work

@testset "ScoreFunction" begin
    Random.seed!(0)
    nnodes = 5; init_node = 1
    problem = myrand(randFODESystem(),nnodes)
    DETsol = FD_L1Solver(problem,init_node;Nt=1000)

    nsims = Int(1e5)
    # testing no save

    # testing save
    Random.seed!(0)
    MCsol = MCSolver(problem,init_node,SaveSamples();nsims=nsims)
    tests_passed = compare(DETsol,MCsol)
    @show test_passed
end

# (test_uT, test_duTdα, test_duTdA) = (10, 9, 10)

Random.seed!(1)
nnodes = 5; init_node = 1
problem = myrand(randFODESystem(),nnodes)
DETsol = FD_L1Solver(problem,init_node;Nt=800)

nsims = Int(1e5)
Random.seed!(1)
MCsol = MCSolver(problem,init_node,SaveSamples();nsims=20000)
tests_passed = compare(DETsol,MCsol)

DETsol.duTdT

using Statistics
var(MCsol.duTdT)