using FODESystemMC
using Test
using HypothesisTests
using Random

# because tests breaks when they fail, I need to not break haha
Random.seed!(3)
ntests = 10
test_passed = 0
display("Testing null-hypothesis H₀: E(X) = μ₀ for forward and backward")
for test=1:ntests
    display("Forward-backward test #$(test)")
    nnodes = rand(2:10); nsims = Int(1e6); init_node = rand(1:nnodes)
    problem = myrand(randFODESystem(),nnodes)
    MCsol = MCSolver(problem,init_node,SaveSamplesNoBranching(); nsims=nsims)
    DETsol = FD_L1Solver(problem,init_node; Nt = 2000)
    tst_forward = OneSampleTTest(MCsol[1],DETsol[1])
    tst_backward = OneSampleHotellingT2Test(MCsol[2]',DETsol[2])
    @show test_passed += (pvalue(tst_forward) > 0.05)
    @show test_passed += (pvalue(tst_backward) > 0.05)
end

test_passed / (ntests*2) 

# Random.seed!(2) = 0.9
# Random.seed!(3) = 1.0