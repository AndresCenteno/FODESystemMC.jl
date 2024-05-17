using Distributed
addprocs(7)
rmprocs(8:15)
workers()
@everywhere begin
using FODESystemMC, Random, Statistics

Random.seed!(0)
nnodes = 10; init_node = 1
problem = myrand(randFODESystem(),nnodes)
end

nsims = Int(1e2)
MCsol = MCSolver(problem,init_node;nsims=nsims)