using FODESystemMC, Random

Random.seed!(1)
problem = myrand(randFODESystem(),100)
MCSolver(problem,1,"none";nsims=Int(1e1))
using BenchmarkTools

Random.seed!(1)
res = @btime MCSolver(problem,1,"none";nsims=Int(1e3)) # 5.967 s (4997059 allocations: 507.83 MiB)

Random.seed!(1)
res2 = @btime MCSolver(problem,1,"matlab";nsims=Int(1e3)) # 1.374 s (10313372 allocations: 1.25 GiB)

Random.seed!(1)
res2 = @btime MCSolver(problem,1,FODESystemMC.StronglyTyped();nsims=Int(1e3)) # 320.554 ms (1916044 allocations: 492.04 MiB)

using Profile
@profview MCSolver(problem,1,"none";nsims=Int(1e3)) # most of the time in ^
@profview MCSolver(problem,1,"matlab";nsims=Int(1e3))
@profview MCSolver(problem,1,FODESystemMC.StronglyTyped();nsims=Int(1e3))
################### check results
res = MCSolver(problem,1)
res2 = MCSolver(problem,1,LaplaceInv())