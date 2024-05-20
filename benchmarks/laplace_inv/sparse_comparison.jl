using FODESystemMC, Random

Random.seed!(1)
problem = myrand(randFODESystem(),100;sparse=true,density=0.1)

using BenchmarkTools

Random.seed!(1)
res = @btime MCSolver(problem,1,"none";nsims=Int(1e3)) # 1.300 s (2660796 allocations: 172.15 MiB)

Random.seed!(1)
res2 = @btime MCSolver(problem,1,"matlab";nsims=Int(1e3)) # 744.459 ms (5484602 allocations: 581.74 MiB)

Random.seed!(1)
res2 = @btime MCSolver(problem,1,FODESystemMC.StronglyTyped();nsims=Int(1e3)) # 188.632 ms (1141447 allocations: 176.75 MiB)

using Profile
@profview MCSolver(problem,1,"none";nsims=Int(1e3)) # most of the time in ^
@profview MCSolver(problem,1,"matlab";nsims=Int(1e3))
@profview MCSolver(problem,1,FODESystemMC.StronglyTyped();nsims=Int(1e3))
################### check results
res = MCSolver(problem,1)
res2 = MCSolver(problem,1,LaplaceInv())