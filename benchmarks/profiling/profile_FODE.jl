using Profile

using FODESystemMC

problem = myrand(randFODESystem(),100)
@profview MCSolver(problem,1;nsims=Int(1e3))