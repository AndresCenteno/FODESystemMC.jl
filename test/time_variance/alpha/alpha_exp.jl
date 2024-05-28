using Pkg
Pkg.activate(".")
using FODESystemMC
using SparseArrays
using Statistics
println(Threads.nthreads())
using DelimitedFiles

nexps = 40
alpha_vec = range(0.1,0.9,9)
Tf = 20

struct nx_sim
    uT_var
    duTdA_var
    duTdα_var
    duTdu0_var
    duTdT_var
    time_mean
    time_quantiles
end

function nx_sim(nx,MCsol::forwback,times::Vector)
    n = size(MCsol.duTdu0,1); nsims = length(MCsol.uT)
    duTdA_var = mean(var(reshape(MCsol.duTdA,n^2,nsims),dims=2))
    duTdα_var = mean(var(MCsol.duTdα,dims=2))
    duTdu0_var = mean(var(MCsol.duTdu0,dims=2))
    duTdT_var = var(MCsol.duTdT)
    [nx;var(MCsol.uT);duTdA_var;duTdα_var;duTdu0_var;duTdT_var;mean(times);quantile(times,0.05);quantile(times,0.95)]
end

simulations = zeros(9,length(alpha_vec),nexps)
nnodes = 5
Threads.@threads for i in 1:nexps
    try
    @show i

    problem = myrand(randFODESystem(),nnodes)
    for j in eachindex(alpha_vec)
        problem_aux = FODESystem(problem.A,problem.u0,alpha_vec[j]*ones(nnodes),Tf)
        MCsol, times = MCSolver(problem,1,SaveSamples();nsims=Int(1e2))
        simulations[:,j,i] = nx_sim(alpha_vec[j],MCsol,times)
    end
    println(simulations[:,:,i])
    catch
    end
    writedlm("test/time_variance/alpha/alpha_Tf.csv",simulations)
end
