using FODESystemMC
using SparseArrays
using Statistics
 
nxvec = [2^i for i=3:9]

function laplacian_1D(nx)
    Δx = 1/nx
    L = spdiagm(-1=>ones(nx-1)/Δx^2,0=>-2*ones(nx)/Δx^2,1=>ones(nx-1)/Δx^2)
    L[1,end] = nx^2
    L[end,1] = nx^2
    return L
end

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

simulations = zeros(9,length(nxvec))

for i in eachindex(nxvec)
    @show nx = nxvec[i]
    problem = FODESystem(laplacian_1D(nx),sin.(2*pi*range(0,1,nx)),ones(nx)*0.8,1.)
    MCsol, times = MCSolver(problem,div(nx,2),SaveSamples();nsims=Int(5e2))
    simulations[:,i] = nx_sim(nx,MCsol,times)
    println(simulations[:,i])
end
using DelimitedFiles
writedlm("test/time_variance/sims.csv",simulations)