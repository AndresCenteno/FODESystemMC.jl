using Pkg
Pkg.activate(".")

using FODESystemMC
using FLoops # to iterate over product
using SparseArrays # to create laplacian
using Statistics
using DataFrames
using CSV

# need to export JULIA_NUM_THREADS=30
println("Using $(Threads.nthreads()) threads")

function laplacian_1D(nx)
    Δx = 1/nx
    L = spdiagm(-1=>ones(nx-1)/Δx^2,0=>-2*ones(nx)/Δx^2,1=>ones(nx-1)/Δx^2)
    L[1,end] = nx^2
    L[end,1] = nx^2
    return L
end

function std_MCsol(MCsol::forwback)
    # this actually should be done only for the non-zero
    # and because it is not random, I don't know which shape does this have
    uT = std(MCsol.uT)
    dA = mean(nonzeros(sparse(dropdims(std(MCsol.duTdA,dims=3),dims=3))))
    du0 = mean(std(MCsol.duTdu0,dims=2))
    dα = mean(std(MCsol.duTdα,dims=2))
    dT = std(MCsol.duTdT)
    uT, dA, du0, dα, dT
end

alpha_vec = [0.1,0.3,0.5,0.7,0.9]
nx_vec = [2^i for i=2:7]
time_vec = [0.5;1;2;4]
sims_per_param = 1000
params = Iterators.product(alpha_vec,nx_vec,time_vec)
params_vec = collect(params)[:]

nx_vec = zeros(length(params_vec))
alpha_vec = zeros(length(params_vec))
T_vec = zeros(length(params_vec))
mean_t = zeros(length(params_vec))
t005 = zeros(length(params_vec))
t095 = zeros(length(params_vec))
uT_std = zeros(length(params_vec))
dA_std = zeros(length(params_vec))
du0_std = zeros(length(params_vec))
dalpha_std = zeros(length(params_vec))
dT_std = zeros(length(params_vec))

Threads.@threads for i in eachindex(params_vec)
    @show Threads.threadid()
    param = params_vec[i]
    @show param;
    alpha = param[1]; nx = param[2]; T = param[3]
    nx_vec[i] = nx; alpha_vec[i] = alpha; T_vec[i] = T
    problem = FODESystem(laplacian_1D(nx),sin.(pi*range(0,1,nx)),ones(nx)*alpha,T)
    MCsol, times = MCSolver(problem,div(nx,2),SaveSamples();nsims=sims_per_param)
    uT_std1, dA_std1, du0_std1, dalpha_std1, dT_std1 = std_MCsol(MCsol)
    uT_std[i] = uT_std1
    dA_std[i] = dA_std1
    du0_std[i] = du0_std1
    dalpha_std[i] = dalpha_std1
    dT_std[i] = dT_std1
    mean_t[i] = mean(times)
    t005[i] = quantile(times,0.05)
    t095[i] = quantile(times,0.95)
end

df = DataFrame(nx = nx_vec, 
               alpha = alpha_vec, 
               T = T_vec,
               mean_t = mean_t, 
               t005 = t005, 
               t095 = t095, 
               uT_std = uT_std, 
               dA_std = dA_std, 
               du0_std = du0_std,
               dalpha_std = dalpha_std, 
               dT_std = dT_std)

CSV.write("test/time_variance/nx/nx_alpha_T_1D_heat_results.csv",df)