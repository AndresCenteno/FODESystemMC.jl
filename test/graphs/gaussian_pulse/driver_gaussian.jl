# HYPERPARAMETERS
EPSILON = sqrt(eps()); NT = 2000; NSIMS = Int(10)
using Pkg; Pkg.activate("../../../.")
using Plots, FODESystemMC, Statistics, Random, StatsBase, DelimitedFiles
include("1D_robin_gaussian.jl")
println(Threads.nthreads())

alphavec = 0.4:0.1:0.9; n = length(alphavec)
det_loss = zeros(n); det_sens = zeros(n)
sto_loss = zeros(3,n); sto_sens = zeros(3,n)

###### true parameters of the problem
# MATRIX
Δx = 0.05; Nt = 20; a1 = -10; a2 = 10; A = mymatrix(Δx,Nt,a1,a2)
# INITIAL CONDITION
u0_vec = myu0(Δx,Nt,0.1,0.025)
# VECTOR OF ALPHAS
true_alpha = 0.7; falpha(α) = α*(sin.(π*Δx.*(1:Nt)) .+1)/4 .+0.5; α = falpha(true_alpha)
# TIME
T = 0.015
# TRUE SOLUTION
problem = FODESystem(A,u0_vec,α,T); init_node = 1
true_sol = L1Solver(problem,init_node;Nt=NT)
################################

Threads.@threads for k in eachindex(alphavec)
    #deterministic stuff
    det_sol = L1Solver(FODESystem(A,u0_vec,falpha(alphavec[k]),T),init_node;Nt=NT); eps_sol = L1Solver(FODESystem(A,u0_vec,falpha(alphavec[k]+EPSILON),T),init_node;Nt=NT)
    det_loss[k] = 0.5*(det_sol-true_sol)^2
    det_sens[k] = (eps_sol - det_sol)*(det_sol-true_sol)/EPSILON
    @show det_loss
end

########## STOCHASTIC
chain_aux = (falpha(0.5) .-0.5)/0.5

function bootstrap_sens(duTdα,uT,true_sol;samples=1000,p=0.05)
    # this function will resample from dudαvec and sto_sol.uT with replacement
    # to create a loss histogram
    n = length(uT)
    loss_quantile = zeros(samples)
    loss_sens = zeros(samples)
    for k = 1:samples
        idx = sample(MersenneTwister(k),1:n,n;replace=true)
        # maybe fix the RNG
        loss_sens[k] = mean(duTdα[idx])*(mean(uT[idx])-true_sol)
        loss_quantile[k] = .5*(mean(uT[idx])-true_sol)^2
    end
    return quantile(loss_quantile,p), quantile(loss_quantile,1-p), 
    quantile(loss_sens,p), quantile(loss_sens,1-p)
end

Threads.@threads for k in eachindex(alphavec)
    @show k
    sto_uT, sto_duTdα, _ = MCSolver(FODESystem(A,u0_vec,falpha(alphavec[k]),T),init_node,chain_aux,QuadLoss();nsims=NSIMS)
    # MC solution
    sto_loss[1,k] = 0.5*(mean(sto_uT)-true_sol)^2
    sto_sens[1,k] = mean(sto_duTdα)*(mean(sto_uT)-true_sol)
    sto_loss[2,k], sto_loss[3,k], sto_sens[2,k], sto_sens[3,k] = bootstrap_sens(sto_duTdα,sto_uT,true_sol)
end
sto_loss[3,:]

writedlm("./alpha/det_loss.csv",det_loss)
writedlm("./alpha/det_sens.csv",det_sens)
writedlm("./alpha/sto_loss.csv",sto_loss)
writedlm("./alpha/sto_sens.csv",sto_sens)

##############################
# belle sourire #            #
# belle sourire #            #
# belle sourire #            #
# belle sourire #            #
# belle sourire #            # 
##############################
