# HYPERPARAMETERS
EPSILON = sqrt(eps()); NT = 2000; NSIMS = Int(5e6)

using Plots, FiniteDifferences, FODESystemMC, Statistics, StatsBase
include("1D_robin_gaussian.jl")

alphavec = 0.4:0.1:0.9; n = length(alphavec)
det_loss = zeros(n); det_sens = zeros(n)
sto_loss = zeros(3,n); sto_sens = zeros(3,n)

###### true parameters of the problem
# MATRIX
Δx = 0.05; Nt = 4; a1 = -10; a2 = 10; A = mymatrix(Δx,Nt,a1,a2)
# INITIAL CONDITION
u0_vec = myu0(Δx,Nt,0.1,0.01)
# VECTOR OF ALPHAS
true_alpha = 0.65; falpha(α) = α*(sin.(π*Δx.*(1:Nt)) .+1)/4 .+0.5; α = falpha(true_alpha)
# TIME
T = 0.015
# TRUE SOLUTION
problem = FODESystem(A,u0_vec,α,T); init_node = 1
true_sol = L1Solver(problem,init_node;Nt=NT)
################################

for k in eachindex(alphavec)
    #deterministic stuff
    det_sol = L1Solver(FODESystem(A,u0_vec,falpha(alphavec[k]),T),init_node;Nt=NT); eps_sol = L1Solver(FODESystem(A,u0_vec,falpha(alphavec[k]+EPSILON),T),init_node;Nt=NT)
    det_loss[k] = 0.5*(det_sol-true_sol)^2
    det_sens[k] = (eps_sol - det_sol)*(det_sol-true_sol)/EPSILON
    @show det_loss
end

p = plot(layout=(1,2))  
plot!(p[1],alphavec,det_loss); plot!(p[2],alphavec,det_sens)

########## STOCHASTIC
function chain_rule(M)
    return (M')*(falpha(0.5) .-0.5)/0.5
end

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

for k in eachindex(alphavec)
    @show k
    #deterministic stuff
    sto_sol, _ = MCSolver(FODESystem(A,u0_vec,falpha(alphavec[k]),T),init_node,SaveSamples();nsims=NSIMS)
    sto_loss[1,k] = 0.5*(mean(sto_sol.uT)-true_sol)^2
    # mean(sto_sol).uT ± 1.96*std(sto_sol.uT)/sqrt(NSIMS) this is our error
    # how does it propagate through the loss f(a) = .5(a - c)^2
    # it is ±b(a-c)
    # need the bootstrappppppp sto_loss[2,k] = 1.96*std(chain_rule(sto_sol.uT))/sqrt(NSIMS)*(mean(sto_sol.uT)-true_sol)
    dudαvec = chain_rule(sto_sol.duTdα)
    sto_sens[1,k] = mean(dudαvec)*(mean(sto_sol.uT)-true_sol)
    sto_loss[2,k], sto_loss[3,k], sto_sens[2,k], sto_sens[3,k] = bootstrap_sens(dudαvec,sto_sol.uT,true_sol)
end

##############################
# belle sourire #            #
# belle sourire #            #
# belle sourire #            #
# belle sourire #            #
# belle sourire #            # 
##############################