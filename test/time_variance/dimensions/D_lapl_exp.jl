using Pkg; Pkg.activate("."); Pkg.instantiate()
using FODESystemMC
using SparseArrays
using Statistics
using LinearAlgebra
 
dimvec = 1:5; nx = 6
using DelimitedFiles

function laplacian_1D(nx)
    Δx = 1/nx
    L = spdiagm(-1=>ones(nx-1)/Δx^2,0=>-2*ones(nx)/Δx^2,1=>ones(nx-1)/Δx^2)
    L[1,end] = nx^2
    L[end,1] = nx^2
    return L
end

function laplacian_matrix(nx,n)
    L = laplacian_1D(nx)
    Ltensor = zeros(nx^n,nx^n)
    for i=1:n
      addL = 1;
      for j=1:n
        if i==j
          addL = kron(addL,L);
        else
          addL = kron(addL,1.0I(nx));
        end
      end
      Ltensor = Ltensor + addL;
    end
    return dropzeros(sparse(Ltensor))
end

function u0(x::Vector) # need to do a proper initial condition here >:(
    res = 1
    for i in eachindex(x)
        res *= sin(x[i]*pi)
    end
    res
end

function eval_u0(d::Int,nx::Int,f::Function)
    deltax = 1/nx
    grid = 0:deltax:1-deltax
    u0_vec = zeros(nx^d)
    mesh = Iterators.product(ntuple(_ -> grid, d)...)
    i = 1
    for point in mesh
        u0_vec[i] = f(collect(point))
        i += 1
    end
    u0_vec
end

# too tired right now, proper mesh in [0,1]^d and evaluate in
# then experiment for random matrices, random initial conditions, but same α, reusing initial
# conditions maybe and matrices for alpha in alphavec

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

simulations = zeros(9,length(dimvec))

Threads.@threads for i in eachindex(dimvec)
    @show d = dimvec[i]
    problem = FODESystem(laplacian_matrix(nx,d),eval_u0(d,nx,u0),ones(nx^d)*0.8,1.)
    MCsol, times = MCSolver(problem,div(nx,2),SaveSamples();nsims=Int(1e3))
    simulations[:,i] = nx_sim(d,MCsol,times)
    println(simulations[:,i])
    writedlm("test/time_variance/augmenting_dims.csv",simulations)
end


  