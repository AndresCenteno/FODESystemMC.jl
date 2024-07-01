using Pkg; Pkg.activate(".")
using FODESystemMC
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using CSV, DataFrames, Statistics

alpha_vec = collect(0.2:0.2:0.8)
d_vec = 1:5
nx = 4
T = 1.

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

params = Iterators.product(alpha_vec,d_vec)
# somewhere to store the benchmarks
n_exp = length(collect(params))
df = DataFrame(
    alpha = zeros(n_exp),
    d = zeros(n_exp),
    mean = zeros(n_exp),
    median = zeros(n_exp),
    std = zeros(n_exp),
    t005 = zeros(n_exp),
    t095 = zeros(n_exp)
)
p = 0.05 # quantile
nsims = 100000
# not gonna do it in parallel to not fuck up the times?
for j_exp in 1:n_exp
    param = collect(params)[j_exp]
    @show param
    alpha = param[1]; d = param[2]
    A = laplacian_matrix(nx,d); Ainv = 1. ./A
    P, Q = MCDecomposition2(A)
    i = div(nx^d,2)
    u0 = ones(nx^d)
    uiT = 0.; duiTdA = zero(A); duiTdu0 = zero(u0); duiTdα = zero(u0); duiTdT = 0.
    α = ones(nx^d)*alpha
    sojo = sojourn(A,α)
    t = @benchmark one_sim!($A,$Ainv,$u0,$α,T,$P,$Q,$sojo,$i,$uiT,$duiTdA,$duiTdu0,$duiTdα,$duiTdT,nsims)
    df[j_exp,:] = (param[1],param[2],mean(t.times),median(t.times),std(t.times),quantile(t.times,p),quantile(t.times,1-p))
end
display("for loop ended")
CSV.write("benchmarks/simulations/varying_nx/benchtimes_d_alpha.csv",df)