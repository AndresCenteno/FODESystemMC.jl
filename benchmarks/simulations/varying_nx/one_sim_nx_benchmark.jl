using Pkg; Pkg.activate(".")
using FODESystemMC
using SparseArrays

using BenchmarkTools
using CSV, DataFrames, Statistics

alpha_vec = collect(0.2:0.2:0.8)
nx_vec = [2^i for i=1:5]

T = 8.

function laplacian_1D(nx)
    Δx = 1/nx
    L = spdiagm(-1=>ones(nx-1)/Δx^2,0=>-2*ones(nx)/Δx^2,1=>ones(nx-1)/Δx^2)
    L[1,end] = nx^2
    L[end,1] = nx^2
    return L
end

params = Iterators.product(alpha_vec,nx_vec)
# somewhere to store the benchmarks
n_exp = length(collect(params))
df = DataFrame(
    alpha = zeros(n_exp),
    Tf = zeros(n_exp),
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
    alpha = param[1]; nx = param[2]
    A = laplacian_1D(nx); Ainv = 1. ./A
    P, Q = MCDecomposition2(A)
    i = div(nx,2)
    Δx = 1/nx; u0 = sin.(collect(0:Δx:1-Δx).*pi)
    uiT = 0.; duiTdA = zero(A); duiTdu0 = zero(u0); duiTdα = zero(u0); duiTdT = 0.
    α = ones(nx)*alpha
    sojo = sojourn(A,α)
    t = @benchmark one_sim!($A,$Ainv,$u0,$α,T,$P,$Q,$sojo,$i,$uiT,$duiTdA,$duiTdu0,$duiTdα,$duiTdT,nsims)
    df[j_exp,:] = (param[1],param[2],mean(t.times),median(t.times),std(t.times),quantile(t.times,p),quantile(t.times,1-p))
end
display("for loop ended")
CSV.write("benchmarks/simulations/varying_nx/benchtimes_nx_alpha.csv",df)