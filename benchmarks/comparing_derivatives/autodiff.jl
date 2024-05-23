using FODESystemMC, Random
using Zygote
include("../../src/aux/mittag_leffler.jl")
Random.seed!(1)
problem = myrand(randFODESystem(),10)
using BenchmarkTools

Random.seed!(1)
res2 = @btime MCSolver(problem,1,FODESystemMC.StronglyTyped();nsims=Int(1e1)) # 290.527 ms (2056926 allocations: 252.70 MiB)

# maybe first co
Random.seed!(1)
res3 = MCSolver(problem,1,Zyg();nsims=Int(1))


function _score_pdf_auto(α,Aii,τ)
    θ = [α; Aii]
    f(θ)=log(-τ^(θ[1]-1)*θ[2]*mittleff(θ[1],θ[1],θ[2]*τ^θ[1]))
    res = Zygote.gradient(f,θ)[1]
    res[1], res[2]
end
@btime _score_pdf_auto(0.1,-0.1,0.2) # 146.316 μs

function _score_cdf_auto(α,Aii,τ)
    θ = [α; Aii]
    f(θ)=log(mittleff(θ[1],θ[2]*τ^θ[1]))
    res = Zygote.gradient(f,θ)[1]
    res[1], res[2]
end

function _score_t_auto(α,Aii,τ)
    f(τ) = τ^(α-1)*Aii*mittleff(α,α,Aii*τ^α)/mittleff(α,Aii*τ^α)
    return Zygote.gradient(f,τ)[1]
end

