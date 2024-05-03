using FODESystemMC
using QuadGK
using MittagLeffler
using HypothesisTests
using Statistics
using ForwardDiff
EPS = sqrt(eps())

p(t,α) = t^(α-1)*mittleff(α,α,-t^α)
function f(t,α) 
    if t < 3; return sin(α*t);
    else # if i dont put this else variance blows up
        return 0.3
    end
end
function dfdα(t,α)
    if t<3;
        return t*cos(α*t);
    else
        return 0
    end
end
g(t,α) = p(t,α)*f(t,α)
Ef(α) = quadgk(t->g(t,α),0,Inf,rtol=1e-8)[1]
dEfdα(α) = (Ef(α+sqrt(eps()))-Ef(α-sqrt(eps())))/(2*sqrt(eps()))

mlf_rng(α,U,V;θ=1) = -abs(θ)^(-1/α)*log.(U).*(sin(α*π)./tan(α*π*V) .-cos(α*π)).^(1/α)
score(t,α) = (log(p(t,α+EPS))-log(p(t,α-EPS)))/(2*EPS)

###
nsims = Int(1e4); α0 = 0.7
U = rand(nsims); V = rand(nsims)
τvec = mlf_rng.(α0,U,V)
forw_sto = f.(τvec,α0)
back_sto = forw_sto.*score.(τvec,α0)+dfdα.(τvec,α0)

OneSampleTTest(forw_sto,Ef(α0))
OneSampleTTest(back_sto,dEfdα(α0)) # this is fucking impossible