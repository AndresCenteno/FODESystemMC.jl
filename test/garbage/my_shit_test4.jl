using FODESystemMC
using HypothesisTests
# 1D score function does it work???
using MittagLeffler
mlf_rng(θ,α) = -abs(θ)^(-1/α)*log(rand())*(sin(α*π)/tan(α*π*rand())-cos(α*π))^(1/α)
T = 0.4

f(α) = mittleff(α,-1*T^α)
d1,d2 = f(0.8), (f(0.8+1e-8)-f(0.8))*1e8
nsims = 1000000
g1,g2=0,0
for sim=1:nsims
    g2 = 
end
