using Random, LinearAlgebra, Distributions
Random.seed!(2024)
# number of states
n = 10
# infinitesimal generator
Q = rand(n,n); Q[1:n+1:end] .= 0; for i=1:n; Q[i,:] ./= sum(Q[i,:]); Q; end; Q[1:n+1:end] .= -1
# embedded Markov chain
P = copy(Q); P[1:n+1:end] .= 0; P = cumsum(P,dims=2)
# rewards
D = -rand(n)
# rates
θ = rand(n)
# final reward
u0 = randn(n)
# initial state
init_state = rand(1:n)
# time
T = 0.5

# deterministic solution
solution = (exp(diagm(θ)*(diagm(D)+Q)*T)*u0)[init_state] # 1.0582327608691373

# stochastic solution
nsims = Int(1e5)
samples = zeros(nsims)

for sim=1:nsims
    i = copy(init_state); t = 0 # reset simulation
    τ = rand(Exponential(1/θ[i])); t += τ
    η = exp(θ[i]*D[i]*τ)
    while t < T
        i = argmax(P[i,:].>rand())
        τ = rand(Exponential(1/θ[i])); t += τ
        η *= exp(θ[i]*D[i]*τ)
    end
    η /= exp(θ[i]*D[i]*(t-T)) # subtract reward that was cutted short by end of simulation
    samples[sim] = η*u0[i]
end

sample_mean = mean(samples)
sample_sem = std(samples)/sqrt(nsims)
println(abs(sample_mean - solution) < 1.96*sample_sem ? "we are so back" : "I should quit the phd") # we are so back in Random.seed!(2024)

# I would like to compute this, of course use other scheme for more precision haha
dsoldθ = zeros(n)
for i=1:n
    θpert = copy(θ); θpert[i] += sqrt(eps(Float64))
    solpert = (exp(diagm(θpert)*(diagm(D)+Q)*T)*u0)[init_state]
    dsoldθ[i] = (solpert-solution)/sqrt(eps(Float64))
end