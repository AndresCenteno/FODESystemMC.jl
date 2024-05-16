using Random, LinearAlgebra, Distributions, ProgressMeter
using StochasticAD
using FODESystemMC
Random.seed!(2024)

nnodes = 10; init_node = rand(1:nnodes); nsims = Int(1e5)
problem = myrand(randFODESystem(),nnodes)

solution = FD_L1Solver(problem,init_node; Nt = 2000)

samples = zeros(nsims);

function get_next_i(i, u)
    return argmax(P[i, :] .> u)
end

# Use https://gaurav-arya.github.io/StochasticAD.jl/stable/devdocs.html#StochasticAD.propagate to handle discrete construct
get_next_i(i::StochasticAD.StochasticTriple, u) = StochasticAD.propagate(get_next_i, i, u)

function primal(α; get_α_bound = (t, i) -> α[i])
    i = copy(init_state); t = 0 # reset simulation

    α = get_α_bound(t, i)
    τ = rand(Exponential(1 / rate_bound))
    t += τ
    η = exp(θ[i]*D[i]*τ)

    while t < T
        true_rate = θ[i]
        # max_rate = true_rate
        @assert true_rate <= rate_bound 
        b = rand(Bernoulli(true_rate / rate_bound))
        # @assert b == 1

        # update i only if b == 1
        i = get_next_i(i, rand()) * b + i * (1 - b)

        # fixed time step
        rate_bound = get_rate_bound(t, i)
        τ = rand(Exponential(1 / rate_bound))
        t += τ

        @assert t isa AbstractFloat # Ensure that t has no dual component when differentiating
        
        # update η
        η *= exp(θ[i]*D[i]*τ)
    end
    η /= exp(θ[i]*D[i]*(t-T)) # subtract reward that was cutted short by end of simulation
    return η*u0[i]
end

primal_unbounded = θ -> primal(θ; get_rate_bound = (t, i) -> θ[i]) # not StochasticAD-compatible
primal_bounded = θ -> primal(θ; get_rate_bound = (t, i) -> 1.0) # is StochasticAD-compatible

@showprogress for sim=1:nsims
    samples[sim] = primal_bounded(θ) 
end

sample_mean = mean(samples) # 1.0582814906302531
sample_sem = std(samples)/sqrt(nsims) # 4.6125383852028576e-5
println(abs(sample_mean - solution) < 1.96*sample_sem) # we are so back in Random.seed!(2024)

## exact diff

# I would like to compute this, of course use other scheme for more precision haha
dsoldθ = zeros(n)
for i=1:n
    θpert = copy(θ); θpert[i] += sqrt(eps(Float64))
    solpert = (exp(diagm(θpert)*(diagm(D)+Q)*T)*u0)[init_state]
    dsoldθ[i] = (solpert-solution)/sqrt(eps(Float64))
end

exact_deriv = dsoldθ

## StochasticAD diff

sample_deriv = @showprogress [derivative_estimate(primal_bounded, θ) for sim in 1:nsims]
sample_deriv_mean = map(i -> mean(d -> d[i], sample_deriv), 1:length(θ))
sample_deriv_sem = map(i -> std(map(d -> d[i], sample_deriv)) / sqrt(nsims), 1:length(θ))

## Compare exact and StochasticAD

using DataFrames
d = DataFrame(; exact_deriv, sample_deriv_mean, sample_deriv_sem)
show(d; display_size=(20, 120))
#=
Row │ exact_deriv  sample_deriv_mean  sample_deriv_sem 
     │ Float64      Float64            Float64          
─────┼──────────────────────────────────────────────────
   1 │ -0.35527          -0.357055          0.00644792
   2 │  0.028692          0.0286772         0.000854623
   3 │ -0.00738235       -0.00777951        0.000659989
   4 │  0.00201763        0.000582919       0.00244559
   5 │  0.00860672        0.00696018        0.000707511
   6 │ -0.0283858        -0.028614          0.000965889
   7 │ -0.0037198        -0.00394399        0.000424466
   8 │ -0.00623535       -0.00546398        0.000493829
   9 │ -0.00565389       -0.00613331        0.000385351
  10 │  0.0187561         0.019082          0.000740168
=#