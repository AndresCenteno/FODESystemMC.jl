using Pkg
Pkg.activate(".")

using DataFrames
using CSV
using Statistics
include("../../src/aux/matlab.jl")
# borrowed from Base/timing.jl
T_cap = 1.
macro elapsed_ns(ex)
    quote
        # Experimental.@force_compile
        local t0 = time_ns()
        $(esc(ex))
        (time_ns() - t0)
    end
end
alpha_vec = 0.1:0.1:0.9
nsims = Int(1e6)
t_pdf = [zeros(nsims) for _ in eachindex(alpha_vec)]
t_cdf = [zeros(nsims) for _ in eachindex(alpha_vec)]

@show Threads.nthreads()

# Threads.@threads for i in eachindex(alpha_vec)
for i in eachindex(alpha_vec)
    @show i, Threads.threadid()
    alpha = alpha_vec[i]
    times_pdf = zeros(nsims)
    times_cdf = zeros(nsims)
    for sim = 1:nsims
        τ = Inf
        while τ > T_cap
            τ = -log(rand())*(sin(alpha*π)/tan(alpha*π*rand())-cos(alpha*π))^(1/alpha)
        end
        times_pdf[sim] = @elapsed_ns log(τ^(alpha-1)*mittleff_matlab(alpha,alpha,-τ^alpha))
        times_cdf[sim] = @elapsed_ns log(mittleff_matlab(alpha,-τ^alpha))
    end
    t_pdf[i] = times_pdf
    t_cdf[i] = times_cdf
end
display("For loop exit")

x = DataFrame(alpha=alpha_vec,mean_t_pdf=mean.(t_pdf),mean_t_cdf=mean.(t_cdf),
qpdf005=quantile.(t_pdf,0.05),qpdf095=quantile.(t_pdf,0.95),stdpdf=std.(t_pdf),
qcdf005=quantile.(t_cdf,0.05),qcdf095=quantile.(t_cdf,0.95),stdcdf=std.(t_cdf))
CSV.write("test/mlf_execution/mean_execution_time_stats_Tcap_$(T_cap).csv",x)