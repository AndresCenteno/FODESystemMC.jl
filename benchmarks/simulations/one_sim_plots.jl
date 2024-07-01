using DataFrames
using CSV
using Plots
using LaTeXStrings
using Statistics

df = CSV.read("benchmarks/simulations/one_sim_benchmarks_more_times.csv",DataFrame)

# que quiero hacer
# plot T final vs T exec for each alpha

begin
    p = plot()
    for df_aux in groupby(df,:alpha)
        plot!(p,log.(df_aux.Tf),log.(df_aux.mean),label="α=$(df_aux.alpha[1]), Δy/Δx≈$(meanslope(log.(df_aux.Tf),log.(df_aux.mean)))")
    end
    xlabel!(L"\log(T_F)"); ylabel!("log(Mean exec. time log(ns))"); title!("Exec. time scaling one random walk")
    display(p)
end
savefig("benchmarks/simulations/logTF_logexec_slopes.png")

function meanslope(x,y)
    round(mean(diff(y)./diff(x)),digits=2)
end