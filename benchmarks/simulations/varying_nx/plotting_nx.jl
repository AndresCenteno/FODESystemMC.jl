using DataFrames
using CSV
using Plots
using LaTeXStrings
using Statistics

df = CSV.read("benchmarks/simulations/varying_nx/benchtimes_d_alpha.csv",DataFrame)

# que quiero hacer
# plot T final vs T exec for each alpha


function meanslope(x,y)
    round(mean(diff(y)./diff(x)),digits=2)
end

begin
    p = plot()
    for df_aux in groupby(df,:alpha)
        plot!(p,df_aux.d,df_aux.mean,label="α=$(df_aux.alpha[1]), Δy/Δx≈$(meanslope(log.(df_aux.d),log.(df_aux.mean)))")
    end
    xlabel!(L"d"); ylabel!("(Mean exec. time log(ns))"); title!("Exec. time scaling one random walk")
    display(p)
end
savefig("benchmarks/simulations/varying_nx/logd_logexec_slopes.png")
