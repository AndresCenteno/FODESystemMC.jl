using DataFrames
using CSV
using Plots
using LaTeXStrings
using Statistics

df = CSV.read("benchmarks/simulations/one_sim_benchmarks_more_times.csv",DataFrame)

# que quiero hacer
# plot T final vs T exec for each alpha

function meanslope(x,y)
    round(mean(diff(y)./diff(x)),digits=2)
end

begin
    p = plot()
    for df_aux in groupby(df,:alpha)
        plot!(p,log.(df_aux.nx),log.(df_aux.mean),label=L"\alpha"*"=$(df_aux.alpha[1]),"*L"\Delta \log(y)/\Delta \log(x)\approx"*"$(meanslope(log.(df_aux.nx),log.(df_aux.mean)))")
    end
    xlabel!(L"T"); ylabel!("Mean exec. time (ms)");
end

begin
    p = plot()
    for df_aux in groupby(df,:alpha)
        plot!(p,(df_aux.nx),(df_aux.mean./1e6),xaxis=:log,yaxis=:log,label=L"\alpha=%$(df_aux.alpha[1]),\Delta \log(y)/\Delta \log(x)\approx %$(meanslope(log.(df_aux.nx),log.(df_aux.mean)))")
    end
    xlabel!(L"T"); ylabel!("Mean execution time (ms)");
end
savefig("benchmarks/simulations/logTF_logexec_slopes_ms.pdf")