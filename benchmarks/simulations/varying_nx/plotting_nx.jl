using DataFrames
using CSV
using Plots
using LaTeXStrings
using Statistics

df = CSV.read("benchmarks/simulations/varying_nx/benchtimes_nx_alpha.csv",DataFrame)

# que quiero hacer
# plot T final vs T exec for each alpha


function meanslope(x,y)
    round(mean(diff(y)./diff(x)),digits=2)
end

begin
    p = plot()
    for df_aux in groupby(df,:alpha)
        plot!(p,df_aux.Tf,df_aux.mean./1e6,xaxis=:log,yaxis=:log,label=L"\alpha=%$(df_aux.alpha[1]),\Delta \log(y)/\Delta\log(x)\approx %$(meanslope(log.(df_aux.Tf),log.(df_aux.mean)))")
    end
    xlabel!(L"n_x"); ylabel!("Mean execution time (ms)");
    display(p)
end
savefig("benchmarks/simulations/varying_nx/lognx_logexec_slopes_ms.pdf")
