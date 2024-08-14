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
        plot!(p,df_aux.d,df_aux.mean./1e6,label=L"\alpha=%$(df_aux.alpha[1]),\Delta \log(y)/\Delta\log(x)\approx %$(meanslope(log.(df_aux.d),log.(df_aux.mean)))")
    end
    xlabel!(L"d"); ylabel!("Mean execution time (ms)");
    display(p)
end
savefig("benchmarks/simulations/varying_nx/logd_logexec_slopes_ms.pdf")
