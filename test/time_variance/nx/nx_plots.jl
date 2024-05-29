using DataFrames
using CSV
using Plots
using LaTeXStrings
using Statistics

df = CSV.read("test/time_variance/nx/nx_alpha_T_1D_heat_results.csv",DataFrame)

grouped_df = groupby(df, :alpha)
Tvec = unique(df[!,:T])
alpha_vec = unique(df[!,:alpha]); n_alpha = length(alpha_vec)
begin
    p = plot(layout=(2,2),size=(1000,1000),plot_title="Mean of 1000 random walks")
    for i in eachindex(Tvec)
        slope = zeros(n_alpha)
        Tf = Tvec[i]
        for j in eachindex(alpha_vec)
            df = grouped_df[j]
            aux = filter(:T => T-> T == Tf, df)
            plot!(p[i],log.(aux[!,:nx]),log.(aux[!,:mean_t]),
            label="α=$(aux[!,:alpha][1])",xlabel=L"\log(n_x)",ylabel="log(Mean execution time (ns))")
            slope[j] = mean_slope(log.(aux[!,:nx]),log.(aux[!,:mean_t]))
        end
        title!(p[i],"Tf = $(Tf), logslope = $(round(mean(slope),digits=2))")
    end
    display(p)
end
savefig("test/time_variance/nx/nx_vs_simtime.png")

T_df = groupby(df, :nx)[end]
begin
    q = plot(layout=(1,2),size=(1000,500))
    for df in groupby(T_df,:alpha)
        plot!(q[1],df[!,:T],log.(df[!,:mean_t]),
        # yerror=(log.(df[!,:t005]),log.(df[!,:t095])),
        label="$(df[1,:alpha])",xlabel=L"T_f",title="1000 runs, 128 space nodes",ylabel="log(Mean execution time (ns))")
        plot!(q[2],df[!,:T],log.(df[!,:mean_t]),
        yerrors=(log.(df[!,:t005]),log.(df[!,:t095])),
        label="$(df[1,:alpha])",xlabel=L"T_f",title="0.05 quantiles added",ylabel="log(Mean execution time (ns))")
    end
    display(q)
end
savefig("test/time_variance/nx/tf_vs_exectime_nx_128.png")

function mean_slope(x,y)
    mean((y[2:end] .- y[1:end-1])./(x[2:end]-x[1:end-1]))
end