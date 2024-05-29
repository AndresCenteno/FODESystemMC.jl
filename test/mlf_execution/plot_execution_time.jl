using DataFrames
using Statistics
using CSV
using Plots

df = CSV.read("test/mlf_execution/mean_execution_time.csv",DataFrame)
df.t_pdf = map(split.(strip.(df.t_pdf, Ref(['[', ']'])), ',')) do nums
    parse.(Float64, nums)
end
df.t_cdf = map(split.(strip.(df.t_cdf, Ref(['[', ']'])), ',')) do nums
    parse.(Float64, nums)
end

begin
    p = plot(layout=(1,2),size=(1000,500),plot_title="1e4 simulations, unit rate")
    quant = 0.05
    plot!(p[1],df.alpha,log.(mean.(df.t_pdf)),
    # yerror=(log.(quantile.(df.t_pdf,quant)),log.(quantile.(df.t_pdf,1-quant))),
    xlabel="α", ylabel = "log(Mean execution time of log(pdf) (ns))")
    plot!(p[2],df.alpha,log.(mean.(df.t_cdf)),
    # yerror=(log.(quantile.(df.t_pdf,quant)),log.(quantile.(df.t_pdf,1-quant))),
    xlabel="α", ylabel = "log(Mean execution time of log(cdf) (ns))")
    display(p)
end