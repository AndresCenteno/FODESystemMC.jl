using CSV
using DataFrames
using Statistics
using Plots
using Plots.PlotMeasures
using LaTeXStrings
df = CSV.read("test/mlf_execution/mean_execution_time_stats_Tcap_1.0.csv",DataFrame)

begin
p = plot(layout=(1,2),size=(1000,500),
plot_title=L"Expected execution time capping at $T=1$")
plot!(p[1],df.alpha,log.(df.mean_t_pdf),yerror=(log.(df.qpdf005),log.(df.qpdf095)),
xlabel="α",ylabel="log(Mean exec. time log(pdf) (ns))",label=nothing,bottom_margin=50px,
left_margin=20px)

plot!(p[2],df.alpha,log.(df.mean_t_cdf),yerror=(log.(df.qcdf005),log.(df.qcdf095)),
xlabel="α",ylabel="log(Mean exec. time log(cdf) (ns))",label=nothing)
end
savefig("test/mlf_execution/execution_time_Tcap_1.0_log.png")


begin
    p = plot(layout=(1,2),size=(1000,500),
    plot_title=L"Expected execution time capping at $T=1$")
    plot!(p[1],df.alpha,(df.mean_t_pdf),yerror=((df.qpdf005),(df.qpdf095)),
    xlabel="α",ylabel="Mean exec. time log(pdf) (ns)",label=nothing,bottom_margin=50px,
    left_margin=20px)
    
    plot!(p[2],df.alpha,(df.mean_t_cdf),yerror=((df.qcdf005),(df.qcdf095)),
    xlabel="α",ylabel="Mean exec. time log(cdf) (ns)",label=nothing)
end

savefig("test/mlf_execution/execution_time_Tcap_1.0.png")
