using Plots, DelimitedFiles, LaTeXStrings, StatsPlots

a11_detloss = readdlm("./a11/det_loss.csv")
a11_detsens = readdlm("./a11/det_sens.csv")
a11_stoloss = readdlm("./a11/sto_loss.csv")
a11_stosens = readdlm("./a11/sto_sens.csv")

alpha_detloss = readdlm("./alpha/det_loss.csv")
alpha_detsens = readdlm("./alpha/det_sens.csv")
alpha_stoloss = readdlm("./alpha/sto_loss.csv")
alpha_stosens = readdlm("./alpha/sto_sens.csv")

A11vec = -395.95238095238085:5.0:-370.95238095238085
alphavec = 0.4:0.1:0.9

p = plot(layout=(2,2),size=(1000,800),left_margin=5Plots.mm)

xaxis!(p[1],L"\alpha")
xaxis!(p[2],L"\alpha")
xaxis!(p[3],L"A_{11}")
xaxis!(p[4],L"A_{11}")

# yaxis!(p[1],L"\frac{1}{2}(u(T;\alpha)-u(T;\alpha_0))^2")
# yaxis!(p[3],L"\frac{1}{2}(u(T;A[1,1])-u(T;A_0[1,1]))^2")
# yaxis!(p[2],L"\frac{du(T;\alpha)}{d\alpha}(u(T;\alpha)-u(T;\alpha_0))")
# yaxis!(p[4],L"\frac{du(T;A[1,1])}{dA[1,1]}(u(T;A[1,1])-u(T;A_0[1,1]))")

plot!(p[1],alphavec,alpha_detloss,label="L1 Scheme")
plot!(p[1],alphavec,alpha_stoloss[1,:],label="Random walks",
yerr=(abs.(alpha_stoloss[1,:]-alpha_stoloss[2,:]),abs.(alpha_stoloss[1,:]-alpha_stoloss[3,:])))
ylabel!(p[1],L"\mathcal{L}(\alpha)")
ylabel!(p[2],L"\frac{d\mathcal{L}(\alpha)}{d\alpha}")
plot!(p[2],alphavec,alpha_detsens,label="Finite differences")
plot!(p[2],alphavec,alpha_stosens[1,:],label="Random walks",
yerr=(abs.(alpha_stosens[1,:]-alpha_stosens[2,:]),abs.(alpha_stosens[1,:]-alpha_stosens[3,:])))
plot!(p[3],A11vec,a11_detloss,label="L1 Scheme")
plot!(p[3],A11vec,a11_stoloss[1,:],label="Random walks",
yerr=(abs.(a11_stoloss[1,:]-a11_stoloss[2,:]),abs.(a11_stoloss[1,:]-a11_stoloss[3,:])))
ylabel!(p[3],L"\mathcal{L}(A_{11})")
ylabel!(p[4],L"\frac{d\mathcal{L}(A_{11})}{dA_{11}}")
plot!(p[4],A11vec,a11_detsens,label="Finite differences")
plot!(p[4],A11vec,a11_stosens[1,:],label="Random walks",
yerr=(abs.(a11_stosens[1,:]-a11_stosens[2,:]),abs.(a11_stosens[1,:]-a11_stosens[3,:])))

savefig("all_graphs_ylabel.pdf")