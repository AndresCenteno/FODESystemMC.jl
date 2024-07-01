using Plots
using LaTeXStrings

p = plot(layout=(1,2),size=(1000,400),xlabel="θ",margin=0.5Plots.cm)
ylabel!(p[1],L"\log(U)^2/θ^2")
ylabel!(p[2],L"1(U>-\log(V)/θ)")
plot!(p[1],θ->log(0.3)^2/θ^2,0.1,2,label=nothing)
u = rand(); v = rand()
function f(u,v,θ)
    return u > -log(v)/θ ? 1 : 0
end
plot!(p[2],θ->f(u,v,θ),0.1,2,label=nothing,color=:red)
savefig("test/graphs/reparam.png")
plot(θ->(log(0.3)^2)/θ^2,0.1,3)