using SpecialFunctions: gamma

α = 0.8
k = 100

function f1(z,k,α)
    s = 1
    for i=1:k
        s += z^k/gamma(α*k+1)
    end
    s
end

function f2(z,k,α)
    s = 1
    t = z
    for i=1:k
        s += t/gamma(α*k+1)
        t *= z
    end
    s
end

using BenchmarkTools

@btime f1(0.5,k,α) # 6.552 μs
@btime f2(0.5,k,α) # 4.506 μs