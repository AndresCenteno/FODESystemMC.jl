# functions
using SpecialFunctions: gamma
function mymatrix(Δx,Nt,a1,a2)
    # A = -L, -L is the Laplacian
    A = zeros(Nt,Nt)
    A[1:Nt+1:end] .= -2
    A[2:Nt+1:end] .= 1
    A[Nt+1:Nt+1:end] .= 1
    A[1,1] = a2/(a1*Δx-a2)
    A[end,end] = -1
    A = A./(Δx^2)
    return A
end

function myu0(Δx,Nt,μ,σ)
    u0 = zeros(Nt)
    for i=1:Nt
        u0[i] = exp(-(i*Δx-μ)^2/(2*σ^2))/sqrt(2*π*σ^2)
    end
    return u0
end