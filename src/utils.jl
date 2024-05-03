using MittagLeffler, ForwardDiff

function MCDecomposition(A::Matrix)
    n = size(A,1); Q = copy(A); Q[1:n+1:end] .= 0; P = copy(abs.(Q)); Q[1:n+1:end] .= -sum(abs.(Q),dims=2)
    # D = diag(A - Q) don't need this as it doesn't appear explicitly
    P = P./ sum(P,dims=2)
    P = cumsum(P,dims=2)
    return Q, P 
end

# random FODESystem that passes all checks
struct randFODESystem end
function myrand(::randFODESystem,n::Int;sign=false,no_potential=false)
    if sign==false
        A = randn(n,n)
    else
        A = rand(n,n)
    end
    A[1:n+1:end] .= 0
    if no_potential==true
        A[1:n+1:end] .= -sum(abs.(A),dims=2)
    else
        A[1:n+1:end] .= -sum(abs.(A),dims=2)*(2+rand())
    end
    u0 = rand(n)
    α = rand(n)*0.4 .+ 0.6 # don't want to kill the computer with Nt really high to solve this
    T = rand()
    return FODESystem(A,u0,α,T)
end

function myrand(S::sojourn,i::Int)
    # not writing check for i in 1:n to avoid overhead  
    ret = -abs(S.diagA[i])^(-1/S.α[i])*log(rand())*(sin(S.α[i]*π)/tan(S.α[i]*π*rand())-cos(S.α[i]*π))^(1/S.α[i])
    #@show ret, S.diagA[i], S.α[i]
    return ret
end

function score(S::sojourn,i::Int,τ::Number;type::Symbol,ϵ=sqrt(eps()))
    if type == :pdf
        # @show "i got in here"
        # in sojourn, rate is positive!
        function f(α) 
            τ^(α-1)*S.diagA[i]*mittleff(α,α,-S.diagA[i]*τ^α)
        end
        return (log(f(S.α[i]+ϵ/2))-log(f(S.α[i]-ϵ/2)))/ϵ
    elseif type == :cdf
        function g(α)
            mittleff(α,-S.diagA[i]*τ^α)
        end
        #@show 
        return (log(g(S.α[i]+ϵ/2))-log(g(S.α[i]-ϵ/2)))/ϵ
    else
        throw("type must be either symbol :pdf or :cdf")
    end
    # return ForwardDiff.derivative(α->log(f(α)),S.α[i]) ForwardDiff not working 
    # return (log(f(S.α[i]+ϵ/2))-log(f(S.α[i]-ϵ/2)))/ϵ # Julia doesn't carry the f outside the ifs
end