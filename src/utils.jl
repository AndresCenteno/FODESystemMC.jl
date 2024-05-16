using MittagLeffler, ForwardDiff

function MCDecomposition(A::Matrix)
    n = size(A,1); Q = copy(A); Q[1:n+1:end] .= 0; P = copy(abs.(Q)); Q[1:n+1:end] .= -sum(abs.(Q),dims=2)
    # D = diag(A - Q) don't need this as it doesn't appear explicitly
    P = P./ sum(P,dims=2)
    P = cumsum(P,dims=2)
    return Q, P 
end

function MCDecomposition2(A::Matrix)
    n = size(A,1); 
    P = copy(abs.(A));
    P[1:n+1:end] .= 0
    Q = sum(P,dims=2);
    P = P./ Q
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

using FiniteDifferences

function score(S::sojourn,i::Int,τ::Number;type::Symbol,p::Int=2)
    if type == :pdf
        function logpdf(α,Aii) 
            return log(-τ^(α-1)*Aii*mittleff(α,α,Aii*τ^α))
        end
        return grad(central_fdm(p,1),α->logpdf(α,S.diagA[i]),S.α[i])[1], grad(central_fdm(p,1),Aii->logpdf(S.α[i],Aii),S.diagA[i])[1]
    elseif type == :cdf
        function logcdf(α,Aii)
            return log(mittleff(α,Aii*τ^α))
        end
        return grad(central_fdm(p,1),α->logcdf(α,S.diagA[i]),S.α[i])[1], grad(central_fdm(p,1),Aii->logcdf(S.α[i],Aii),S.diagA[i])[1]
    else
        throw("type must be either symbol :pdf or :cdf")
    end
end