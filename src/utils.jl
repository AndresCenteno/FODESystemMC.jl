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
    elseif type == :finaltime
        return -(-τ^(S.α[i]-1)*S.diagA[i]*mittleff(S.α[i],S.α[i],S.diagA[i]*τ^S.α[i]))
    else
        throw("type must be either symbol :pdf or :cdf")
    end
end

function compare(forwback_det::forwback,forwback_sto::forwback)
    # want to implement a method for when fields of forwback have an aditional dimension
    if size(forwback_sto.duTdu0,2) == 1
        return _compare(forwback_det,forwback_sto)
    end
    # do hotelling tests
    return _test(forwback_det,forwback_sto)
end

function _compare(forwback_det::forwback,forwback_sto::forwback)
    # want to implement a method for when fields of forwback have an aditional dimension
    rel_err(a,b) = norm(a .- b)/norm(a)
    err_vec = zeros(length(fieldnames(forwback)))
    for i=1:length(err_vec)
        err_vec[i] = rel_err(getfield(forwback_det,i),getfield(forwback_sto,i))
    end
    err_vec
end

using HypothesisTests

function _test(forwback_det::forwback,forwback_sto::forwback;p=0.05)
    @show "Hypothesis testing p = $(p)"
    nsims = length(forwback_sto.uT); n = size(forwback_sto.duTdu0,1)
    test_passed = zeros(Bool,length(fieldnames(forwback)))
    test_passed[1] = pvalue(OneSampleTTest(forwback_sto.uT,forwback_det.uT)) > p
    test_passed[2] = pvalue(OneSampleHotellingT2Test(
        reshape(forwback_sto.duTdA,n^2,nsims)',reshape(forwback_det.duTdA,n^2))) > p
    test_passed[3] = pvalue(OneSampleHotellingT2Test(forwback_sto.duTdu0',forwback_det.duTdu0)) > p
    test_passed[4] = pvalue(OneSampleHotellingT2Test(forwback_sto.duTdα',forwback_det.duTdα)) > p
    test_passed[5] = pvalue(OneSampleTTest(forwback_sto.duTdT,forwback_det.duTdT)) > p
    test_passed
end