struct StronglyTyped end

function MCSolver(problem::FODESystem,init_state::Int,::StronglyTyped;nsims::Int=Int(1e6))
    @unpack A, u0, α, T = problem
    if !(init_state in 1:length(u0))
        throw(DomainError(init_state,"Choose a node within 1:$(length(u0))"))
    end
    if nsims <= 0
        throw(DomainError(nsims,"Number of simulation must be a positive integer"))
    end
    Q, P = MCDecomposition2(A)
    Ainv = 1 ./A
    sojo = sojourn(A,α)
    res = @distributed _add_forwback for sim=1:nsims
        i = copy(init_state) # state
        L = 1; score_α = zero(α); score_A = zero(A)
        τ = myrand(sojo,i); t = τ
        s1, s2 = _score_pdf(α[i],A[i,i],τ)
        score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time
        score_A[i,i] -= Ainv[i,i]
        while t < T
            # I forgot the score for the probability transition lmao
            k = findfirst(rand().<P[i,:]) # new state
            L *= sign(A[i,k])*Q[i]/(-A[i,i]) # update mass + sign, just need a Q vector, not a matrix
            # matrix is for transition probabilities
            score_A[i,k] += Ainv[i,k]
            score_A[i,[1:i-1;i+1:end]] .-= sign.(A[i,[1:i-1;i+1:end]])./Q[i]
            score_A[i,k] += Ainv[i,k]
            i = k # update state
            τ = myrand(sojo,i); t += τ
            s1, s2 = _score_pdf(α[i],A[i,i],τ)
            score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time
            score_A[i,i] -= Ainv[i,i]
        end
        score_α[i] -= s1; score_A[i,i] -= s2; # score for sojourn time
        score_A[i,i] += Ainv[i,i]
        s1, s2 =  _score_cdf(α[i],A[i,i],τ)
        score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time

        res = L*u0[i]/nsims
        forwback(res,score_A*res,L/nsims,score_α*res,_score_t(α[i],A[i,i],τ)*res)
    end
    return res
end

function _score_pdf(α,Aii,τ; ϵ = sqrt(eps()))
    f(α,Aii)=log(-τ^(α-1)*Aii*mittleff_matlab(α,α,Aii*τ^α))
    sol = f(α,Aii)
    return (f(α+ϵ,Aii)-sol)/ϵ, (f(α,Aii+ϵ)-sol)/ϵ
end

function _score_cdf(α,Aii,τ; ϵ = sqrt(eps()))
    f(α,Aii)=log(mittleff_matlab(α,Aii*τ^α))
    sol = f(α,Aii)
    return (f(α+ϵ,Aii)-sol)/ϵ, (f(α,Aii+ϵ)-sol)/ϵ
end

function _score_t(α,Aii,τ;ϵ = sqrt(eps()))
    f(τ) = τ^(α-1)*Aii*mittleff_matlab(α,α,Aii*τ^α)/mittleff_matlab(α,Aii*τ^α)
    return (f(τ+ϵ)-f(τ-ϵ))/(2*ϵ)
end

