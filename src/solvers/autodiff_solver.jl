using Zygote

struct Zyg end

function MCSolver(problem::FODESystem,init_state::Int,::Zyg;nsims::Int=Int(1e6))
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
        display("entered here")
        i = copy(init_state) # state
        L = 1; score_α = zero(α); score_A = zero(A)
        τ = myrand(sojo,i); t = τ
        s1, s2 = _score_pdf_auto(α[i],A[i,i],τ)
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
            s1, s2 = _score_pdf_auto(α[i],A[i,i],τ)
            score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time
            score_A[i,i] -= Ainv[i,i]
            display("inside while")
        end
        score_α[i] -= s1; score_A[i,i] -= s2; # score for sojourn time
        score_A[i,i] += Ainv[i,i]
        s1, s2 =  _score_cdf_auto(α[i],A[i,i],τ)
        score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time

        res = L*u0[i]/nsims
        forwback(res,score_A*res,L/nsims,score_α*res,_score_t_auto(α[i],A[i,i],τ)*res)
    end
    return res
end

function _score_pdf_auto(α,Aii,τ)
    θ = [α; Aii]
    f(θ)=log(-τ^(θ[1]-1)*θ[2]*mittleff(θ[1],θ[1],θ[2]*τ^θ[1]))
    res = Zygote.gradient(f,θ)[1]
    res[1], res[2]
end

function _score_cdf_auto(α,Aii,τ)
    θ = [α; Aii]
    f(θ)=log(mittleff(θ[1],θ[2]*τ^θ[1]))
    res = Zygote.gradient(f,θ)[1]
    res[1], res[2]
end

function _score_t_auto(α,Aii,τ)
    f(τ) = τ^(α-1)*Aii*mittleff(α,α,Aii*τ^α)/mittleff(α,Aii*τ^α)
    return Zygote.gradient(f,τ)[1]
end

