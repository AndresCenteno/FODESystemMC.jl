struct QuadLoss end
function MCSolver(problem::FODESystem,init_state::Int,chain_aux::Vector,::QuadLoss;nsims::Int=Int(1e6))
    @unpack A, u0, α, T = problem
    Q, P = MCDecomposition2(A)
    Ainv = 1 ./A
    sojo = sojourn(A,α)
    uiT = zeros(nsims);
    # personalized loss
    duiTdα = zeros(nsims);
    duiTdA = zeros(nsims)
    for sim=1:nsims
        i = copy(init_state) # state
        L = 1;
        score_α = zero(α);
        score_A = 0
        τ = myrand(sojo,i); t = τ
        s1, s2 = _score_pdf(α[i],A[i,i],τ)
        score_α[i] += s1; # score for sojourn time
        if i==1; score_A += s2; score_A -= Ainv[i,i]; end
        while t < T
            k = findfirst(rand().<P[i,:]) # new state
            L *= sign(A[i,k])*Q[i]/(-A[i,i]) # update mass + sign, just need a Q vector, not a matrix
            # matrix is for transition probabilities
            i = k # update state
            τ = myrand(sojo,i); t += τ
            s1, s2 = _score_pdf(α[i],A[i,i],τ)
            score_α[i] += s1;
            if i==1; score_A += s2; end # score for sojourn time
            if i==1; score_A -= Ainv[i,i]; end
        end
        score_α[i] -= s1;
        if i==1
            score_A -= s2; # score for sojourn time
            score_A += Ainv[i,i]
        end
        s1, s2 =  _score_cdf(α[i],A[i,i],T-(t-τ))
        score_α[i] += s1;
        if i==1; score_A += s2; end

        res = L*u0[i]
        uiT[sim] = res;
        duiTdα[sim] = (chain_aux'*score_α)*res
        duiTdA[sim] = score_A*res;
    end
    return uiT, duiTdα, duiTdA
end