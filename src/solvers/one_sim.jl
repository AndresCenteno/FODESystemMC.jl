function one_sim!(A,Ainv,u0::Vector,α::Vector,T::Float64,
    P::Matrix,Q,sojo::sojourn,i::Int,uiT,duiTdA,duiTdu0,duiTdα,duiTdT,nsims)
    score_α = zero(α); score_A = zero(A)
    L = 1;
    τ = myrand(sojo,i); t = τ
    s1, s2 = _score_pdf(α[i],A[i,i],τ)
    score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time
    score_A[i,i] -= Ainv[i,i]
    while t < T
        k = findfirst(rand().<P[i,:]) # new state
        L *= sign(A[i,k])*Q[i]/(-A[i,i]) # update mass + sign, just need a Q vector, not a matrix
        # matrix is for transition probabilities
        score_A[i,k] += Ainv[i,k]
        i = k # update state
        τ = myrand(sojo,i); t += τ
        s1, s2 = _score_pdf(α[i],A[i,i],τ)
        score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time
        score_A[i,i] -= Ainv[i,i]
    end
    score_α[i] -= s1; score_A[i,i] -= s2; # score for sojourn time
    score_A[i,i] += Ainv[i,i]
    s1, s2 =  _score_cdf(α[i],A[i,i],T-(t-τ))
    score_α[i] += s1; score_A[i,i] += s2;

    # updating stuff
    res = L*u0[i]
    uiT += res/nsims
    duiTdA .+= score_A*res/nsims
    duiTdu0[i] += L/nsims
    duiTdα .+= score_α*res/nsims
    k = findfirst(rand().<P[i,:])
    duiTdT += L*(sign(A[i,k])*Q[i]/(-A[i,i])*u0[k]-u0[i])*_score_t(α[i],A[i,i],T-(t-τ))/mittleff_matlab(α[i],A[i,i]*(T-(t-τ))^α[i])/nsims
end