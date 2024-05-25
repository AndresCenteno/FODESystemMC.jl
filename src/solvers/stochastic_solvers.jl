struct SaveSamples end
function MCSolver(problem::FODESystem,init_state::Int,::SaveSamples;nsims::Int=Int(1e6))
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
    uiT = zeros(nsims); duiTdα = zeros(length(α),nsims); duiTdA = zeros(size(A,1),size(A,2),nsims)
    duiTdu0 = zeros(length(α),nsims); duiTdT = zeros(nsims)
    timevec = zeros(nsims)
    for sim=1:nsims
        t1 = time_ns()
        i = copy(init_state) # state
        L = 1; score_α = zero(α); score_A = zero(A);
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

        res = L*u0[i]
        uiT[sim] = res; duiTdα[:,sim] .= score_α*res; duiTdA[:,:,sim] .= score_A*res;
        duiTdu0[i,sim] = L;

        # SPA FOR TIME 
        # (next_value-standing_value)*infinitesimal_probability_of_jumping_to_next
        k = findfirst(rand().<P[i,:])
        duiTdT[sim] = L*(sign(A[i,k])*Q[i]/(-A[i,i])*u0[k]-u0[i])*_score_t(α[i],A[i,i],T-(t-τ))/mittleff_matlab(α[i],A[i,i]*(T-(t-τ))^α[i])
        timevec[sim] = time_ns() - t1
    end
    return forwback(uiT, duiTdA, duiTdu0, duiTdα, duiTdT), timevec
end

#### GARBAGE

# no method Solver
function MCSolver(problem::FODESystem,init_state::Int,method;nsims::Int=Int(1e6))
    if isnothing(method); method = "none"; end
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
    # forward
    uiT = zero(u0[1]);
    # sensitivities
    duiTdα = zero(α); duiTdA = zero(A); duiTdu0 = zero(u0); duiTdT = zero(T)
    res = @distributed _add_forwback for sim=1:nsims
        i = copy(init_state) # state
        L = 1; score_α = zero(α); score_A = zero(A)
        τ = myrand(sojo,i); t = τ
        s1, s2 = score(sojo,i,τ,method,type=:pdf)
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
            s1, s2 = score(sojo,i,τ,method,type=:pdf)
            score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time
            score_A[i,i] -= Ainv[i,i]
        end
        score_α[i] -= s1; score_A[i,i] -= s2; # score for sojourn time
        score_A[i,i] += Ainv[i,i]
        s1, s2 =  _score_cdf(α[i],A[i,i],T-(t-τ))
        score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time

        res = L*u0[i]/nsims
        forwback(res,score_A*res,L/nsims,score_α*res,score(sojo,i,T-(t-τ),method,type=:finaltime)*res)
    end
    return res
end
