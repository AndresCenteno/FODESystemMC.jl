using SharedArrays

function _add_forwback(s1::forwback,s2::forwback)
    # reduce function, just fieldwise addition
    forwback(s1.uT+s2.uT,s1.duTdA+s2.duTdA,s1.duTdu0+s2.duTdu0,s1.duTdα+s2.duTdα,s1.duTdT+s2.duTdT)
end

struct SharedSamples end

function MCSolver(problem::FODESystem,init_state::Int,::SharedSamples;nsims::Int=Int(1e6))
    # preprocessing
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

    # allocation
    uiT = SharedArray{Float64}(nsims)
    duiTdα = SharedArray{Float64}(length(α),nsims)
    duiTdA = SharedArray{Float64}(size(A,1),size(A,2),nsims)
    duiTdu0 = SharedArray{Float64}(length(α),nsims)
    duiTdT = SharedArray{Float64}(nsims)

    @distributed for sim=1:nsims # no reduction
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
    end
    return forwback(uiT, duiTdA, duiTdu0, duiTdα, duiTdT)
end