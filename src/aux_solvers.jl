struct SaveSamples end

function MCSolver(problem::FODESystem,init_state::Int,::SaveSamples;nsims::Int=Int(1e6))
    @unpack A, u0, α, T = problem
    if !(init_state in 1:length(u0))
        throw(DomainError(init_state,"Choose a node within 1:$(length(u0))"))
    end
    if nsims <= 0
        throw(DomainError(nsims,"Number of simulation must be a positive integer"))
    end
    Q, P = MCDecomposition(A)
    sojo = sojourn(A,α)
    uT_at_i = zeros(nsims)
    n = length(u0);
    duTdαj_at_i = zeros(n,nsims)
    for sim=1:nsims
        # directly update duTdαj_at_i
        L = 1 # signed mass
        i = copy(init_state) # state
        τ = myrand(sojo,i); t = τ
        duTdαj_at_i[i,sim] += score(sojo,i,τ,type=:pdf)
        while t<T
            #TODO
            k = findfirst(rand().<P[i,:]) # new state
            L *= sign(A[i,k])*Q[i,i]/A[i,i] # update mass + sign
            i = k # update state
            τ = myrand(sojo,i); t += τ
            duTdαj_at_i[i,sim] += score(sojo,i,τ,type=:pdf)
        end
        # last score is of survival
        survival = T - (t - τ)
        duTdαj_at_i[i,sim] -= score(sojo,i,survival,type=:pdf)
        duTdαj_at_i[i,sim] += score(sojo,i,survival,type=:cdf)
        # survival!!
        uT_at_i[sim] = L*u0[i]
        duTdαj_at_i[:,sim] *= L*u0[i]
    end
    return uT_at_i, duTdαj_at_i
end