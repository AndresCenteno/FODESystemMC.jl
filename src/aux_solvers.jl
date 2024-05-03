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
        if t > T
            duTdαj_at_i[i,sim] += score(sojo,i,T,type=:cdf)
        else
            duTdαj_at_i[i,sim] += score(sojo,i,τ,type=:pdf)
            while t<T
                k = findfirst(rand().<P[i,:]) # new state
                L *= sign(A[i,k])*Q[i,i]/A[i,i] # update mass + sign
                i = k # update state
                τ = myrand(sojo,i); t += τ
                if t > T
                    duTdαj_at_i[i,sim] += score(sojo,i,T-(t-τ),type=:cdf)
                else
                    duTdαj_at_i[i,sim] += score(sojo,i,τ,type=:pdf)
                end
            end
        end
        uT_at_i[sim] = L*u0[i]
        duTdαj_at_i[:,sim] *= L*u0[i]
    end
    return uT_at_i, duTdαj_at_i
end

struct SaveSamplesNoBranching end
# this should actually be equivalent
function MCSolver(problem::FODESystem,init_state::Int,::SaveSamplesNoBranching;nsims::Int=Int(1e6))
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
        while t < T
            k = findfirst(rand().<P[i,:]) # new state
            L *= sign(A[i,k])*Q[i,i]/A[i,i] # update mass + sign
            i = k # update state
            τ = myrand(sojo,i); t += τ
            duTdαj_at_i[i,sim] += score(sojo,i,τ,type=:pdf)
        end
        duTdαj_at_i[i,sim] -= score(sojo,i,τ,type=:pdf)
        duTdαj_at_i[i,sim] += score(sojo,i,T-(t-τ),type=:cdf)
        uT_at_i[sim] = L*u0[i]
        duTdαj_at_i[:,sim] *= L*u0[i]
    end
    return uT_at_i, duTdαj_at_i
end