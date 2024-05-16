struct SaveSamplesNoBranching end

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
    duTdAij_at_i = zeros(n,n,nsims)
    for sim=1:nsims
        # directly update duTdαj_at_i
        L = 1 # signed mass
        i = copy(init_state) # state
        τ = myrand(sojo,i); t = τ
        duTdαj_at_i[i,sim] += score(sojo,i,τ,type=:pdf)
        # should do conditional Monte-Carlo here
        while t < T
            k = findfirst(rand().<P[i,:]) # new state
            L *= sign(A[i,k])*Q[i,i]/A[i,i] # update mass + sign
            duTdAij_at_i[i,k,sim] += 1
            i = k # update state
            τ = myrand(sojo,i); t += τ
            duTdαj_at_i[i,sim] += score(sojo,i,τ,type=:pdf)
        end
        duTdαj_at_i[i,sim] -= score(sojo,i,τ,type=:pdf)
        duTdαj_at_i[i,sim] += score(sojo,i,T-(t-τ),type=:cdf)
        uT_at_i[sim] = L*u0[i]
        duTdαj_at_i[:,sim] *= L*u0[i]
        duTdAij_at_i[:,:,sim] *= 1/(L*u0[i])
    end
    return uT_at_i, duTdαj_at_i, duTdAij_at_i
end

############### WORKING HERE AS IT WILL BE QUICKER

struct NoSave end
function MCSolver(problem::FODESystem,init_state::Int,::NoSave;nsims::Int=Int(1e4))
    @unpack A, u0, α, T = problem
    if !(init_state in 1:length(u0))
        throw(DomainError(init_state,"Choose a node within 1:$(length(u0))"))
    end
    if nsims <= 0
        throw(DomainError(nsims,"Number of simulation must be a positive integer"))
    end
    Q, P = MCDecomposition(A)
    sojo = sojourn(A,α)
    uT_at_i = 0
    n = length(u0);
    duTdαj = zeros(n)
    duTdAij = zeros(n,n)
    for sim=1:nsims
        duTdαj_at_i = zeros(n)
        duTdAij_at_i = zeros(n,n)
        # directly update duTdαj_at_i
        L = 1 # signed mass
        i = copy(init_state) # state
        τ = myrand(sojo,i); t = τ
        duTdαj_at_i[i] += score(sojo,i,τ,type=:pdf)
        duTdAij_at_i[i,i] += diag_score(A[i,i],τ,α[i],type=:pdf)
        # should do conditional Monte-Carlo here
        while t < T
            k = findfirst(rand().<P[i,:]) # new state
            L *= sign(A[i,k])*Q[i,i]/A[i,i] # update mass + sign
            duTdAij_at_i[i,:] += matrix_score(A[i,:],i,k)
            i = k # update state
            τ = myrand(sojo,i); t += τ
            duTdαj_at_i[i] += score(sojo,i,τ,type=:pdf)
            duTdAij_at_i[i,i] += diag_score(A[i,i],τ,α[i],type=:pdf)
        end
        duTdαj_at_i[i] -= score(sojo,i,τ,type=:pdf)
        duTdαj_at_i[i] += score(sojo,i,T-(t-τ),type=:cdf)
        uT_at_i += L*u0[i]/nsims
        duTdαj .+= duTdαj_at_i*L*u0[i]/nsims
        duTdAij_at_i[i,i] -= diag_score(A[i,i],τ,α[i],type=:pdf)
        duTdAij_at_i[i,i] += diag_score(A[i,i],τ,α[i],type=:cdf)
        duTdAij .+= duTdAij_at_i*(L*u0[i])/nsims
    end
    return uT_at_i, duTdαj, duTdAij
end

function matrix_score(Ai::Vector,i::Int,k::Int)
    # we jump from i to k, updatze score for matrix
    # except for i :D
    score_mat = -sign.(Ai)/sum(abs.(Ai)).*(-abs.(Ai))/Ai[i].-sign(Ai[k])/Ai[i]
    score_mat[k] += 1/Ai[k]
    score_mat[i] = 0
    score_mat
end

function diag_score(Aii::Float64,τ::Float64,α::Float64;type)
    if type == :pdf
        return -1/Aii + grad(central_fdm(2,1),Aii->log(-τ^(α-1)*Aii*mittleff(α,α,Aii*τ^α)),Aii)[1]
    elseif type == :cdf
        return grad(central_fdm(2,1),Aii->log(mittleff(α,Aii*τ^α)),Aii)[1]
    end
end

# no method Solver
function MCSolver(problem::FODESystem,init_state::Int;nsims::Int=Int(1e6))
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
    uiT = zero(u0[1]); duiTdα = zero(α); duiTdA = zero(A)
    for sim=1:nsims
        i = copy(init_state) # state
        L = 1; score_α = zero(α); score_A = zero(A)
        τ = myrand(sojo,i); t = τ
        s1, s2 = score(sojo,i,τ,type=:pdf)
        score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time
        score_A[i,i] -= Ainv[i,i]
        while t < T
            k = findfirst(rand().<P[i,:]) # new state
            L *= sign(A[i,k])*Q[i]/(-A[i,i]) # update mass + sign, just need a Q vector, not a matrix
            # matrix is for transition probabilities
            score_A[i,k] += Ainv[i,k]
            i = k # update state
            τ = myrand(sojo,i); t += τ
            s1, s2 = score(sojo,i,τ,type=:pdf)
            score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time
            score_A[i,i] -= Ainv[i,i]
        end
        score_α[i] -= s1; score_A[i,i] -= s2; # score for sojourn time
        score_A[i,i] += Ainv[i,i]
        s1, s2 = score(sojo,i,τ,type=:cdf)
        score_α[i] += s1; score_A[i,i] += s2; # score for sojourn time

        res = L*u0[i]/nsims
        uiT += res; duiTdα += score_α*res; duiTdA += score_A*res;
    end
    return uiT, duiTdα, duiTdA
end