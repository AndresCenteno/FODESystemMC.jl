using SpecialFunctions: gamma

############STOCHASTIC SOLVER

function MCSolver(problem::FODESystem,init_state::Int;nsims::Int=Int(1e6))
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
    for sim=1:nsims
        L = 1 # signed mass
        i = init_state # state
        #@show i
        t = myrand(sojo,i)
        while t < T
            #we are not getting inside??
            k = findfirst(rand().<P[i,:]) # new state
            L *= sign(A[i,k])*Q[i,i]/A[i,i] # update mass + sign
            i = k # update state (I MADE A HUGE MISTAKE HERE BEFORE, I HAD k=i instead of i=k)
            #@show i
            t += myrand(sojo,i)
        end
        uT_at_i += L*u0[i]/nsims
    end
    return uT_at_i
end

#############DETERMINISTIC SOLVER
function L1Solver(problem::FODESystem; Nt::Int=1000)
    @unpack A, u0, α, T = problem
    Δt = T/Nt
    uT = zeros(length(u0),Nt+1); uT[:,1] = u0
    nnodes = length(u0)
    for n=1:Nt
        uT[:,n+1] = Δt*A*uT[:,n]
        for node = 1:nnodes
            for k=1:n-1
                uT[node,n+1] -= L1_weights(α[node],k,Δt*n,Δt)*(uT[node,k+1]-uT[node,k])
            end
            uT[node,n+1] = uT[node,n+1]/L1_weights(α[node],n,Δt*n,Δt)+uT[node,n]
        end
    end
    return uT
end

function L1_weights(α,k,t,Δt)
    # for equispaced nodes in [0,t], 1 is for u[Δt]-u[0]
    return (-(t-Δt*k)^(1-α)+(t-Δt*(k-1))^(1-α))/gamma(1-α)/(1-α)
end