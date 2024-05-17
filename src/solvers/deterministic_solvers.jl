using SpecialFunctions: gamma

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

function L1Solver(A::Matrix,u0::Vector,α::Vector,T::Float64; Nt::Int=1000)
    return L1Solver(FODESystem(A,u0,α,T);Nt=Nt)
end

function L1Solver(A::Matrix,u0::Vector,α::Vector,T::Float64,init_node::Int; Nt::Int=1000)
    return L1Solver(FODESystem(A,u0,α,T);Nt=Nt)[init_node,end]
end


function L1Solver(problem::FODESystem,init_node::Int; Nt::Int=1000)
    return L1Solver(problem;Nt=Nt)[init_node,end]
end

function L1_weights(α,k,t,Δt)
    # for equispaced nodes in [0,t], 1 is for u[Δt]-u[0]
    return (-(t-Δt*k)^(1-α)+(t-Δt*(k-1))^(1-α))/gamma(1-α)/(1-α)
end

function FD_L1Solver(problem::FODESystem,init_state::Int;Nt::Int=1000,p::Int=2)
    uTD = L1Solver(problem,init_state;Nt=Nt)
    f(A,u0,α,T) = L1Solver(A,u0,α,T,init_state;Nt=Nt)
    duTdAD = grad(central_fdm(p, 1), x->f(x,problem.u0,problem.α,problem.T), problem.A)[1]
    duTdu0D = grad(central_fdm(p, 1), x->f(problem.A,x,problem.α,problem.T), problem.u0)[1]
    duTdαD = grad(central_fdm(p, 1), x->f(problem.A,problem.u0,x,problem.T), problem.α)[1]
    duTdT = grad(central_fdm(p, 1), x->f(problem.A,problem.u0,problem.α,x), problem.T)[1]

    return forwback(uTD,duTdAD,duTdu0D,duTdαD,duTdT)
end