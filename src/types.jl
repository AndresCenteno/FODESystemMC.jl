struct FODESystem
    A::Matrix
    u0::Vector
    α::Vector
    T::Number
    function FODESystem(A::Matrix,u0::Vector,α::Vector,T::Number)
        # checks that result in error, type checks are done by Julia because I typed the parameters in the definition
        if size(A,1) != size(A,2)
            throw(DomainError(A,"Matrix needs to be squared."))
        end
        if size(A,1) != length(u0)
            throw(DomainError(u0,"Initial condition needs to be same size as matrix."))
        end
        if any(α.<0.1) || any(α.>=1)
            throw(DomainError(α,"Fractional coefficients must lie between (0.1,1), higher than 0.1 due to stability."))
        end
        if any(diag(A).>=0)
            throw(DomainError(diag(A),"Diagonal elements of matrix must be negative."))
        end
        if T <= 0
            throw(DomainError(T,"Time must be positive."))
        end
        # check that results in warning
        A2 = copy(A); n = size(A,1); A2[1:n+1:end] .= 0
        if any(sum(A2,dims=2)./diag(A) .>= 1)
            warn("Matrix is not diagonally dominant, solution might blow up.") 
        end
        new(A,u0,α,T)
    end
end

struct sojourn
    diagA::Vector # of rates
    α::Vector # of fractional coefficients
    # not writing checks because this should be checked first at FODESystem
    function sojourn(A::Matrix,α::Vector)
        diagA = copy(diag(A))
        new(diagA,α)
    end
end