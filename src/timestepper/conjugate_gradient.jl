const CG_TOL_DEFAULT :: Float64 = 10^-10 #10^-5

struct ConjugateGradientSolver! 
    p :: Array{ComplexF64,2}
    q :: Array{ComplexF64,2}
    LHS! :: LHS!
    tol :: Float64

    function ConjugateGradientSolver!(p,q,LHS!,cg_tol = CG_TOL_DEFAULT)
        cg_tol > 0 || throw(DomainError(cg_tol,"tolerance must be positive"))
        new(p,q,LHS!,cg_tol)
    end
end

function (cg::ConjugateGradientSolver!)(fi,r)
    # extract variables from cg
    p,q,LHS! = cg.p, cg.q, cg.LHS!

    # termination condition
    condition = cg.tol*real(dot(r,r))

    # setup
    LHS!(q,fi)
    @. r -= q
    @. p = r

    for _ = 1:200
        LHS!(q,p)
        s = dot(r,r)
        α = s/dot(p,q)
        @. fi += α*p
        @. r -= α*q
        t = dot(r,r)
        if real(t) < condition; break; end
        @. p = r + t/s*p
    end
    return nothing
end