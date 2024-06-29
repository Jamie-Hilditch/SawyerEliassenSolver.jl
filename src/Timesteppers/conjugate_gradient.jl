const CG_TOL_DEFAULT::Float64 = 10^-10 #10^-5
const MAX_ITERATIONS_DEFAULT::Int = 200

struct ConjugateGradientSolver{T}
    domain::Domain{T}
    p::FSVariable{T}
    q::FSVariable{T}
    z::FSVariable{T}
    max_iterations::Int
    tol::T

    function ConjugateGradientSolver(
        domain::Domain{T}, max_iterations::Int, cg_tol::T
    ) where {T}
        cg_tol > 0 || throw(DomainError(cg_tol, "tolerance must be positive"))
        return new{T}(
            domain,
            FSVariable(domain),
            FSVariable(domain),
            FSVariable(domain),
            max_iterations,
            cg_tol,
        )
    end
end

function ConjugateGradientSolver(domain::Domain{T}, ::Nothing, cg_tol::T) where {T}
    return ConjugateGradientSolver(domain, MAX_ITERATIONS_DEFAULT, cg_tol)
end

function ConjugateGradientSolver(
    domain::Domain{T}, max_iterations::Int, ::Nothing
) where {T}
    return ConjugateGradientSolver(domain, max_iterations, convert(T, CG_TOL_DEFAULT))
end

function ConjugateGradientSolver(domain::Domain{T}, ::Nothing, ::Nothing) where {T}
    return ConjugateGradientSolver(
        domain, MAX_ITERATIONS_DEFAULT, convert(T, CG_TOL_DEFAULT)
    )
end

Domains.get_domain(cgs::ConjugateGradientSolver) = cgs.domain

@inline function solve_implicit_equation!(
    cgs::ConjugateGradientSolver{T},
    𝓛ᴵ!::ImplicitSawyerEliassenOperator!,
    x::FSVariable{T},
    b::FSVariable{T},
    𝓟::AbstractPreconditioner{T},
) where {T}
    # extract variables from cg and A
    (; domain, p, q, z, max_iterations, tol) = cgs
    aᵢᵢ, h = 𝓛ᴵ!.aᵢᵢ, 𝓛ᴵ!.h

    @boundscheck consistent_domains(domain, 𝓛ᴵ!, x, b, 𝓟) || throw(
        ArgumentError("`cgs`, `𝓛ᴵ!`, `x`, `b` and `𝓟` must have the same domain.")
    )

    # termination condition
    condition = tol * real(dot(b, b))
    r = b # use the input array to store the residuals

    # setup
    @inbounds 𝓛ᴵ!(q, x)
    @inbounds @. r -= q # r₀ = b - Ax₀
    @inbounds apply_preconditioner!(𝓟, z, r, aᵢᵢ, h) # Mz₀ = r₀
    @inbounds @. p = z # p₀ = z₀

    for _ in 1:max_iterations
        @inbounds 𝓛ᴵ!(q, p) # qₖ = Apₖ
        s = r ⋅ z # s = rₖᵀzₖ
        α = s / (p ⋅ q) # αₖ = rₖᵀzₖ / pₖᵀApₖ = s / pₖᵀqₖ
        @inbounds @. x += α * p # xₖ₊₁ = xₖ + αₖpₖ
        @inbounds @. r -= α * q # rₖ₊₁ = rₖ - αₖApₖ = rₖ - αₖqₖ
        if real(r ⋅ r) < condition # rₖ₊₁ᵀrₖ₊₁ < tol * bᵀb
            return nothing
        end
        @inbounds apply_preconditioner!(𝓟, z, r, aᵢᵢ, h) # Mzₖ₊₁ = rₖ₊₁
        β = (r ⋅ z) / s # βₖ = rₖ₊₁ᵀzₖ₊₁ / rₖᵀzₖ
        @inbounds @. p = z + β * p # pₖ₊₁ = zₖ₊₁ + βₖpₖ
    end
    @warn "Conjugate gradient solver did not terminate after $(max_iterations) iterations."
    return nothing
end
