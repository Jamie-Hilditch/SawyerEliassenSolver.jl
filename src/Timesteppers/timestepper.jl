const DEFAULT_C::Float64 = 17 / 14

struct DIRKNCoefficients{T}
    a₁₁::T
    a₂₁::T
    a₂₂::T
    b₁::T
    b₂::T
    b₁ᵗ::T
    b₂ᵗ::T
    c₁::T
    c₂::T

    function DIRKNCoefficients(c::T) where {T}
        denom = 4 * (3 * c^2 - 3 * c + 1)

        a₁₁ = c^2 / 2
        a₂₁ = -2 * (9 * c^4 - 9 * c^3 + 3 * c - 1) / (9 * (2 * c - 1)^2)
        a₂₂ = c^2 / 2
        b₁ = (1 - c) / denom
        b₂ = (3 * c - 1) * (2 * c - 1) / denom
        b₁ᵗ = 1 / denom
        b₂ᵗ = 3 * (2 * c - 1)^2 / denom
        c₁ = c
        c₂ = (3 * c - 2) / (3 * (2 * c - 1))
        return new{T}(a₁₁, a₂₁, a₂₂, b₁, b₂, b₁ᵗ, b₂ᵗ, c₁, c₂)
    end
end

DIRKNCoefficients(::Type{T}) where {T} = DIRKNCoefficients(convert(T, DEFAULT_C))

"""Auxillary variables for computing ζⁿ⁺¹ and ζₜⁿ⁺¹"""
struct AuxillaryVariables{T}
    ζⁿ⁺ᶜ¹::FSVariable{T}
    ζⁿ⁺ᶜ²::FSVariable{T}
    tmp::FSVariable{T}
    rhs::FSVariable{T}
end

function AuxillaryVariables(domain::Domain{T}) where {T}
    ζⁿ⁺ᶜ¹ = FSVariable(domain)
    ζⁿ⁺ᶜ² = FSVariable(domain)
    tmp = FSVariable(domain)
    rhs = FSVariable(domain)
    return AuxillaryVariables{T}(ζⁿ⁺ᶜ¹, ζⁿ⁺ᶜ², tmp, rhs)
end

struct WorkingVariables{T}
    FC::FCVariable{T}
    XC::XCVariable{T}
    XS::XSVariable{T}
    XZ::XZVariable{T}
end

function WorkingVariables(
    domain::Domain{T},
    fc_array::Matrix{Complex{T}},
    xsc_array::Matrix{T},
    xz_array::Matrix{T},
) where {T}
    size(fc_array) == size(domain.spectral) &&
        size(xsc_array) == size(xz_array) == size(domain.grid) ||
        throw(ArgumentError("arrays are not compatible size with `domain`"))
    return WorkingVariables{T}(
        FCVariable(domain, fc_array),
        XCVariable(domain, xsc_array),
        XSVariable(domain, xsc_array),
        XZVariable(domain, xz_array),
    )
end

"""$(TYPEDEF)
Object that stores all the variables and operators required to advance a problem one timestep.
"""
struct Timestepper{T}
    problem::Problem{T}
    h::T
    𝓒::DIRKNCoefficients{T}
    auxillary_variables::AuxillaryVariables{T}
    working_variables::WorkingVariables{T}
    𝓛!::SawyerEliassenOperator!{T}
    𝓛ᴵ!::ImplicitSawyerEliassenOperator!{T}
    cgs::ConjugateGradientSolver{T}
    𝓟::AbstractPreconditioner{T}

    function Timestepper(
        problem::Problem{T},
        h::T,
        𝓒::DIRKNCoefficients{T},
        cgs::ConjugateGradientSolver{T},
        𝓟::AbstractPreconditioner{T},
    ) where {T}
        consistent_domains(problem, cgs, 𝓟) ||
            throw(ArgumentError("`problem`, `cgs` and `𝓟` must have the same domain."))

        domain = get_domain(problem)
        auxillary_variables = AuxillaryVariables(domain)

        # create working arrays for working_variable and 𝓛!
        fsc_array = zeros(Complex{T}, size(domain.spectral))
        xsc_array = zeros(T, size(domain.grid))
        xz_array = zeros(T, size(domain.grid))
        𝓛ζ_array = zeros(T, size(domain.grid))

        working_variables = WorkingVariables(domain, fsc_array, xsc_array, xz_array)
        𝓛! = SawyerEliassenOperator!(problem, fsc_array, xsc_array, xz_array, 𝓛ζ_array)
        𝓛ᴵ! = ImplicitSawyerEliassenOperator!(𝓒.a₁₁, h, 𝓛!)

        return new{T}(
            problem, h, 𝓒, auxillary_variables, working_variables, 𝓛!, 𝓛ᴵ!, cgs, 𝓟
        )
    end
end

function Timestepper(
    problem::Problem{T},
    h::T;
    c=nothing,
    cg_max_iterations=nothing,
    cg_tol=nothing,
    𝓟=nothing,
) where {T}
    domain = get_domain(problem)
    𝓒 = isnothing(c) ? DIRKNCoefficients(T) : DIRKNCoefficients(c)
    cgs = ConjugateGradientSolver(domain, cg_max_iterations, cg_tol)
    if isnothing(𝓟)
        𝓟 = IdentityPreconditioner(domain)
    end
    return Timestepper(problem, h, 𝓒, cgs, 𝓟)
end

function Base.show(io::IO, ::MIME"text/plain", ts::Timestepper)
    return print(
        io,
        "Problem:\n",
        "  ├───────── problem: $(summary(ts.problem))\n",
        "  ├──────── timestep: h = $(sfmt(ts.h))\n",
        "  ├─────────────── 𝓒: $(summary(ts.𝓒))\n",
        "  ├───────────── cgs: $(summary(ts.cgs))\n",
        "  └─────────────── 𝓟: $(summary(ts.𝓟))\n",
    )
end

function Base.summary(io::IO, ts::Timestepper)
    return print(io, "Timestepper with timestep $(sfmt(ts.h))")
end
