"""Abstract supertype for conjugate gradient preconditioners."""
abstract type AbstractPreconditioner{T} end

@inline function apply_preconditioner!(
    P::AbstractPreconditioner{T}, out::FSVariable{T}, in::FSVariable{T}, aᵢᵢ::T, h::T
) where {T}
    @boundscheck consistent_domains(P, out, in) ||
        throw(ArgumentError("`P`, `out` and `in` must have the same domains."))
    @inbounds solve_preconditioner_equation!(P, out, in, aᵢᵢ, h)
end

Domains.get_domain(P::AbstractPreconditioner) = P.domain

"""
`solve_preconditioner_equation(P::AbstractPreconditioner{T}, z::FSVariable{T}, r::FSVariable{T}, aᵢᵢ::T, h::T)::Nothing where {T}`

Function that solves the preconditioner equation `Mz = r` where `M ≈ I + aᵢᵢ h² 𝓛`.
User defined preconditioners must implement this function with the given signature.
It may be assumed that `P`, `z` and `r` have the same domain.
The output `z` should be updated in-place and the input `r` must not be modified.
"""
function solve_preconditioner_equation! end

struct IdentityPreconditioner{T} <: AbstractPreconditioner{T}
    domain::Domain{T}
end

@inline function solve_preconditioner_equation!(
    ::IdentityPreconditioner{T}, out::FSVariable{T}, in::FSVariable{T}, aᵢᵢ::T, h::T
) where {T}
    @inbounds out .= in
    return nothing
end

"""
    $(TYPEDEF)

A preconditioner that approximates `𝓛ζ` in spectral space by
``𝓛ζ ≈ (ω₀² k_x² + ω₁² k_z^2) / (k_x² + k_z²) ζ``

# Fields
$(TYPEDFIELDS)
"""
struct DiagonalQuadraticPreconditioner{T} <: AbstractPreconditioner{T}
    domain::Domain{T}
    ω₀²::T
    ω₁²::T
end

@inline function solve_preconditioner_equation!(
    P::DiagonalQuadraticPreconditioner{T},
    out::FSVariable{T},
    in::FSVariable{T},
    aᵢᵢ::T,
    h::T,
) where {T}
    (; domain, ω₀², ω₁²) = P
    CNX = domain.spectral.CNX
    CNZ = domain.spectral.CNZ
    kx = domain.spectral.kx
    kz = domain.kz
    @inbounds @. out[1:CNX, 1:CNZ] =
        in[1:CNX, 1:CNZ] / (1 + aᵢᵢ * h^2 * (ω₀² * kx^2 + ω₁² * kz^2) / (kx^2 + kz^2))
    @inbounds @. out[(CNX + 1):end, :] = 0
    @inbounds @. out[1:CNX, (CNZ + 1):end] = 0
    return nothing
end
