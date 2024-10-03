
"""Abstract super type for all forcings."""
abstract type AbstractForcingFunction{T} end

"""Type representing no forcing."""
struct NoForcing{T} <: AbstractForcingFunction{T}
    domain::Domain{T}
end

"""
    $(TYPEDEF)

Apply pointwise forcing to the Sawyer-Eliassen equation in physical space.
`func` is a function that computes the forcing at a given position and time.
`func` should have signature `func(x::T,z::T,t::T,params::P)::T where {T<:Union{Float32,Float64}}`.

# Fields
$(TYPEDFIELDS)
"""
struct PointwisePhysicalForcing{T,P} <: AbstractForcingFunction{T}
    domain::Domain{T}
    func::Function
    params::P
end

"""
    $(TYPEDEF)

Apply global forcing to the Sawyer-Eliassen equation in physical space.
`func` is a function that updates the entire forcing variable in-place at time `t`.
`func` should have signature `func(F::XZVariable{T},t::T,params::P)::Nothing where {T<:Union{Float32,Float64}}`.

# Fields
$(TYPEDFIELDS)
"""
struct GlobalPhysicalForcing{T,P} <: AbstractForcingFunction{T}
    domain::Domain{T}
    func::Function
    params::P
end

"""
    $(TYPEDEF)

Apply pointwise forcing to the Sawyer-Eliassen equation in spectral space.
`func` is a function that computes the forcing at given wavenumbers and time.
`func` should have signature `func(kx::T,kz::T,t::T,params::P)::complex{T} where {T<:Union{Float32,Float64}}`.

# Fields
$(TYPEDFIELDS)
"""
struct PointwiseSpectralForcing{T,P} <: AbstractForcingFunction{T}
    domain::Domain{T}
    func::Function
    params::P
end

"""
    $(TYPEDEF)

Apply global forcing to the Sawyer-Eliassen equation in spectral space.
`func` is a function that updates the entire forcing variable in-place at time `t`.
`func` should have signature `func(F::FSVariable{T},t::T,params::P)::Nothing where {T<:Union{Float32,Float64}}`.
Note that the `FSVariable{T}` has entries of type `Complex{T}`.

# Fields
$(TYPEDFIELDS)
"""
struct GlobalSpectralForcing{T,P} <: AbstractForcingFunction{T}
    domain::Domain{T}
    func::Function
    params::P
end

"""Union of types representing forcing in physical (XZ) space."""
PhysicalForcing{T} =
    Union{PointwisePhysicalForcing{T,P},GlobalPhysicalForcing{T,Q}} where {T,P,Q}
"""Union of types representing forcing in spectral (FS) space."""
SpectralForcing{T} =
    Union{PointwiseSpectralForcing{T,P},GlobalSpectralForcing{T,Q}} where {T,P,Q}

Domains.get_domain(F::AbstractForcingFunction) = F.domain

@inline function evaluate_physical_forcing!(
    ::NoForcing{T}, XZ::XZVariable{T}, t::T
) where {T}
    @. XZ = 0
    return nothing
end

@inline function evaluate_physical_forcing!(
    F::PointwisePhysicalForcing{T,P}, XZ::XZVariable{T}, t::T
) where {T,P}
    @boundscheck consistent_domains(F, XZ) ||
        throw(ArgumentError("`F` and `XZ` must have the same domain."))
    x, z = gridpoints(F)
    @inbounds @. XZ = F.func(x, z, tuple(t), tuple(F.params))
    return nothing
end

@inline function evaluate_physical_forcing!(
    F::GlobalPhysicalForcing{T,P}, XZ::XZVariable{T}, t::T
) where {T,P}
    @boundscheck consistent_domains(F, XZ) ||
        throw(ArgumentError("`F` and `XZ` must have the same domain."))
    F.func(F, XZ, t, F.params)
    return nothing
end



@inline function evaluate_ζ_forcing!(
    ::NoForcing, out::FSVariable{T}, t::T, ::XSVariable, ::XZVariable
) where {T}
    @. out = 0
    return nothing
end

@inline function evaluate_ζ_forcing!(
    F::PointwiseSpectralForcing{T,P},
    out::FSVariable{T},
    t::T,
    ::XSVariable{T},
    ::XZVariable{T},
) where {T,P}
    @boundscheck consistent_domains(F, out) ||
        throw(ArgumentError("`F` and `out` must have the same domain."))
    kx, kz = wavenumbers_full(F)
    @inbounds @. F.func(kx, kz, tuple(t), tuple(params))
    return nothing
end

@inline function evaluate_ζ_forcing!(
    F::GlobalSpectralForcing{T,P}, out::FSVariable{T}, t::T, ::XSVariable, ::XZVariable
) where {T,P}
    @boundscheck consistent_domains(F, out) ||
        throw(ArgumentError("`F` and `out` must have the same domain."))
    F.func(out, t, F.params)
    return nothing
end

@inline function evaluate_ζ_forcing!(
    F::PhysicalForcing{T}, out::FSVariable{T}, t::T, XS::XSVariable, XZ::XZVariable
) where {T}
    @boundscheck consistent_domains(F, out, XS, XZ) ||
        throw(ArgumentError("`F`, `out`, `XS` and`XZ` must have the same domain."))
    @inbounds evaluate_physical_forcing!(F, XZ, t)
    Tˢ!(XS, XZ)
    Tᴴ!(out, XS)
    return nothing
end
