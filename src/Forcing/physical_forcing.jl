"""
    $(TYPEDEF)

Apply pointwise forcing to the Sawyer-Eliassen equation in physical space.
`func` is a function that computes the forcing at a given position and time.
`func` should have signature `func(x::T,z::T,t::T,params::P)::T where {T<:Union{Float32,Float64}}`
unless `params` is `nothing`, in which case `func` should have signature
`func(x::T,z::T,t::T)::T where {T<:Union{Float32,Float64}}`.

# Fields
$(TYPEDFIELDS)
"""
struct PointwisePhysicalForcing{T,P} <: AbstractForcing{T}
    domain::Domain{T}
    func::Function
    params::P
end

"""
    $(TYPEDEF)

Pointwise physical forcing with no parameters.
`func` should have signature `func(x::T,z::T,t::T)::T where {T<:Union{Float32,Float64}}`.
"""
function PointwisePhysicalForcing(domain::Domain{T}, func::Function) where {T}
    return PointwisePhysicalForcing{T,Nothing}(domain, func, nothing)
end

"""
    $(TYPEDEF)

Apply global forcing to the Sawyer-Eliassen equation in physical space.
`func` is a function that updates the entire forcing variable in-place at time `t`.
`func` should have signature `func(F::XZVariable{T},t::T,params::P)::Nothing where {T<:Union{Float32,Float64}}`.
unless `params` is `nothing`, in which case `func` should have signature
`func(F::XZVariable{T},t::T)::T where {T<:Union{Float32,Float64}}`.

# Fields
$(TYPEDFIELDS)
"""
struct GlobalPhysicalForcing{T,P} <: AbstractForcing{T}
    domain::Domain{T}
    func::Function
    params::P
end

"""
    $(TYPEDEF)

Global physical forcing with no parameters.
`func` should have signature `func(F::XZVariable{T},t::T)::T where {T<:Union{Float32,Float64}}`.
"""
function GlobalPhysicalForcing(domain::Domain{T}, func::Function) where {T}
    return GlobalPhysicalForcing{T,Nothing}(domain, func, nothing)
end

"""Union of types representing forcing in physical (XZ) space."""
PhysicalForcing{T} =
    Union{PointwisePhysicalForcing{T,P},GlobalPhysicalForcing{T,Q}} where {T,P,Q}

####################
# Evaluate forcing #
####################

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

@inline function evaluate_physical_forcing!(
    F::PhysicalForcing{T}, out::XZVariable{T}, t::T
) where {T}
    @boundscheck consistent_domains(F, out) ||
        throw(ArgumentError("`F` and `out` must have the same domain."))
    return _evaluate_physical_forcing!(F, out, t)
end

@inline function _evaluate_physical_forcing!(
    F::PointwisePhysicalForcing{T,P}, out::XZVariable{T}, t::T
) where {T,P}
    x, z = gridpoints(F)
    @inbounds out .= F.func.(x, z, t, tuple(F.params))
    return nothing
end

@inline function _evaluate_physical_forcing!(
    F::PointwisePhysicalForcing{T,Nothing}, out::XZVariable{T}, t::T
) where {T}
    x, z = gridpoints(F)
    @inbounds out .= F.func.(x, z, t)
    return nothing
end

@inline function _evaluate_physical_forcing!(
    F::GlobalPhysicalForcing{T,P}, out::XZVariable{T}, t::T
) where {T,P}
    F.func(F, out, t, F.params)
    return nothing
end

@inline function _evaluate_physical_forcing!(
    F::GlobalPhysicalForcing{T,Nothing}, out::XZVariable{T}, t::T
) where {T}
    F.func(F, out, t)
    return nothing
end
