"""
    $(TYPEDEF)

Apply pointwise forcing to the Sawyer-Eliassen equation in spectral space.
`func` is a function that computes the forcing at given wavenumbers and time.
`func` should have signature `func(kx::T,kz::T,t::T,params::P)::complex{T} where {T<:Union{Float32,Float64}}`.
unless `params` is `nothing`, in which case `func` should have signature
`func(kx::T,kz::T,t::T)::T where {T<:Union{Float32,Float64}}`

# Fields
$(TYPEDFIELDS)
"""
struct PointwiseSpectralForcing{T,P} <: AbstractForcing{T}
    domain::Domain{T}
    func::Function
    params::P
end

"""
    $(TYPEDEF)

Pointwise spectral forcing with no parameters.
`func` should have signature `func(kx::T,kz::T,t::T)::T where {T<:Union{Float32,Float64}}`.
"""
function PointwiseSpectralForcing(domain::Domain{T}, func::Function) where {T}
    return PointwiseSpectralForcing{T,Nothing}(domain, func, nothing)
end

"""
    $(TYPEDEF)

Apply global forcing to the Sawyer-Eliassen equation in spectral space.
`func` is a function that updates the entire forcing variable in-place at time `t`.
`func` should have signature `func(F::FSVariable{T},t::T,params::P)::Nothing where {T<:Union{Float32,Float64}}`
unless `params` is `nothing`, in which case `func` should have signature
`func(F::FSVariable{T},t::T)::T where {T<:Union{Float32,Float64}}`
Note that the `FSVariable{T}` has entries of type `Complex{T}`.

# Fields
$(TYPEDFIELDS)
"""
struct GlobalSpectralForcing{T,P} <: AbstractForcing{T}
    domain::Domain{T}
    func::Function
    params::P
end

"""
    $(TYPEDEF)

Global spectral forcing with no parameters.
`func` should have signature `func(F::FSVariable{T},t::T)::T where {T<:Union{Float32,Float64}}`.
"""
function GlobalSpectralForcing(domain::Domain{T}, func::Function) where {T}
    return GlobalSpectralForcing{T,Nothing}(domain, func, nothing)
end

"""Union of types representing forcing in spectral (FS) space."""
SpectralForcing{T} = Union{PointwiseSpectralForcing{T,P},GlobalSpectralForcing{T,Q}} where {T,P,Q}

####################
# Evaluate forcing #
####################

@inline function evaluate_ζ_forcing!(
    F::SpectralForcing{T},
    out::FSVariable{T},
    t::T,
    XS::XSVariable{T},
    XZ::XZVariable{T}
) where {T}
    @boundscheck consistent_domains(F, out) ||
    throw(ArgumentError("`F` and `out` must have the same domain."))
    _evaluate_ζ_forcing!(F, out, t, XS, XZ)
end

@inline function _evaluate_ζ_forcing!(
    F::PointwiseSpectralForcing{T,P},
    out::FSVariable{T},
    t::T,
    ::XSVariable{T},
    ::XZVariable{T},
) where {T,P}
    kx, kz = wavenumbers_full(F)
    @inbounds out .= F.func.(kx, kz, t, tuple(F.params))
    return nothing
end

@inline function _evaluate_ζ_forcing!(
    F::PointwiseSpectralForcing{T,Nothing},
    out::FSVariable{T},
    t::T,
    ::XSVariable{T},
    ::XZVariable{T},
) where {T}
    kx, kz = wavenumbers_full(F)
    @inbounds out .= F.func.(kx, kz, t)
    return nothing
end

@inline function _evaluate_ζ_forcing!(
    F::GlobalSpectralForcing{T,P}, out::FSVariable{T}, t::T, ::XSVariable, ::XZVariable
) where {T,P}
    F.func(out, t, F.params)
    return nothing
end

@inline function _evaluate_ζ_forcing!(
    F::GlobalSpectralForcing{T,Nothing}, out::FSVariable{T}, t::T, ::XSVariable, ::XZVariable
) where {T}
    F.func(out, t)
    return nothing
end
