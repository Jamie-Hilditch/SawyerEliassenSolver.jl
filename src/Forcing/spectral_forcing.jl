"""
    $(TYPEDEF)

Apply pointwise forcing to the Sawyer-Eliassen equation in spectral space.
`func` is a function that computes the forcing at given wavenumbers and time.
`func` should have signature `func(kx::T,kz::T,t::T,params::P)::complex{T} where {T<:Union{Float32,Float64}}`.

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

Apply global forcing to the Sawyer-Eliassen equation in spectral space.
`func` is a function that updates the entire forcing variable in-place at time `t`.
`func` should have signature `func(F::FSVariable{T},t::T,params::P)::Nothing where {T<:Union{Float32,Float64}}`.
Note that the `FSVariable{T}` has entries of type `Complex{T}`.

# Fields
$(TYPEDFIELDS)
"""
struct GlobalSpectralForcing{T,P} <: AbstractForcing{T}
    domain::Domain{T}
    func::Function
    params::P
end

"""Union of types representing forcing in spectral (FS) space."""
SpectralForcing{T} =
    Union{PointwiseSpectralForcing{T,P},GlobalSpectralForcing{T,Q}} where {T,P,Q}

####################
# Evaluate forcing #
####################

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
    @inbounds out .= F.func.(kx, kz, tuple(t), tuple(F.params))
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
