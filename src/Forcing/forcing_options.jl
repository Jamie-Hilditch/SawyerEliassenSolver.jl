
abstract type ForcingFunction end

"""
    $(TYPEDEF)

Apply pointwise forcing to the Sawyer-Eliassen equation in physical space.
`func` is a function that computes the forcing at a given position and time.
`func` should have signature `func(x::T,z::T,t::T,params::P)::T where {T<:Union{Float32,Float64}}`.

# Fields
$(TYPEDFIELDS)
"""
struct PointwisePhysicalForcing{T,P} <: ForcingFunction
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
struct GlobalPhysicalForcing{T,P} <: ForcingFunction
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
struct PointwiseSpectralForcing{T,P} <: ForcingFunction
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
struct GlobalSpectralForcing{T,P} <: ForcingFunction
    domain::Domain{T}
    func::Function
    params::P
end

PhysicalForcing{T} = Union{PointwisePhysicalForcing{T,P},GlobalPhysicalForcing{T,Q}} where {T,P,Q}
SpectralForcing{T} = Union{PointwiseSpectralForcing{T,P},GlobalSpectralForcing{T,Q}} where {T,P,Q}
