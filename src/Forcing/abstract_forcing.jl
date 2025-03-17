
"""Abstract super type for all forcings."""
abstract type AbstractForcing{T} end

Domains.get_domain(F::AbstractForcing) = F.domain

# show methods

function Base.show(io::IO, ::MIME"text/plain", forcing::AbstractForcing)
    forcing_type = typeof(forcing)
    forcing_type_name = nameof(forcing_type)
    forcing_type_string = string(forcing_type_name)
    return print(
        io,
        "$(forcing_type_string):\n",
        "  ├──── domain: $(summary(get_domain(forcing)))\n",
        "  └──── params: $(forcing.params)\n",
    )
end

"""Type representing no forcing."""
struct NoForcing{T} <: AbstractForcing{T}
    domain::Domain{T}
end

function Base.show(io::IO, ::MIME"text/plain", forcing::NoForcing)
    return print(io, "NoForcing:\n", "  └──── domain: $(summary(get_domain(forcing)))\n")
end

@inline function evaluate_physical_forcing!(
    ::NoForcing{T}, XZ::XZVariable{T}, t::T
) where {T}
    @. XZ = 0
    return nothing
end

@inline function evaluate_ζ_forcing!(
    ::NoForcing{T}, out::FSVariable{T}, t::T, ::XSVariable{T}, ::XZVariable{T}
) where {T}
    @. out = 0
    return nothing
end
