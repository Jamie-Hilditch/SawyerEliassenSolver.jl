"""$(TYPEDEF)

The domain, both physical and spectral, on which a problem is defined. This struct also stores
the FFTW plans for transforming between the physical and spectral representations of variables.
"""
struct Domain{T<:SingleOrDouble}
    grid::Grid{T}
    spectral::Spectral{T}
    transforms::Transforms{T}
end

"""$(TYPEDSIGNATURES)
Construct a domain from a grid. To dealias the products with the background flow a number of
the higher wavenumbers can be zeroed out. The number of retained wavenumbers are
`CNX = grid.NX ÷ 2 - dealias_x` and
`CNZ = grid.NZ - 1 - dealias_z`
Note that we already drop the Nyquist wavenumber in ``x`` as well as the zeroth (barotropic)
cosine mode / highest sine mode.
"""
function Domain(grid::Grid{T}; dealias_x=0, dealias_z=0) where {T}
    dealias_x >= 0 || throw(DomainError(dealias_x, "dealias_x must be non negative"))
    dealias_z >= 0 || throw(DomainError(dealias_z, "dealias_z must be non negative"))

    CNX = grid.NX ÷ 2 - dealias_x
    CNZ = grid.NZ - 1 - dealias_z
    CNX > 0 || throw(DomainError(dealias_x, "dealias_x must be smaller than NX ÷ 2"))
    CNZ > 0 || throw(DomainError(dealias_z, "dealias_z must be smaller than NZ - 1"))

    spectral = Spectral(grid, CNX, CNZ)
    transforms = Transforms(grid)
    return Domain{T}(grid, spectral, transforms)
end

Base.eltype(::Domain{T}) where {T} = T

# forward some functions onto the grid
Base.size(domain::Domain) = size(domain.grid)
"""$(TYPEDSIGNATURES)"""
xgridpoints(domain::Domain) = xgridpoints(domain.grid)
"""$(TYPEDSIGNATURES)"""
zgridpoints(domain::Domain) = zgridpoints(domain.grid)
"""$(TYPEDSIGNATURES)"""
gridpoints(domain::Domain) = gridpoints(domain.grid)

# and forward some on to the spectral domain
"""$(TYPEDSIGNATURES)
The ``x`` wavenumbers ``k_x``.
"""
xwavenumbers(domain::Domain) = xwavenumbers(domain.spectral)
"""$(TYPEDSIGNATURES)
The ``z`` wavenumbers ``k_z``.
"""
zwavenumbers(domain::Domain) = zwavenumbers(domain.spectral)
"""$(TYPEDSIGNATURES)
``k_x`` and ``k_z`` reshaped to size `(CNX,1)` and `(1,CNZ)` ready for broadcasting.
"""
wavenumbers(domain::Domain) = wavenumbers(domain.spectral)

"""Get the domain this object is defined over."""
function get_domain end

get_domain(domain::Domain) = domain

"""$(TYPEDSIGNATURES)
Check that these variables have the same domain.
"""
@inline consistent_domains(A...) = isempty(A) || mapreduce(get_domain, ==, A)

xgridpoints(x) = xgridpoints(get_domain(x))
zgridpoints(x) = zgridpoints(get_domain(x))
gridpoints(x) = gridpoints(get_domain(x))
xwavenumbers(x) = xwavenumbers(get_domain(x))
zwavenumbers(x) = zwavenumbers(get_domain(x))
wavenumbers(x) = wavenumbers(get_domain(x))
