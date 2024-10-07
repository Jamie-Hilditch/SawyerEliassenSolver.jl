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

# Examples
Build an domain with no dealiasing.
```jldoctest
julia> grid = Grid(16,8,1,1);

julia> domain = Domain(grid)
Domain:
  ├─────────── grid: Grid with eltype Float64 and size (16, 8)
  ├─────── spectral: Spectral domain of size (9, 8) and spectral resolution (8, 7)
  └───── transforms: FFTW transforms: rfft, type II DST and type II DCT.
```
A domain with dealising in the vertical.
```jldoctest
julia> grid = Grid(16,8,1,1);

julia> domain = Domain(grid, dealias_z=2)
Domain:
  ├─────────── grid: Grid with eltype Float64 and size (16, 8)
  ├─────── spectral: Spectral domain of size (9, 8) and spectral resolution (8, 5)
  └───── transforms: FFTW transforms: rfft, type II DST and type II DCT.
```
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
"""$(TYPEDSIGNATURES)
The ``x`` wavenumbers ``k_x``.
"""
xwavenumbers_full(domain::Domain) = xwavenumbers_full(domain.spectral)
"""$(TYPEDSIGNATURES)
The ``z`` wavenumbers ``k_z``.
"""
zwavenumbers_full(domain::Domain) = zwavenumbers_full(domain.spectral)
"""$(TYPEDSIGNATURES)
``k_x`` and ``k_z`` reshaped to size `(CNX,1)` and `(1,CNZ)` ready for broadcasting.
"""
wavenumbers_full(domain::Domain) = wavenumbers_full(domain.spectral)
"""$(TYPEDSIGNATURES)"""
resolved_fourier_indices(domain::Domain) = resolved_fourier_indices(domain.spectral)
"""$(TYPEDSIGNATURES)"""
resolved_sine_indices(domain::Domain) = resolved_sine_indices(domain.spectral)
"""$(TYPEDSIGNATURES)"""
resolved_cosine_indices(domain::Domain) = resolved_cosine_indices(domain.spectral)

function Base.show(io::IO, ::MIME"text/plain", domain::Domain)
    return print(
        io,
        "Domain:\n",
        "  ├─────────── grid: $(summary(domain.grid))\n",
        "  ├─────── spectral: $(summary(domain.spectral))\n",
        "  └───── transforms: $(summary(domain.transforms))\n",
    )
end

function Base.summary(io::IO, domain::Domain)
    return print(
        io, "Domain with eltype $(eltype(domain)) and physical size $(size(domain))"
    )
end

"""Get the domain this object is defined over."""
function get_domain end

get_domain(domain::Domain) = domain

"""$(TYPEDSIGNATURES)
Check that these variables have the same domain.
"""
@inline function consistent_domains(a, B...)
    domain_of_a = get_domain(a)
    all(b -> get_domain(b) === domain_of_a, B)
end
# could check that `a` actually has a domain but for now let the upstream function trigger a
# method error on `get_domain`
@inline consistent_domains(a) = true
@inline consistent_domains() = true

# forward some functions onto domain
xgridpoints(x) = xgridpoints(get_domain(x))
zgridpoints(x) = zgridpoints(get_domain(x))
gridpoints(x) = gridpoints(get_domain(x))
xwavenumbers(x) = xwavenumbers(get_domain(x))
zwavenumbers(x) = zwavenumbers(get_domain(x))
wavenumbers(x) = wavenumbers(get_domain(x))
xwavenumbers_full(x) = xwavenumbers_full(get_domain(x))
zwavenumbers_full(x) = zwavenumbers_full(get_domain(x))
wavenumbers_full(x) = wavenumbers_full(get_domain(x))
resolved_fourier_indices(x) = resolved_fourier_indices(get_domain(x))
resolved_sine_indices(x) = resolved_sine_indices(get_domain(x))
resolved_cosine_indices(x) = resolved_cosine_indices(get_domain(x))
