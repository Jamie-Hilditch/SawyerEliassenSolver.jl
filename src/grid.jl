"""
    $(TYPEDEF)

The physical grid. ``x`` is uniformly spaced and periodic with the left boundary at
``x[1]`` and the right boundary at ``x[NX + 1] = x[1]`` as is standard for DFTs. ``z`` is a
standard DCT/DST type II grid -- uniformly spaced interior points with the boundaries at
``z[1/2]`` and ``z[NZ + 1/2]``.

# Fields
$(TYPEDFIELDS)
"""
struct Grid
    """Number of grid points in ``x``"""
        NX::Int
    """Number of grid points in ``z``"""
        NZ::Int
    """Lower and upper ``x`` bounds: ``[x0,x1)``"""
        x_bounds::NTuple{2,Float64}
    """Lower and upper ``z`` bounds: ``[z0,z1]``"""
        z_bounds::NTuple{2,Float64}
    """Vector of ``x`` gridpoints"""
        x::Vector{Float64}
    """Vector of ``z`` gridpoints"""
        z::Vector{Float64}

    @doc """$(TYPEDSIGNATURES)
    Create a physical [`Grid`](@ref) with `(NX,NZ)` points. The lateral bounds, `[x0,x1)` are
    set with a 2-tuple `x_bounds` or a width `LX` in which case `x0 = 0, x1 = LX`. The
    vertical bounds `[z0,z1]` are set by the 2-tuple `z_bounds` or a height `LZ` in which
    case `z0 = -LZ, z1 = 0`.

    # Examples
    ```jldoctest
    julia> grid = Grid(256,128,(-1000,1000),(-1,0))
    Grid:
      ├── NX: 256
      ├── NZ: 128
      ├─── x: [-1000,1000)
      └─── z: [-1,0]
    ```
    """
    function Grid(NX, NZ, x_bounds, z_bounds)
        NX > 0 || throw(DomainError(NX, "NX must be positive"))
        NZ > 0 || throw(DomainError(NZ, "NZ must be positive"))
        x_bounds[1] != x_bounds[2] ||
            throw(ArgumentError("Length (x) of the domain cannot be 0"))
        z_bounds[1] != z_bounds[2] ||
            throw(ArgumentError("Length (z) of the domain cannot be 0"))

        x = LinRange(x_bounds[1], x_bounds[2], NX + 1)[1:(end - 1)]
        tmp = LinRange(z_bounds[1], z_bounds[2], NZ + 1)
        z = (tmp[1:(end - 1)] + tmp[2:end]) / 2

        return new(NX, NZ, x_bounds, z_bounds, x, z)
    end
end

"""$(TYPEDSIGNATURES)
See [`Grid(NX,NZ,x_bounds,z_bounds)`](@ref)

# Examples
```jldoctest
julia> grid = Grid(256,128,(-1000,1000),1)
Grid:
  ├── NX: 256
  ├── NZ: 128
  ├─── x: [-1000,1000)
  └─── z: [-1,0]
```
"""
Grid(NX, NZ, LX::Real, z_bounds) = Grid(NX, NZ, (0, LX), z_bounds)

"""$(TYPEDSIGNATURES)
See [`Grid(NX,NZ,x_bounds,z_bounds)`](@ref)

# Examples
```jldoctest
julia> grid = Grid(256,128,2000,(-1,0))
Grid:
  ├── NX: 256
  ├── NZ: 128
  ├─── x: [0,2000)
  └─── z: [-1,0]
```
"""
Grid(NX, NZ, x_bounds, LZ::Real) = Grid(NX, NZ, x_bounds, (-LZ, 0))

"""$(TYPEDSIGNATURES)
See [`Grid(NX,NZ,x_bounds,z_bounds)`](@ref)

# Examples
```jldoctest
julia> grid = Grid(256,128,2000,1)
Grid:
  ├── NX: 256
  ├── NZ: 128
  ├─── x: [0,2000)
  └─── z: [-1,0]
```
"""
Grid(NX, NZ, LX::Real, LZ::Real) = Grid(NX, NZ, (0, LX), (-LZ, 0))

"""$(TYPEDSIGNATURES)"""
Base.size(grid::Grid) = (grid.NX, grid.NZ)

"""$(TYPEDSIGNATURES)"""
function Base.show(io::IO, ::MIME"text/plain", grid::Grid)
    return print(
        io,
        "Grid:\n",
        "  ├── NX: $(grid.NX)\n",
        "  ├── NZ: $(grid.NZ)\n",
        "  ├─── x: [$(sfmt(grid.x_bounds[1])),$(sfmt(grid.x_bounds[2])))\n",
        "  └─── z: [$(sfmt(grid.z_bounds[1])),$(sfmt(grid.z_bounds[2]))]\n",
    )
end

"""$(TYPEDSIGNATURES)"""
Base.summary(io::IO, grid::Grid) = print(io, "Grid of size ($(grid.NX),$(grid.NZ))")
