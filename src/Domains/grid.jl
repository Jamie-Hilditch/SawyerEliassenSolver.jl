"""$(TYPEDEF)

The physical grid.
``x`` is uniformly spaced and periodic with the left boundary at ``x[1]`` and the right boundary at ``x[NX + 1] = x[1]`` as is standard for DFTs.
``z`` is a standard DCT/DST type II grid -- uniformly spaced interior points with the boundaries at ``z[1/2]`` and ``z[NZ + 1/2]``.
The float type of the grid determines the float point precision of the simulation and since we use FFTW is restriced to `Float32` or `Float64`.

# Fields
$(TYPEDFIELDS)
"""
struct Grid{T<:SingleOrDouble}
    """Number of grid points in ``x``"""
    NX::Int
    """Number of grid points in ``z``"""
    NZ::Int
    """Lower and upper ``x`` bounds: ``[x₀,x₁)``"""
    x_bounds::NTuple{2,T}
    """Lower and upper ``z`` bounds: ``[z₀,z₁]``"""
    z_bounds::NTuple{2,T}

    @doc """$(TYPEDSIGNATURES)
    Create a physical [`Grid`](@ref) with float type `FT` and `(NX,NZ)` points.
    The lateral bounds, `[x₀,x₁)` are set with a 2-tuple `x_bounds` or a single number `LX` which is equivalent to passing `(0, LX)`.
    Similarly, the vertical bounds `[z₀,z₁]` are set by the 2-tuple `z_bounds` or a height `LZ` equivalent to `(-LZ, 0)`.

    # Examples
    ```jldoctest
    julia> grid = Grid(256,128,(-1000,1000),(-1,1))
    Grid{Float64}:
      ├── NX: 256
      ├── NZ: 128
      ├─── x: [-1000,1000)
      └─── z: [-1,1]
    ```
    ```jldoctest
    julia> grid = Grid(Float32,256,128,(-1000,1000),10)
    Grid{Float32}:
      ├── NX: 256
      ├── NZ: 128
      ├─── x: [-1000,1000)
      └─── z: [-10,0]
    ```
    """
    function Grid(FT, NX, NZ, x_bounds, z_bounds)
        NX > 0 || throw(DomainError(NX, "NX must be positive"))
        NZ > 0 || throw(DomainError(NZ, "NZ must be positive"))
        x_bounds[1] != x_bounds[2] ||
            throw(ArgumentError("Length (x) of the domain cannot be 0"))
        z_bounds[1] != z_bounds[2] ||
            throw(ArgumentError("Length (z) of the domain cannot be 0"))

        return new{FT}(NX, NZ, x_bounds, z_bounds)
    end
end

function Grid(FT, NX, NZ, LX::Real, z_bounds::Tuple{Real,Real})
    return Grid(FT, NX, NZ, (zero(FT), LX), z_bounds)
end
function Grid(FT, NX, NZ, x_bounds::Tuple{Real,Real}, LZ::Real)
    return Grid(FT, NX, NZ, x_bounds, (-LZ, zero(FT)))
end
Grid(FT, NX, NZ, LX::Real, LZ::Real) = Grid(FT, NX, NZ, (zero(FT), LX), (-LZ, zero(FT)))
"""$(TYPEDSIGNATURES)
4 argument constructor for  `Grid` uses the default float type of `Float64`.
"""
Grid(NX, NZ, x_bounds, z_bounds) = Grid(Float64, NX, NZ, x_bounds, z_bounds)

# define some standard functions from base for the Grid struct
Base.size(grid::Grid) = (grid.NX, grid.NZ)
Base.eltype(::Grid{T}) where {T} = T
function Base.summary(io::IO, grid::Grid)
    return print(io, "Grid with eltype $(eltype(grid)) and size $(size(grid))")
end

function Base.show(io::IO, ::MIME"text/plain", grid::Grid{T}) where {T}
    return print(
        io,
        "Grid{$(T)}:\n",
        "  ├── NX: $(grid.NX)\n",
        "  ├── NZ: $(grid.NZ)\n",
        "  ├─── x: [$(sfmt(grid.x_bounds[1])),$(sfmt(grid.x_bounds[2])))\n",
        "  └─── z: [$(sfmt(grid.z_bounds[1])),$(sfmt(grid.z_bounds[2]))]\n",
    )
end

# Gridpoints

"""$(TYPEDSIGNATURES)
`LinRange` of the ``x`` gridpoints.
"""
@inline xgridpoints(grid::Grid) = LinRange(grid.x_bounds[1], grid.x_bounds[2], grid.NX + 1)[1:(end - 1)]

"""$(TYPEDSIGNATURES)
`LinRange` of the ``z`` gridpoints.
"""
@inline function zgridpoints(grid::Grid)
    z_faces = LinRange(grid.z_bounds[1], grid.z_bounds[2], grid.NZ + 1)
    return z_faces[1:(end - 1)] .+ step(z_faces) / 2
end

"""$(TYPEDSIGNATURES)
Get both the ``x`` and ``z`` gridpoints, reshaped to size `(NX,1)` and `(1,NZ)` ready for broadcasting.
"""
@inline function gridpoints(grid::Grid)
    x = xgridpoints(grid)
    z = zgridpoints(grid)
    return reshape(x, :, 1), reshape(z, 1, :)
end

"""$(TYPEDSIGNATURES)"""
@inline xstepsize(grid::Grid) = (grid.x_bounds[2] - grid.x_bounds[1]) / grid.NX
"""$(TYPEDSIGNATURES)"""
@inline zstepsize(grid::Grid) = (grid.z_bounds[2] - grid.z_bounds[1]) / (grid.NZ + 1)
"""$(TYPEDSIGNATURES)"""
@inline stepsize(grid::Grid) = xstepsize(grid), zstepsize(grid)
