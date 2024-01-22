"""Physical grid"""
struct Grid 
    NX :: Int 
    NZ :: Int 
    x_bounds :: NTuple{2,Float64}
    z_bounds :: NTuple{2,Float64}
    x :: Array{Float64,1}
    z :: Array{Float64,1}

    function Grid(NX,NZ,x_bounds,z_bounds)
        NX > 0 || throw(DomainError(NX,"NX must be positive"))
        NZ > 0 || throw(DomainError(NZ,"NZ must be positive"))
        x_bounds[1] != x_bounds[2] || throw(ArgumentError("Length (x) of the domain cannot be 0"))
        z_bounds[1] != z_bounds[2] || throw(ArgumentError("Length (z) of the domain cannot be 0"))
        
        x = LinRange(x_bounds[1],x_bounds[2],NX+1)[1:end-1]
        tmp = LinRange(z_bounds[1],z_bounds[2],NZ+1)
        z = (tmp[1:end-1] + tmp[2:end])/2 

        new(NX, NZ, x_bounds, z_bounds, x, z)
    end
end

Grid(NX,NZ,LX :: Real,z_bounds) = Grid(NX,NZ,(0,LX),z_bounds)
Grid(NX,NZ,x_bounds, LZ :: Real) = Grid(NX,NZ,x_bounds,(-LZ,0))
Grid(NX,NZ,LX :: Real, LZ :: Real) = Grid(NX,NZ,(0,LX), (-LZ,0))

Base.size(grid :: Grid) = (grid.NX,grid.NZ)

Base.show(io::IO, ::MIME"text/plain", grid::Grid) = print(io,
    "Grid:\n",
    "  ├── NX: $(grid.NX)\n",
    "  ├── NZ: $(grid.NZ)\n",
    "  ├─── x: [$(sfmt(grid.x_bounds[1])),$(sfmt(grid.x_bounds[2])))\n",
    "  └─── z: [$(sfmt(grid.z_bounds[1])),$(sfmt(grid.z_bounds[2]))]\n",
)

Base.summary(io::IO, grid::Grid) = print(io,"Grid of size ($(grid.NX),$(grid.NZ))")