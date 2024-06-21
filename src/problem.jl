"""Main problem struct"""

include("problem/spectral_domain.jl")
include("problem/transforms.jl")
include("problem/scratch.jl")

"""
    $(TYPEDEF)

Components of a Sawyer-Eliassen problem.

# Fields
$(TYPEDFIELDS)
"""
struct Problem
    grid::Grid
    background_flow::BackgroundFlow
    """Wavenumbers and metadata for spectral domain"""
    spectral_domain::SpectralDomain
    """Transforms between physical and spectral space"""
    transforms::Transforms
    """Temporary working arrays"""
    scratch::Scratch

    @doc """
        Problem(
            grid::Grid,
            background_flow::BackgroundFlow;
            dealias_x::Int=0,
            dealias_z::Int=0
        ) -> Problem

    Set-up a Sawyer-Eliassen problem. The Sawyer-Eliassen equation involves a product
    between the background gradients and the solution. For smooth background fields the
    aliasing is confined to just the highest wavenumbers. Dealiasing should be applied
    equal to the bandwidth of the background flow.

    Keyword arguments
    ===================
      - `dealias_x`: the number (default 0) of high-wavenumbers that are zero-ed out
      - `dealias_z`: the number (default 0) of high-modes that are zero-ed out

    # Examples
    ```jldoctest
    grid = Grid(256,128,2000,1)
    bg = BackgroundFlow(grid)
    problem = Problem(grid,bg,dealias_x=10,dealias_z=10)

    # output
    Problem:
      ├───────────── grid: Grid of size (256,128)
      ├── background_flow: BackgroundFlow: f = 1, Vx,Bx,Bz = 256×128 Matrix{Float64}
      ├── spectral_domain: SawyerEliassenSolver.SpectralDomain
      ├─────── transforms: SawyerEliassenSolver.Transforms
      └────────── scratch: SawyerEliassenSolver.Scratch
    ```
    """
    function Problem(
        grid::Grid, background_flow::BackgroundFlow; dealias_x::Int=0, dealias_z::Int=0
    )
        check_consistent_grid(grid, background_flow)
        scratch = Scratch(grid.NX, grid.NZ)
        return new(
            grid,
            background_flow,
            SpectralDomain(grid, dealias_x, dealias_z),
            Transforms(scratch.rtmp1),
            scratch,
        )
    end
end

"""$(TYPEDSIGNATURES)"""
function Base.show(io::IO, ::MIME"text/plain", problem::Problem)
    return print(
        io,
        "Problem:\n",
        "  ├───────────── grid: $(summary(problem.grid))\n",
        "  ├── background_flow: $(summary(problem.background_flow))\n",
        "  ├── spectral_domain: $(summary(problem.spectral_domain))\n",
        "  ├─────── transforms: $(summary(problem.transforms))\n",
        "  └────────── scratch: $(summary(problem.scratch))\n",
    )
end

"""$(TYPEDSIGNATURES)"""
Base.size(problem::Problem) = size(problem.grid)

"""$(TYPEDSIGNATURES)

Construct the timestepper given a timestep `h`. Optionally modify the tolerance for the
conjugate gradient solver.

# Examples
```
grid = Grid(256,128,2000,1)
bg = BackgroundFlow(grid)
problem = Problem(grid,bg,dealias_x=10,dealias_z=10)
timestepper, state = setup_simulation(problem,0.1)

# output
[...]
```
"""
function setup_simulation(problem::Problem, h, cg_tol::Union{Float64,Nothing}=nothing)
    timestepper = Timestepper!(problem, h, cg_tol)
    return timestepper, timestepper.state
end
