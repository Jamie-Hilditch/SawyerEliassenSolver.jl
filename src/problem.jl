"""Main problem struct"""

include("problem/spectral_domain.jl")
include("problem/transforms.jl")
include("problem/scratch.jl")

struct Problem
    grid::Grid
    background_flow::BackgroundFlow
    spectral_domain::SpectralDomain
    transforms::Transforms
    scratch::Scratch

    function Problem(grid, background_flow; dealias_x=0, dealias_z=0)
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

Base.size(problem::Problem) = size(problem.grid)

function setup_simulation(problem::Problem, h, cg_tol::Union{Float64,Nothing}=nothing)
    timestepper = Timestepper!(problem, h, cg_tol)
    return timestepper, timestepper.state
end
