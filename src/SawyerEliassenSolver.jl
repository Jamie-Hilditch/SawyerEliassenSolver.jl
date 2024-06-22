"""
`SawyerEliassenSolver.jl` -- Solve the Sawyer-Eliassen equation using a pseudo-spectral
discretisation and 3rd order accurate implicit timestepping for arbitrary balanced
background flows.

# Exports

"""
module SawyerEliassenSolver

using DocStringExtensions
using LinearAlgebra
using Printf

using FFTW: FFTW
FFTW.set_num_threads(Threads.nthreads())

export
    # structs
    Grid,
    Domain,
    XZVariable,
    XSVariable,
    XCVariable,
    FZVariable,
    FSVariable,
    FCVariable
    # BackgroundFlow,
    # Problem,
    # Timestepper!,
    # State,

    # functions
    # setup_simulation,

    # submodules
    # Tools

include("display.jl")

include("domain/grid.jl")
include("domain/spectral.jl")
include("domain/transforms.jl")
include("domain/domain.jl")

include("variables/variables.jl")
include("variables/transforms.jl")

# include("background.jl")
# include("problem.jl")
# include("timestepper.jl")

# submodules
# include("Tools/Tools.jl")

end
