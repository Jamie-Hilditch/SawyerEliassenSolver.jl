module SawyerEliassenSolver

using LinearAlgebra 
using Printf 

import FFTW
FFTW.set_num_threads(Threads.nthreads())


include("display.jl")
include("grid.jl")
include("background.jl")
include("problem.jl")
include("timestepper.jl")

# submodules
include("Tools/Tools.jl")

export 
    # structs
    Grid,
    BackgroundFlow,
    Problem,
    Timestepper!,
    State,

    # functions
    setup_simulation,

    # submodules
    Tools

end
