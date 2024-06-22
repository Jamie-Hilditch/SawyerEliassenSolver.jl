"""
`SawyerEliassenSolver.jl` -- Solve the Sawyer-Eliassen equation using a pseudo-spectral
discretisation and 3rd order accurate implicit timestepping for arbitrary balanced
background flows.

# Exports

"""
module SawyerEliassenSolver

using Printf
using Reexport

using FFTW: FFTW
FFTW.set_num_threads(Threads.nthreads())

include("Utils/Utils.jl")
include("Domains/Domains.jl")
include("Variables/Variables.jl")

@reexport using .Domains
@reexport using .Variables

# include("background.jl")
# include("problem.jl")
# include("timestepper.jl")

# submodules
# include("Tools/Tools.jl")

end
