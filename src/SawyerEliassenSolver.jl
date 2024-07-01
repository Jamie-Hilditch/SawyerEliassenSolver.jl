module SawyerEliassenSolver

using DocStringExtensions
using Printf
using Reexport

using FFTW: FFTW
FFTW.set_num_threads(Threads.nthreads())

include("Utils/Utils.jl")
include("Domains/Domains.jl")
include("Variables/Variables.jl")
include("Forcing/Forcing.jl")
include("Problems/Problems.jl")
include("Timesteppers/Timesteppers.jl")

@reexport using .Domains
@reexport using .Variables
@reexport using .Forcing
@reexport using .Problems
@reexport using .Timesteppers

# add a module docstring.
# N.B. to use DocStringExtensions .e. $(EXPORTS) we need to define the docstring inside the module.

"""
`SawyerEliassenSolver.jl` -- Solve the Sawyer-Eliassen equation using a pseudo-spectral
discretisation and 3rd order accurate implicit timestepping for arbitrary balanced
background flows.

# Exports
$(EXPORTS)
"""
SawyerEliassenSolver
end
