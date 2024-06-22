module Domains

using DocStringExtensions
using FFTW

using SawyerEliassenSolver.Utils

export
    Grid,
    Domain

include("grid.jl")
include("spectral.jl")
include("transforms.jl")
include("domain.jl")

end
