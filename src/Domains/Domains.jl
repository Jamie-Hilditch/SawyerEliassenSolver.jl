module Domains

using DocStringExtensions
using FFTW

using SawyerEliassenSolver.Utils

export
    Grid,
    Domain

# Since we use FFTW we are limited to Float32 or Float64
SingleOrDouble = Union{Float32, Float64}

include("grid.jl")
include("spectral.jl")
include("transforms.jl")
include("domain.jl")

end
