module Domains

using DocStringExtensions
using FFTW

using SawyerEliassenSolver.Utils

export Grid,
    xgridpoints,
    zgridpoints,
    gridpoints,
    xstepsize,
    zstepsize,
    stepsize,
    xwavenumbers,
    zwavenumbers,
    wavenumbers,
    xwavenumbers_full,
    zwavenumbers_full,
    wavenumbers_full,
    resolved_fourier_indices,
    resolved_sine_indices,
    resolved_cosine_indices,
    Domain,
    get_domain,
    consistent_domains

# Since we use FFTW we are limited to Float32 or Float64
SingleOrDouble = Union{Float32,Float64}

include("grid.jl")
include("spectral.jl")
include("transforms.jl")
include("domain.jl")

end
