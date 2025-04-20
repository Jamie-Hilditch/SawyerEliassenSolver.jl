module Forcing

using DocStringExtensions

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables

using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!

export AbstractForcing,
    NoForcing,
    PhysicalForcing,
    SpectralForcing,
    PointwisePhysicalForcing,
    GlobalPhysicalForcing,
    PointwiseSpectralForcing,
    GlobalSpectralForcing

include("abstract_forcing.jl")
include("physical_forcing.jl")
include("spectral_forcing.jl")

"""
Module for implementing different RHS forcings.

# Exports
$(EXPORTS)
"""
Forcing
end
