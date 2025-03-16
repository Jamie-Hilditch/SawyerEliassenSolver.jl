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

include("forcing_options.jl")

"""
Module for implementing different RHS forcings.

# Exports
$(EXPORTS)
"""
Forcing
end
