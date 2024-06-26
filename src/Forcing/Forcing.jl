module Forcing

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables

export AbstractForcingFunction,
    PhysicalForcing,
    SpectralForcing,
    PointwisePhysicalForcing,
    GlobalPhysicalForcing,
    PointwiseSpectralForcing,
    GlobalSpectralForcing

include("forcing_options.jl")
end
