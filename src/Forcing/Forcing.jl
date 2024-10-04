module Forcing

using DocStringExtensions

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables

using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!

export AbstractForcingFunction,
    NoForcing,
    PhysicalForcing,
    SpectralForcing,
    PointwisePhysicalForcing,
    GlobalPhysicalForcing,
    PointwiseSpectralForcing,
    GlobalSpectralForcing

include("forcing_options.jl")
end
