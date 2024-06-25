module Forcing

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables

export PointwisePhysicalForcing,
    GlobalPhysicalForcing, PointwiseSpectralForcing, GlobalSpectralForcing

    include("forcing_options.jl")
end
