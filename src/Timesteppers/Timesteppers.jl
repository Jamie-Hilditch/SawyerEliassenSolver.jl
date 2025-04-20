module Timesteppers

using DocStringExtensions
using LinearAlgebra

using SawyerEliassenSolver.Utils
using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Forcing
using SawyerEliassenSolver.Problems

using Base: @propagate_inbounds
using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!
using SawyerEliassenSolver.Problems: update_clock!
using SawyerEliassenSolver.Forcing: evaluate_ζ_forcing!, evaluate_physical_forcing!

export AbstractPreconditioner,
    IdentityPreconditioner, DiagonalQuadraticPreconditioner, Timestepper, advance!

include("linear_operators.jl")
include("preconditioners.jl")
include("conjugate_gradient.jl")
include("timestepper.jl")
include("advance.jl")
end
