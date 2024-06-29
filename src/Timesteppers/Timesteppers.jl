module Timesteppers

using DocStringExtensions
using LinearAlgebra

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Forcing
using SawyerEliassenSolver.Problems

using Base: @propagate_inbounds
using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!

export AbstractPreconditioner

include("linear_operators.jl")
include("preconditioners.jl")
include("conjugate_gradient.jl")
include("timestepper.jl")
end
