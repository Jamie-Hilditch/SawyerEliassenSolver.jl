module Timestepper

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Forcing
using SawyerEliassenSolver.Problems
using Base: @propagate_inbounds

export AbstractPreconditioner

include("linear_operators.jl")
include("preconditioners.jl")
include("conjugate_gradient.jl")
end
