module Problems

using DocStringExtensions

using SawyerEliassenSolver.Utils
using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Forcing

export BackgroundFlow, State, Clock, Problem

include("background.jl")
include("state.jl")
include("problem.jl")

end
