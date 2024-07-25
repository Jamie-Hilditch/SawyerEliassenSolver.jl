module Output

using DocStringExtensions

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Problems

using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!

export compute_u_from_Ψ!,
    compute_w_from_Ψ!,
    compute_u_and_w!

include("derived_variables.jl")
end
