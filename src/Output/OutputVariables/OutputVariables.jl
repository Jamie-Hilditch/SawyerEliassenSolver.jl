module OutputVariables

using DocStringExtensions

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Problems
using SawyerEliassenSolver.Output

using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!

export u,
    ∂u∂x,
    ∂²u∂x²,
    ∂u∂z,
    ∂²u∂z²,
    ∂²u∂x∂z,
    w,
    ∂w∂x,
    ∂²w∂x²,
    ∂w∂z,
    ∂²w∂z²,
    ∂²w∂x∂z


include("u_and_derivatives.jl")
include("w_and_derivatives.jl")

"""
Implementation of common output variables for the Sawyer-Eliassen solver.

# Exports
$(EXPORTS)
"""
OutputVariables

end  # module OutputVariables
