"""Implementations of common output variables.

$(Exports)

"""
module OutputVariables

using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Problems
using SawyerEliassenSolver.Output

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

end  # module OutputVariables
