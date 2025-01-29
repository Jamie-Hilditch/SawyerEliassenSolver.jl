module Problems

using DocStringExtensions
using Base: @propagate_inbounds

using SawyerEliassenSolver.Utils
using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Forcing

using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!

export BackgroundFlow,
    State,
    Clock,
    Scratch,
    Problem,
    set_ζ!,
    set_v!,
    set_b!,
    compute_ζₜ!,
    set_vb!,
    get_time,
    get_iteration,
    get_f,
    get_Vx,
    get_Bx,
    get_Bz,
    get_problem

include("background.jl")
include("state.jl")
include("problem.jl")
include("initial_conditions.jl")

end
