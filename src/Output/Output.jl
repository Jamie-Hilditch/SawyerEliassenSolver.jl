module Output

using DocStringExtensions
using HDF5

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Problems

using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!

export OutputVariable,
    OutputWriter,
    compute_u_from_Ψ!,
    compute_w_from_Ψ!,
    compute_u_and_w!

include("derived_variables.jl")
include("output_variables.jl")
include("output_writer.jl")
end
