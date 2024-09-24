module Output

using DocStringExtensions
using HDF5

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Problems

using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!

export OutputVariable,
    OutputWriter,
    add_output_variables!,
    OutputVariables

include("output_variables.jl")
include("output_writer.jl")
include("OutputVariables/OutputVariables.jl")

end  # module Output
