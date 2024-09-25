module Output

using DocStringExtensions
using HDF5

using SawyerEliassenSolver.Domains
using SawyerEliassenSolver.Variables
using SawyerEliassenSolver.Problems

using SawyerEliassenSolver.Variables: Tᴴ!, Tˢ!, Tᶜ!

export OutputVariable,
    compute!,
    OutputWriter,
    add_output_variables!,
    write!,
    write_attributes!,
    OutputVariables

include("output_variable.jl")
include("output_writer.jl")
include("OutputVariables/OutputVariables.jl")

end  # module Output
