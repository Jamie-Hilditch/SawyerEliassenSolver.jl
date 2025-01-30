"""
"""
module Output

using Base: @propagate_inbounds

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
    OutputVariables,
    ConstantOutputVariable,
    write_constant_array!,
    write_f!,
    write_Vx!,
    write_Bx!,
    write_Bz!,
    write_background_flow!,
    compute_background_buoyancy,
    write_background_buoyancy!

include("output_variable.jl")
include("output_writer.jl")
include("OutputVariables/OutputVariables.jl")
include("constant_output_variable.jl")

end  # module Output
