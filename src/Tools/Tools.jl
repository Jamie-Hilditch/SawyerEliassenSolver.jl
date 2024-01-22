"""convenience functions"""

module Tools

using HDF5
using LinearAlgebra

using ..SawyerEliassenSolver

include("initial_conditions.jl")
include("output_writer.jl")
include("wall_time_stop.jl")

export OutputWriter,
    VerticalSliceWriter, write_attribute, WallTimeStop, ics_from_u!, set_ψt!

end
