"""convenience functions"""

module Tools

using HDF5
using LinearAlgebra

import ..SawyerEliassen as SE

include("initial_conditions.jl")
include("output_writer.jl")
include("wall_time_stop.jl")

export
    OutputWriter,
    VerticalSliceWriter,
    write_attribute,
    WallTimeStop,
    ics_from_u!,
    set_ψt!

end