"""$(TYPEDEF)

A struct to create output variables. `func` is a function with signature
`func(problem::Problem, output::AbstractArray{T,N}, args...) -> AbstractArray` that computes the output field.
`args` are the arguments to `func` and `dims` are the dimensions of the output field.

# Fields
$(FIELDS)
"""
struct OutputVariable{A,T,N}
    func::Function
    output_array::AbstractArray{T,N}
    args::A
    dims::NTuple{N,Union{Symbol,Int}}
end

Base.size(output_variable::OutputVariable) = size(output_variable.output_array)
Base.eltype(output_variable::OutputVariable) = eltype(output_variable.output_array)
Base.ndims(output_variable::OutputVariable) = ndims(output_variable.output_array)

"""$(TYPEDSIGNATURES)

Check if the dimensions of the output variable exist and have the same length as the dimension coordinates.
"""
function validate_dimensions(output_variable::OutputVariable{A,T,N}, dimension_coordinates::NamedTuple{Vector{T}}) where {A,T,N}
    # Check if the dimensions of the output variable exist and have the correct length
    for (dim, dim_name) in enumerate(output_variable.dims)
        if !haskey(possible_dims, dim_name)
            throw(ArgumentError("Dimension $dim_name on `output_variable` not found in `dimension_coordinates`"))
        end
        dimension_length = length(dimension_coordinates[dim_name])
        if dimension_length != size(output_variable.output_array, dim)
            throw(ArgumentError("Dimension length mismatch for dimension $dim_name: expected $(size(output_variable.output, dim)), got $dimension_length"))
        end
    end
    return nothing
end

function create_output_variable!(parent::HDF5.Group, path::String, output_variable::OutputVariable)
    # Create the dataspace
    dimensions_sizes = size(output_variable)
    dims = tuple(dimensions_sizes...,0)
    maxdims = tuple(dimensions_sizes...,-1)
    dspace = dataspace(dims,maxdims)
    # Create a dataset for the output variable
    create_dataset(parent, path, eltype(output_variable), dspace)
end
