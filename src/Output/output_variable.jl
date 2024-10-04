"""$(TYPEDEF)

A struct to create output variables.

# Fields
$(FIELDS)

`func` is a function with signature
    `func(problem::Problem, output_array::AbstractArray{T,N}, args...)`
that computes the output field in-place in `output_array`.
`args` are the arguments to `func`.
`dimension_labels` is a tuple of dimension labels. When added to an [`OutputWriter`](@ref)
these labels are used to attach dimension scales (coordinates) to the output variable.

!!! note
    Dimension labels are not required to be associated with a dimension scale.
    But if they are then the size of the dimension must match.
"""
struct OutputVariable{T,N,O,A} <: AbstractArray{O,N}
    problem::Problem{T}
    func::Function
    dimension_labels::NTuple{N,Union{Symbol,Nothing}}
    output_array::AbstractArray{O,N}
    args::A
end

"""$(TYPEDSIGNATURES)

Create an output variable with a function that takes no additional arguments.
"""
function OutputVariable(
    problem::Problem{T},
    func::Function,
    dimension_labels::NTuple{N,Union{Symbol,Nothing}},
    output_array::AbstractArray{O,N},
) where {T,N,O}
    return OutputVariable{T,N,O,Tuple{}}(
        problem, func, dimension_labels, output_array, tuple()
    )
end

"""$(TYPEDSIGNATURES)

Create an output variable creating the output array of shape `sizes` and type `type`.
"""
function OutputVariable(
    problem::Problem{T},
    func::Function,
    dimension_labels::NTuple{N,Union{Symbol,Nothing}},
    sizes::NTuple{N,Integer},
    type::Type,
    args,
) where {T,N}
    output_array = zeros(type, sizes)
    return OutputVariable{T,N,type,typeof(args)}(
        problem, func, dimension_labels, output_array, args
    )
end

"""$(TYPEDSIGNATURES)"""
function OutputVariable(
    problem::Problem{T},
    func::Function,
    dimension_labels::NTuple{N,Union{Symbol,Nothing}},
    sizes::NTuple{N,Integer},
    type::Type,
) where {T,N}
    output_array = zeros(type, sizes)
    return OutputVariable{T,N,type,Tuple{}}(
        problem, func, dimension_labels, output_array, tuple()
    )
end

"""$(TYPEDSIGNATURES)

Create an output variable specifying the dimension labels and sizes in a named tuple.
"""
function OutputVariable(
    problem::Problem{T},
    func::Function,
    dimensions::NamedTuple{S,NTuple{N,Integer}},
    type::Type,
    args,
) where {T,S,N}
    dimension_labels = keys(dimensions)
    sizes = values(dimensions)
    return OutputVariable(problem, func, dimension_labels, sizes, type, args)
end

"""$(TYPEDSIGNATURES)"""
function OutputVariable(
    problem::Problem{T},
    func::Function,
    dimensions::NamedTuple{S,NTuple{N,Integer}},
    type::Type,
) where {T,S,N}
    dimension_labels = keys(dimensions)
    sizes = values(dimensions)
    return OutputVariable(problem, func, dimension_labels, sizes, type)
end

Base.size(output_variable::OutputVariable) = size(output_variable.output_array)
Base.eltype(output_variable::OutputVariable) = eltype(output_variable.output_array)
Base.ndims(output_variable::OutputVariable) = ndims(output_variable.output_array)

function create_output_variable!(
    h5::HDF5.File, path::String, output_variable::OutputVariable
)
    # Create the dataspace
    dimensions_sizes = size(output_variable)
    dims = tuple(dimensions_sizes..., 0)
    maxdims = tuple(dimensions_sizes..., -1)
    dspace = dataspace(dims, maxdims)
    chunks = tuple(dimensions_sizes..., 1)
    N = ndims(output_variable) + 1 # total number of dimensions including time

    # Create a dataset for the output variable
    var_dset = create_dataset(h5, path, eltype(output_variable), dspace; chunk=chunks)

    # label dimensions and attach dimension scales
    # Note the low-level API is just a thin wrapper around the C API
    # so we need to use zero-indexed row major indexing i.e. i -> N - i
    for (i, dim) in enumerate(output_variable.dimension_labels)
        # get dimension name as a string
        dim_name = isnothing(dim) ? "$(path)_$(i)" : String(dim)

        # label the dimension
        label = Vector{UInt8}(dim_name)
        HDF5.API.h5ds_set_label(var_dset.id, N - i, label)

        # if there is no dimension scale associated with this dimension then we are done
        if !haskey(h5, dim_name)
            continue
        end

        # open the dimension scale dataset
        dim_dset = open_dataset(h5, dim_name)

        # check it is actually a scale
        if !HDF5.API.h5ds_is_scale(dim_dset.id)
            @warn "Not attaching dimension $dim_name to $path as it is not a dimension scale"
            continue
        end

        # enforce that the length matches the dimension size
        if length(dim_dset) != dimensions_sizes[i]
            @warn "Not attaching dimension $dim_name as it has length $(length(dim_dset)) but dimension $i on $path has length $(dimensions_sizes[i])"
            continue
        end

        HDF5.API.h5ds_attach_scale(var_dset.id, dim_dset.id, N - i)
    end

    # finally attach the time dimension scale
    time_dset = open_dataset(h5, "time")
    HDF5.API.h5ds_attach_scale(var_dset.id, time_dset.id, 0)
    return HDF5.API.h5ds_set_label(var_dset.id, 0, Vector{UInt8}("time"))
end

"""$(TYPEDSIGNATURES)

Compute the output field in-place in `output_variable.output_array`.
"""
function compute!(output_variable::OutputVariable)
    return output_variable.func(
        output_variable.problem, output_variable.output_array, output_variable.args...
    )
end

# methods for reading data from an output variables
# Base.IndexStyle(::OutputVariable) = IndexLinear()
@propagate_inbounds Base.getindex(v::OutputVariable, i::Int) = getindex(v.output_array, i)
@propagate_inbounds Base.getindex(v::OutputVariable, I::Vararg{Int,N}) where {N} =
    getindex(v.output_array, I...)
@propagate_inbounds Base.getindex(v::OutputVariable, I...) = getindex(v.output_array, I...)
