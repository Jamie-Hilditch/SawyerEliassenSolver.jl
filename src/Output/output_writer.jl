function initialise_output_file(
    filepath::String, ::Type{T}, dimensions::NamedTuple{S,<:NTuple{N,AbstractVector{T}}}, overwrite::Bool
) where {S,N,T}
    # Check if the file already exists
    if !overwrite && isfile(filepath)
        throw(
            ArgumentError(
                "File $(filepath) already exists. Set `overwrite` to true to overwrite the file.",
            ),
        )
    end

    # dimensions we always save


    # open the file
    h5open(filepath, "w") do h5
        # Create dimension scales
        grp = create_group(h5, "dimensions")
        for (dim_symbol, coords) in pairs(dimensions)
            dim_name = String(dim_symbol)
            # first create the dataset
            write_dataset(grp, dim_name, collect(coords))
            # now to convert to a dimension scale
            dset = open_dataset(grp, dim_name)
            HDF5.API.h5ds_set_scale(dset.id, dim_name)
        end

        # create the time dimension
        dset = create_dataset(grp, "time", T, dataspace(tuple(0),tuple(-1)), chunk=(1,))
        HDF5.API.h5ds_set_scale(dset.id, "time")

        # Create variables group
        create_group(h5, "variables")

    end
end

"""$(TYPEDEF)

Struct for writing output to an HDF5 file.

# Fields
$(FIELDS)
"""
struct OutputWriter{T,D}
    filepath::String
    problem::Problem{T}
    dimension_coordinates::D
    output_variables::Dict{String,OutputVariable}

    function OutputWriter(
        filepath::String,
        problem::Problem{T},
        dimension_coordinates::NamedTuple{S,<:NTuple{N,AbstractVector{T}}};
        overwrite::Bool=false,
    ) where {T,S,N}
        # add `x` and `z` to dimensions if not already defined
        default_dimensions = (
            x=xgridpoints(problem.domain),
            z=zgridpoints(problem.domain),
        )
        dimensions = merge(default_dimensions, dimension_coordinates)

        # initialise the output file
        initialise_output_file(filepath, T, dimensions, overwrite)

        # initialise an empty dict for the variables
        output_variables = Dict{String,OutputVariable}()
        return new{T,typeof(dimensions)}(filepath, problem, dimensions, output_variables)
    end
end

function OutputWriter(filepath::String, problem::Problem{T}; overwrite::Bool=false) where {T}
    return OutputWriter(filepath, problem, NamedTuple(); overwrite=overwrite)
end

function Base.show(io::IO, ::MIME"text/plain", ow::OutputWriter)
    dimensions_string = join(keys(ow.dimension_coordinates),", ")
    variables_string = join(keys(ow.output_variables),", ")
return print(
    io,
    "OutputWriter:\n",
    "  ├─────── filepath: $(ow.filepath)\n",
    "  ├──────── problem: $(summary(ow.problem))\n",
    "  ├───── dimensions: $(dimensions_string)\n",
    "  └────── variables: $(variables_string)\n",
)
end

function Base.summary(io::IO, ow::OutputWriter)
    return print(
        io, "OutputWriter to $(ow.filepath) for variables $(keys(ow.output_variables)...)"
    )
end

"""$(TYPEDSIGNATURES)

Add `OutputVariable`s to the `OutputWriter`.
Excepts keyword arguments where the key is the name of the variable and the value an `OutputVariable`.
"""
function add_output_variables!(output_writer::OutputWriter; kwargs...)
    # open the file
    h5open(output_writer.filepath, "a") do h5
        # open the dimensions and variables group
        dim_grp = open_group(h5, "dimensions")
        var_grp = open_group(h5, "variables")
        # Iterate over the keyword arguments
        for (path_symbol, output_variable) in kwargs

            path = String(path_symbol)

            # Create a dataset for the output variable
            var_dset = create_output_variable!(var_grp, path, output_variable)

            # attach dimension scales
            for (i,dim) in enumerate(output_variable.dims)
                if isnothing(dim)
                    continue
                end
                dim_dset = open_dataset(dim_grp, String(dim))
                HDF5.API.h5ds_attach_scale(var_dset.id, dim_dset.id, i)
            end
            time_dset = open_dataset(dim_grp, "time")
            HDF5.API.h5ds_attach_scale(var_dset.id, time_dset.id, ndims(output_variable) + 1)

            # Store the output variable in the output writer
            output_writer.output_variables[path] = output_variable


        end
    end
end
