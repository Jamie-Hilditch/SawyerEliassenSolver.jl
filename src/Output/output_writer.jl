function initialise_output_file(
    ::Type{T}, filepath::String, dimensions::NamedTuple, overwrite::Bool
) where {T}
    # Check if the file already exists
    if !overwrite && isfile(filepath)
        throw(
            ArgumentError(
                "File $(filepath) already exists. Set `overwrite` to true to overwrite the file.",
            ),
        )
    end

    # open the file
    h5open(filepath, "w") do h5
        # Create dimension coordinates
        for (dim_symbol, coords) in pairs(dimensions)
            dim_name = String(dim_symbol)
            # first create the dataset
            write_dataset(h5, dim_name, collect(coords))
            # now to convert to a dimension scale
            dset = open_dataset(h5, dim_name)
            HDF5.API.h5ds_set_scale(dset.id, dim_name)
        end

        # create the time dimension
        dset = create_dataset(h5, "time", T, dataspace(tuple(0),tuple(-1)), chunk=(1,))
        HDF5.API.h5ds_set_scale(dset.id, "time")
    end
end

"""$(TYPEDEF)

Struct for writing output to an HDF5 file.
"""
struct OutputWriter{T,D}
    problem::Problem{T}
    filepath::String
    coordinates::D
    _output_variables::Dict{String,OutputVariable}

    @doc """$(TYPEDSIGNATURES)"""
    function OutputWriter(
        problem::Problem{T},
        filepath::String,
        coordinates::NamedTuple;
        overwrite::Bool=false,
    ) where {T}
        # add `x` and `z` to dimensions if not already defined
        default_coordinates = (
            x=xgridpoints(problem.domain),
            z=zgridpoints(problem.domain),
        )
        all_coordinates = merge(default_coordinates, coordinates)

        # initialise the output file
        initialise_output_file(T, filepath, all_coordinates, overwrite)

        # initialise an empty dict for the variables
        output_variables = Dict{String,OutputVariable}()
        return new{T,typeof(all_coordinates)}(problem, filepath, all_coordinates, output_variables)
    end
end

"""$(TYPEDSIGNATURES)"""
function OutputWriter(problem::Problem{T}, filepath::String; overwrite::Bool=false) where {T}
    return OutputWriter(problem, filepath, NamedTuple(); overwrite=overwrite)
end

function Base.show(io::IO, ::MIME"text/plain", ow::OutputWriter)
    coordinates_string = join(keys(ow.coordinates),", ")
    variables_string = join(keys(ow._output_variables),", ")
return print(
    io,
    "OutputWriter:\n",
    "  ├─────── filepath: $(ow.filepath)\n",
    "  ├──────── problem: $(summary(ow.problem))\n",
    "  ├──── coordinates: $(coordinates_string)\n",
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
    h5open(output_writer.filepath, "r+") do h5

        # Iterate over the keyword arguments
        for (path_symbol, output_variable) in kwargs

            path = String(path_symbol)

            # Check if the variable already exists
            if haskey(h5, path)
                throw(ArgumentError("$path already exists in the output file. Variables can be deleted using `delete!(output_writer, variable_name)`."))
            end

            # Create a dataset for the output variable
            create_output_variable!(h5, path, output_variable)

            # Store the output variable in the output writer
            output_writer._output_variables[path] = output_variable


        end
    end
end

function Base.delete!(output_writer::OutputWriter, variable_name::String)
    # check if the key is a dimension
    if haskey(output_writer.coordinates, Symbol(variable_name))
        throw(ArgumentError("Cannot delete a coordinate from the output file."))
    end
    # delete the variable from the output file
    h5open(output_writer.filepath, "a") do h5
        delete_object(h5, variable_name)
    end
    # delete the variable from the output writer dict
    delete!(output_writer._output_variables, variable_name)
end

Base.delete!(output_writer::OutputWriter, variable_name::Symbol) = delete!(output_writer, String(variable_name))
