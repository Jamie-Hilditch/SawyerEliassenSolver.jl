
# define a Named Tuple of functions for computing different variables
# functions are defined in compute_variables.jl and should return the array (typical problem.scratch.rtmp1) in which the variable is stored
include("compute_variables.jl")
const Variables = (
    u = compute_u!,
    w = compute_w!,
    v = get_v,
    b = get_b,
    ux = compute_ux!,
    uz = compute_uz!,
    wx = compute_wx!,
    vz = compute_vz!,
    bz = compute_bz!,
    ψ = compute_ψ!
)

function extend_dims(dset :: HDF5.Dataset)
    dims,_ = HDF5.get_extent_dims(dset)
    new_dims = (dims[1:end-1]...,dims[end] + 1)
end

const HDF5_OBJECT = Union{HDF5.File,HDF5.Group}

abstract type AbstractWriter end

mutable struct OutputWriter <: AbstractWriter
    const fid :: HDF5.File
    const object_handle :: HDF5_OBJECT
    const problem :: SE.Problem
    const state :: SE.State
    const variables :: NTuple{N,Symbol} where N

    function OutputWriter(filepath :: String, problem :: SE.Problem, state :: SE.State, variables :: Symbol ...; group :: Union{String,Nothing} = nothing, overwrite :: Bool = false)
        size(problem) == size(state) || throw(ArgumentError("problem and state don't have the same size"))
        
        # open the group
        fid = h5open(filepath,"cw")
        if group === nothing 
            object_handle = fid 
        else
            haskey(fid,group) || create_group(fid,group)
            object_handle = fid[group]
        end

        # clear the group
        for key ∈ keys(object_handle)
            overwrite ? delete_object(object_handle,key) : throw(ErrorException("$(key) already exists"))
        end

        # clear attributes 
        for key ∈ keys(HDF5.attrs(object_handle))
            overwrite ? HDF5.delete_attribute(object_handle,key) : throw(ErrorException("$(key) already exists"))
        end

        # create the datasets
        # directly save x and z
        write(object_handle,"x",problem.grid.x)
        write(object_handle,"z",problem.grid.z)
        write(object_handle,"Vx",problem.background_flow.Vx)
        write(object_handle,"Bx",problem.background_flow.Bx)
        write(object_handle,"Bz",problem.background_flow.Bz)

        # time
        dtype = HDF5.datatype(Float64)
        dspace = dataspace((0,),max_dims=(-1,))
        create_dataset(object_handle,"t",dtype,dspace,chunk = (1,))

        for var ∈ variables 
            haskey(Variables,var) || throw(ArgumentError("variable $(var) has not been implemented"))
            dspace = dataspace((size(problem)...,0),max_dims= (size(problem)...,-1))
            create_dataset(object_handle,String(var),dtype,dspace,chunk = (size(problem)...,1))
        end


        output_writer = new(fid,object_handle,problem,state,variables)
        finalizer(close,output_writer)
    end
end 

function (ow :: OutputWriter)()
    # save time
    dset = ow.object_handle["t"]
    HDF5.set_extent_dims(dset,(HDF5.get_extent_dims(dset)[1][1] + 1,))
    dset[end] = ow.state.clock.t

    for var ∈ ow.variables
        data = get(Variables,var,missing_variable)(ow.problem,ow.state)

        dset = ow.object_handle[String(var)]
        HDF5.set_extent_dims(dset,extend_dims(dset))
        dset[:,:,end] = data
    end
end

mutable struct VerticalSliceWriter <: AbstractWriter
    const fid :: HDF5.File
    const object_handle :: HDF5_OBJECT
    const problem :: SE.Problem
    const state :: SE.State
    const slices :: NTuple{M,Pair{Int,String}} where M
    const variables :: NTuple{N,Symbol} where N
    

    function VerticalSliceWriter(filepath :: String, problem :: SE.Problem, state :: SE.State,
        slices :: NTuple{N,Pair{Int,String}}, variables :: Symbol ...;
        group :: Union{String,Nothing} = nothing, overwrite :: Bool = false) where N
        size(problem) == size(state) || throw(ArgumentError("problem and state don't have the same size"))
        
        # open root 
        fid = h5open(filepath,"cw")
        if group === nothing 
            object_handle = fid 
        else
            haskey(fid,group) || create_group(fid,group)
            object_handle = fid[group]
        end

        # clear the group
        for key ∈ keys(object_handle)
            overwrite ? delete_object(object_handle,key) : throw(ErrorException("$(key) already exists"))
        end

        # clear attributes 
        for key ∈ keys(HDF5.attrs(object_handle))
            overwrite ? HDF5.delete_attribute(object_handle,key) : throw(ErrorException("$(key) already exists"))
        end

        # directly save x and z
        haskey(object_handle,"Problem") || create_group(object_handle,"Problem")
        problem_group = object_handle["Problem"]
        write(problem_group,"x",problem.grid.x)
        write(problem_group,"z",problem.grid.z)
        write(problem_group,"Vx",problem.background_flow.Vx)
        write(problem_group,"Bx",problem.background_flow.Bx)
        write(problem_group,"Bz",problem.background_flow.Bz)

        # create dataset for time
        dtype = HDF5.datatype(Float64)
        dspace = dataspace((0,),max_dims=(-1,))
        create_dataset(problem_group,"t",dtype,dspace,chunk = (1,))

        # now loop through slices 
        for (idx,group) in slices

            group != "Problem" || throw(ArgumentError("'Problem' is a reserved group name"))
            haskey(object_handle,group) || create_group(object_handle,group)
            gid = object_handle[group]

            # clear the group
            for key ∈ keys(gid)
                overwrite ? delete_object(gid,key) : throw(ErrorException("$(key) already exists"))
            end

            # clear attributes 
            for key ∈ keys(HDF5.attrs(gid))
                overwrite ? HDF5.delete_attribute(gid,key) : throw(ErrorException("$(key) already exists"))
            end

            # write x attribute
            # this also checks that idx is inbounds
            HDF5.write_attribute(gid,"x",problem.grid.x[idx])

            for var ∈ variables 
                haskey(Variables,var) || throw(ArgumentError("variable $(var) has not been implemented"))
                dspace = dataspace((size(problem.grid.z)...,0),max_dims= (size(problem.grid.z)...,-1))
                create_dataset(gid,String(var),dtype,dspace,chunk = (size(problem.grid.z)...,1))
            end

        end

        vertical_slice_writer = new(fid,object_handle,problem,state,slices,variables)
        finalizer(close,vertical_slice_writer)
    end
end 

function (vsw :: VerticalSliceWriter)()
    # save time
    dset = vsw.object_handle["Problem"]["t"]
    HDF5.set_extent_dims(dset,(HDF5.get_extent_dims(dset)[1][1] + 1,))
    dset[end] = vsw.state.clock.t

    for var ∈ vsw.variables
        data = get(Variables,var,missing_variable)(vsw.problem,vsw.state)

        for (idx,group) in vsw.slices
            dset = vsw.object_handle[group][String(var)]
            HDF5.set_extent_dims(dset,extend_dims(dset))
            @inbounds dset[:,end] = data[idx,:]
        end
    end
end

Base.close(aw :: AbstractWriter) = close(aw.fid)

"""
Thin wrapper around HDF5.write_attribute
"""
function write_attribute(ow :: OutputWriter,attributes :: Pair{String,T}...) where T
    for attribute in attributes
        HDF5.write_attribute(ow.object_handle, attribute.first, attribute.second)
    end
    return nothing 
end

function write_attribute(vsw :: VerticalSliceWriter,attributes :: Pair{String,T}...) where T
    for attribute in attributes
        HDF5.write_attribute(vsw.object_handle, attribute.first, attribute.second)
    end
    return nothing 
end