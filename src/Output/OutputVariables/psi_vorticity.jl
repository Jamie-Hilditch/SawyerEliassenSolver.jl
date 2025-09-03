function compute_Ψ!(problem::Problem{T}, Ψ::XZVariable{T}) where {T}
    @boundscheck consistent_domains(Ψ, problem) ||
        throw(ArgumentError("`Ψ `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ_FS, Ψ_XS = problem.scratch.FS_tmp, problem.scratch.XS_tmp

    # first compute the streamfunction Ψ in FS space from the vorticity ζ
    solve_poisson!(Ψ_FS, problem.state.ζ)

    # transform to XZ
    Tᴴ!(Ψ_XS, Ψ_FS)
    Tˢ!(Ψ, Ψ_XS)

    return nothing
end

function compute_ζ!(problem::Problem{T}, ζ::XZVariable{T}) where {T}
    @boundscheck consistent_domains(ζ, problem) ||
        throw(ArgumentError("`ζ` and `problem` must have the same domain."))
    # use the scratch variables for storage
    ζ_XS = problem.scratch.XS_tmp

    # transform to XZ
    Tᴴ!(ζ_XS, problem.state.ζ)
    Tˢ!(ζ, ζ_XS)

    return nothing
end

# define unsafe versions of the functions with no bounds checking
@inline _compute_Ψ!(problem, Ψ) = @inbounds compute_Ψ!(problem, Ψ)
@inline _compute_ζ!(problem, ζ) = @inbounds compute_ζ!(problem, ζ)

# define the output variables
"""$(TYPEDSIGNATURES)

The streamfunction Ψ.
"""
Ψ(problem::Problem) = OutputVariable(problem, _compute_Ψ!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)

The vorticity ζ.
"""
ζ(problem::Problem) = OutputVariable(problem, _compute_ζ!, (:x, :z), problem.scratch.XZ_tmp)
