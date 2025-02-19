function compute_ψ!(problem::Problem{T}, ψ::XZVariable{T}) where {T}
    @boundscheck consistent_domains(ψ, problem) ||
        throw(ArgumentError("`ψ` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ_FS, ψ_XS = problem.scratch.FS_tmp, problem.scratch.XS_tmp

    # first compute the streamfunction Ψ in FS space from the vorticity ζ
    ∇⁻²!(Ψ_FS, problem.state.ζ)

    # transform to XZ
    Tᴴ!(ψ_XS, ψ_FS)
    Tˢ!(Ψ, ψ_XS)

    return nothing
end

function compute_ζ!(problem::Problem{T}, ζ::XZVariable{T}) where {T}
    @boundscheck consistent_domains(ζ, problem) ||
        throw(ArgumentError("`ζ` and `problem` must have the same domain."))
    # use the scratch variables for storage
    ζ_FS, ζ_XS = problem.scratch.FS_tmp, problem.scratch.XS_tmp

    # transform to XZ
    Tᴴ!(ζ_XS, problem.state.ζ)
    Tˢ!(ζ, ζ_XS)

    return nothing
end

# define unsafe versions of the functions with no bounds checking
@inline _compute_ψ!(problem, ψ) = @inbounds compute_ψ!(problem, ψ)
@inline _compute_ζ!(problem, ζ) = @inbounds compute_ζ!(problem, ζ)

# define the output variables
"""$(TYPEDSIGNATURES)"""
ψ(problem::Problem) = OutputVariable(problem, _compute_ψ!, (:x, :z), problem.scratch.XZ_tmp)
ζ(problem::Problem) = OutputVariable(problem, _compute_ζ!, (:x, :z), problem.scratch.XZ_tmp)
