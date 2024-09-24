
function compute_u!(problem::Problem{T}, u::XZVariable{T}) where {T}
    @boundscheck consistent_domains(u, problem) || throw(ArgumentError("`u` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, u_FC, u_XC = problem.scratch.FS_tmp, problem.scratch.FC_tmp, problem.scratch.XC_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute -u = ψz in FC space
    @inbounds ∂z!(u_FC, Ψ)

    # transform to XZ
    Tᴴ!(u_XC, u_FC)
    Tᶜ!(u, u_XC)

    # negate to get u
    u *= -1
    return nothing
end

function compute_∂u∂x!(problem::Problem{T}, ∂u∂x::XZVariable{T}) where {T}
    @boundscheck consistent_domains(∂u∂x, problem) || throw(ArgumentError("`∂u∂x` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, u_FC, ∂u∂x_XC = problem.scratch.FS_tmp, problem.scratch.FC_tmp, problem.scratch.XC_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute -u = ψz in FC space
    @inbounds ∂z!(u_FC, Ψ)

    # compute the x derivative in place
    @inbounds ∂x!(u_FC)

    # transform to XZ
    Tᴴ!(∂u∂x_XC, u_FC)
    Tᶜ!(∂u∂x, ∂u∂x_XC)

    # negate to get ∂u/∂x
    ∂u∂x *= -1
    return nothing
end

function compute_∂²u∂x²!(problem::Problem{T}, ∂²u∂x²::XZVariable{T}) where {T}
    @boundscheck consistent_domains(∂²u∂x², problem) || throw(ArgumentError("`∂²u∂x²` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, u_FC, ∂²u∂x²_XC = problem.scratch.FS_tmp, problem.scratch.FC_tmp, problem.scratch.XC_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute -u = ψz in FC space
    @inbounds ∂z!(u_FC, Ψ)

    # compute the x derivatives in place
    @inbounds ∂x!(u_FC,2)

    # transform to XZ
    Tᴴ!(∂²u∂x²_XC, u_FC)
    Tᶜ!(∂u∂x, ∂²u∂x²_XC)

    # negate to get ∂²u/∂x²
    ∂²u∂x² *= -1
    return nothing
end

function compute_∂u∂z!(problem::Problem{T}, ∂u∂z::XZVariable{T}) where {T}
    @boundscheck consistent_domains(∂u∂z, problem) || throw(ArgumentError("`∂u∂z` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, ∂u∂z_XS = problem.scratch.FS_tmp, problem.scratch.XS_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute -∂u/∂z = ∂²ψ/∂z² in FS space
    @inbounds ∂z²!(Ψ)

    # transform to XZ
    Tᴴ!(∂u∂z_XS, ψ)
    Tˢ!(∂u∂z, ∂u∂z_XS)

    # negate to get ∂u/∂z
    ∂u∂z *= -1
    return nothing
end

function compute_∂²u∂z²!(problem::Problem{T}, ∂²u∂z²::XZVariable{T}) where {T}
    @boundscheck consistent_domains(u, problem) || throw(ArgumentError("`∂²u∂z²` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, u_FC, u_XC = problem.scratch.FS_tmp, problem.scratch.FC_tmp, problem.scratch.XC_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute -u = ψz in FC space
    @inbounds ∂z!(u_FC, Ψ)

    # compute -∂²u/∂z²
    @inbounds ∂z²!(u_FC)

    # transform to XZ
    Tᴴ!(∂²u∂z²_XC, u_FC)
    Tᶜ!(∂²u∂z², ∂²u∂z²_XC)

    # negate to get u
    ∂²u∂z² *= -1
    return nothing
end

function compute_∂²u∂x∂z!(problem::Problem{T}, ∂²u∂x∂z::XZVariable{T}) where {T}
    @boundscheck consistent_domains(∂²u∂x∂z, problem) || throw(ArgumentError("`∂²u∂x∂z` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, ∂²u∂x∂z_XS = problem.scratch.FS_tmp, problem.scratch.XS_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute -∂u/∂z = ∂²ψ/∂z² in FS space
    @inbounds ∂z²!(Ψ)

    # compute the x derivative in place
    @inbounds ∂x!(Ψ)

    # transform to XZ
    Tᴴ!(∂²u∂x∂z_XS, ψ)
    Tˢ!(∂²u∂x∂z, ∂²u∂x∂z_XS)

    # negate to get ∂²u/∂x∂z
    ∂²u∂x∂z *= -1
    return nothing
end


# define unsafe versions of the functions with no bounds checking
@inline _compute_u!(problem, u) = @inbounds compute_u!(problem, u)
@inline _compute_∂u∂x!(problem, ∂u∂x) = @inbounds compute_∂u∂x!(problem, ∂u∂x)
@inline _compute_∂²u∂x²!(problem, ∂²u∂x²) = @inbounds compute_∂²u∂x²!(problem, ∂²u∂x²)
@inline _compute_∂u∂z!(problem, ∂u∂z) = @inbounds compute_∂u∂z!(problem, ∂u∂z)
@inline _compute_∂²u∂z²!(problem, ∂²u∂z²) = @inbounds compute_∂²u∂z²!(problem, ∂²u∂z²)
@inline _compute_∂²u∂x∂z!(problem, ∂²u∂x∂z) = @inbounds compute_∂²u∂x∂z!(problem, ∂²u∂x∂z)

# define the output variables
"""$(TYPEDSIGNATURES)"""
u(problem::Problem) = OutputVariable(_compute_u!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂u∂x(problem::Problem) = OutputVariable(_compute_∂u∂x!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂²u∂x²(problem::Problem) = OutputVariable(_compute_∂²u∂x²!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂u∂z(problem::Problem) = OutputVariable(_compute_∂u∂z!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂²u∂z²(problem::Problem) = OutputVariable(_compute_∂²u∂z²!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂²u∂x∂z(problem::Problem) = OutputVariable(_compute_∂²u∂x∂z!, (:x, :z), problem.scratch.XZ_tmp)
