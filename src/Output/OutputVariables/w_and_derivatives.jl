function compute_w!(problem::Problem{T}, w::XZVariable{T}) where {T}
    @boundscheck consistent_domains(w, problem) || throw(ArgumentError("`w` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, w_XS = problem.scratch.FS_tmp, problem.scratch.XS_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute w = ψx in FS space - this is an inplace operation
    @inbounds ∂x!(ψ)

    # transform to XZ
    Tᴴ!(w_XS, ψ)
    Tˢ!(w, w_XS)
    return nothing
end

function compute_∂w∂x!(problem::Problem{T}, ∂w∂x::XZVariable{T}) where {T}
    @boundscheck consistent_domains(∂w∂x, problem) || throw(ArgumentError("`∂w∂x` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, ∂w∂x_XS = problem.scratch.FS_tmp, problem.scratch.XS_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute wx = ψxx in FS space - this is an inplace operation
    @inbounds ∂x!(ψ,2)

    # transform to XZ
    Tᴴ!(∂w∂x_XS, ψ)
    Tˢ!(∂w∂x, ∂w∂x_XS)
    return nothing
end

function compute_∂²w∂x²!(problem::Problem{T}, ∂²w∂x²::XZVariable{T}) where {T}
    @boundscheck consistent_domains(∂²w∂x², problem) || throw(ArgumentError("`∂²w∂x²` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, ∂²w∂x²_XS = problem.scratch.FS_tmp, problem.scratch.XS_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute wxx = ψxxx in FS space - this is an inplace operation
    @inbounds ∂x!(ψ,3)

    # transform to XZ
    Tᴴ!(∂²w∂x²_XS, ψ)
    Tˢ!(∂²w∂x², ∂²w∂x²_XS)
    return nothing
end

function compute_∂w∂z!(problem::Problem{T}, ∂w∂z::XZVariable{T}) where {T}
    @boundscheck consistent_domains(∂w∂z, problem) || throw(ArgumentError("`∂w∂z` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, ∂w∂z_FC, ∂w∂z_XC = problem.scratch.FS_tmp, problem.scratch.FC_tmp, problem.scratch.XC_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute w = ψx in FS space - this is an inplace operation
    @inbounds ∂x!(ψ)

    # compute z derivative
    @inbounds ∂z!(∂w∂z_FC, ψ)

    # transform to XZ
    Tᴴ!(∂w∂z_XC, ∂w∂z_FC)
    Tᶜ!(∂w∂z, ∂w∂z_XC)
    return nothing
end

function compute_∂²w∂z²!(problem::Problem{T}, ∂²w∂z²::XZVariable{T}) where {T}
    @boundscheck consistent_domains(∂²w∂z², problem) || throw(ArgumentError("`∂²w∂z²` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, ∂²w∂z²_XS = problem.scratch.FS_tmp, problem.scratch.XS_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute w = ψx in FS space - this is an inplace operation
    @inbounds ∂x!(ψ)

    # compute second z derivative
    @inbounds ∂z²!(ψ)

    # transform to XZ
    Tᴴ!(∂²w∂z²_XS, ψ)
    Tˢ!(∂²w∂z², ∂²w∂z²_XS)
    return nothing
end

function compute_∂²w∂x∂z!(problem::Problem{T}, ∂²w∂x∂z::XZVariable{T}) where {T}
    @boundscheck consistent_domains(∂²w∂x∂z, problem) || throw(ArgumentError("`∂²w∂x∂z` and `problem` must have the same domain."))
    # use the scratch variables for storage
    Ψ, ∂²w∂x∂z_FC, ∂²w∂x∂z_XC = problem.scratch.FS_tmp, problem.scratch.FC_tmp, problem.scratch.XC_tmp

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, problem.state.ζ)

    # now compute wx = ψxx in FS space - this is an inplace operation
    @inbounds ∂x!(ψ,2)

    # compute z derivative
    @inbounds ∂z!(∂²w∂x∂z_FC, ψ)

    # transform to XZ
    Tᴴ!(∂²w∂x∂z_XC, ∂²w∂x∂z_FC)
    Tᶜ!(∂²w∂x∂z, ∂²w∂x∂z_XC)
    return nothing
end

# functions without boundschecks
@inline _compute_w!(problem, w) = @inbounds compute_w!(problem, w)
@inline _compute_∂w∂x!(problem, ∂w∂x) = @inbounds compute_∂w∂x!(problem, ∂w∂x)
@inline _compute_∂²w∂x²!(problem, ∂²w∂x²) = @inbounds compute_∂²w∂x²!(problem, ∂²w∂x²)
@inline _compute_∂w∂z!(problem, ∂w∂z) = @inbounds compute_∂w∂z!(problem, ∂w∂z)
@inline _compute_∂²w∂z²!(problem, ∂²w∂z²) = @inbounds compute_∂²w∂z²!(problem, ∂²w∂z²)
@inline _compute_∂²w∂x∂z!(problem, ∂²w∂x∂z) = @inbounds compute_∂²w∂x∂z!(problem, ∂²w∂x∂z)

"""$(TYPEDSIGNATURES)"""
w(problem::Problem) = OutputVariable(_compute_w!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂w∂x(problem::Problem) = OutputVariable(_compute_∂w∂x!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂²w∂x²(problem::Problem) = OutputVariable(_compute_∂²w∂x²!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂w∂z(problem::Problem) = OutputVariable(_compute_∂w∂z!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂²w∂z²(problem::Problem) = OutputVariable(_compute_∂²w∂z²!, (:x, :z), problem.scratch.XZ_tmp)
"""$(TYPEDSIGNATURES)"""
∂²w∂x∂z(problem::Problem) = OutputVariable(_compute_∂²w∂x∂z!, (:x, :z), problem.scratch.XZ_tmp)
