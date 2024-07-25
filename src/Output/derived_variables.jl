"""$(TYPEDSIGNATURES)"""
function compute_u_from_Ψ!(
    u::XZVariable{T}, Ψ::FSVariable{T}, FC_tmp::FCVariable{T}, XC_tmp::XCVariable{T}
) where {T}
    @boundscheck consistent_domains(u, Ψ, FC_tmp, XC_tmp) || throw(
        ArgumentError("`u`, `FC_tmp`, `XC_tmp` and `Ψ` must have the same domain.")
    )
    @inbounds ∂z!(FC_tmp, Ψ)
    Tᴴ!(XC_tmp, FC_tmp)
    Tᶜ!(u, XC_tmp)
    u *= -1
    return nothing
end

"""$(TYPEDSIGNATURES)"""
function compute_w_from_Ψ!(
    w::XZVariable{T}, Ψ::FSVariable{T}, FS_tmp::FSVariable{T}, XS_tmp::XSVariable{T}
) where {T}
    @boundscheck consistent_domains(w, Ψ, FS_tmp, XS_tmp) || throw(
        ArgumentError("`w`, `FS_tmp`, `XS_tmp` and `Ψ` must have the same domain.")
    )
    @inbounds ∂x!(FS_tmp, Ψ)
    Tᴴ!(XS_tmp, FS_tmp)
    Tˢ!(w, XS_tmp)
    return nothing
end

"""$(TYPEDSIGNATURES)"""
function compute_u_and_w!(scratch::Scratch{T}, ζ::FSVariable{T}) where {T}
    @boundscheck consistent_domains(ζ, scratch) ||
        throw(ArgumentError("`ζ` and `scratch` must have the same domain."))
    # use the scratch variables for storage
    Ψ, u, w = scratch.FS_tmp, scratch.XZ_tmp, scratch.XZ_tmp2

    # first compute the streamfunction Ψ from the vorticity ζ
    ∇⁻²!(Ψ, ζ)
    # next get u
    compute_u_from_Ψ!(u, Ψ, scratch.FC_tmp, scratch.XC_tmp)
    # finally get w
    # we reuse the Ψ variable for w in FS which is equivalent to taking the x derivative inplace
    compute_w_from_Ψ!(w, Ψ, Ψ, scratch.XS_tmp)

    return u, w
end

"""$(TYPEDSIGNATURES)"""
compute_u_and_w!(problem::Problem{T}) where {T} = compute_u_and_w!(problem.scratch, problem.state.ζ)
