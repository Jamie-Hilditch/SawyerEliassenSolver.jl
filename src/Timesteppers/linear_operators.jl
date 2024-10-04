"""$(TYPEDSIGNATURES)

The Sawyer-Eliassen operator 𝓛 is defined as `ζₜₜ = -𝓛ζ + F`
"""
@inline function 𝓛!(problem::Problem{T}, out::FSVariable{T}, in::FSVariable{T}) where {T}

    # if the domains are the same then everything will be inbounds
    @boundscheck consistent_domains(problem, out, in) ||
        throw(ArgumentError("Domains of `problem`, `out` and `in` must match."))

    # unpack working arrays and the background flow
    (; FS_tmp, FC_tmp, XS_tmp, XC_tmp, XZ_tmp, XZ_tmp2) = problem.scratch
    (; f, Vx, Bx, Bz) = problem.background

    # first compute ψ in fs with the inverse Laplacian, we can store this in the output array
    @inbounds ∇⁻²!(out, in)

    # now we build 𝓛ζ in xz
    # use XZ_tmp2 to store 𝓛ζ
    𝓛ζ = XZ_tmp2

    # first term is Bz * ψxx
    @inbounds ∂x!(FS_tmp, out, 2) # ψxx in fs
    Tᴴ!(XS_tmp, FS_tmp) # ψxx in xs
    Tˢ!(XZ_tmp, XS_tmp) # ψxx in xz
    @inbounds @. 𝓛ζ = Bz * XZ_tmp

    # second term is  -2 * Bx * ψxz
    @inbounds ∂z!(FC_tmp, out) # ψz in fc
    ∂x!(FC_tmp) # ψxz in fc
    Tᴴ!(XC_tmp, FC_tmp) # ψxz in xc
    Tᶜ!(XZ_tmp, XC_tmp) # ψxz in xz
    @inbounds @. 𝓛ζ -= 2 * Bx * XZ_tmp

    # third term is f * (f + Vx) * ψzz
    ∂z²!(out) # ψzz in fs, we don't need ψ again so do this inplace
    Tᴴ!(XS_tmp, out) # ψzz in xs
    Tˢ!(XZ_tmp, XS_tmp) # ψzz in xz
    @inbounds @. 𝓛ζ += f * (f + Vx) * XZ_tmp

    # finally transform 𝓛ζ to fs
    Tˢ!(XS_tmp, 𝓛ζ)
    Tᴴ!(out, XS_tmp)

    return nothing
end

"""$(TYPEDSIGNATURES)

𝓛ᴵ = 1 + aᵢᵢ h² 𝓛
"""
@propagate_inbounds function 𝓛ᴵ!(
    problem::Problem{T}, out::FSVariable{T}, in::FSVariable{T}, aᵢᵢh²::T
) where {T}
    𝓛!(problem, out, in)
    @inbounds out .*= aᵢᵢh²
    @inbounds out .+= in
    return nothing
end
