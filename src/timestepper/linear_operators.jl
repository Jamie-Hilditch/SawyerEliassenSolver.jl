"""ζₜₜ = -𝓛ζ + F"""
struct 𝓛!{T}
    problem::Problem{T}
    fs_tmp::FSVariable{T}
    xs_tmp::XSVariable{T}
    xz_tmp::XZVariable{T}
    fc_tmp::FCVariable{T}
    xc_tmp::XCVariable{T}
    𝓛ζ::XZVariable{T}
end

"""Construct 𝓛! given working arrays.
We only need 4 as we are never working in sine and cosine space at the same time.
"""
function 𝓛!{T}(
    problem::Problem{T},
    fsc_array::Matrix{Complex{T}},
    xsc_array::Matrix{T},
    xz_array::Matrix{T},
    𝓛ζ_array::Matrix{T},
) where {T}
    domain = problem.domain
    return new{T}(
        problem,
        FSVariable(domain, fsc_array),
        XSVariable(domain, xsc_array),
        XZVariable(domain, xz_array),
        FCVariable(domain, fsc_array),
        XCVariable(domain, xsc_array),
        XZVariable(domain, 𝓛ζ_array),
    )
end

@inline function (L::𝓛!{T})(out::FSVariable{T}, in::FSVariable{T}) where {T}

    # unpack working arrays and the background flow
    (; problem, fs_tmp, xs_tmp, xz_tmp, fc_tmp, xc_tmp, 𝓛ζ) = 𝓛!
    (; f, Vx, Bx, Bz) = problem.background

    # if the domains are the same then everything will be inbounds
    @boundscheck consistent_domains(𝓛!, out, in) ||
        throw(ArgumentError("Domains of `𝓛!`, `out` and `in` must match."))

    # first compute ψ in xs with the inverse Laplacian, we can store this in the output array
    @inbounds ∇⁻²!(out, in)

    # first term is - Bz * ψxx
    @inbounds ∂x!(fs_tmp, out, 2) # ψxx in fs
    Tᴴ!(xs_tmp, fs_tmp) # ψxx in xs
    Τˢ!(xz_tmp, xs_tmp) # ψxx in xz
    @inbounds @. 𝓛ζ = -Bz * xz_tmp

    # second term is  2 * Bx * ψxz
    @inbounds ∂z!(fc_tmp, out) # ψz in fc
    ∂x!(fc_tmp) # ψxz in fc
    Tᴴ!(xc_tmp, fc_tmp) # ψxz in xc
    Tᶜ!(xz_tmp, xc_tmp) # ψxz in xz
    @inbounds @. 𝓛ζ += 2 * Bx * xz_tmp

    # third term is - f * (f + Vx) * ψzz
    ∂z²!(out) # ψzz in fs, we don't need ψ again so do this inplace
    Tᴴ!(xs_tmp, out) # ψzz in xs
    Tˢ!(xz_tmp, xs_tmp) # ψzz in xz
    @inbounds @. 𝓛ζ -= f * (f + Vx) * xz_tmp

    # now transform back to xs
    Tˢ!(xs_tmp, 𝓛ζ)
    Tᴴ!(out, xs_tmp)

    return nothing
end

@inline Domains.get_domain(L::𝓛!) = get_domain(L.problem)

"""𝓛♯ = 1 + aᵢᵢ h² 𝓛"""
struct 𝓛♯!{T}
    aᵢᵢ::T
    h::T
    L::𝓛!{T}
end

@propagate_inbounds function (L♯::𝓛♯!{T})(out::XSVariable, in::XSVariable)
    (; aᵢᵢ, h, L) = L♯
    return 1 + aᵢᵢ * h^2 * L(out, in)
end

@inline Domains.get_domain(A::𝓛♯!) = get_domain(A.L)

abstract type Preconditioner end
