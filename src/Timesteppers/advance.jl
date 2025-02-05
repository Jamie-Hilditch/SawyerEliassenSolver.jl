@inline function get_ζ_forcing!(problem::Problem{T}, out::FSVariable{T}, t::T) where {T}
    @boundscheck consistent_domains(problem, out) ||
        throw(ArgumentError("`problem` and `out` must have the same domain."))
    ζ_forcing = problem.ζ_forcing
    XS_tmp = problem.scratch.XS_tmp
    XZ_tmp = problem.scratch.XZ_tmp
    evaluate_ζ_forcing!(ζ_forcing, out, t, XS_tmp, XZ_tmp)
    return nothing
end

function advect_background_flow!(problem::Problem{T}, Ψ::FSVariable{T}) where {T}
    @boundscheck consistent_domains(problem, Ψ) ||
        throw(ArgumentError("Domains of `problem` and `Ψ` must match."))

    # unpack variables
    v, b = problem.state.v, problem.state.b
    (; f, Vx, Bx, Bz) = problem.background
    (; FS_tmp, FC_tmp, XS_tmp, XC_tmp, XZ_tmp, XZ_tmp2) = problem.scratch

    # form -U in physical space
    @inbounds ∂z!(FC_tmp, Ψ)
    Tᴴ!(XC_tmp, FC_tmp)
    Tᶜ!(XZ_tmp, XC_tmp)

    # update v and b with -U *(f + Vx) and -U * Bx
    @inbounds @. v += XZ_tmp * (f + Vx)
    @inbounds @. b += XZ_tmp * Bx

    # form W in physical space
    @inbounds ∂x!(FS_tmp, Ψ)
    Tᴴ!(XS_tmp, FS_tmp)
    Tˢ!(XZ_tmp, XS_tmp)

    # update v and b with - W * Bx / f and - W * Bz
    @inbounds @. v -= XZ_tmp * Bx / f
    @inbounds @. b -= XZ_tmp * Bz
    return nothing
end

"""$(TYPEDSIGNATURES)

Add forcing for v and b for one of the DIRKN substeps.
"""
function add_physical_forcing!(problem::Problem{T}, h::T, cᵢ::T, bᵢᵗ::T) where {T}
    # unpack variables
    v, b, clock = problem.state.v, problem.state.b, problem.state.clock
    v_forcing, b_forcing = problem.v_forcing, problem.b_forcing
    XZ_tmp = problem.scratch.XZ_tmp

    # add on the forcing terms for v
    evaluate_physical_forcing!(v_forcing, XZ_tmp, clock.t + cᵢ * h)
    @inbounds @. v += bᵢᵗ * h * XZ_tmp

    # and for b
    evaluate_physical_forcing!(b_forcing, XZ_tmp, clock.t + cᵢ * h)
    @inbounds @. b += bᵢᵗ * h * XZ_tmp

    return nothing
end

"""
    $(TYPEDSIGNATURES)

Advance the timestepper one timestep.
"""
function advance!(ts::Timestepper)
    # unpack structs
    (; problem, h, 𝓒, auxillary_variables, cgs, 𝓟) = ts
    (; a₁₁, a₂₁, a₂₂, b₁, b₂, b₁ᵗ, b₂ᵗ, c₁, c₂) = 𝓒
    (; ζ, ζₜ, v, b, clock) = problem.state
    (; ζⁿ⁺ᶜ¹, ζⁿ⁺ᶜ², tmp, rhs) = auxillary_variables

    # get forcing at n + c₁
    @inbounds get_ζ_forcing!(problem, tmp, clock.t + c₁ * h)
    # construct rhs of implicit equation for ζⁿ⁺ᶜ¹
    @inbounds @. rhs = ζ + c₁ * h * ζₜ + a₁₁ * h^2 * tmp
    @inbounds solve_implicit_equation!(cgs, ζⁿ⁺ᶜ¹, rhs, 𝓟)

    # start constructing the rhs of implicit equation at ζⁿ⁺ᶜ²
    # include ζⁿ, ζₜⁿ and Fⁿ⁺ᶜ¹ terms
    @inbounds @. rhs = ζ + c₂ * h * ζₜ + a₂₁ * h^2 * tmp

    # we no longer need ζⁿ and ζₜⁿ so we can start forming ζⁿ⁺¹ and ζₜⁿ⁺¹
    # add on the ζⁿ and Fⁿ⁺ᶜ¹ terms
    @inbounds @. ζ += h * ζₜ + b₁ * h^2 * tmp
    @inbounds @. ζₜ += b₁ᵗ * h * tmp

    # now we are done with Fⁿ⁺ᶜ¹ and can use tmp for 𝓛ζⁿ⁺ᶜ¹
    @inbounds 𝓛!(problem, tmp, ζⁿ⁺ᶜ¹)
    # add 𝓛ζⁿ⁺ᶜ¹ term to rhs and ζⁿ⁺¹, ζₜⁿ⁺¹
    @inbounds @. rhs -= a₂₁ * h^2 * tmp
    @inbounds @. ζ -= b₁ * h^2 * tmp
    @inbounds @. ζₜ -= b₁ᵗ * h * tmp

    # now we are done with 𝓛ζⁿ⁺ᶜ¹ and can use tmp for Fⁿ⁺ᶜ²
    @inbounds get_ζ_forcing!(problem, tmp, clock.t + c₂ * h)
    # add the Fⁿ⁺ᶜ² terms to rhs, ζⁿ⁺¹ and ζₜⁿ⁺¹
    @inbounds @. rhs += a₂₂ * h^2 * tmp
    @inbounds @. ζ += b₂ * h^2 * tmp
    @inbounds @. ζₜ += b₂ᵗ * h * tmp

    # we have fully formed the rhs of the implicit equation for ζⁿ⁺ᶜ² so we solve
    @inbounds solve_implicit_equation!(cgs, ζⁿ⁺ᶜ², rhs, 𝓟)

    # now compute 𝓛ζⁿ⁺ᶜ² and add those terms to ζⁿ⁺¹ and ζₜⁿ⁺¹
    @inbounds 𝓛!(problem, tmp, ζⁿ⁺ᶜ²)
    @inbounds @. ζ -= b₂ * h^2 * tmp
    @inbounds @. ζₜ -= b₂ᵗ * h * tmp

    # this concludes the computation of ζⁿ⁺¹ and ζₜⁿ⁺¹

    # now we can do the advection of the background flow to update v and b
    # for this we require Ψ = h * (b₁ᵗ * Ψⁿ⁺ᶜ¹ + b₂ᵗ * Ψⁿ⁺ᶜ²)
    # first we form Ζ = h * (b₁ᵗ * ζⁿ⁺ᶜ¹ + b₂ᵗ * ζⁿ⁺ᶜ²) and then we take the inverse Laplacian
    # we'll alias the ζⁿ⁺ᶜ¹ and ζⁿ⁺ᶜ² variables for these computations
    Z, Ψ = ζⁿ⁺ᶜ¹, ζⁿ⁺ᶜ²
    @inbounds @. Z *= h * b₁ᵗ
    @inbounds @. Z += h * b₂ᵗ * ζⁿ⁺ᶜ²
    @inbounds ∇⁻²!(Ψ, Z)

    @inbounds advect_background_flow!(problem, Ψ)

    # finally add on the forcing terms for v and b
    add_physical_forcing!(problem, h, c₁, b₁ᵗ)
    add_physical_forcing!(problem, h, c₂, b₂ᵗ)

    # finally update the clock
    return update_clock!(clock, h)
end

"""
    $(TYPEDSIGNATURES)

Advance the timestepper `n` timesteps.
"""
function advance!(ts::Timestepper, n::Integer)
    for _ in 1:n
        advance!(ts)
    end
end
