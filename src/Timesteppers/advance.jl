"""
    $(TYPEDSIGNATURES)

Advance the timestepper one timestep.
"""
function advance!(ts::Timestepper)
    # unpack structs
    (; problem, h, 𝓒, auxillary_variables, working_variables, 𝓛, 𝓛ᴵ!, cgs, 𝓟) = ts
    (; a₁₁, a₂₁, a₂₂, b₁, b₂, b₁ᵗ, b₂ᵗ, c₁, c₂) = 𝓒
    (; ζ, ζₜ, v, b, clock) = problem.state
    (; f, Vx, Bx, Bz) = problem.background
    ζ_forcing, v_forcing, b_forcing = problem.ζ_forcing,
    problem.v_forcing,
    problem.b_forcing
    (; ζⁿ⁺ᶜ¹, ζⁿ⁺ᶜ², rhs) = auxillary_variables
    (; FC, XC, XS, XZ) = working_variables

    # get forcing at n + c₁
    @inbounds evaluate_ζ_forcing!(ζ_forcing, tmp, clock.t + c₁ * h, XS, XZ)
    # construct rhs of implicit equation for ζⁿ⁺ᶜ¹
    @inbounds @. rhs = ζ + c₁ * h * ζₜ + a₁₁ * h^2 * tmp
    solve_implicit_equation(cgs, 𝓛ᴵ!, ζⁿ⁺ᶜ¹, rhs, 𝓟)

    # start constructing the rhs of implicit equation at ζⁿ⁺ᶜ²
    # include ζⁿ, ζₜⁿ and Fⁿ⁺ᶜ¹ terms
    @inbounds @. rhs = ζ + c₂ * h * ζₜ + a₂₁ * h^2 * tmp

    # we no longer need ζⁿ and ζₜⁿ so we can start forming ζⁿ⁺¹ and ζₜⁿ⁺¹
    # add on the ζⁿ and Fⁿ⁺ᶜ¹ terms
    @inbounds @. ζ += h * ζₜ + b₁ * h^2 * tmp
    @inbounds @. ζₜ += b₁ᵗ * h * tmp

    # now we are done with Fⁿ⁺ᶜ¹ and can use tmp for 𝓛ζⁿ⁺ᶜ¹
    @inbounds 𝓛(tmp, ζⁿ⁺ᶜ¹)
    # add 𝓛ζⁿ⁺ᶜ¹ term to rhs and ζⁿ⁺¹, ζₜⁿ⁺¹
    @inbounds @. rhs -= a₂₁ * h^2 * tmp
    @inbounds @. ζ -= b₁ * h^2 * tmp
    @inbounds @. ζₜ -= b₁ᵗ * h * tmp

    # now we are done with 𝓛ζⁿ⁺ᶜ¹ and can use tmp for Fⁿ⁺ᶜ²
    @inbounds evaluate_ζ_forcing!(ζ_forcing, tmp, clock.t + c₂ * h, XS, XZ)
    # add the Fⁿ⁺ᶜ² terms to rhs, ζⁿ⁺¹ and ζₜⁿ⁺¹
    @inbounds @. rhs += a₂₂ * h^2 * tmp
    @inbounds @. ζ += b₂ * h^2 * tmp
    @inbounds @. ζₜ += b₂ᵗ * h * tmp

    # we have fully formed the rhs of the implicit equation for ζⁿ⁺ᶜ² so we solve
    solve_implicit_equation(cgs, 𝓛ᴵ!, ζⁿ⁺ᶜ², rhs, 𝓟)

    # now compute 𝓛ζⁿ⁺ᶜ² and add those terms to ζⁿ⁺¹ and ζₜⁿ⁺¹
    @inbounds 𝓛(tmp, ζⁿ⁺ᶜ²)
    @inbounds @. ζ -= b₂ * h^2 * tmp
    @inbounds @. ζₜ -= b₂ᵗ * h * tmp

    # this concludes the computation of ζⁿ⁺¹ and ζₜⁿ⁺¹
    # now we can do the advection of the background flow to update v and b

    # first we need to form U = h * (b₁ᵗ * uⁿ⁺ᶜ¹ + b₂ᵗ * uⁿ⁺ᶜ²) and W = h * (b₁ᵗ * wⁿ⁺ᶜ¹ + b₂ᵗ * wⁿ⁺ᶜ²)
    # first form Ζ = h * (b₁ᵗ * ζⁿ⁺ᶜ¹ + b₂ᵗ * ζⁿ⁺ᶜ²)
    # we'll alias the ζⁿ⁺ᶜ¹ variable for this
    Z = ζⁿ⁺ᶜ¹
    @inbounds @. Z *= h * b₁ᵗ
    @inbounds @. Z += h * b₁ᵗ * ζⁿ⁺ᶜ²

    # now compute Ψ =  h * (b₁ᵗ * ζⁿ⁺ᶜ¹ + b₂ᵗ * ζⁿ⁺ᶜ²) aliasing the ζⁿ⁺ᶜ² variable
    Ψ = ζⁿ⁺ᶜ²
    @inbounds ∇⁻²!(Ψ, Z)

    # form -U in physical space
    @inbounds ∂z!(FC, ψ)
    Tᴴ!(XC, FC)
    Tᶜ!(XZ, XC)

    # update v and b with -U *(f + Vx) and -U * Bx
    @inbounds @. v += XZ * (f + Vx)
    @inbounds @. b += XZ * Bx

    # form W in physical space
    @inbounds ∂x!(tmp, ψ)
    Tᴴ!(XS, tmp)
    Tˢ!(XZ, XS)

    # update v and b with - W * Bx / f and - W * Bz
    @inbounds @. v -= XZ * Bx / f
    @inbounds @. b -= XZ * Bz

    # now we need to add on the forcing terms for v
    evaluate_physical_forcing(v_forcing, XZ, clock.t + c₁ * h)
    @inbounds @. v += b₁ᵗ * h * XZ
    evaluate_physical_forcing(v_forcing, XZ, clock.t + c₂ * h)
    @inbounds @. v += b₂ᵗ * h * XZ

    # and for b
    evaluate_physical_forcing(b_forcing, XZ, clock.t + c₁ * h)
    @inbounds @. b += b₁ᵗ * h * XZ
    evaluate_physical_forcing(b_forcing, XZ, clock.t + c₂ * h)
    @inbounds @. b += b₂ᵗ * h * XZ

    # finally update the clock
    return update_clock!(clock, h)
end
