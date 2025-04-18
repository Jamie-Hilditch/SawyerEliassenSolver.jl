"""$(TYPEDSIGNATURES)
Set ζ from a variable in spectral space.
"""
function set_ζ!(state::State, ζ::FSVariable)
    @boundscheck consistent_domains(state, ζ)
    @inbounds state.ζ .= ζ
    return nothing
end

"""$(TYPEDSIGNATURES)
Set ζ from a variable in XS space.
"""
function set_ζ!(state::State, ζ::XSVariable)
    @boundscheck consistent_domains(state, ζ)
    @inbounds horizontal_transform!(state.ζ, ζ)
    return nothing
end

"""$(TYPEDSIGNATURES)
Set ζ from a variable in physical space.
"""
function set_ζ!(state::State, ζ::XZVariable)
    @boundscheck consistent_domains(state, ζ)
    ζˣˢ = sine_transform(ζ)
    @inbounds set_ζ!(state, ζˣˢ)
    return nothing
end

"""$(TYPEDSIGNATURES)
Set ζ from a function. `func` should specify `ζ = func(x,z)` pointwise.
"""
function set_ζ!(state::State, func::Function)
    domain = get_domain(state)
    ζ = XZVariable(domain)
    x, z = gridpoints(domain)
    @inbounds @. ζ = func(x, z)
    @inbounds set_ζ!(state, ζ)
    return nothing
end

function _add_u_to_ζ!(state::State, u::FCVariable)
    @boundscheck consistent_domains(state, u) ||
        throw(ArgumentError("`state` and `u` must have the same domain."))
    uz = ∂z(u)
    state.ζ .-= uz
    return nothing
end

function _add_u_to_ζ!(state::State, u::XCVariable)
    @boundscheck consistent_domains(state, u) ||
        throw(ArgumentError("`state` and `u` must have the same domain."))
    uᶠᶜ = horizontal_transform(u)
    @inbounds _add_u_to_ζ!(state, uᶠᶜ)
    return nothing
end

function _add_u_to_ζ!(state::State, u::XZVariable)
    @boundscheck consistent_domains(state, u) ||
        throw(ArgumentError("`state` and `u` must have the same domain."))
    uˣᶜ = cosine_transform(u)
    @inbounds _add_u_to_ζ!(state, uˣᶜ)
    return nothing
end

function _add_u_to_ζ!(state::State, func::Function)
    domain = get_domain(state)
    u = XZVariable(domain)
    x, z = gridpoints(domain)
    @inbounds @. u = func(x, z)
    @inbounds _add_u_to_ζ!(state, u)
    return nothing
end

function _add_w_to_ζ!(state::State, w::FSVariable)
    @boundscheck consistent_domains(state, w) ||
        throw(ArgumentError("`state` and `w` must have the same domain."))
    wx = ∂x(w)
    state.ζ .+= wx
    return nothing
end

function _add_w_to_ζ!(state::State, w::XSVariable)
    @boundscheck consistent_domains(state, w) ||
        throw(ArgumentError("`state` and `u` must have the same domain."))
    wᶠˢ = horizontal_transform(w)
    @inbounds _add_w_to_ζ!(state, wᶠˢ)
    return nothing
end

function _add_w_to_ζ!(state::State, w::XZVariable)
    @boundscheck consistent_domains(state, w) ||
        throw(ArgumentError("`state` and `u` must have the same domain."))
    wˣˢ = sine_transform(w)
    @inbounds _add_w_to_ζ!(state, wˣˢ)
    return nothing
end

function _add_w_to_ζ!(state::State, func::Function)
    domain = get_domain(state)
    w = XZVariable(domain)
    x, z = gridpoints(domain)
    @inbounds @. w = func(x, z)
    @inbounds _add_w_to_ζ!(state, w)
    return nothing
end

"""$(TYPEDSIGNATURES)
Set ζ from u and w. `u` and `w` may be a `XZVariable`, `XCVariable`/ `XSVariable`, `FCVariable` / `FSVariable`
or a function specifying `u(x,z)` / `w(x,z)` pointwise. If not specified `u` and `w` default to 0.
"""
function set_ζ!(state::State; u=nothing, w=nothing)
    state.ζ .= 0
    isnothing(u) || _add_u_to_ζ!(state, u)
    isnothing(w) || _add_w_to_ζ!(state, w)
    return nothing
end

"""$(TYPEDSIGNATURES)"""
@propagate_inbounds set_ζ!(problem::Problem, ζ) = set_ζ!(problem.state, ζ)
"""$(TYPEDSIGNATURES)"""
@propagate_inbounds set_ζ!(problem::Problem; u=nothing, w=nothing) =
    set_ζ!(problem.state; u=u, w=w)

"""$(TYPEDSIGNATURES)
Set `v` from a variable.
"""
function set_v!(state::State, v::XZVariable)
    @boundscheck consistent_domains(state, v) ||
        throw(ArgumentError("`state` and `v` must have the same domain."))
    state.v .= v
    return nothing
end

"""$(TYPEDSIGNATURES)
Set `v` from a function `v(x,z)`.
"""
function set_v!(state::State, func::Function)
    domain = get_domain(state)
    x, z = gridpoints(domain)
    @. state.v = func(x, z)
    return nothing
end

"""$(TYPEDSIGNATURES)"""
set_v!(problem::Problem, v) = set_v!(problem.state, v)

"""$(TYPEDSIGNATURES)
Set `b` from a variable.
"""
function set_b!(state::State, b::XZVariable)
    @boundscheck consistent_domains(state, b) ||
        throw(ArgumentError("`state` and `b` must have the same domain."))
    state.b .= b
    return nothing
end

"""$(TYPEDSIGNATURES)
Set `b` from a function `b(x,z)`.
"""
function set_b!(state::State, func::Function)
    domain = get_domain(state)
    x, z = gridpoints(domain)
    @. state.b = func(x, z)
    return nothing
end

"""$(TYPEDSIGNATURES)"""
set_b!(problem::Problem, b) = set_b!(problem.state, b)

"""$(TYPEDSIGNATURES)
Compute and set ζₜ by projecting `b` onto sine space, `v` onto cosine space and using `ζₜ = bx - fvz`
"""
function compute_ζₜ!(problem::Problem)
    b = problem.state.b
    v = problem.state.v
    ζₜ = problem.state.ζₜ
    f = problem.background.f

    bˣˢ = sine_transform(b)
    Tᴴ!(ζₜ, bˣˢ)
    ∂x!(ζₜ)

    vˣᶜ = cosine_transform(v)
    vᶠᶜ = horizontal_transform(vˣᶜ)
    vzᶠˢ = ∂z(vᶠᶜ)
    @inbounds @. ζₜ -= f * vzᶠˢ
    return nothing
end

"""$(TYPEDSIGNATURES)
Set `v` and / or `b` and then compute and set `ζₜ`
"""
function set_vb!(problem::Problem; v=nothing, b=nothing)
    isnothing(v) || set_v!(problem, v)
    isnothing(b) || set_b!(problem, b)
    compute_ζₜ!(problem)
    return nothing
end
