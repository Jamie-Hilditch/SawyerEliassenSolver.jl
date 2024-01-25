"""Diagonally Implicit Runge-Kutta Nystrom (DIRKN) scheme"""

include("timestepper/sawyer_eliassen_rhs.jl")
include("timestepper/update_v.jl")
include("timestepper/update_b.jl")
include("timestepper/conjugate_gradient.jl")
include("timestepper/state.jl")

# define the constants
const c::Float64 = 17 / 14
const denom = 4 * (3 * c^2 - 3 * c + 1)

const b11 = (1 - c) / denom
const b12 = (3 * c - 1) * (2 * c - 1) / denom
const b21 = 1 / denom
const b22 = 3 * (2 * c - 1)^2 / denom
const c1 = c
const c2 = (3 * c - 2) / (3 * (2 * c - 1))
const a11 = c^2 / 2
const a22 = c^2 / 2
const a21 = -2 * (9 * c^4 - 9 * c^3 + 3 * c - 1) / (9 * (2 * c - 1)^2)

"""
    $(TYPEDEF)

Timesteps the state using a 2-stage 3rd order accurate Diagonally Implicit Runge-Kutta
Nyström (DIRKN) scheme.

# Fields
$(TYPEDFIELDS)
"""
struct Timestepper!
    state::State
    h::Float64
    f1::Array{ComplexF64,2}
    f2::Array{ComplexF64,2}
    r::Array{ComplexF64,2}
    L!::SawyerEliassenRHS!
    update_v!::UpdateV!
    update_b!::UpdateB!
    cgs!::ConjugateGradientSolver!
end

function Timestepper!(problem::Problem, h::Float64, cg_tol::Union{Float64,Nothing}=nothing)
    NX, NZ = problem.grid.NX, problem.grid.NZ
    L! = SawyerEliassenRHS!(problem)
    lhs! = LHS!(h, a11, L!)

    if cg_tol === nothing
        cgs! = ConjugateGradientSolver!(
            zeros(ComplexF64, NX ÷ 2 + 1, NZ), zeros(ComplexF64, NX ÷ 2 + 1, NZ), lhs!
        )
    else
        cgs! = ConjugateGradientSolver!(
            zeros(ComplexF64, NX ÷ 2 + 1, NZ),
            zeros(ComplexF64, NX ÷ 2 + 1, NZ),
            lhs!,
            cg_tol,
        )
    end

    return Timestepper!(
        State(problem),
        h,
        zeros(ComplexF64, NX ÷ 2 + 1, NZ),
        zeros(ComplexF64, NX ÷ 2 + 1, NZ),
        zeros(ComplexF64, NX ÷ 2 + 1, NZ),
        L!,
        UpdateV!(problem),
        UpdateB!(problem),
        cgs!,
    )
end

function Timestepper!(problem::Problem, h, cg_tol)
    return Timestepper!(problem, convert(Float64, h), cg_tol)
end

"""$(TYPEDSIGNATURES)

Advance the state one timestep using the DIRKN scheme

# Examples

```
[...]
timestepper, state = setup_simulation(problem,0.1)
timestepper()
```
"""
function (t::Timestepper!)()
    # alias state arrays
    ψ = t.state.ψ
    ψt = t.state.ψt
    v = t.state.v
    b = t.state.b

    # construct RHS of f1 equation, storing in r
    @. t.r = ψ + c1 * t.h * ψt
    t.L!(t.r)

    # solve for f1
    t.cgs!(t.f1, t.r)

    # construct RHS of f2 equation, storing in r
    @. t.r = ψ + c2 * t.h * ψt + t.h^2 * a21 * t.f1
    t.L!(t.r)

    # solve for f2
    t.cgs!(t.f2, t.r)

    # construct the input for update_v and update_b, storing in r
    # @. t.r = t.h*b21*(ψ + t.h*c1*ψt + t.h^2*a11*t.f1) + t.h*b22*(ψ + t.h*c2*ψt + t.h^2*a21*t.f1 + t.h^2*a22*t.f2)
    @. t.r = t.h * (b21 + b22) * ψ
    @. t.r += t.h^2 * (b21 * c1 + b22 * c2) * ψt
    @. t.r += t.h^3 * (b21 * a11 + b22 * a21) * t.f1
    @. t.r += t.h^3 * b22 * a22 * t.f2

    # update the state
    @. ψ += t.h * ψt + t.h^2 * (b11 * t.f1 + b12 * t.f2)
    @. ψt += t.h * (b21 * t.f1 + b22 * t.f2)
    t.update_v!(v, t.r)
    t.update_b!(b, t.r)
    return update_clock!(t.state.clock, t.h)
end

"""$(TYPEDSIGNATURES)

Advance the state `n` timesteps

# Examples

```
[...]
timestepper, state = setup_simulation(problem,0.1)
timestepper(10)
```
"""
function (t::Timestepper!)(n::Int)
    for _ in 1:n
        t()
    end
end
