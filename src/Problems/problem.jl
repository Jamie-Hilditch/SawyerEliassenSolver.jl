"""
    $(TYPEDEF)

Scratch space for temporary variables used in the Sawyer-Eliassen problem.
These variables have two primary uses in the code:

    1) intermediate terms in the computation of the Sawyer-Eliassen operator 𝓛
    2) intermediate terms in the advection of the background flow

Between timesteps these variables are available for other purposes e.g. setting the initial
conditions, computing output.

# Fields
$(TYPEDFIELDS)
"""
struct Scratch{T}
    FS_tmp::FSVariable{T}
    FC_tmp::FCVariable{T}
    XS_tmp::XSVariable{T}
    XC_tmp::XCVariable{T}
    XZ_tmp::XZVariable{T}
    XZ_tmp2::XZVariable{T}
end

function Scratch(domain::Domain{T}) where {T}
    return Scratch{T}(
        FSVariable(domain),
        FCVariable(domain),
        XSVariable(domain),
        XCVariable(domain),
        XZVariable(domain),
        XZVariable(domain),
    )
end

"""
    $(TYPEDEF)

A struct representing a Sawyer-Eliassen problem.

# Fields
$(TYPEDFIELDS)
"""
struct Problem{T,F,G,H}
    domain::Domain{T}
    background::BackgroundFlow{T}
    ζ_forcing::F
    v_forcing::Union{PhysicalForcing{T},NoForcing{T}}
    b_forcing::Union{PhysicalForcing{T},NoForcing{T}}
    state::State{T}
    scratch::Scratch{T}

    function Problem(
        domain::Domain{T},
        background::BackgroundFlow{T},
        ζ_forcing::AbstractForcingFunction{T},
        v_forcing::Union{PhysicalForcing{T},NoForcing{T}},
        b_forcing::Union{PhysicalForcing{T},NoForcing{T}},
    ) where {T}
        # first validate that the background and forcing match domain
        size(domain.grid) == size(background) ||
        throw(ArgumentError("`domain` and `background` must be the same size."))

        consistent_domains(domain, ζ_forcing, v_forcing, b_forcing) ||
        throw(ArgumentError("`domain` must match forcing domains"))

        return new{T,typeof(ζ_forcing),typeof(v_forcing),typeof(b_forcing)}(
            domain,
            background,
            ζ_forcing,
            v_forcing,
            b_forcing,
            State(domain),
            Scratch(domain),
        )
    end
end

function Problem(
    domain::Domain{T},
    background::BackgroundFlow{T};
    ζ_forcing=nothing,
    v_forcing=nothing,
    b_forcing=nothing,
) where {T}
    ζ_forcing = isnothing(ζ_forcing) ? NoForcing(domain) : ζ_forcing
    v_forcing = isnothing(v_forcing) ? NoForcing(domain) : v_forcing
    b_forcing = isnothing(b_forcing) ? NoForcing(domain) : b_forcing
    return Problem(domain, background, ζ_forcing, v_forcing, b_forcing)
end

"""$(TYPEDSIGNATURES)"""
function Base.show(io::IO, ::MIME"text/plain", problem::Problem)
    return print(
        io,
        "Problem:\n",
        "  ├─────────── domain: $(summary(problem.domain))\n",
        "  ├─────── background: $(summary(problem.background))\n",
        "  ├──────── ζ_forcing: $(summary(problem.ζ_forcing))\n",
        "  ├──────── v_forcing: $(summary(problem.v_forcing))\n",
        "  ├──────── b_forcing: $(summary(problem.b_forcing))\n",
        "  └──────────── state: $(summary(problem.state))\n",
    )
end

"""$(TYPEDSIGNATURES)"""
get_ζ_forcing(problem::Problem) = problem.ζ_forcing

"""$(TYPEDSIGNATURES)"""
get_scratch(problem::Problem) = problem.scratch

"""$(TYPEDSIGNATURES)"""
get_problem(problem::Problem) = problem

# get the domain of a problem and use this as the fallback for get_domain
Domains.get_domain(problem::Problem) = problem.domain
Domains.get_domain(a::Any) = get_domain(get_problem(a))

"""$(TYPEDSIGNATURES)"""
get_time(problem::Problem) = get_time(problem.state)
"""$(TYPEDSIGNATURES)"""
get_iteration(problem::Problem) = get_iteration(problem.state)

# forward some stuff onto background flow
"""$(TYPEDSIGNATURES)"""
get_f(problem::Problem) = get_f(problem.background)
"""$(TYPEDSIGNATURES)"""
get_Vx(problem::Problem) = get_Vx(problem.background)
"""$(TYPEDSIGNATURES)"""
get_Bx(problem::Problem) = get_Bx(problem.background)
"""$(TYPEDSIGNATURES)"""
get_Bz(problem::Problem) = get_Bz(problem.background)
