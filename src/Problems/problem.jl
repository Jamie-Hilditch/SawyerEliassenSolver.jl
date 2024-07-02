"""
    $(TYPEDEF)

A struct representing a Sawyer-Eliassen problem.

# Fields
$(TYPEDFIELDS)
"""
struct Problem{T}
    domain::Domain{T}
    background::BackgroundFlow{T}
    ζ_forcing::AbstractForcingFunction{T}
    v_forcing::Union{PhysicalForcing{T},NoForcing{T}}
    b_forcing::Union{PhysicalForcing{T},NoForcing{T}}
    state::State{T}
    function Problem(
        domain::Domain{T},
        background::BackgroundFlow{T};
        ζ_forcing=nothing,
        v_forcing=nothing,
        b_forcing=nothing,
    ) where {T}
        # first validate that the background and forcing match domain
        size(domain.grid) == size(background) ||
            throw(ArgumentError("`domain` and `background` must be the same size."))

        ζ_forcing = isnothing(ζ_forcing) ? NoForcing(domain) : ζ_forcing
        v_forcing = isnothing(v_forcing) ? NoForcing(domain) : v_forcing
        b_forcing = isnothing(b_forcing) ? NoForcing(domain) : b_forcing

        consistent_domains(domain, ζ_forcing, v_forcing, b_forcing) ||
            throw(ArgumentError("`domain` must match forcing domains"))
        # now create the state
        return new{T}(domain, background, ζ_forcing, v_forcing, b_forcing, State(domain))
    end
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

Domains.get_domain(problem::Problem) = problem.domain
