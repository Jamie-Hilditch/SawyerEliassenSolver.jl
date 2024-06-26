
struct Problem{T}
    domain::Domain{T}
    background::BackgroundFlow{T}
    ζ_forcing::Union{AbstractForcingFunction{T},Nothing}
    v_forcing::Union{PhysicalForcing{T},Nothing}
    b_forcing::Union{PhysicalForcing{T},Nothing}
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
        isnothing(ζ_forcing) ||
            ζ_forcing.domain == domain ||
            throw(ArgumentError("`ζ_forcing` has incompatible domain."))
        isnothing(v_forcing) ||
            v_forcing.domain == domain ||
            throw(ArgumentError("`v_forcing` has incompatible domain."))
        isnothing(b_forcing) ||
            b_forcing.domain == domain ||
            throw(ArgumentError("`b_forcing` has incompatible domain."))
        # now create the state
        return new{T}(domain, background, ζ_forcing, v_forcing, b_forcing, State(domain))
    end
end

"""$(TYPEDSIGNATURES)"""
function Base.show(io::IO, ::MIME"text/plain", problem::Problem)
    return print(
        io,
        "Problem:\n",
        "  ├─────────── domain: $(summary(problem.grid))\n",
        "  ├── background_flow: $(summary(problem.background_flow))\n",
        "  ├──────── ζ_forcing: $(summary(problem.ζ_forcing))\n",
        "  ├──────── v_forcing: $(summary(problem.v_forcing))\n",
        "  ├──────── b_forcing: $(summary(problem.b_forcing))\n",
        "  └──────────── state: $(summary(problem.state))\n",
    )
end
