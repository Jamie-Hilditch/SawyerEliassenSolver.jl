"""
    $(TYPEDEF)

Current simulation time and iteration.

$(TYPEDFIELDS)
"""
mutable struct Clock{T}
    t::T
    iteration::Int
end

"""Default initialiser for Clock"""
Clock{T}() where {T} = Clock{T}(0, 0)

"""
    $(TYPEDEF)

State variables of the problem.

$(TYPEDFIELDS)
"""
struct State{T}
    ζ::FSVariable{T}
    ζₜ::FSVariable{T}
    v::XZVariable{T}
    b::XZVariable{T}
    clock::Clock{T}

    function State(domain::Domain{T}, clock::Clock{T}) where {T}
        return new{T}(
            FSVariable(domain),
            FSVariable(domain),
            XZVariable(domain),
            XZVariable(domain),
            clock,
        )
    end
end

State(domain::Domain{T}) where {T} = State(domain, Clock{T}())

@inline function update_clock!(clock::Clock{T}, h::T) where {T}
    clock.t += h
    clock.iteration += 1
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", state::State)
    return print(
        io,
        "State:\n",
        "  ├────── ζ: $(summary(state.ζ)))\n",
        "  │          $(sprint(show,state.ζ, context=:limit => true))\n",
        "  ├───── ζₜ: $(summary(state.ζt))\n",
        "  │          $(sprint(show,state.ζt, context=:limit => true))\n",
        "  ├────── v: $(summary(state.v))\n",
        "  │          $(sprint(show,state.v, context=:limit => true))\n",
        "  ├────── b: $(summary(state.b))\n",
        "  │          $(sprint(show,state.b, context=:limit => true))\n",
        "  └── clock: $(summary(state.clock))\n",
    )
end

function Base.show(io::IO, ::MIME"text/plain", clock::Clock)
    return print(
        io,
        "Clock:n",
        "  ├────────── t: $(sfmt(clock.t))\n",
        "  └── iteration: $(clock.iteration)\n",
    )
end

function Base.summary(io::IO, state::State)
    return print(
        io,
        "ζ,ζₜ = $(summary(state.ζ)), v,b = $(summary(state.v)), clock = $(summary(state.clock))",
    )
end

function Base.summary(io::IO, clock::Clock)
    return print(io, "Clock(t = $(sfmt(clock.t)), iteration = $(clock.iteration))")
end

Domains.get_domain(state::State) = get_domain(state.ζ)
