mutable struct Clock
    t :: Float64 
    iteration :: Int
end 

"""Default initialiser for Clock"""
function Clock()
    Clock(0,0)
end

struct State
    ψ :: Matrix{ComplexF64}
    ψt :: Matrix{ComplexF64}
    v :: Matrix{Float64}
    b :: Matrix{Float64}
    clock :: Clock

    function State(problem :: Problem,clock :: Clock)
        new(
            zeros(ComplexF64,size(problem.spectral_domain)),
            zeros(ComplexF64,size(problem.spectral_domain)),
            zeros(Float64,size(problem.grid)),
            zeros(Float64,size(problem.grid)),
            clock
        )
    end
end

function State(problem :: Problem)
    State(problem, Clock())
end

function update_clock!(clock :: Clock, h :: Float64)
    clock.t += h 
    clock.iteration += 1 
    return nothing
end

Base.size(state :: State) = size(state.v)

Base.show(io::IO, ::MIME"text/plain", state::State) = print(io,
    "State:\n",
    "  ├────── ψ: $(summary(state.ψ)))\n",
    "  │          $(sprint(show,state.ψ, context=:limit => true))\n",
    "  ├──────ψt: $(summary(state.ψt))\n",
    "  │          $(sprint(show,state.ψt, context=:limit => true))\n",
    "  ├────── v: $(summary(state.v))\n",
    "  │          $(sprint(show,state.v, context=:limit => true))\n",
    "  ├────── b: $(summary(state.b))\n",
    "  │          $(sprint(show,state.b, context=:limit => true))\n",
    "  └── clock: $(summary(state.clock))\n"
)

Base.show(io::IO, ::MIME"text/plain", clock::Clock) = print(io,
    "Clock:n",
    "  ├────────── t: $(sfmt(clock.t))\n",
    "  └── iteration: $(clock.iteration)\n"
)

Base.summary(io::IO, state::State) = print(io,"State:  ψ,ψt = $(summary(state.ψ)), v,b = $(summary(state.v)), clock = $(summary(state.clock))")

Base.summary(io::IO, clock::Clock) = print(io,"Clock(t = $(sfmt(clock.t)), iteration = $(clock.iteration))")