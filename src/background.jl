"""
    $(TYPEDEF)

The balanced background gradients. The background flow is in thermal wind balance and
therefore assumed to be steady.

# Fields
$(TYPEDFIELDS)
"""
struct BackgroundFlow{T}
    """Coriolis frequency ``f``"""
    f::T
    raw"""Lateral shear ``\partial V/\partial x``"""
    Vx::Matrix{T}
    raw"""Horizontal buoyancy gradient ``\partial B/\partial x = f\partial V/\partial z``"""
    Bx::Matrix{T}
    raw"""Vertical buoyancy gradient ``\partial B/\partial z``"""
    Bz::Matrix{T}

    @doc """$(TYPEDSIGNATURES)

    Create a new [`BackgroundFlow`](@ref) given a [`Grid`](@ref) and Coriolis frequency
    ``f`` which defaults to 1. Background gradients are initialised to 0.

    # Examples
    ```jldoctest
    grid = Grid(256,128,2000.0,1.0)
    bg = BackgroundFlow(grid)

    # output
    BackgroundFlow:
      ├─── f: 1
      ├── Vx: 256×128 Matrix{Float64}
      │       [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
      ├── Bx: 256×128 Matrix{Float64}
      │       [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
      └── Bz: 256×128 Matrix{Float64}
              [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
    ```
    """
    function BackgroundFlow(grid::Grid{T}, f::T=one(T)) where {T}
        NX, NZ = grid.NX, grid.NZ
        return new{T}(f, zeros(T, NX, NZ), zeros(T, NX, NZ), zeros(T, NX, NZ))
    end
end

function check_consistent_grid(grid::Grid, bg::BackgroundFlow)
    if size(bg.Vx) == size(bg.Bx) == size(bg.Bz) == size(grid)
        return nothing
    else
        error("grid and background flow do not have compatible sizes")
    end
end

"""$(TYPEDSIGNATURES)"""
function Base.show(io::IO, ::MIME"text/plain", bg::BackgroundFlow)
    return print(
        io,
        "BackgroundFlow:\n",
        "  ├─── f: $(sfmt(bg.f))\n",
        "  ├── Vx: $(summary(bg.Vx))\n",
        "  │       $(sprint(show,bg.Vx, context=:limit => true))\n",
        "  ├── Bx: $(summary(bg.Bx))\n",
        "  │       $(sprint(show,bg.Bx, context=:limit => true))\n",
        "  └── Bz: $(summary(bg.Bz))\n",
        "          $(sprint(show,bg.Bz, context=:limit => true))\n",
    )
end

"""$(TYPEDSIGNATURES)"""
function Base.summary(io::IO, bg::BackgroundFlow)
    return print(io, "BackgroundFlow: f = $(sfmt(bg.f)), Vx,Bx,Bz = $(summary(bg.Vx))")
end
