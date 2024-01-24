"""
    $(TYPEDEF)

The balanced background gradients. The background flow is in thermal wind balance and
therefore assumed to be steady.

# Fields
$(TYPEDFIELDS)
"""
struct BackgroundFlow
    """Coriolis frequency ``f``"""
    f::Float64
    raw"""Lateral shear ``\partial V/\partial x``"""
    Vx::Array{Float64,2}
    raw"""Horizontal buoyancy gradient ``\partial B/\partial x = f\partial V/\partial z``"""
    Bx::Array{Float64,2}
    raw"""Vertical buoyancy gradient ``\partial B/\partial z``"""
    Bz::Array{Float64,2}

    @doc """$(TYPEDSIGNATURES)

    Create a new [`BackgroundFlow`](@ref) given a [`Grid`](@ref) and Coriolis frequency
    ``f`` which defaults to 1. Background gradients are initialised to 0.

    # Examples
    ```jldoctest
    grid = Grid(256,128,2000,1)
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
    function BackgroundFlow(grid::Grid, f::Float64=1.0)
        NX, NZ = grid.NX, grid.NZ
        return new(
            f, zeros(Float64, NX, NZ), zeros(Float64, NX, NZ), zeros(Float64, NX, NZ)
        )
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
