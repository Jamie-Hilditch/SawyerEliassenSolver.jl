struct BackgroundFlow
    f :: Float64
    Vx :: Array{Float64,2}
    Bx :: Array{Float64,2}
    Bz :: Array{Float64,2}

    function BackgroundFlow(grid :: Grid,f = 1.)
        NX, NZ = grid.NX, grid.NZ
        new(
            f,
            zeros(Float64,NX,NZ),
            zeros(Float64,NX,NZ),
            zeros(Float64,NX,NZ)
        )
    end
end

function check_consistent_grid(grid:: Grid, bg :: BackgroundFlow)
    if size(bg.Vx) == size(bg.Bx) == size(bg.Bz) == size(grid)
        return nothing
    else
        error("grid and background flow do not have compatible sizes")

    end
end

Base.show(io::IO, ::MIME"text/plain", bg::BackgroundFlow) = print(io,
    "BackgroundFlow:\n",
    "  ├─── f: $(sfmt(bg.f))\n",
    "  ├── Vx: $(summary(bg.Vx))\n",
    "  │       $(sprint(show,bg.Vx, context=:limit => true))\n",
    "  ├── Bx: $(summary(bg.Bx))\n",
    "  │       $(sprint(show,bg.Bx, context=:limit => true))\n",
    "  └── Bz: $(summary(bg.Bz))\n",
    "          $(sprint(show,bg.Bz, context=:limit => true))\n"
)


Base.summary(io::IO, bg::BackgroundFlow) = print(io,"BackgroundFlow: f = $(sfmt(bg.f)), Vx,Bx,Bz = $(summary(bg.Vx))")