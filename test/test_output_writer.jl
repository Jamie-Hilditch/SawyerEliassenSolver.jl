NX, NZ = 4, 4
LX, LZ = 1, 1

function create_problem(::FT)
    grid = Grid(FT, NZ, NX, LX, LZ)
    domain = Domain(grid)
    bg = BackgroundFlow(grid)
    return Problem(domain, bg)
end
