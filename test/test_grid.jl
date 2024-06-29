
function test_grid_eltype()
    NX, NZ = 4, 4
    LX, LZ = 1, 1

    @test_throws TypeError Grid(Float16, NX, NZ, LX, LZ)
    @test eltype(Grid(Float32, NX, NZ, LX, LZ)) == Float32
    @test eltype(Grid(Float64, NX, NZ, LX, LZ)) == Float64
    @test_throws TypeError Grid(Int, NX, NZ, LX, LZ)
end

function test_grid_methods()
    NX, NZ = 4, 4
    LX, LZ = 1, 1
    x_bounds = (-1, 1)
    z_bounds = (-1, 1)

    @test Grid(NX, NZ, LX, LZ) isa Grid
    @test Grid(NX, NZ, x_bounds, LZ) isa Grid
    @test Grid(NX, NZ, LX, z_bounds) isa Grid
    @test Grid(NX, NZ, x_bounds, z_bounds) isa Grid

    @test Grid(Float32, NX, NZ, LX, LZ) isa Grid
    @test Grid(Float32, NX, NZ, x_bounds, LZ) isa Grid
    @test Grid(Float32, NX, NZ, LX, z_bounds) isa Grid
    @test Grid(Float32, NX, NZ, x_bounds, z_bounds) isa Grid
end

function test_gridpoints(FT)
    NX, NZ = 4, 4
    LX, LZ = 1, 1
    grid = Grid(FT, NX, NZ, LX, LZ)
    @test xgridpoints(grid) == FT[0, 0.25, 0.5, 0.75]
    @test zgridpoints(grid) == FT[-0.875, -0.625, -0.375, -0.125]
end

@testset "Grid eltype" begin
    test_grid_eltype()
end
@testset "Grid methods" begin
    test_grid_methods()
end

@testset "Gridpoints" for FT in float_types
    test_gridpoints(FT)
end
