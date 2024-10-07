function test_z_derivatives_method_types(FT)
    grid = Grid(FT, 8, 8, 1, 1)
    domain = Domain(grid)

    @test ∂z(XZVariable(domain)) isa XZVariable{FT}
    @test ∂z(XSVariable(domain)) isa XCVariable{FT}
    @test ∂z(XCVariable(domain)) isa XSVariable{FT}
    @test ∂z(FZVariable(domain)) isa FZVariable{FT}
    @test ∂z(FSVariable(domain)) isa FCVariable{FT}
    @test ∂z(FCVariable(domain)) isa FSVariable{FT}
end

function test_z_second_derivatives_method_types(FT)
    grid = Grid(FT, 8, 8, 1, 1)
    domain = Domain(grid)

    @test_throws MethodError ∂z²(XZVariable(domain))
    @test ∂z²(XSVariable(domain)) isa XSVariable{FT}
    @test ∂z²(XCVariable(domain)) isa XCVariable{FT}
    @test_throws MethodError ∂z²(FZVariable(domain))
    @test ∂z²(FSVariable(domain)) isa FSVariable{FT}
    @test ∂z²(FCVariable(domain)) isa FCVariable{FT}
end

function test_z_integral_method_types(FT)
    grid = Grid(FT, 8, 8, 1, 1)
    domain = Domain(grid)

    @test_throws MethodError ∫dz(XZVariable(domain))
    @test ∫dz(XSVariable(domain)) isa XCVariable{FT}
    @test ∫dz(XCVariable(domain)) isa XSVariable{FT}
    @test_throws MethodError ∫dz(FZVariable(domain))
    @test ∫dz(FSVariable(domain)) isa FCVariable{FT}
    @test ∫dz(FCVariable(domain)) isa FSVariable{FT}
end

function test_z_second_integral_method_types(FT)
    grid = Grid(FT, 8, 8, 1, 1)
    domain = Domain(grid)

    @test_throws MethodError ∫dz²(XZVariable(domain))
    @test ∫dz²(XSVariable(domain)) isa XSVariable{FT}
    @test ∫dz²(XCVariable(domain)) isa XCVariable{FT}
    @test_throws MethodError ∫dz²(FZVariable(domain))
    @test ∫dz²(FSVariable(domain)) isa FSVariable{FT}
    @test ∫dz²(FCVariable(domain)) isa FCVariable{FT}
end

@testset "z derivatives method types" for FT in FLOAT_TYPES
    @info "\tTesting z derivative methods with $(FT)"
    test_z_derivatives_method_types(FT)
    @info "\tTesting second z derivative methods with $(FT)"
    test_z_second_derivatives_method_types(FT)
    @info "\tTesting z integral methods with $(FT)"
    test_z_integral_method_types(FT)
    @info "\tTesting second z integral methods with $(FT)"
    test_z_second_integral_method_types(FT)
end
