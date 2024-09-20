function test_inplace_transforms_from_variable(variable_type, domain, valid_outputs)
    @info "\tTesting inplace transforms of $(variable_type)"
    variables = (
        xz=XZVariable(domain),
        xs=XSVariable(domain),
        xc=XCVariable(domain),
        fz=FZVariable(domain),
        fs=FSVariable(domain),
        fc=FCVariable(domain),
    )
    for (output_space, output_variable) in pairs(variables)
        variable = variable_type(domain)
        if output_space in valid_outputs
            @test transform!(output_variable, variable) isa Nothing
        else
            @test_throws MethodError transform!(output_variable, variable)
        end
    end
    return nothing
end

function test_inplace_transforms(FT)
    NX, NZ = 4, 4
    LX, LZ = one(FT), one(FT)
    grid = Grid(NZ, NX, LX, LZ)
    domain = Domain(grid)

    test_inplace_transforms_from_variable(XZVariable, domain, [:xs, :xc, :fz])
    test_inplace_transforms_from_variable(XSVariable, domain, [:xz, :fs])
    test_inplace_transforms_from_variable(XCVariable, domain, [:xz, :fc])
    test_inplace_transforms_from_variable(FZVariable, domain, [:xz])
    test_inplace_transforms_from_variable(FSVariable, domain, [:xs])
    test_inplace_transforms_from_variable(FCVariable, domain, [:xc])
    return nothing
end

function test_horizontal_transforms(FT)
    NX, NZ = 4, 4
    LX, LZ = one(FT), one(FT)
    grid = Grid(NZ, NX, LX, LZ)
    domain = Domain(grid)

    @test horizontal_transform(XZVariable(domain)) isa FZVariable
    @test horizontal_transform(XSVariable(domain)) isa FSVariable
    @test horizontal_transform(XCVariable(domain)) isa FCVariable
    @test horizontal_transform(FZVariable(domain)) isa XZVariable
    @test horizontal_transform(FSVariable(domain)) isa XSVariable
    @test horizontal_transform(FCVariable(domain)) isa XCVariable
end

function test_sine_transforms(FT)
    NX, NZ = 4, 4
    LX, LZ = one(FT), one(FT)
    grid = Grid(NZ, NX, LX, LZ)
    domain = Domain(grid)

    @test sine_transform(XZVariable(domain)) isa XSVariable
    @test sine_transform(XSVariable(domain)) isa XZVariable
    @test_throws MethodError sine_transform(XCVariable(domain))
    @test_throws MethodError sine_transform(FZVariable(domain))
    @test_throws MethodError sine_transform(FSVariable(domain))
    @test_throws MethodError sine_transform(FCVariable(domain))
end

function test_cosine_transforms(FT)
    NX, NZ = 4, 4
    LX, LZ = one(FT), one(FT)
    grid = Grid(NZ, NX, LX, LZ)
    domain = Domain(grid)

    @test cosine_transform(XZVariable(domain)) isa XCVariable
    @test_throws MethodError cosine_transform(XSVariable(domain))
    @test cosine_transform(XCVariable(domain)) isa XZVariable
    @test_throws MethodError cosine_transform(FZVariable(domain))
    @test_throws MethodError cosine_transform(FSVariable(domain))
    @test_throws MethodError cosine_transform(FCVariable(domain))
end

@testset "Inplace transforms" begin
    for FT in FLOAT_TYPES
        @info "Testing inplace transforms with $(FT)"
        test_inplace_transforms(FT)
    end
end

@testset "Horizontal transforms" begin
    for FT in FLOAT_TYPES
        @info "Testing horizontal transforms with $(FT)"
        test_horizontal_transforms(FT)
    end
end

@testset "Sine transforms" begin
    for FT in FLOAT_TYPES
        @info "Testing sine transforms with $(FT)"
        test_sine_transforms(FT)
    end
end

@testset "Cosine transforms" for FT in FLOAT_TYPES
    @info "Testing cosine transforms with $(FT)"
    test_cosine_transforms(FT)
end
