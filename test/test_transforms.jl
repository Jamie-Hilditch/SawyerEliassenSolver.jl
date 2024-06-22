NX, NZ = 4, 4
SX = NX ÷ 2 + 1

function test_inplace_transforms_from_variable(variable_type, domain, valid_outputs)
    @info "Testing inplace transforms of $(variable_type)"
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

@testset "Inplace transforms" begin
    for FT in float_types
        test_inplace_transforms(FT)
    end
end
