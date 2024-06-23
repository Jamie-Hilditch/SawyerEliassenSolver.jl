NX, NZ = 4, 4
SX = NX ÷ 2 + 1
real_variable_types = [XZVariable, XSVariable, XCVariable]
complex_variable_types = [FZVariable, FSVariable, FCVariable]

function test_data_element_type(FT)
    LX, LZ = one(FT), one(FT)
    grid = Grid(NZ, NX, LX, LZ)
    domain = Domain(grid)
    complex_types = map(T -> Complex{T}, float_types)
    test_data_types = (float_types..., complex_types..., Integer, Float16, Float64)

    for data_FT in test_data_types
        for variable_type in real_variable_types
            data = zeros(data_FT, NX, NZ)
            if data_FT == FT
                @test variable_type(domain, data) isa variable_type
            else
                @test_throws Union{MethodError,ArgumentError} variable_type(domain, data)
            end
        end

        for variable_type in complex_variable_types
            data = zeros(data_FT, SX, NZ)
            if data_FT == Complex{FT}
                @test variable_type(domain, data) isa variable_type
            else
                @test_throws Union{MethodError,ArgumentError} variable_type(domain, data)
            end
        end
    end
end

function test_domain_and_data_are_compatible_size(FT)
    LX, LZ = one(FT), one(FT)
    grid = Grid(NZ, NX, LX, LZ)
    domain = Domain(grid)
    correct_real_array = zeros(FT, NX, NZ)
    correct_complex_array = zeros(Complex{FT}, SX, NZ)
    incorrect_real_array = zeros(FT, SX, NZ)
    incorrect_complex_array = zeros(Complex{FT}, NX, NZ)

    # test constructers that create the data array
    for variable_type in real_variable_types
        var = variable_type(domain)
        @test size(var.data) == size(domain.grid)
    end

    for variable_type in complex_variable_types
        var = variable_type(domain)
        @test size(var.data) == size(domain.spectral)
    end

    # test when we pass in an array of the correct size
    for variable_type in real_variable_types
        var = variable_type(domain, correct_real_array)
        @test size(var.data) == size(domain.grid)
    end

    for variable_type in complex_variable_types
        var = variable_type(domain, correct_complex_array)
        @test size(var.data) == size(domain.spectral)
    end

    # test we error when we pass in an array of the incorrect size
    for variable_type in real_variable_types
        @test_throws ArgumentError variable_type(domain, incorrect_real_array)
    end

    for variable_type in complex_variable_types
        @test_throws ArgumentError variable_type(domain, incorrect_complex_array)
    end
end

@testset "Variable data types" begin
    for FT in float_types
        test_data_element_type(FT)
    end
end

@testset "Variable sizes" begin
    for FT in float_types
        test_domain_and_data_are_compatible_size(FT)
    end
end
