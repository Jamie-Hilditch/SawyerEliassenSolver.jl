const real_variable_types = (XZVariable, XSVariable, XCVariable)
const complex_variable_types = (FZVariable, FSVariable, FCVariable)

function test_data_element_type(domain)
    FT = eltype(domain)
    NX, NZ = size(domain)
    SX, _ = size(domain.spectral)
    complex_types = map(T -> Complex{T}, FLOAT_TYPES)
    test_data_types = (FLOAT_TYPES..., complex_types..., Integer, Float16, Float64)

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

function test_domain_and_data_are_compatible_size(domain)
    FT = eltype(domain)
    NX, NZ = size(domain)
    SX, _ = size(domain.spectral)
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
    for FT in FLOAT_TYPES
        grid = Grid(FT, 4, 4, 1, 1)
        domain = Domain(grid)
        test_data_element_type(domain)
    end
end

@testset "Variable sizes" begin
    for FT in FLOAT_TYPES
        grid = Grid(FT, 4, 4, 1, 1)
        domain = Domain(grid)
        test_domain_and_data_are_compatible_size(domain)
    end
end
