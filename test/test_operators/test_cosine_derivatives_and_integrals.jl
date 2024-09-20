@inline trial_function_for_cosine_operators(z) = cos(π * z) - cos(5 * π * z)
@inline derivative_of_trial_function_for_cosine_operators(z) = -π * sin(π * z) + 5 * π * sin(5 * π * z)
@inline second_derivative_of_trial_function_for_cosine_operators(z) =
    -π^2 * cos(π * z) + 25 * π^2 * cos(5 * π * z)
@inline integral_of_trial_function_for_cosine_operators(z) = 1 / π * sin(π * z) - 1 / (5 * π) * sin(5 * π * z)
@inline second_integral_of_trial_function_for_cosine_operators(z) =
    -1 / π^2 * cos(π * z) + 1 / (25 * π^2) * cos(5 * π * z)

function setup_input_XC_for_cosine_operators(domain)
    z_array = zgridpoints(domain)'
    input_XZ = XZVariable(domain)
    @. input_XZ = trial_function_for_cosine_operators(z_array)
    return cosine_transform(input_XZ)
end

function test_cosine_derivative(FT)
    grid = Grid(FT, 4, 16, 1, 1)
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    @. expected_output = derivative_of_trial_function_for_cosine_operators(z_array)

    # out-of-place
    input_XC = setup_input_XC_for_cosine_operators(domain)
    output_XS = XSVariable(domain)
    ∂z!(output_XS, input_XC)
    @test expected_output ≈ sine_transform(output_XS)

    # new_output
    input_XC = setup_input_XC_for_cosine_operators(domain)
    output_XS = ∂z(input_XC)
    @test expected_output ≈ sine_transform(output_XS)
end

function test_cosine_second_derivative(FT)
    grid = Grid(FT, 4, 16, 1, 1)
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    @. expected_output = second_derivative_of_trial_function_for_cosine_operators(z_array)

    # out-of-place
    input_XC = setup_input_XC_for_cosine_operators(domain)
    output_XC = XCVariable(domain)
    ∂z²!(output_XC, input_XC)
    @test expected_output ≈ cosine_transform(output_XC)

    # new_output
    input_XC = setup_input_XC_for_cosine_operators(domain)
    output_XC = ∂z²(input_XC)
    @test expected_output ≈ cosine_transform(output_XC)
end

function test_cosine_integral(FT)
    grid = Grid(FT, 4, 16, 1, 1)
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    @. expected_output = integral_of_trial_function_for_cosine_operators(z_array)

    # out-of-place
    input_XC = setup_input_XC_for_cosine_operators(domain)
    output_XS = XSVariable(domain)
    ∫dz!(output_XS, input_XC)
    @test expected_output ≈ sine_transform(output_XS)

    # new_output
    input_XC = setup_input_XC_for_cosine_operators(domain)
    output_XS = ∫dz(input_XC)
    @test expected_output ≈ sine_transform(output_XS)
end

function test_cosine_second_integral(FT)
    grid = Grid(FT, 4, 16, 1, 1)
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    @. expected_output = second_integral_of_trial_function_for_cosine_operators(z_array)

    # out-of-place
    input_XC = setup_input_XC_for_cosine_operators(domain)
    output_XC = XCVariable(domain)
    ∫dz²!(output_XC, input_XC)
    @test expected_output ≈ cosine_transform(output_XC)

    # new_output
    input_XC = setup_input_XC_for_cosine_operators(domain)
    output_XC = ∫dz²(input_XC)
    @test expected_output ≈ cosine_transform(output_XC)
end

@testset "z - cosine derivatives" for FT in FLOAT_TYPES
    @info "\tTesting z derivatives in cosine space for $(FT)"
    test_cosine_derivative(FT)
end

@testset "z - cosine 2nd derivatives" for FT in FLOAT_TYPES
    @info "\tTesting z 2nd derivatives in sine space for $(FT)"
    test_cosine_second_derivative(FT)
end

@testset "z - cosine integral" for FT in FLOAT_TYPES
    @info "\tTesting z integrals in sine space for $(FT)"
    test_cosine_integral(FT)
end

@testset "z - cosine 2nd integral" for FT in FLOAT_TYPES
    @info "\tTesting z 2nd integrals in sine space for $(FT)"
    test_cosine_second_integral(FT)
end
