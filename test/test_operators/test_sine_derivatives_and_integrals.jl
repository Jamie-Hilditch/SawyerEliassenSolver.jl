@inline trial_function_for_sine_operators(z) = sin(π * z) - sin(5 * π * z)
@inline derivative_of_trial_function_for_sine_operators(z) = π * cos(π * z) - 5 * π * cos(5 * π * z)
@inline second_derivative_of_trial_function_for_sine_operators(z) =
    -π^2 * sin(π * z) + 25 * π^2 * sin(5 * π * z)
@inline integral_of_trial_function_for_sine_operators(z) = -1 / π * cos(π * z) + 1 / (5 * π) * cos(5 * π * z)
@inline second_integral_of_trial_function_for_sine_operators(z) =
    -1 / π^2 * sin(π * z) + 1 / (25 * π^2) * sin(5 * π * z)

function setup_input_XS_for_sine_operators(domain)
    z_array = zgridpoints(domain)'
    input_XZ = XZVariable(domain)
    @. input_XZ = trial_function_for_sine_operators(z_array)
    return sine_transform(input_XZ)
end

function test_sine_derivative(FT)
    grid = Grid(FT, 4, 16, 1, 1)
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    @. expected_output = derivative_of_trial_function_for_sine_operators(z_array)

    # out-of-place
    input_XS = setup_input_XS_for_sine_operators(domain)
    output_XC = XCVariable(domain)
    ∂z!(output_XC, input_XS)
    @test expected_output ≈ cosine_transform(output_XC)

    # new_output
    input_XS = setup_input_XS_for_sine_operators(domain)
    output_XC = ∂z(input_XS)
    @test expected_output ≈ cosine_transform(output_XC)
end

function test_sine_second_derivative(FT)
    grid = Grid(FT, 4, 16, 1, 1)
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    @. expected_output = second_derivative_of_trial_function_for_sine_operators(z_array)

    # out-of-place
    input_XS = setup_input_XS_for_sine_operators(domain)
    output_XS = XSVariable(domain)
    ∂z²!(output_XS, input_XS)
    @test expected_output ≈ sine_transform(output_XS)

    # new_output
    input_XS = setup_input_XS_for_sine_operators(domain)
    output_XS = ∂z²(input_XS)
    @test expected_output ≈ sine_transform(output_XS)
end

function test_sine_integral(FT)
    grid = Grid(FT, 4, 16, 1, 1)
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    @. expected_output = integral_of_trial_function_for_sine_operators(z_array)

    # out-of-place
    input_XS = setup_input_XS_for_sine_operators(domain)
    output_XC = XCVariable(domain)
    ∫dz!(output_XC, input_XS)
    @test expected_output ≈ cosine_transform(output_XC)

    # new_output
    input_XS = setup_input_XS_for_sine_operators(domain)
    output_XC = ∫dz(input_XS)
    @test expected_output ≈ cosine_transform(output_XC)
end

function test_sine_second_integral(FT)
    grid = Grid(FT, 4, 16, 1, 1)
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    @. expected_output = second_integral_of_trial_function_for_sine_operators(z_array)

    # out-of-place
    input_XS = setup_input_XS_for_sine_operators(domain)
    output_XS = XSVariable(domain)
    ∫dz²!(output_XS, input_XS)
    @test expected_output ≈ sine_transform(output_XS)

    # new_output
    input_XS = setup_input_XS_for_sine_operators(domain)
    output_XS = ∫dz²(input_XS)
    @test expected_output ≈ sine_transform(output_XS)
end



@testset "z - sine derivatives" for FT in FLOAT_TYPES
    @info "\tTesting z derivatives in sine space for $(FT)"
    test_sine_derivative(FT)
end

@testset "z - sine 2nd derivatives" for FT in FLOAT_TYPES
    @info "\tTesting z 2nd derivatives in sine space for $(FT)"
    test_sine_second_derivative(FT)
end

@testset "z - sine integral" for FT in FLOAT_TYPES
    @info "\tTesting z integrals in sine space for $(FT)"
    test_sine_integral(FT)
end

@testset "z - sine 2nd integral" for FT in FLOAT_TYPES
    @info "\tTesting z 2nd integrals in sine space for $(FT)"
    test_sine_second_integral(FT)
end
