@inline trial_function_for_physical_z_operators(z) = z^3
@inline derivative_of_trial_function_for_physical_z_operators(z) = 3 * z^2

function setup_input_XZ_for_physical_z_operators(domain)
    z_array = zgridpoints(domain)'
    input_XZ = XZVariable(domain)
    @. input_XZ = trial_function_for_physical_z_operators(z_array)
    return input_XZ
end

function test_physical_z_derivative(FT)
    grid = Grid(FT, 4, 128, 1, 1)
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    compare_to_expected_output = isapprox(expected_output; atol=1e-2, rtol=1e-2)
    @. expected_output = derivative_of_trial_function_for_physical_z_operators(z_array)

    # out-of-place
    input_XZ = setup_input_XZ_for_physical_z_operators(domain)
    output_XZ = XZVariable(domain)
    ∂z!(output_XZ, input_XZ)
    @test compare_to_expected_output(output_XZ)

    # new_output
    input_XZ = setup_input_XZ_for_physical_z_operators(domain)
    output_XZ = ∂z(input_XZ)
    @test compare_to_expected_output(output_XZ)
end

@testset "physical z derivative" for FT in FLOAT_TYPES
    @info "\tTesting physical z derivative with $(FT)"
    test_physical_z_derivative(FT)
end
