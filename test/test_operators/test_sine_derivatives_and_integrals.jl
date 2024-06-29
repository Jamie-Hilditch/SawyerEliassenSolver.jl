
NX = 4
NZ = 16

@inline trial_function(z) = sin(π * z) - sin(5 * π * z)
@inline derivative_of_trial_function(z) = π * cos(π * z) - 5 * π * cos(5 * π * z)
@inline second_derivative_of_trial_function(z) =
    -π^2 * sin(π * z) + 25 * π^2 * sin(5 * π * z)
@inline integral_of_trial_function(z) = -1 / π * cos(π * z) + 1 / (5 * π) * cos(5 * π * z)
@inline second_integral_of_trial_function(z) =
    -1 / π^2 * sin(π * z) + 1 / (25 * π^2) * cos(5 * π * z)

function test_z_derivatives_method_types()
    grid = Grid(NX, NZ, one(FT), one(FT))
    domain = Domain(grid)

    @test_throws MethodError ∂z(XZVariable(domain))
    @test ∂z(XSVariable(domain)) isa XCVariable
    @test ∂z(XCVariable(domain)) isa XSVariable
    @test_throws MethodError ∂z(FZVariable(domain))
    @test ∂z(FSVariable(domain)) isa FCVariable
    @test ∂z(FCVariable(domain)) isa FSVariable
end

function setup_input_XS(domain)
    z_array = zgridpoints(domain)'
    input_XZ = XZVariable(domain)
    @. input_XZ = trial_function(z_array)
    return sine_transform(input_XZ)
end

function test_sine_derivative(FT)
    grid = Grid(NX, NZ, one(FT), one(FT))
    domain = Domain(grid)
    z_array = zgridpoints(domain)'

    # expected output in physical x and z space
    expected_output = XZVariable(domain)
    @. expected_output = derivative_of_trial_function(z_array)

    # out-of-place
    input_XS = setup_input_XS(domain)
    output_XC = XCVariable(domain)
    ∂z!(output_XC, input_XS)
    @test expected_output ≈ cosine_transform(output_XC)

    # new_output
    input_XS = setup_input_XS(domain)
    output_XC = ∂z(input_XS)
    @test expected_output ≈ cosine_transform(output_XC)
end

@testset "z - sine derivatives" begin
    for FT in float_types
        @info "Testing z derivatives in sine space for $(FT)"
        test_sine_derivative(FT)
    end
end
