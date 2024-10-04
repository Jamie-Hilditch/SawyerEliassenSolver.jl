
NX = 16
NZ = 4

@inline trial_function_for_x_operators(x) = sin(π * x) + cos(2 * π * x)
@inline derivative_of_trial_function_for_x_operators(x) = π * cos(π * x) - 2 * π * sin(2 * π * x)
@inline second_derivative_of_trial_function_for_x_operators(x) =
    -π^2 * sin(π * x) - 4 * π^2 * cos(2 * π * x)
@inline integral_of_trial_function_for_x_operators(x) = -1 / π * cos(π * x) + 1 / (2 * π) * sin(2 * π * x)
@inline second_integral_of_trial_function_for_x_operators(x) =
    -1 / π^2 * sin(π * x) - 1 / (4 * π^2) * cos(2 * π * x)

function test_physical_x_derivatives(FT)
    x_bounds = (-one(FT), one(FT))
    LZ = one(FT)
    grid = Grid(NX, NZ, x_bounds, LZ)
    domain = Domain(grid)

    for x_variable_type in (XZVariable, XSVariable, XCVariable)
        x_array = xgridpoints(grid)
        expected_output = x_variable_type(domain)
        @. expected_output = derivative_of_trial_function_for_x_operators(x_array)

        # out-of-place
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        output = x_variable_type(domain)
        ∂x!(output, x_variable)
        @test output ≈ expected_output

        # in-place
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        ∂x!(x_variable)
        @test x_variable ≈ expected_output

        # new output
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        @test ∂x(x_variable) ≈ expected_output

        # out-of-place with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        output = x_variable_type(domain)
        tmp = horizontal_counterpart(x_variable)(domain)
        ∂x!(output, x_variable, tmp)
        @test output ≈ expected_output

        # in-place with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        tmp = horizontal_counterpart(x_variable)(domain)
        ∂x!(x_variable, tmp)
        @test x_variable ≈ expected_output

        # new output with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        tmp = horizontal_counterpart(x_variable)(domain)
        @test ∂x(x_variable, tmp) ≈ expected_output
    end
end

function test_physical_x_second_derivatives(FT)
    x_bounds = (-one(FT), one(FT))
    LZ = one(FT)
    grid = Grid(NX, NZ, x_bounds, LZ)
    domain = Domain(grid)

    for x_variable_type in (XZVariable, XSVariable, XCVariable)
        x_array = xgridpoints(grid)
        expected_output = x_variable_type(domain)
        @. expected_output = second_derivative_of_trial_function_for_x_operators(x_array)

        # out-of-place
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        output = x_variable_type(domain)
        ∂x!(output, x_variable, 2)
        @test output ≈ expected_output

        # in-place
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        ∂x!(x_variable, 2)
        @test x_variable ≈ expected_output

        # new output
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        @test ∂x(x_variable, 2) ≈ expected_output

        # out-of-place with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        output = x_variable_type(domain)
        tmp = horizontal_counterpart(x_variable)(domain)
        ∂x!(output, x_variable, tmp, 2)
        @test output ≈ expected_output

        # in-place with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        tmp = horizontal_counterpart(x_variable)(domain)
        ∂x!(x_variable, tmp, 2)
        @test x_variable ≈ expected_output

        # new output with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        tmp = horizontal_counterpart(x_variable)(domain)
        @test ∂x(x_variable, tmp, 2) ≈ expected_output
    end
end

function test_physical_x_integrals(FT)
    x_bounds = (-one(FT), one(FT))
    LZ = one(FT)
    grid = Grid(NX, NZ, x_bounds, LZ)
    domain = Domain(grid)

    for x_variable_type in (XZVariable, XSVariable, XCVariable)
        x_array = xgridpoints(domain)
        expected_output = x_variable_type(domain)
        @. expected_output = integral_of_trial_function_for_x_operators(x_array)

        # out-of-place
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        output = x_variable_type(domain)
        ∫dx!(output, x_variable)
        @test output ≈ expected_output

        # in-place
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        ∫dx!(x_variable)
        @test x_variable ≈ expected_output

        # new output
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        @test ∫dx(x_variable) ≈ expected_output

        # out-of-place with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        output = x_variable_type(domain)
        tmp = horizontal_counterpart(x_variable)(domain)
        ∫dx!(output, x_variable, tmp)
        @test output ≈ expected_output

        # in-place with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        tmp = horizontal_counterpart(x_variable)(domain)
        ∫dx!(x_variable, tmp)
        @test x_variable ≈ expected_output

        # new output with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        tmp = horizontal_counterpart(x_variable)(domain)
        @test ∫dx(x_variable, tmp) ≈ expected_output
    end
end

function test_physical_x_second_integrals(FT)
    x_bounds = (-one(FT), one(FT))
    LZ = one(FT)
    grid = Grid(NX, NZ, x_bounds, LZ)
    domain = Domain(grid)

    for x_variable_type in (XZVariable, XSVariable, XCVariable)
        x_array = xgridpoints(grid)
        expected_output = x_variable_type(domain)
        @. expected_output = second_integral_of_trial_function_for_x_operators(x_array)

        # out-of-place
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        output = x_variable_type(domain)
        ∫dx!(output, x_variable, 2)
        @test output ≈ expected_output

        # in-place
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        ∫dx!(x_variable, 2)
        @test x_variable ≈ expected_output

        # new output
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        @test ∫dx(x_variable, 2) ≈ expected_output

        # out-of-place with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        output = x_variable_type(domain)
        tmp = horizontal_counterpart(x_variable)(domain)
        ∫dx!(output, x_variable, tmp, 2)
        @test output ≈ expected_output

        # in-place with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        tmp = horizontal_counterpart(x_variable)(domain)
        ∫dx!(x_variable, tmp, 2)
        @test x_variable ≈ expected_output

        # new output with tmp
        x_variable = x_variable_type(domain)
        @. x_variable = trial_function_for_x_operators(x_array)
        tmp = horizontal_counterpart(x_variable)(domain)
        @test ∫dx(x_variable, tmp, 2) ≈ expected_output
    end
end

@testset "X derivatives and integrals in physical space" begin
    for FT in FLOAT_TYPES
        @info "\tTesting x derivatives in physical space with $(FT)"
        test_physical_x_derivatives(FT)
        @info "\tTesting second x derivatives in physical space with $(FT)"
        test_physical_x_second_derivatives(FT)
        @info "\tTesting x integrals in physical space with $(FT)"
        test_physical_x_integrals(FT)
        @info "\tTesting second x integrals in physical space with $(FT)"
        test_physical_x_second_integrals(FT)
    end
end
