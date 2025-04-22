function dummy_function_for_output_variables(problem, output_array)
    output_array .= 0
    return nothing
end

function dummy_function_for_output_variables_with_args(problem, output_array, a, b)
    output_array .= a + b
    return nothing
end

function create_problem_for_testing_output_variables(FT)
    grid = Grid(FT, 4, 4, 1, 1)
    domain = Domain(grid)
    bg = BackgroundFlow(grid)
    return Problem(domain, bg)
end

function test_valid_output_variable_constructors(problem::Problem)
    # default
    ov1 = OutputVariable(
        problem,
        dummy_function_for_output_variables,
        (:x, :z),
        zeros(eltype(problem), 4, 4),
        tuple(),
    )
    @test ov1 isa OutputVariable

    # no args
    ov2 = OutputVariable(
        problem,
        dummy_function_for_output_variables,
        (:x, :z),
        zeros(eltype(problem), 4, 4),
    )
    @test ov2 isa OutputVariable

    # build output array from sizes and type
    ov3 =  OutputVariable(
        problem,
        dummy_function_for_output_variables,
        (:x, :z),
        (4, 4),
        eltype(problem),
        tuple(),
    )
    @test ov3 isa OutputVariable

    # build output array from sizes and type, no args
    ov4 = OutputVariable(
        problem,
        dummy_function_for_output_variables,
        (:x, :z),
        (4, 4),
        eltype(problem),
    )
    @test ov4 isa OutputVariable

    # build output array from named tuple of sizes and type
    ov5 = OutputVariable(
        problem,
        dummy_function_for_output_variables,
        (; x=4, z=4),
        eltype(problem),
        tuple(),
    )
    @test ov5 isa OutputVariable

    # build output array from named tuple of sizes and type, no args
    ov6 = OutputVariable(
        problem,
        dummy_function_for_output_variables,
        (; x=4, z=4),
        eltype(problem),
    )
    @test ov6 isa OutputVariable

    # test all constructors give the same type
    @test allequal(typeof, [ov1, ov2, ov3, ov4, ov5, ov6])
end

function test_inconsistent_dimensions_and_array_sizes(problem::Problem)
    @test_throws MethodError OutputVariable(
        problem,
        dummy_function_for_output_variables,
        (:x, :z),
        zeros(eltype(problem), 4, 4, 4),
    )
    @test_throws MethodError OutputVariable(
        problem,
        dummy_function_for_output_variables,
        (:x, :z),
        (4, 4, 4),
    )
end

function test_scalar_output_variable(problem::Problem)
    @test OutputVariable(
        problem,
        dummy_function_for_output_variables,
        tuple(),
        tuple(),
        eltype(problem),
    ) isa OutputVariable
    @test OutputVariable(
        problem,
        dummy_function_for_output_variables,
        tuple(),
        Array{eltype(problem), 0}(undef),
    ) isa OutputVariable
end


@testset "OutputVariable Constructor" for FT in FLOAT_TYPES
    @info "\tTesting OutputVariable constructors for $FT"
    problem = create_problem_for_testing_output_writer(FT)
    test_valid_output_variable_constructors(problem)
    test_inconsistent_dimensions_and_array_sizes(problem)
    test_scalar_output_variable(problem)
end
