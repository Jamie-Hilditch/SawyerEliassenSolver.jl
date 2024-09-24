const DUMMY_COORDINATES_FOR_CONSTRUCTOR = [1.0,2.0,3.0,4.0,5.0]
const DIMENSION_OPTIONS_FOR_CONSTRUCTOR = (
    x=DUMMY_COORDINATES_FOR_CONSTRUCTOR,
    z=DUMMY_COORDINATES_FOR_CONSTRUCTOR,
    p=DUMMY_COORDINATES_FOR_CONSTRUCTOR,
    q=DUMMY_COORDINATES_FOR_CONSTRUCTOR,
)

function create_problem_for_testing_output_writer(FT)
    grid = Grid(FT, 4, 4, 1, 1)
    domain = Domain(grid)
    bg = BackgroundFlow(grid)
    return Problem(domain, bg)
end

function create_dimensions_for_testing_output_writer(FT)
    coordinates = Vector{FT}([1.0,2.0,3.0,4.0,5.0])
    dimension_names = (:x, :z, :p, :q)
    combinations = powerset(dimension_names)
    return map(names -> (; (names .=> tuple(coordinates))...), combinations)
end

function test_output_writer_constructor_no_dims(FT)
    problem = create_problem_for_testing_output_writer(FT)
    filepath = joinpath(TEST_OUTPUT_DIRECTORY, "ow_no_dims.h5")

    try
        # create a file
        OutputWriter(problem, filepath)
        @test isfile(filepath)

        # check we can overwrite it
        OutputWriter(problem, filepath, overwrite=true)
        @test isfile(filepath)

        # check we can't overwrite it
        @test_throws ArgumentError OutputWriter(problem, filepath)
    finally
        rm(filepath, force=true)
    end
end

function test_output_writer_constructor_with_dims(FT, dims)
    problem = create_problem_for_testing_output_writer(FT)
    filepath = joinpath(TEST_OUTPUT_DIRECTORY, "ow_with_dims.h5")

    try
        # create a file
        OutputWriter(problem, filepath, dims)
        @test isfile(filepath)

        # check we can overwrite it
        OutputWriter(problem, filepath, dims, overwrite=true)
        @test isfile(filepath)

        # check we can't overwrite it
        @test_throws ArgumentError OutputWriter(problem, filepath, dims)
    finally
        rm(filepath, force=true)
    end
end

@testset "OutputWriter Constructor" for FT in FLOAT_TYPES
    @testset "No dimensions" begin
        @info "\tTesting OutputWriter constructor with no dimensions for $FT"
        test_output_writer_constructor_no_dims(FT)
    end

    @testset "With dimensions" for dims in create_dimensions_for_testing_output_writer(FT)
        @info "\tTesting OutputWriter constructor with dimensions $(keys(dims)) for $FT"
        test_output_writer_constructor_with_dims(FT, dims)
    end
end
