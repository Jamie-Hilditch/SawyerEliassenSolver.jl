using Aqua
using Combinatorics
using JET
using Test

using SawyerEliassenSolver

const FLOAT_TYPES = (Float32, Float64)
const TEST_OUTPUT_DIRECTORY = joinpath(@__DIR__, "_test_output")

@testset "SawyerEliassenSolver.jl" begin
    @info "Creating clean test output directory at $TEST_OUTPUT_DIRECTORY"
    if isdir(TEST_OUTPUT_DIRECTORY)
        rm(TEST_OUTPUT_DIRECTORY; force=true, recursive=true)
    end
    mkdir(TEST_OUTPUT_DIRECTORY)

    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SawyerEliassenSolver)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(SawyerEliassenSolver; target_defined_modules=true)
    end
    @testset "Grid" begin
        @info "Testing grids ..."
        include("test_grid.jl")
    end
    @testset "Variables" begin
        @info "Testing variables ..."
        include("test_variables.jl")
    end
    @testset "Transforms" begin
        @info "Testing transforms ..."
        include("test_transforms.jl")
    end
    @testset "Operators" begin
        @info "Testing operators ..."
        @testset "x operators" begin
            include("test_operators/test_x_derivatives_and_integrals.jl")
        end
        @testset "z operators" begin
            include("test_operators/test_z_derivatives_and_integrals.jl")
            @testset "sine operators" begin
                include("test_operators/test_sine_derivatives_and_integrals.jl")
            end
            @testset "cosine operators" begin
                include("test_operators/test_cosine_derivatives_and_integrals.jl")
            end
            @testset "physical z operators" begin
                include("test_operators/test_physical_z_operators.jl")
            end
        end
    end

    @testset "Forcings" begin
        @info "Testing forcings ..."
        include("test_forcings.jl")
    end

    @testset "OutputWriter" begin
        @info "Testing OutputWriter ..."
        include("test_output_writer.jl")
    end

    @info "Cleaning up test output directory at $TEST_OUTPUT_DIRECTORY"
    rm(TEST_OUTPUT_DIRECTORY; force=true, recursive=true)
end
