using SawyerEliassenSolver
using Test
using Aqua
using JET

float_types = (Float32, Float64)

@testset "SawyerEliassenSolver.jl" begin
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
        @testset "sine operators" begin
            include("test_operators/test_sine_derivatives_and_integrals.jl")
        end
    end
end
