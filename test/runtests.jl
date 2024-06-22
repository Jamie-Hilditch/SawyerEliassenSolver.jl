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
    @testset "Variables" begin
        @info "Testing variables ..."
        include("test_variables.jl")
    end
end
