using SawyerEliassenSolver
using Test
using Aqua
using JET

@testset "SawyerEliassenSolver.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(SawyerEliassenSolver)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(SawyerEliassenSolver; target_defined_modules=true)
    end
    # Write your tests here.
end
