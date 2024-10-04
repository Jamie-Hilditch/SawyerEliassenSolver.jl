"""$(TYPEDSIGNATURES)"""
b(problem::Problem) = OutputVariable(problem, (p,a) -> nothing, (:x,:z), problem.state.b)
