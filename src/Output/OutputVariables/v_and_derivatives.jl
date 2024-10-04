"""$(TYPEDSIGNATURES)"""
v(problem::Problem) = OutputVariable(problem, (p,a) -> nothing, (:x,:z), problem.state.v)
