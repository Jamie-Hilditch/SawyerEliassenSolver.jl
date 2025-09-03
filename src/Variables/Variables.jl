module Variables

using Base: @propagate_inbounds
using DocStringExtensions
using LinearAlgebra

using SawyerEliassenSolver.Domains

export XZVariable,
    XSVariable,
    XCVariable,
    FZVariable,
    FSVariable,
    FCVariable,
    horizontal_counterpart,
    transform!,
    horizontal_transform,
    horizontal_transform!,
    sine_transform,
    sine_transform!,
    cosine_transform,
    cosine_transform!,
    resolved,
    ∂x!,
    ∂x,
    ∫dx!,
    ∫dx,
    ∂z!,
    ∂z,
    ∂z²!,
    ∂z²,
    ∫dz!,
    ∫dz,
    ∫dz²!,
    ∫dz²,
    solve_poisson!


include("variable_constructors.jl")
include("transforms.jl")
include("operators.jl")

"""
Module for different variable representations, transforms between representations and operators that act on variables.

# Exports
$(EXPORTS)

Most of the operators have many methods signatures to allow for different combinations of input and output variables.
!!! note
    To see all the method signatures for a particular function e.g. [`∂x!`](@ref), run
    ```julia-repl
    julia> ?
    help?> ∂x!
    ```
"""
Variables
end
