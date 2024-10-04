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
    ∂z!,
    ∂z,
    ∂z²!,
    ∂z²,
    ∫dx!,
    ∫dx,
    ∫dz!,
    ∫dz,
    ∫dz²!,
    ∫dz²,
    ∇⁻²!

include("variable_constructors.jl")
include("transforms.jl")
include("operators.jl")

"""
Module for different variable representations, transforms between representations and operators that act on variables.

# Exports
$(EXPORTS)
"""
Variables
end
