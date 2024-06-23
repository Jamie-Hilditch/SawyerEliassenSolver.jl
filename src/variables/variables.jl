module Variables

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
    sine_transform,
    cosine_transform,
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
    ∫dz²

include("variable_constructors.jl")
include("transforms.jl")
include("operators.jl")
end
