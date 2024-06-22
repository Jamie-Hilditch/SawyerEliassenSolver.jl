module Variables

using DocStringExtensions
using LinearAlgebra

using SawyerEliassenSolver.Domains

export
    XZVariable,
    XSVariable,
    XCVariable,
    FZVariable,
    FSVariable,
    FCVariable,
    transform!,
    horizontal_transform,
    sine_transform,
    cosine_transform

include("variable_constructors.jl")
include("transforms.jl")
end
