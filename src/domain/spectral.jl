"""Size of spectral domain and wavenumbers
"""
struct Spectral{T <: AbstractFloat}
    SX::Int
    NZ::Int
    CNX::Int
    CNZ::Int
    kx::Array{T,2}
    kz::Array{T,2}
end

function Spectral(grid::Grid{T}, CNX, CNZ) where {T}
    SX = grid.NX ÷ 2 + 1
    LX = grid.x_bounds[2] - grid.x_bounds[1]
    LZ = grid.z_bounds[2] - grid.z_bounds[1]
    kx = 2 * π / LX * UnitRange{T}(0, CNX - 1)
    kz = π / LZ * UnitRange{T}(1, CNZ)
    return Spectral{T}(SX, grid.NZ, CNX, CNZ, reshape(kx, :, 1), reshape(kz, 1, :))
end

"""$(TYPEDSIGNATURES)"""
Base.eltype(::Spectral{T}) where{T} = T

"""$(TYPEDSIGNATURES)"""
Base.size(spectral::Spectral) = (spectral.SX, spectral.NZ)
