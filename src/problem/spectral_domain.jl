"""Size of spectral domain and wavenumbers
"""
struct SpectralDomain
    SX::Int
    NZ::Int
    CNX::Int
    CNZ::Int
    kx::Array{Float64,2}
    kz::Array{Float64,2}

    function SpectralDomain(grid::Grid, dealias_x::Int, dealias_z::Int)
        dealias_x >= 0 || throw(DomainError(dealias_x, "dealias_x must be positive"))
        dealias_z >= 0 || throw(DomainError(dealias_z, "dealias_z must be positive"))

        CNX = grid.NX//2 - dealias_x
        CNZ = grid.NZ - 1 - dealias_z
        CNX > 0 || throw(DomainError(dealias_x, "dealias_x must be smaller than NX ÷ 2"))
        CNZ > 0 || throw(DomainError(dealias_z, "dealias_z must be smaller than NZ - 1"))

        LX = grid.x_bounds[2] - grid.x_bounds[1]
        LZ = grid.z_bounds[2] - grid.z_bounds[1]
        kx = 2 * π / LX * UnitRange(0, CNX - 1)
        kz = π / LZ * UnitRange(1, CNZ)
        return new(grid.NX//2 + 1, grid.NZ, CNX, CNZ, reshape(kx, :, 1), reshape(kz, 1, :))
    end
end

Base.size(spectral_domain::SpectralDomain) = (spectral_domain.SX, spectral_domain.NZ)
