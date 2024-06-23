struct Domain{T<:SingleOrDouble}
    grid::Grid{T}
    spectral::Spectral{T}
    transforms::Transforms{T}
end

function Domain(grid::Grid{T}, dealias_x::Int=0, dealias_z::Int=0) where {T}
    dealias_x >= 0 || throw(DomainError(dealias_x, "dealias_x must be non negative"))
    dealias_z >= 0 || throw(DomainError(dealias_z, "dealias_z must be non negative"))

    CNX = grid.NX ÷ 2 - dealias_x
    CNZ = grid.NZ - 1 - dealias_z
    CNX > 0 || throw(DomainError(dealias_x, "dealias_x must be smaller than NX ÷ 2"))
    CNZ > 0 || throw(DomainError(dealias_z, "dealias_z must be smaller than NZ - 1"))

    spectral = Spectral(grid, CNX, CNZ)
    transforms = Transforms(grid)
    return Domain{T}(grid, spectral, transforms)
end

Base.size(domain::Domain) = size(domain.grid)
