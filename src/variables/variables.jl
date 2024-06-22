"""
    $(TYPEDEF)

Abstract supertype for all variables with eltype `T`.
"""
abstract type AbstractVariable{T} <: AbstractMatrix{T} end

Base.size(v::AbstractVariable) = size(v.data)
Base.getindex(v::AbstractVariable, i) = getindex(v.data, i)
Base.getindex(v::AbstractVariable, I::NTuple{2,Int}) = getindex(v.data, I)

struct XZVariable{T<:AbstractFloat} <: AbstractVariable{T}
    domain::Domain{T}
    data::Matrix{T}
end

struct FZVariable{T<:AbstractFloat} <: AbstractVariable{Complex{T}}
    domain::Domain{T}
    data::Matrix{Complex{T}}
end

struct XSVariable{T<:AbstractFloat} <: AbstractVariable{T}
    domain::Domain{T}
    data::Matrix{T}
end

struct FSVariable{T<:AbstractFloat} <: AbstractVariable{Complex{T}}
    domain::Domain{T}
    data::Matrix{Complex{T}}
end

struct XCVariable{T<:AbstractFloat} <: AbstractVariable{T}
    domain::Domain{T}
    data::Matrix{T}
end

struct FCVariable{T<:AbstractFloat} <: AbstractVariable{Complex{T}}
    domain::Domain{T}
    data::Matrix{Complex{T}}
end

XVariable = Union{XZVariable,XSVariable,XCVariable}
FVariable = Union{FZVariable,FSVariable,FCVariable}
ZVariable = Union{XZVariable,FZVariable}
SVariable = Union{XSVariable,FSVariable}
CVariable = Union{XCVariable,FCVariable}

# constructors
XZVariable(domain::Domain{T}) where {T} = XZVariable(domain, zeros(T, size(domain)))
function FZVariable(domain::Domain{T}) where {T}
    return FZVariable(domain, zeros(Complex{T}, size(domain.spectral)))
end
XSVariable(domain::Domain{T}) where {T} = XSVariable(domain, zeros(T, size(domain)))
function FSVariable(domain::Domain{T}) where {T}
    return FSVariable(domain, zeros(Complex{T}, size(domain.spectral)))
end
XCVariable(domain::Domain{T}) where {T} = XCVariable(domain, zeros(T, size(domain)))
function FCVariable(domain::Domain{T}) where {T}
    return FCVariable(domain, zeros(Complex{T}, size(domain.spectral)))
end
