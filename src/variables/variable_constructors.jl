"""
    $(TYPEDEF)

Abstract supertype for all variables with eltype `T`.
"""
abstract type AbstractVariable{T} <: AbstractMatrix{T} end

Base.size(v::AbstractVariable) = size(v.data)
Base.getindex(v::AbstractVariable, i) = getindex(v.data, i)
Base.getindex(v::AbstractVariable, I::Vararg{Int, N}) where {N} = getindex(v.data, I)
Base.setindex!(v::AbstractVariable, value, i) = setindex!(v.data, value, i)
Base.setindex!(v::AbstractVariable, value, I...) = setindex!(v.data, value, I...)

struct XZVariable{T<:AbstractFloat} <: AbstractVariable{T}
    domain::Domain{T}
    data::Matrix{T}
    function XZVariable(domain::Domain{T},data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain,data)
    end
end

struct XSVariable{T<:AbstractFloat} <: AbstractVariable{T}
    domain::Domain{T}
    data::Matrix{T}
    function XSVariable(domain::Domain{T},data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain,data)
    end
end

struct XCVariable{T<:AbstractFloat} <: AbstractVariable{T}
    domain::Domain{T}
    data::Matrix{T}
    function XCVariable(domain::Domain{T},data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain,data)
    end
end

struct FZVariable{T<:AbstractFloat} <: AbstractVariable{Complex{T}}
    domain::Domain{T}
    data::Matrix{Complex{T}}
    function FZVariable(domain::Domain{T},data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain,data)
    end
end

struct FSVariable{T<:AbstractFloat} <: AbstractVariable{Complex{T}}
    domain::Domain{T}
    data::Matrix{Complex{T}}
    function FSVariable(domain::Domain{T},data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain,data)
    end
end

struct FCVariable{T<:AbstractFloat} <: AbstractVariable{Complex{T}}
    domain::Domain{T}
    data::Matrix{Complex{T}}
    function FCVariable(domain::Domain{T},data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain,data)
    end
end

XVariable{T} = Union{XZVariable{T},XSVariable{T},XCVariable{T}} where T
FVariable{T} = Union{FZVariable{T},FSVariable{T},FCVariable{T}} where T

# validate arrays are correct size
function validate_variable_array(domain::Domain{T}, A::Matrix{T}) where {T<:AbstractFloat}
    size(domain.grid) == size(A) || throw(
        ArgumentError(
            "Array $(size(A)) must be the same size as `domain.grid` $(size(domain.grid))",
        ),
    )
    return nothing
end
function validate_variable_array(
    domain::Domain{T}, A::Matrix{Complex{T}}
) where {T<:AbstractFloat}
    size(domain.spectral) == size(A) || throw(
        ArgumentError(
            "Array $(size(A)) must be the same size as `domain.spectral` $(size(domain.spectral))",
        ),
    )
    return nothing
end

# constructors
XZVariable(domain::Domain{T}) where {T} = XZVariable(domain, zeros(T, size(domain)))
XSVariable(domain::Domain{T}) where {T} = XSVariable(domain, zeros(T, size(domain)))
XCVariable(domain::Domain{T}) where {T} = XCVariable(domain, zeros(T, size(domain)))

function FZVariable(domain::Domain{T}) where {T}
    return FZVariable(domain, zeros(Complex{T}, size(domain.spectral)))
end
function FSVariable(domain::Domain{T}) where {T}
    return FSVariable(domain, zeros(Complex{T}, size(domain.spectral)))
end
function FCVariable(domain::Domain{T}) where {T}
    return FCVariable(domain, zeros(Complex{T}, size(domain.spectral)))
end
