"""
    $(TYPEDEF)

Abstract supertype for all variables with eltype `T`.
"""
abstract type AbstractVariable{T} <: AbstractMatrix{T} end

# Variables are just wrappers around a matrix with some domain metadata
Base.parent(v::AbstractVariable) = v.data
Base.IndexStyle(::AbstractVariable) = IndexLinear()
Base.size(v::AbstractVariable) = size(parent(v))

@propagate_inbounds Base.getindex(v::AbstractVariable, i::Int) = getindex(v.data, i)
# @propagate_inbounds Base.getindex(v::AbstractVariable, I::Vararg{Int,N}) where {N} =
#     getindex(v.data, I...)
# @propagate_inbounds Base.getindex(v::AbstractVariable, I...) = getindex(v.data, I...)
@propagate_inbounds Base.setindex!(v::AbstractVariable, value, i::Int) = setindex!(
    v.data, value, i
)
# @propagate_inbounds Base.setindex!(v::AbstractVariable, value, I::Vararg{Int,N}) where {N} =
#     setindex!(v.data, value, I...)
# @propagate_inbounds Base.setindex!(v::AbstractVariable, value, I...) =
#     setindex!(v.data, value, I...)

# forward strided arrays methods onto underlying matrix
Base.strides(v::AbstractVariable) = strides(parent(v))
function Base.unsafe_convert(::Type{Ptr{T}}, v::AbstractVariable{T}) where {T}
    return Base.unsafe_convert(Type{Ptr{T}}, parent(v))
end
Base.elsize(::AbstractVariable{T}) where {T} = Base.elsize(Matrix{T})
Base.eltype(::AbstractVariable{T}) where {T} = T
Base.pointer(v::AbstractVariable) = pointer(parent(v))

"""A variable expressed in physical space in the horizontal and vertical."""
struct XZVariable{T} <: AbstractVariable{T}
    domain::Domain{T}
    data::Matrix{T}
    function XZVariable(domain::Domain{T}, data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain, data)
    end
end

"""A variable expressed in physical space in the horizontal and sine space in the vertical."""
struct XSVariable{T} <: AbstractVariable{T}
    domain::Domain{T}
    data::Matrix{T}
    function XSVariable(domain::Domain{T}, data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain, data)
    end
end

"""A variable expressed in physical space in the horizontal and cosine space in the vertical."""
struct XCVariable{T} <: AbstractVariable{T}
    domain::Domain{T}
    data::Matrix{T}
    function XCVariable(domain::Domain{T}, data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain, data)
    end
end

"""A variable expressed in Fourier space in the horizontal and physical space in the vertical."""
struct FZVariable{T} <: AbstractVariable{Complex{T}}
    domain::Domain{T}
    data::Matrix{Complex{T}}
    function FZVariable(domain::Domain{T}, data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain, data)
    end
end

"""A variable expressed in Fourier space in the horizontal and sine space in the vertical."""
struct FSVariable{T} <: AbstractVariable{Complex{T}}
    domain::Domain{T}
    data::Matrix{Complex{T}}
    function FSVariable(domain::Domain{T}, data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain, data)
    end
end

"""A variable expressed in Fourier space in the horizontal and cosine space in the vertical."""
struct FCVariable{T} <: AbstractVariable{Complex{T}}
    domain::Domain{T}
    data::Matrix{Complex{T}}
    function FCVariable(domain::Domain{T}, data) where {T}
        validate_variable_array(domain, data)
        return new{T}(domain, data)
    end
end

"""Union of types for variables expressed in physical space in the horizontal."""
XVariable{T} = Union{XZVariable{T},XSVariable{T},XCVariable{T}} where {T}
"""Union of types for variables expressed in Fourier space in the horizontal."""
FVariable{T} = Union{FZVariable{T},FSVariable{T},FCVariable{T}} where {T}

# validate arrays are correct size
function validate_variable_array(domain::Domain{T}, A::Matrix{T}) where {T}
    size(domain.grid) == size(A) || throw(
        ArgumentError(
            "Array $(size(A)) must be the same size as `domain.grid` $(size(domain.grid))",
        ),
    )
    return nothing
end
function validate_variable_array(domain::Domain{T}, A::Matrix{Complex{T}}) where {T}
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

XZVariable(v::AbstractVariable) = XZVariable(v.domain)
XSVariable(v::AbstractVariable) = XSVariable(v.domain)
XCVariable(v::AbstractVariable) = XCVariable(v.domain)
FZVariable(v::AbstractVariable) = FZVariable(v.domain)
FSVariable(v::AbstractVariable) = FSVariable(v.domain)
FCVariable(v::AbstractVariable) = FCVariable(v.domain)

Base.similar(v::XZVariable) = XZVariable(v.domain)
Base.similar(v::XSVariable) = XSVariable(v.domain)
Base.similar(v::XCVariable) = XCVariable(v.domain)
Base.similar(v::FZVariable) = FZVariable(v.domain)
Base.similar(v::FSVariable) = FSVariable(v.domain)
Base.similar(v::FCVariable) = FCVariable(v.domain)

# do we actually need to deinfe these?
# Base.copyto!(x::AbstractVariable, y) = copyto!(parent(x), y)
# Base.copyto!(x, y::AbstractVariable) = copyto!(x,parent(y))
# Base.copyto!(x::AbstractVariable, y::AbstractVariable) = copyto!(parent(x), parent(y))
# Base.fill!(v::AbstractVariable, x) = fill!(parent(v), x)

# comparison operators
function Base.isapprox(x::V, y::V; kws...) where {V<:AbstractVariable}
    return x.domain == y.domain && isapprox(parent(x), parent(y); kws...)
end
function Base.:(==)(x::V, y::V) where {V<:AbstractVariable}
    return x.domain == y.domain && parent(x) == parent(y)
end

"""Return the fourier (physical) counterpart to a physical (fourier) variable."""
horizontal_counterpart(::XZVariable) = FZVariable
horizontal_counterpart(::XSVariable) = FSVariable
horizontal_counterpart(::XCVariable) = FCVariable
horizontal_counterpart(::FZVariable) = XZVariable
horizontal_counterpart(::FSVariable) = XSVariable
horizontal_counterpart(::FCVariable) = XCVariable

@inline Domains.get_domain(v::AbstractVariable) = v.domain

"""Return a view into the resolved wavenumbers."""
function resolved end
@inline resolved(v::XZVariable) = view(v.data, :, :)
@inline resolved(v::FZVariable) = view(v.data, resolved_fourier_indices(v), :)
@inline resolved(v::XSVariable) = view(v.data, :, resolved_sine_indices(v))
@inline resolved(v::FSVariable) = view(
    v.data, resolved_fourier_indices(v), resolved_sine_indices(v)
)
@inline resolved(v::XCVariable) = view(v.data, :, resolved_cosine_indices(v))
@inline resolved(v::FCVariable) = view(
    v.data, resolved_fourier_indices(v), resolved_cosine_indices(v)
)
