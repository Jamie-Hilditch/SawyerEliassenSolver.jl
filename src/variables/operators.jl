"""Document me!"""
function ∂x! end
"""Document me!"""
function ∂x end
"""Document me!"""
function ∫dx! end
"""Document me!"""
function ∫dx end
"""Document me!"""
function ∂z! end
"""Document me!"""
function ∂z end
"""Document me!"""
function ∂z²! end
"""Document me!"""
function ∂z² end
"""Document me!"""
function ∫dz! end
"""Document me!"""
function ∫dz end
"""Document me!"""
function ∫dz²! end
"""Document me!"""
function ∫dz² end

# out-of-place nth x derivative in Fourier space
@inline function ∂x!(out::V, in::V, n::Int) where {V<:FVariable}
    CNX = in.domain.spectral.CNX
    kx = in.domain.spectral.kx
    @. out[1:CNX, :] = (1im * kx)^n * in[1:CNX, :]
    @. out[(CNX + 1, :)] = 0
    return nothing
end

# new output nth derivative in Fourier space
@inline function ∂x(v::FVariable, n::Int)
    out = similar(v)
    ∂x!(out, in, n)
    return out
end

# out-of-place nth x derivative in physical space
@inline function ∂x!(out::V, in::V, n::Int) where {V<:XVariable}
    v˜ = horizontal_transform(in)
    ∂x!(v˜, n)
    transform!(out, v˜)
    return nothing
end

# new output nth derivative in physical space
@inline function ∂x(v::XVariable, n::Int)
    v˜ = horizontal_transform(v)
    ∂x!(v˜, n)
    return horizontal_transform(v˜)
end

# inplace nth x derivatives in both Fourier and physical space
@inline ∂x!(v::AbstractVariable, n::Int) = ∂x!(v, v, n)

# single derivative
@inline ∂x!(out::V, in::V) where {V<:AbstractVariable} = ∂x!(out, in, 1)
@inline ∂x(v::AbstractVariable) = ∂x(v, 1)
@inline ∂x!(v::AbstractVariable) = ∂x!(v, v, 1)

# nth x integral
@inline ∫dx!(out::V, in::V, n::Int) where {V<:AbstractVariable} = ∂x!(out, in, -n)
@inline ∫dx(v::AbstractVariable, n::Int) = ∂x(v, -n)
@inline ∫dx!(v::AbstractVariable, n::Int) = ∂x!(v, -n)

# single x integral
@inline ∫dx!(out::V, in::V) where {V<:AbstractVariable} = ∂x!(out, in, -1)
@inline ∫dx(v::AbstractVariable) = ∂x(v, -1)
@inline ∫dx!(v::AbstractVariable) = ∂x!(v, -1)

###################
## z derivatives ##
###################

@inline function ∂z!(out::XCVariable, in::XSVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1] = 0
    @. out[:, 2:(CNZ + 1)] = kz * in[:, 1:CNZ]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

@inline function ∂z!(out::XSVariable, in::XCVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = -kz * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

@inline function ∂z!(out::FCVariable, in::FSVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1] = 0
    @. out[:, 2:(CNZ + 1)] = kz * in[:, 1:CNZ]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

@inline function ∂z!(out::FSVariable, in::FCVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = -kz * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

@inline function ∂z(v::XSVariable)
    out = XCVariable(v.domain)
    ∂z!(out, v)
    return out
end

@inline function ∂z(v::XCVariable)
    out = XSVariable(v.domain)
    ∂z!(out, v)
    return out
end


@inline function ∂z(v::FSVariable)
    out = FCVariable(v.domain)
    ∂z!(out, v)
    return out
end

@inline function ∂z(v::FCVariable)
    out = FSVariable(v.domain)
    ∂z!(out, v)
    return out
end

@inline function ∂z²!(out::T, in::T) where {T <: Union{XSVariable,FSVariable}}
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = -kz^2 * in[:, 1:CNZ]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

@inline function ∂z²!(out::T, in::T) where {T <: Union{XCVariable,FCVariable}}
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 0] = 0
    @. out[:, 2:(CNZ + 1)] = -kz^2 * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

@inline ∂z²!(v::AbstractVariable) = ∂z²!(v,v)

@inline function ∫dz!(out::XCVariable,in::XSVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1] = 0
    @. out[:, 2:(CNZ + 1)] = -1 / kz * in[:, 1:CNZ]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

@inline function ∫dz!(out::XSVariable,in::XCVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = 1 / kz * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

@inline function ∫dz!(out::FCVariable, in::FSVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1] = 0
    @. out[:, 2:(CNZ + 1)] = -1 / kz * in[:, 1:CNZ]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

@inline function ∫dz!(out::FSVariable, in::FCVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = 1 / kz * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

@inline function ∫dz(v::XSVariable)
    out = XCVariable(v.domain)
    ∫dz!(out, v)
    return out
end

@inline function ∫dz(v::XCVariable)
    out = XSVariable(v.domain)
    ∫dz!(out, v)
    return out
end


@inline function ∫dz(v::FSVariable)
    out = FCVariable(v.domain)
    ∫dz!(out, v)
    return out
end

@inline function ∫dz(v::FCVariable)
    out = FSVariable(v.domain)
    ∫dz!(out, v)
    return out
end

@inline function ∫dz²!(out::T, in::T) where {T <: Union{XSVariable,FSVariable}}
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = -1 / kz^2 * in[:, 1:CNZ]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

@inline function ∫dz²!(out::T, in::T) where {T <: Union{XCVariable,FCVariable}}
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 0] = 0
    @. out[:, 2:(CNZ + 1)] = -1 / kz^2 * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

@inline ∫dz²!(v::AbstractVariable) = ∫dz²!(v,v)
