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

###################
## X derivatives ##
###################

# out-of-place nth x derivative in Fourier space
@inline function ∂x!(out::V, in::V, n::Int) where {V<:FVariable}
    CNX = in.domain.spectral.CNX
    kx = in.domain.spectral.kx
    @. out[1, :] = 0
    @. out[2:CNX, :] = (1im * kx[2:CNX])^n * in[2:CNX, :]
    @. out[(CNX + 1):end, :] = 0
    return nothing
end

# new output nth derivative in Fourier space
@inline function ∂x(v::FVariable, n::Int)
    out = similar(v)
    ∂x!(out, v, n)
    return out
end

# out-of-place nth x derivative in physical space
@inline function ∂x!(out::V, in::V, n::Int) where {V<:XVariable}
    v˜ = horizontal_transform(in)
    ∂x!(v˜, n)
    horizontal_transform!(out, v˜)
    return nothing
end

# new output nth derivative in physical space
@inline function ∂x(v::XVariable, n::Int)
    v˜ = horizontal_transform(v)
    ∂x!(v˜, n)
    return horizontal_transform(v˜)
end

# out-of-place nth x derivative in physical space with scratch variable
@inline function ∂x!(out::V, in::V, tmp::FVariable, n::Int) where {V<:XVariable}
    horizontal_transform!(tmp, in) # n.b. this errors if in and tmp aren't compatible
    ∂x!(tmp, n)
    horizontal_transform!(out, tmp)
    return nothing
end

# new output nth x derivative in physical space with scratch variable
@inline function ∂x(v::XVariable, tmp::FVariable, n::Int)
    horizontal_transform!(tmp, v) # n.b. this errors if v and tmp aren't compatible
    ∂x!(tmp, n)
    return horizontal_transform(tmp)
end

# inplace nth x derivatives in both Fourier and physical space and physical space with scratch
@inline ∂x!(v::AbstractVariable, n::Int) = ∂x!(v, v, n)
@inline ∂x!(v::XVariable, tmp::FVariable, n::Int) = ∂x!(v, v, tmp, n)

# single derivative
@inline ∂x!(out::V, in::V) where {V<:AbstractVariable} = ∂x!(out, in, 1)
@inline ∂x(v::AbstractVariable) = ∂x(v, 1)
@inline ∂x!(out::V, in::V, tmp::FVariable) where {V<:XVariable} = ∂x!(out, in, tmp, 1)
@inline ∂x!(v::AbstractVariable) = ∂x!(v, v, 1)
@inline ∂x(v::XVariable, tmp::FVariable) = ∂x(v, tmp, 1)
@inline ∂x!(v::XVariable, tmp::FVariable) = ∂x!(v, v, tmp, 1)

#################
## X integrals ##
#################

# nth x integral
@inline ∫dx!(out::V, in::V, n::Int) where {V<:AbstractVariable} = ∂x!(out, in, -n)
@inline ∫dx!(out::V, in::V, tmp::FVariable, n::Int) where {V<:XVariable} =
    ∂x!(out, in, tmp, -n)
@inline ∫dx(v::AbstractVariable, n::Int) = ∂x(v, -n)
@inline ∫dx!(v::AbstractVariable, n::Int) = ∂x!(v, -n)
@inline ∫dx(v::XVariable, tmp::FVariable, n::Int) = ∂x(v, tmp, -n)
@inline ∫dx!(v::XVariable, tmp::FVariable, n::Int) = ∂x!(v, tmp, -n)

# single x integral
@inline ∫dx!(out::V, in::V) where {V<:AbstractVariable} = ∂x!(out, in, -1)
@inline ∫dx!(out::V, in::V, tmp::FVariable) where {V<:XVariable} = ∂x!(out, in, tmp, -1)
@inline ∫dx(v::AbstractVariable) = ∂x(v, -1)
@inline ∫dx!(v::AbstractVariable) = ∂x!(v, -1)
@inline ∫dx(v::XVariable, tmp::FVariable) = ∂x(v, tmp, -1)
@inline ∫dx!(v::XVariable, tmp::FVariable) = ∂x!(v, tmp, -1)

###################
## z derivatives ##
###################

# out-of-place sine to cosine in physical X space
@inline function ∂z!(out::XCVariable, in::XSVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1] = 0
    @. out[:, 2:(CNZ + 1)] = kz * in[:, 1:CNZ]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

# out-of-place cosine to sine in physical X space
@inline function ∂z!(out::XSVariable, in::XCVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = -kz * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

# out-of-place sine to cosine in Fourier space
@inline function ∂z!(out::FCVariable, in::FSVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1] = 0
    @. out[:, 2:(CNZ + 1)] = kz * in[:, 1:CNZ]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

# out-of-place cosine to sine in Fourier space
@inline function ∂z!(out::FSVariable, in::FCVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = -kz * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

# new output sine to cosine in physical X space
@inline function ∂z(v::XSVariable)
    out = XCVariable(v.domain)
    ∂z!(out, v)
    return out
end

# new output cosine to sine in physical X space
@inline function ∂z(v::XCVariable)
    out = XSVariable(v.domain)
    ∂z!(out, v)
    return out
end

# new output sine to cosine in Fourier space
@inline function ∂z(v::FSVariable)
    out = FCVariable(v.domain)
    ∂z!(out, v)
    return out
end

# new output cosine to sine in Fourier space
@inline function ∂z(v::FCVariable)
    out = FSVariable(v.domain)
    ∂z!(out, v)
    return out
end

#######################
## 2nd Z derivatives ##
#######################

# out-of-place sine to sine
@inline function ∂z²!(out::T, in::T) where {T<:Union{XSVariable,FSVariable}}
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = -kz^2 * in[:, 1:CNZ]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

# out-of-place cosine to cosine
@inline function ∂z²!(out::T, in::T) where {T<:Union{XCVariable,FCVariable}}
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 0] = 0
    @. out[:, 2:(CNZ + 1)] = -kz^2 * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

# in-place 2nd z derivative
@inline ∂z²!(v::AbstractVariable) = ∂z²!(v, v)

# new output 2nd z derivative
@inline function ∂z²(v::T) where {T<:Union{XSVariable,FSVariable,XCVariable,FCVariable}}
    out = similar(v)
    ∂z²!(out, v)
    return out
end

#################
## Z integrals ##
#################

# out-of-place sine to cosine in physical X space
@inline function ∫dz!(out::XCVariable, in::XSVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1] = 0
    @. out[:, 2:(CNZ + 1)] = -1 / kz * in[:, 1:CNZ]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

# out-of-place cosine to sine in physical X space
@inline function ∫dz!(out::XSVariable, in::XCVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = 1 / kz * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

# out-of-place sine to cosine in Fourier space
@inline function ∫dz!(out::FCVariable, in::FSVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1] = 0
    @. out[:, 2:(CNZ + 1)] = -1 / kz * in[:, 1:CNZ]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

# out-of-place cosine to sine in Fourier space
@inline function ∫dz!(out::FSVariable, in::FCVariable)
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = 1 / kz * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

# new output sine to cosine in physical X space
@inline function ∫dz(v::XSVariable)
    out = XCVariable(v.domain)
    ∫dz!(out, v)
    return out
end

# new output cosine to sine in physical X space
@inline function ∫dz(v::XCVariable)
    out = XSVariable(v.domain)
    ∫dz!(out, v)
    return out
end

# new output sine to cosine in Fourier space
@inline function ∫dz(v::FSVariable)
    out = FCVariable(v.domain)
    ∫dz!(out, v)
    return out
end

# new output cosine to sine in Fourier space
@inline function ∫dz(v::FCVariable)
    out = FSVariable(v.domain)
    ∫dz!(out, v)
    return out
end

#####################
## 2nd Z integrals ##
#####################

# out-of-place sine to sine
@inline function ∫dz²!(out::T, in::T) where {T<:Union{XSVariable,FSVariable}}
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 1:CNZ] = -1 / kz^2 * in[:, 1:CNZ]
    @. out[:, (CNZ + 1):end] = 0
    return nothing
end

# out-of-place cosine to cosine
@inline function ∫dz²!(out::T, in::T) where {T<:Union{XCVariable,FCVariable}}
    CNZ = in.domain.spectral.CNZ
    kz = in.domain.spectral.kz
    @. out[:, 0] = 0
    @. out[:, 2:(CNZ + 1)] = -1 / kz^2 * in[:, 2:(CNZ + 1)]
    @. out[:, (CNZ + 2):end] = 0
    return nothing
end

# in-place 2nd Z integral
@inline ∫dz²!(v::AbstractVariable) = ∫dz²!(v, v)

# new output 2nd z integral
@inline function ∫dz²(v::T) where {T<:Union{XSVariable,FSVariable,XCVariable,FCVariable}}
    out = similar(v)
    ∫dz²!(out, v)
    return out
end
