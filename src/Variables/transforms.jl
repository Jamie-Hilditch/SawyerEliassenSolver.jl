## Unsafe horizontal transforms for internal use. ##

"""Unsafe horizontal transformation from physical to Fourier space."""
@inline function Tᴴ!(out::FVariable{T}, in::XVariable{T}) where {T}
    fourier = in.domain.transforms.fourier
    mul!(out.data, fourier, in.data)
    # zero out the high wavenumbers for dealiasing
    CNX = out.domain.spectral.CNX
    @inbounds out[(CNX + 1):end, :] .= 0
    return out
end

"""Unsafe horizontal transformation from Fourier to physical space."""
@inline function Tᴴ!(out::XVariable{T}, in::FVariable{T}) where {T}
    inverse_fourier = in.domain.transforms.inverse_fourier
    mul!(out.data, inverse_fourier, in.data)
    return out
end

## Unsafe sine transform for internal use. ##

"""Unsafe sine transformation from physical to spectral space."""
@inline function Tˢ!(out::XSVariable{T}, in::XZVariable{T}) where {T}
    sine = in.domain.transforms.sine
    mul!(out.data, sine, in.data)
    # zero out the high wavenumbers for dealiasing
    CNZ = out.domain.spectral.CNZ
    @inbounds out[:, (CNZ + 1):end] .= 0
    return out
end

"""Unsafe sine transformation from spectral to physical space."""
@inline function Tˢ!(out::XZVariable{T}, in::XSVariable{T}) where {T}
    inverse_sine = in.domain.transforms.inverse_sine
    mul!(out.data, inverse_sine, in.data)
    return out
end

## Unsafe cosine transform for internal use. ##

"""Unsafe cosine transformation from physical to spectral space."""
@inline function Tᶜ!(out::XCVariable{T}, in::XZVariable{T}) where {T}
    cosine = in.domain.transforms.cosine
    mul!(out.data, cosine, in.data)
    # zero out the high wavenumbers for dealiasing
    CNZ = out.domain.spectral.CNZ
    @inbounds out[:, (CNZ + 1):end] .= 0
    return out
end

"""Unsafe cosine transformation from spectral to physical space."""
@inline function Tᶜ!(out::XZVariable{T}, in::XCVariable{T}) where {T}
    inverse_cosine = in.domain.transforms.inverse_cosine
    mul!(out.data, inverse_cosine, in.data)
    return out
end

## all valid combinations of inplace transformations ##
# horizontal
@inline _transform!(out::FZVariable{T}, in::XZVariable{T}) where {T} = Tᴴ!(out, in)
@inline _transform!(out::FSVariable{T}, in::XSVariable{T}) where {T} = Tᴴ!(out, in)
@inline _transform!(out::FCVariable{T}, in::XCVariable{T}) where {T} = Tᴴ!(out, in)
@inline _transform!(out::XZVariable{T}, in::FZVariable{T}) where {T} = Tᴴ!(out, in)
@inline _transform!(out::XSVariable{T}, in::FSVariable{T}) where {T} = Tᴴ!(out, in)
@inline _transform!(out::XCVariable{T}, in::FCVariable{T}) where {T} = Tᴴ!(out, in)
# sine
@inline _transform!(out::XSVariable{T}, in::XZVariable{T}) where {T} = Tˢ!(out, in)
@inline _transform!(out::XZVariable{T}, in::XSVariable{T}) where {T} = Tˢ!(out, in)
# cosine
@inline _transform!(out::XCVariable{T}, in::XZVariable{T}) where {T} = Tᶜ!(out, in)
@inline _transform!(out::XZVariable{T}, in::XCVariable{T}) where {T} = Tᶜ!(out, in)

######################
## Public functions ##
######################

"""Safe transforms with domain validation for public use."""
@inline function transform!(out::AbstractVariable, in::AbstractVariable)
    consistent_domains(out, in) || error("in and out must have the same domain")
    _transform!(out, in)
    return nothing
end

# horizontal transformations creating the output variable
"""Transform the variable from physical space to Fourier space or vice verse in the horizontal."""
function horizontal_transform end
horizontal_transform(in::XZVariable) = Tᴴ!(FZVariable(in.domain), in)
horizontal_transform(in::XSVariable) = Tᴴ!(FSVariable(in.domain), in)
horizontal_transform(in::XCVariable) = Tᴴ!(FCVariable(in.domain), in)
horizontal_transform(in::FZVariable) = Tᴴ!(XZVariable(in.domain), in)
horizontal_transform(in::FSVariable) = Tᴴ!(XSVariable(in.domain), in)
horizontal_transform(in::FCVariable) = Tᴴ!(XCVariable(in.domain), in)
"""$(TYPEDSIGNATURES)"""
@propagate_inbounds horizontal_transform!(out::FVariable, in::XVariable) = transform!(
    out, in
)
"""$(TYPEDSIGNATURES)"""
@propagate_inbounds horizontal_transform!(out::XVariable, in::FVariable) = transform!(
    out, in
)

# sine transforms creating the output variable
"""Transform the variable from physical space to sine space or vice verse in the vertical."""
function sine_transform end
sine_transform(in::XZVariable) = Tˢ!(XSVariable(in.domain), in)
sine_transform(in::XSVariable) = Tˢ!(XZVariable(in.domain), in)
"""$(TYPEDSIGNATURES)"""
@propagate_inbounds sine_transform!(out::XSVariable, in::XZVariable) = transform!(out, in)
"""$(TYPEDSIGNATURES)"""
@propagate_inbounds sine_transform!(out::XZVariable, in::XSVariable) = transform!(out, in)

# cosine transforms creating the output variable
"""Transform the variable from physical space to cosine space or vice verse in the vertical."""
function cosine_transform end
cosine_transform(in::XZVariable) = Tᶜ!(XCVariable(in.domain), in)
cosine_transform(in::XCVariable) = Tᶜ!(XZVariable(in.domain), in)
"""$(TYPEDSIGNATURES)"""
@propagate_inbounds cosine_transform!(out::XCVariable, in::XZVariable) = transform!(out, in)
"""$(TYPEDSIGNATURES)"""
@propagate_inbounds cosine_transform!(out::XZVariable, in::XCVariable) = transform!(out, in)
