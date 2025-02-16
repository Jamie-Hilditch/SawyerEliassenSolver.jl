struct Transforms{T<:SingleOrDouble}
    fourier
    sine
    cosine
    inverse_fourier
    inverse_sine
    inverse_cosine
end

function Transforms(real_variable::Array{T,2}) where {T<:SingleOrDouble}

    # plan the forward transforms
    fourier = FFTW.plan_rfft(real_variable, 1; flags=FFTW.PATIENT)
    sine = FFTW.plan_r2r(real_variable, FFTW.RODFT10, 2; flags=FFTW.PATIENT)
    cosine = FFTW.plan_r2r(real_variable, FFTW.REDFT10, 2; flags=FFTW.PATIENT)
    # plan the inverse transforms
    inverse_fourier = inv(fourier)
    inverse_sine = inv(sine)
    inverse_cosine = inv(cosine)

    return Transforms{T}(fourier, sine, cosine, inverse_fourier, inverse_sine, inverse_cosine)
end

Transforms(grid::Grid{T}) where {T} = Transforms(zeros(T, size(grid)))

function Base.show(io::IO, ::MIME"text/plain", transforms::Transforms)
    return print(
        io,
        "Transforms:\n",
        "  ├── fourier: $(summary(transforms.fourier))\n",
        "  ├───── sine: $(summary(transforms.sine))\n",
        "  └─── cosine: $(summary(transforms.cosine))\n",
    )
end

function Base.summary(io::IO, ::Transforms)
    return print(io, "FFTW transforms: rfft, type II DST and type II DCT.")
end
