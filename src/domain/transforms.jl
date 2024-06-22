struct Transforms{T<:AbstractFloat}
    fourier
    sine
    cosine
end

function Transforms(real_variable::Array{T,2}) where {T<:AbstractFloat}

    # plan the forward transforms
    fourier = FFTW.plan_rfft(real_variable, 1; flags=FFTW.PATIENT)
    sine = FFTW.plan_r2r(real_variable, FFTW.RODFT10, 2; flags=FFTW.PATIENT)
    cosine = FFTW.plan_r2r(real_variable, FFTW.REDFT10, 2; flags=FFTW.PATIENT)
    # plan the inverse transforms
    inv(fourier)
    inv(sine)
    inv(cosine)

    return Transforms{T}(fourier, sine, cosine)
end

Transforms(grid::Grid) = Transforms(zeros(size(grid)))

function Base.show(io::IO, ::MIME"text/plain", transforms::Transforms)
    return print(
        io,
        "Transforms:\n",
        "  ├── fourier: $(summary(transforms.fourier))\n",
        "  ├───── sine: $(summary(transforms.sine))\n",
        "  └─── cosine: $(summary(transforms.cosine))\n",
    )
end
