function missing_variable(problem::SES.Problem, state::SES.State)
    problem.scratch.rtmp1 .= NaN
    return problem.scratch.rtmp1
end

function compute_u!(problem::SES.Problem, state::SES.State)
    ctmp = problem.scratch.ctmp
    out = problem.scratch.rtmp1
    rtmp2 = problem.scratch.rtmp2

    CNX = problem.spectral_domain.CNX
    CNZ = problem.spectral_domain.CNZ
    kz = problem.spectral_domain.kz

    @inbounds @. ctmp[1:CNX, 2:(CNZ + 1)] = -kz * state.ψ[1:CNX, 1:CNZ]
    @inbounds @. ctmp[(CNX + 1):end, :] = 0
    @inbounds @. ctmp[1:CNX, 1] = 0 # constant 0th cosine mode 
    @inbounds @. ctmp[1:CNX, (CNZ + 2):end] = 0
    ldiv!(rtmp2, problem.transforms.fourier, ctmp)
    ldiv!(out, problem.transforms.cosine, rtmp2)
    return out
end

function compute_w!(problem::SES.Problem, state::SES.State)
    ctmp = problem.scratch.ctmp
    out = problem.scratch.rtmp1
    rtmp2 = problem.scratch.rtmp2

    CNX = problem.spectral_domain.CNX
    CNZ = problem.spectral_domain.CNZ
    kx = problem.spectral_domain.kx

    @inbounds @. ctmp[1:CNX, 1:CNZ] = kx * im * state.ψ[1:CNX, 1:CNZ]
    @inbounds @. ctmp[(CNX + 1):end, :] = 0
    @inbounds @. ctmp[1:CNX, (CNZ + 1):end] = 0
    ldiv!(rtmp2, problem.transforms.fourier, ctmp)
    ldiv!(out, problem.transforms.sine, rtmp2)
    return out
end

get_v(problem::SES.Problem, state::SES.State) = state.v
get_b(problem::SES.Problem, state::SES.State) = state.b

function compute_ux!(problem::SES.Problem, state::SES.State)
    ctmp = problem.scratch.ctmp
    out = problem.scratch.rtmp1
    rtmp2 = problem.scratch.rtmp2

    CNX = problem.spectral_domain.CNX
    CNZ = problem.spectral_domain.CNZ
    kx = problem.spectral_domain.kx
    kz = problem.spectral_domain.kz

    @inbounds @. ctmp[1:CNX, 2:(CNZ + 1)] = -kx * im * kz * state.ψ[1:CNX, 1:CNZ]
    @inbounds @. ctmp[(CNX + 1):end, :] = 0
    @inbounds @. ctmp[1:CNX, 1] = 0 # constant 0th cosine mode 
    @inbounds @. ctmp[1:CNX, (CNZ + 2):end] = 0
    ldiv!(rtmp2, problem.transforms.fourier, ctmp)
    ldiv!(out, problem.transforms.cosine, rtmp2)
    return out
end

function compute_uz!(problem::SES.Problem, state::SES.State)
    ctmp = problem.scratch.ctmp
    out = problem.scratch.rtmp1
    rtmp2 = problem.scratch.rtmp2

    CNX = problem.spectral_domain.CNX
    CNZ = problem.spectral_domain.CNZ
    kz = problem.spectral_domain.kz

    @inbounds @. ctmp[1:CNX, 1:CNZ] = kz^2 * state.ψ[1:CNX, 1:CNZ]
    @inbounds @. ctmp[(CNX + 1):end, :] = 0
    @inbounds @. ctmp[1:CNX, (CNZ + 1):end] = 0
    ldiv!(rtmp2, problem.transforms.fourier, ctmp)
    ldiv!(out, problem.transforms.sine, rtmp2)
    return out
end

function compute_wx!(problem::SES.Problem, state::SES.State)
    ctmp = problem.scratch.ctmp
    out = problem.scratch.rtmp1
    rtmp2 = problem.scratch.rtmp2

    CNX = problem.spectral_domain.CNX
    CNZ = problem.spectral_domain.CNZ
    kx = problem.spectral_domain.kx

    @inbounds @. ctmp[1:CNX, 1:CNZ] = -kx^2 * state.ψ[1:CNX, 1:CNZ]
    @inbounds @. ctmp[(CNX + 1):end, :] = 0
    @inbounds @. ctmp[1:CNX, (CNZ + 1):end] = 0
    ldiv!(rtmp2, problem.transforms.fourier, ctmp)
    ldiv!(out, problem.transforms.sine, rtmp2)
    return out
end

function _differentiate_z!(out::Matrix{T}, in::Matrix{T}, dz::Float64) where {T<:Number}
    size(in)[2] >= 5 ||
        throw(ArgumentError("z grid is too short to apply five point stencil"))

    @inbounds @. out[:, 1] =
        -25 / 12 * in[:, 1] + 4 * in[:, 2] - 3 * in[:, 3] + 4 / 3 * in[:, 4] -
        1 / 4 * in[:, 5]
    @inbounds @. out[:, 2] =
        -1 / 4 * in[:, 1] - 5 / 6 * in[:, 2] + 3 / 2 * in[:, 3] - 1 / 2 * in[:, 4] +
        1 / 12 * in[:, 5]
    @inbounds @. out[:, 3:(end - 2)] =
        1 / 12 * in[:, 1:(end - 4)] - 2 / 3 * in[:, 2:(end - 3)] +
        2 / 3 * in[:, 4:(end - 1)] - 1 / 12 * in[:, 5:end]
    @inbounds @. out[:, end - 1] =
        -1 / 12 * in[:, end - 4] + 1 / 2 * in[:, end - 3] - 3 / 2 * in[:, end - 2] +
        5 / 6 * in[:, end - 1] +
        1 / 4 * in[:, end]
    @inbounds @. out[:, end] =
        1 / 4 * in[:, end - 4] - 4 / 3 * in[:, end - 3] + 3 * in[:, end - 2] -
        4 * in[:, end - 1] + 25 / 12 * in[:, end]
    @. out /= dz

    return nothing
end

function compute_vz!(problem::SES.Problem, state::SES.State)
    out = problem.scratch.rtmp1
    dz = problem.grid.z[2] - problem.grid.z[1]
    _differentiate_z!(out, state.v, dz)
    return out
end

function compute_bz!(problem::SES.Problem, state::SES.State)
    out = problem.scratch.rtmp1
    dz = problem.grid.z[2] - problem.grid.z[1]
    _differentiate_z!(out, state.b, dz)
    return out
end

function compute_ψ!(problem::SES.Problem, state::SES.State)
    ctmp = problem.scratch.ctmp
    out = problem.scratch.rtmp1
    rtmp2 = problem.scratch.rtmp2

    CNX = problem.spectral_domain.CNX
    CNZ = problem.spectral_domain.CNZ

    @inbounds @. ctmp[1:CNX, 2:(CNZ + 1)] = state.ψ[1:CNX, 1:CNZ]
    @inbounds @. ctmp[(CNX + 1):end, :] = 0
    @inbounds @. ctmp[1:CNX, 1] = 0 # constant 0th cosine mode 
    @inbounds @. ctmp[1:CNX, (CNZ + 2):end] = 0
    ldiv!(rtmp2, problem.transforms.fourier, ctmp)
    ldiv!(out, problem.transforms.cosine, rtmp2)
    return out
end
