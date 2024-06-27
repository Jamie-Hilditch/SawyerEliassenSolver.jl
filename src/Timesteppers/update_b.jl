struct UpdateB!
    spectral_domain::SpectralDomain
    background_flow::BackgroundFlow
    transforms::Transforms
    scratch::Scratch

    function UpdateB!(problem::Problem)
        return new(
            problem.spectral_domain,
            problem.background_flow,
            problem.transforms,
            problem.scratch,
        )
    end
end

"""Time evolution of buoyancy perturbations

    arguments:
        b: real array
        in: complex array (not modified) containing h*(b1'*ψ^n+c1 + b2'*ψ^n+c2)
    """

function (bt::UpdateB!)(b, in)
    # extract variables out of SE
    CNX = bt.spectral_domain.CNX
    CNZ = bt.spectral_domain.CNZ
    kx = bt.spectral_domain.kx
    kz = bt.spectral_domain.kz

    bg = bt.background_flow

    # alias the temporary working arrays
    ctmp = bt.scratch.ctmp
    rtmp1 = bt.scratch.rtmp1
    rtmp2 = bt.scratch.rtmp2

    # begin algorithm

    # compute -u
    @inbounds @. ctmp[1:CNX, 2:(CNZ + 1)] = kz * in[1:CNX, 1:CNZ]
    @inbounds @. ctmp[(CNX + 1):end, :] = 0
    @inbounds @. ctmp[1:CNX, 1] = 0 # constant 0th cosine mode
    @inbounds @. ctmp[1:CNX, (CNZ + 2):end] = 0
    ldiv!(rtmp1, bt.transforms.fourier, ctmp)
    ldiv!(rtmp2, bt.transforms.cosine, rtmp1)

    # b = b -u*Bx
    @. b += rtmp2 * bg.Bx

    # compute w
    @inbounds @. ctmp[1:CNX, 1:CNZ] = kx * im * in[1:CNX, 1:CNZ]
    @inbounds @. ctmp[(CNX + 1):end, :] = 0
    @inbounds @. ctmp[1:CNX, (CNZ + 1):end] = 0
    # transform to physical space
    ldiv!(rtmp1, bt.transforms.fourier, ctmp)
    ldiv!(rtmp2, bt.transforms.sine, rtmp1)

    # b = b - w*Bz
    @inbounds @. b -= rtmp2 * bg.Bz

    return nothing
end
