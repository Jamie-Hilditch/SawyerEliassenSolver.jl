struct UpdateV!
    spectral_domain :: SpectralDomain
    background_flow :: BackgroundFlow 
    transforms :: Transforms
    scratch :: Scratch

    function UpdateV!(problem :: Problem)
        new(
            problem.spectral_domain,
            problem.background_flow,
            problem.transforms,
            problem.scratch
        )
    end
end

"""Time evolution of along-front velocity perturbations

    arguments:
        v: real array 
        in: complex array (not modified) containing h*(b1'*ψ^n+c1 + b2'*ψ^n+c2)
    """

function (vt::UpdateV!)(v,in)
    # extract variables out of SE
    CNX = vt.spectral_domain.CNX
    CNZ = vt.spectral_domain.CNZ
    kx = vt.spectral_domain.kx 
    kz = vt.spectral_domain.kz

    bg = vt.background_flow
    f = bg.f

    # alias the temporary working arrays
    ctmp = vt.scratch.ctmp
    rtmp1 = vt.scratch.rtmp1
    rtmp2 = vt.scratch.rtmp2

    # begin algorithm

    # compute -u
    @inbounds @. ctmp[1:CNX,2:CNZ+1] = kz*in[1:CNX,1:CNZ]
    @inbounds @. ctmp[CNX+1:end,:] = 0
    @inbounds @. ctmp[1:CNX,1] = 0 # constant 0th cosine mode 
    @inbounds @. ctmp[1:CNX,CNZ + 2:end] = 0
    ldiv!(rtmp1,vt.transforms.fourier,ctmp)
    ldiv!(rtmp2,vt.transforms.cosine,rtmp1)

    # v = v -u*(f + V_x)
    @. v += rtmp2*(f + bg.Vx)

    # compute w 
    @inbounds @. ctmp[1:CNX,1:CNZ] = kx*im*in[1:CNX,1:CNZ]
    @inbounds @. ctmp[CNX+1:end,:] = 0 
    @inbounds @. ctmp[1:CNX,CNZ+1:end] = 0
    # transform to physcial space
    ldiv!(rtmp1,vt.transforms.fourier,ctmp) 
    ldiv!(rtmp2,vt.transforms.sine,rtmp1)

    # v = v - w*Bx/f
    @. v -= rtmp2*bg.Bx/f

    return nothing
end