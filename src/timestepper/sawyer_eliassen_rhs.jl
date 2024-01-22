struct SawyerEliassenRHS!
    spectral_domain :: SpectralDomain
    background_flow :: BackgroundFlow 
    transforms :: Transforms
    scratch :: Scratch

    function SawyerEliassenRHS!(problem :: Problem)
        new(
            problem.spectral_domain,
            problem.background_flow,
            problem.transforms,
            problem.scratch
        )
    end
end

"""Linear operator: RHS of the Sawyer-Eliassen equation

    arguments:
        out: complex array
        in: complex array (not modified unless in = out)
"""
function (SE::SawyerEliassenRHS!)(out,in)
    # extract variables out of SE
    CNX = SE.spectral_domain.CNX
    CNZ = SE.spectral_domain.CNZ
    kx = SE.spectral_domain.kx 
    kz = SE.spectral_domain.kz

    bg = SE.background_flow
    f = bg.f

    # alias the scratch arrays
    ctmp = SE.scratch.ctmp
    tmp1 = SE.scratch.rtmp1 # store intermediate stage of transforms
    tmp2 = SE.scratch.rtmp2 # store derivatives of ψ
    ζtt = SE.scratch.rtmp3

    # now start the algorithm

    # since we overwrite out we can use it as temporary work space
    # compute -f(f+Vx)*ψzz accounting for dealiasing
    @inbounds @. ctmp[1:CNX,1:CNZ] = -kz^2*in[1:CNX,1:CNZ]
    @inbounds @. ctmp[CNX+1:end,:] = 0 
    @inbounds @. ctmp[1:CNX,CNZ+1:end] = 0
    # transform to physical space
    ldiv!(tmp1,SE.transforms.fourier,ctmp) 
    ldiv!(tmp2,SE.transforms.sine,tmp1)
    # start constructing ζtt 
    @. ζtt = -f*(f + bg.Vx)*tmp2

    # now 2*Bx*ψxz which requires a cosine transform and a shift in the coefficients
    @inbounds @. ctmp[1:CNX,2:CNZ+1] = 1*im*kx*kz*in[1:CNX,1:CNZ]
    @inbounds @. ctmp[CNX+1:end,:] = 0
    @inbounds @. ctmp[1:CNX,1] = 0 # constant 0th cosine mode 
    @inbounds @. ctmp[1:CNX,CNZ + 2:end] = 0
    # transform to physical space 
    ldiv!(tmp1,SE.transforms.fourier,ctmp) 
    ldiv!(tmp2,SE.transforms.cosine,tmp1)
    # add to ζtt 
    @. ζtt += 2*bg.Bx*tmp2
    
    # now - Bz*ψxx 
    @inbounds @. ctmp[1:CNX,1:CNZ] = -kx^2*in[1:CNX,1:CNZ]
    @inbounds @. ctmp[CNX+1:end,:] = 0 
    @inbounds @. ctmp[1:CNX,CNZ+1:end] = 0
    # transform to physcial space
    ldiv!(tmp1,SE.transforms.fourier,ctmp) 
    ldiv!(tmp2,SE.transforms.sine,tmp1)
    # start constructing ζtt 
    @. ζtt -= bg.Bz*tmp2

    # transform back to spectral space 
    mul!(tmp1,SE.transforms.sine,ζtt)
    mul!(ctmp,SE.transforms.fourier,tmp1)       

    # invert the Laplacian to get the output
    @inbounds @. out[1:CNX,1:CNZ] = -1*ctmp[1:CNX,1:CNZ]/(kx^2 + kz^2)
    @inbounds @. out[CNX+1:end,:] = 0
    @inbounds @. out[1:CNX,CNZ+1:end] = 0

    return nothing
end

function (SE::SawyerEliassenRHS!)(inplace)
    (SE::SawyerEliassenRHS!)(inplace,inplace)
end

struct LHS!
    h :: Float64
    a :: Float64
    L! :: SawyerEliassenRHS!
end

function (lhs::LHS!)(out,in)
    lhs.L!(out,in)
    @. out = in - lhs.a*lhs.h^2*out
    return nothing
end

