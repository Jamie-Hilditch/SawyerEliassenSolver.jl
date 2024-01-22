function ics_from_u!(state :: SES.State, problem :: SES.Problem, u :: Union{Float64,Vector{Float64},Matrix{Float64}})
    # copy u into the temporary storage broadcasting as necessary
    rtmp1 = problem.scratch.rtmp1
    rtmp2 = problem.scratch.rtmp2
    ctmp = problem.scratch.ctmp 
    CNX = problem.spectral_domain.CNX
    CNZ = problem.spectral_domain.CNZ
    kz = problem.spectral_domain.kz
    
    rtmp1 .= u
    mul!(rtmp2,problem.transforms.cosine,rtmp1)
    mul!(ctmp,problem.transforms.fourier,rtmp2)
    @inbounds @. state.ψ[1:CNX,1:CNZ] = -ctmp[1:CNX,2:CNZ+1]/kz
    @inbounds @. state.ψ[CNX+1:end,:] = 0 
    @inbounds @. state.ψ[1:CNX,CNZ + 1:end] = 0
    return nothing 
end

"""
Compute ψt in spectral space from b and v in physical space
"""
function set_ψt!(state :: SES.State, problem :: SES.Problem)
        
    CNX = problem.spectral_domain.CNX
    CNZ = problem.spectral_domain.CNZ
    kx = problem.spectral_domain.kx
    kz = problem.spectral_domain.kz
    rtmp2 = problem.scratch.rtmp2
    ctmp = problem.scratch.ctmp 
    ψt = state.ψt

    # compute bx and put in ψt
    mul!(rtmp2,problem.transforms.sine,state.b)
    mul!(ctmp,problem.transforms.fourier,rtmp2)
    @inbounds @. ψt[1:CNX,1:CNZ] = 1*im*kx*ctmp[1:CNX,1:CNZ]

    # compute vz0
    vz = problem.scratch.rtmp1
    dz = problem.grid.z[2] - problem.grid.z[1]
    _differentiate_z!(vz,state.v,dz)

    # thermal wind imbalance 
    f = problem.background_flow.f
    mul!(rtmp2,problem.transforms.sine,vz)
    mul!(ctmp,problem.transforms.fourier,rtmp2)
    @inbounds @. ψt[1:CNX,1:CNZ] -= f*ctmp[1:CNX,1:CNZ]

    # invert laplacian
    @inbounds @. ψt[1:CNX,1:CNZ] *= -1/(kx^2 + kz^2)
    @inbounds @. ψt[CNX+1:end,:] = 0
    @inbounds @. ψt[1:CNX,CNZ+1:end] = 0

    return nothing

end


