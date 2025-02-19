"""$(TYPEDSIGNATURES)

Integrate periodic gradients.
"""
function integrate_periodic_gradients(
    Bx::XZVariable, Bz::XZVariable; out::XZVariable=XZVariable(Bx), FZ_working_space::FZVariable=FZVariable(Bx)
)
    @boundscheck consistent_domains(Bx, Bz, out)
    domain = get_domain(Bx)
    Bx_FZ = FZ_working_space

    # first get horizontal buoyancy gradient in FZ space
    horizontal_transform!(Bx_FZ,Bx)

    # save the mean component and integrate the rest (note the integration sets the mean to 0)
    Bx_mean = Bx_FZ[1, :] ./ size(Bx, 1)
    ∫dx!(Bx_FZ)

    # transform back to physical space and add on linear part
    horizontal_transform!(out, Bx_FZ)
    out .+= xgridpoints(domain) * Bx_mean'

    # Now we need to add on the mean z dependence
    Bz_mean = mean(Bz; dims=1)

    # integrate Bz_mean
    # use trapezoid rule but don't subtract half the first value
    # this is equivalent to setting the horizontal mean component to zero on the bottom boundary
    B_mean = (cumsum(Bz_mean; dims=2) .- Bz_mean ./ 2) .* zstepsize(domain)
    out .+= B_mean

    return out
end
