"""$(TYPEDSIGNATURES)

Integrate periodic gradients.
"""
function integrate_periodic_gradients(
    x_derivative::XZVariable, z_derivative::XZVariable; out::XZVariable=XZVariable(x_derivative), FZ_working_space::FZVariable=FZVariable(x_derivative)
)
    @boundscheck consistent_domains(x_derivative, z_derivative, out)
    domain = get_domain(x_derivative)
    x_derivative_FZ = FZ_working_space

    # first get horizontal buoyancy gradient in FZ space
    horizontal_transform!(x_derivative_FZ,x_derivative)

    # save the mean component and integrate the rest (note the integration sets the mean to 0)
    x_derivative_mean = x_derivative_FZ[1, :] ./ size(x_derivative, 1)
    ∫dx!(x_derivative_FZ)

    # transform back to physical space and add on linear part
    horizontal_transform!(out, x_derivative_FZ)
    out .+= xgridpoints(domain) * x_derivative_mean'

    # Now we need to add on the mean z dependence
    z_derivative_mean = mean(z_derivative; dims=1)

    # integrate z_derivative_mean
    # use trapezoid rule but don't subtract half the first value
    # this is equivalent to setting the horizontal mean component to zero on the bottom boundary
    mean = (cumsum(z_derivative_mean; dims=2) .- z_derivative_mean ./ 2) .* zstepsize(domain)
    out .+= mean

    return out
end
