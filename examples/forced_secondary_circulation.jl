# # [Forced secondary circulation](@id forced_secondary_circulation_example)
#
# This example demonstrates how to apply RHS forcing and how to create output variables.
#
# With the flow initially at rest we will apply forcing to the RHS of the vorticity equation and generate a secondary circulation.
# The forcing function will take the form
# ```math
# \mathcal{F}(x,z,t) = \mathcal{S}\left(\frac{x}{L},\frac{z}{H}\right)\sigma\left(\frac{t}{\tau}\right)
# ```
# where ``L``, ``H`` and ``\tau`` are horizontal, vertical and temporal scales respectively.
# The spatial structure is given by a Gaussian,
# ```math
# \mathcal{S}(x',z') = \exp\left(-\frac{x'^2}{2} + \frac{z'^2}{2}\right),
# ```
# and the temporal part is given by (minus) the derivative of a Gaussian,
# ```math
# \sigma(t') = t'\exp(-t'^2 / 2).
# ```
#
# ### Load packages

using CairoMakie
using HDF5
using Printf

using SawyerEliassenSolver

# ## Grid, background flow and domain
#
# We non-dimensionalise by the half-height of the domain and `f`.
# The background flow has uniform gradients and thus we require no dealiasing.

const NX, NZ = 4096, 256
const LX::Float64 = 1000
const M²::Float64 = 100
const N²::Float64 = 1e5

grid = Grid(NX, NZ, (-LX / 2, LX / 2), (-1, 1))

background_flow = BackgroundFlow(grid)
background_flow.Bx .= M²
background_flow.Bz .= N²

domain = Domain(grid);

# ## Forcing
# ``\zeta`` forcing can be applied in four different ways:
# - [`PointwisePhysicalForcing`](@ref)
# - [`GlobalPhysicalForcing`](@ref)
# - [`PointwiseSpectralForcing`](@ref)
# - [`GlobalSpectralForcing`](@ref)
# Here we implement Pointwise Physical Forcing where we specify the forcing as a function of `x`, `z` and `t`.
# The function can also optionally accepted parameters as the final arguments.
# We specify the spatial and temporal parts of the RHS forcing as non-dimensional functions
# and then provide the length and time scales as parameters.

@inline forcing_spatial(x, z) = exp(-(x^2 + z^2) / 2)
@inline forcing_temporal(t) = t * exp(-t^2 / 2)
@inline forcing(x, z, t, p) = forcing_spatial(x / p.L, z / p.H) * forcing_temporal(t / p.τ)

# We now construct the RHS forcing for the domain by passing in the function and parameters.
# We choose a timescale of 3 inertial periods.

const parameters = (; L=20.0, H=0.05, τ=6π)
physical_forcing = PointwisePhysicalForcing(domain, forcing, parameters)

# ### Time dependence of the forcing

fig = fig = Figure(; size=(1200, 400))
Label(fig[1, :], "Forcing time dependence"; tellwidth=false)
ax = Axis(
    fig[2, :];
    xlabel=L"t \ [\text{inertial periods}]",
    ylabel=L"\sigma\left(\frac{t}{\tau}\right)",
)
time = 0.0:0.1:10.0
lines!(ax, time, forcing_temporal.(time / 3))
xlims!(ax, 0, 10)
ylims!(ax, 0, nothing)

current_figure() # hide
fig

# ## Setting up the problem
# Now we pass in the forcing as a keyword argument when constructing the problem.

problem = Problem(domain, background_flow; ζ_forcing=physical_forcing)

# We can also pass in [`PointwisePhysicalForcing`](@ref) or [`PointwiseSpectralForcing`](@ref)
# as keyword arguments `v_forcing` and `b_forcing` for `v` and `b` respectively.
#
# We don't need to set any initial conditions but we'll use a preconditioner.

preconditioner = DiagonalQuadraticPreconditioner(domain, 1.0, N²)
ts = Timestepper(problem, 2π / 50, preconditioner)

# ## Output
# Let's save `u` and `w` to visualise the overturning circulation that we're forcing.

output_writer = OutputWriter(problem, "forced_secondary_circulation.h5"; overwrite=true)
add_output_variables!(
    output_writer; u=OutputVariables.u(problem), w=OutputVariables.w(problem)
)
write_attributes!(
    output_writer; f=1.0, M²=M², N²=N², L=parameters.L, H=parameters.H, tau=parameters.τ
)
write!(output_writer)

# Run the simulation.

for ii in 1:100
    advance!(ts, 5)
    write!(output_writer)
end

# ## Movie

output = h5open("forced_secondary_circulation.h5", "r") do h5
    (
        u=read_dataset(h5, "u"),
        w=read_dataset(h5, "w"),
        time=read_dataset(h5, "time"),
        x=read_dataset(h5, "x"),
        z=read_dataset(h5, "z"),
    )
end;

n = Observable(1)
uₙ = @lift output[:u][:, :, $n]
wₙ = @lift output[:w][:, :, $n]
title = @lift @sprintf "t = %.2f inertial periods" output[:time][$n] / 2π

fig = Figure(; size=(1200, 400))
Label(fig[1, 1:2], title; tellwidth=false)
ax_u = Axis(fig[2, 1]; ylabel="z")
ax_w = Axis(fig[2, 2]; xlabel="x", ylabel="z")
linkaxes!(ax_u, ax_w)

cf_u = heatmap!(
    ax_u, output[:x], output[:z], uₙ; colormap=:balance, colorrange=(-0.01, 0.01)
)
Colorbar(fig[3, 1], cf_u; vertical=false, label=L"u", labelpadding=10)

cf_w = heatmap!(
    ax_w, output[:x], output[:z], wₙ; colormap=:balance, colorrange=(-1e-4, 1e-4)
)
Colorbar(fig[3, 2], cf_w; vertical=false, label=L"w", labelpadding=10)

record(fig, "forced_secondary_circulation.mp4", 1:length(output[:time]); framerate=10) do i
    n[] = i
end

# ![](forced_secondary_circulation.mp4)
