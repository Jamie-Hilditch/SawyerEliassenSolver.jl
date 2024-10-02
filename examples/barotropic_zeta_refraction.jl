# # [Barotropic zeta refraction](@id barotropic_zeta_refraction_example)
#
# This example demonstrates how to apply dealiasing, explores different variable representations
# and uses a preconditioner to speed up the iterative solver.
#
# ### Grid
# The first step is to load the packages and then create a grid. Here the grid is non-dimensionalised
# by the half-width of the domain.

using CairoMakie
using HDF5
using Printf

using SawyerEliassenSolver

using SpecialFunctions: erf

#-

NX, NZ = 128, 512
const LX, LZ = 2, 0.1
grid = Grid(NX, NZ, (-LX / 2, LX / 2), LZ)

# ### Background Flow
# The background flow is barotropic with constant stratification ``N^2 = 200^2f^2`` but spatially varying
# vorticity. We are modelling a cyclonic Gaussian filament in an otherwise anticyclonic flow.
# Non-dimensionalising by ``f`` the lateral shear is given by
# ```math
# V_x = Ro\left(-1 + 2L^{-1}\mathrm{e}^{-\pi x^2/L^2}\right)
# ```
# where ``Ro = 0.5`` is the Rossby number of the anticyclonic region and ``L = 0.2`` is the half-width
# of the Gaussian filament.

const N²::Float64 = 200^2
const Ro = 0.5
const L = 0.2
Vx(x) = Ro * (-1 + 2 / L * exp(-π * x^2 / L^2));

# We create a background flow with the default ``f = 1``.

xgrid, zgrid = gridpoints(grid)
background_flow = BackgroundFlow(grid)
background_flow.Vx .= Vx.(xgrid);
background_flow.Bz .= N²;

# ## Dealiasing
# The background flow is spatially varying and so it is a good idea to dealias by zeroing out
# the highest wavenumbers. This is done when creating the [`Domain`](@ref). Here we drop the
# highest 10 wavenumbers in the horizontal.

const dealias = 10
domain = Domain(grid; dealias_x=dealias)

# ### Variables
# Let's visualise the background flow in both physical and Fourier space. We can do this by creating an
# [`XZVariable`](@ref) and then transforming in the horizontal.

Vx_physical = XZVariable(domain)
Vx_physical .= background_flow.Vx
Vx_fourier = horizontal_transform(Vx_physical);
summary(Vx_fourier)

# Notice that `Vx_fourier` is an [`FZVariable`](@ref). Variables can be any of six types
# determined by their horizontal and vertical representations. The six types are
# * [`XZVariable`](@ref) - physical in horizontal and vertical
# * [`FZVariable`](@ref) - Fourier in horizontal, physical in vertical
# * [`XSVariable`](@ref) - physical in horizontal, sine in vertical
# * [`FSVariable`](@ref) - Fourier in horizontal, sine in vertical
# * [`XCVariable`](@ref) - physical in horizontal, cosine in vertical
# * [`FCVariable`](@ref) - Fourier in horizontal, cosine in vertical
# Furthermore, the functions [`horizontal_transform`](@ref), [`sine_transform`](@ref),
# [`cosine_transform`](@ref), and the corresponding inplace versions [`horizontal_transform!`](@ref),
# [`sine_transform!`](@ref), [`cosine_transform!`](@ref), can be used to transform between these representations.

# Let's plot the background flow as a function of ``x`` and its power spectrum as a function of horizontal wavenumber ``k_x``.
# We also indicate where the 10th wavenumber is.

kx = xwavenumbers_full(domain)

fig = Figure(; size=(1200, 600))
Label(fig[1, 1], "Gaussian vorticity filament"; tellwidth=false)
ax_physical = Axis(fig[2, 1]; xlabel=L"x", ylabel=L"V_x")
lines!(ax_physical, xgrid[:,1], Vx_physical[:,1]; color=:black)
xlims!(ax_physical, -1, 1)

ax_spectral = Axis(fig[3, 1]; xlabel=L"k_x", ylabel=L"\Phi_{V_x}")
scatterlines!(ax_spectral, kx, abs2.(Vx_fourier[:, 1]); color=:black)
vlines!(ax_spectral, kx[dealias]; color=:red)
xlims!(ax_spectral, 0, 200)

current_figure() # hide
fig

# We see that all the power of the background flow is contained in the first 10 wavenumbers
# and that our dealiasing will be sufficient.

# ## Initial Conditions
# We initialise a slab of across-front velocity
# ```math
# u_0 = \frac{1}{2}\left(1 + \mathrm{erf}\left(3 + \frac{2 z}{H}\right)\right)
# ```
# where ``H = 2\times10^{-3}``.

const H = 2e-3
u₀(x, z) = 0.5 * (1 + erf(3 + 2 * z / H));

#-

problem = Problem(domain, background_flow)
set_ζ!(problem; u=u₀)

# ## Preconditioning
# We can speed up the iterative solver by using a preconditioner. Here we use a [`DiagonalQuadraticPreconditioner`](@ref)
# This preconditioner approximates the Sawyer-Eliassen equation in spectral space by
# ```math
# \mathcal{L}\zeta \approx \left(\omega_0^2 k_z^2 + \omega_1^2 k_x^2\right) / (k_x^2 + k_z^2) \zeta
# ```
# Here we set ``\omega_0^2 = 1`` and ``\omega_1^2 = N²``.

ω₀² = 1.0
ω₁² = N²
preconditioner = DiagonalQuadraticPreconditioner(domain, ω₀², ω₁²)

# Then we create a timestepper with timestep ``\Delta t = 2\pi / f / 50``.

ts = Timestepper(problem, 2π / 50, preconditioner)

# ## Output
# Save ``u``, ``v``, ``w``, ``b`` and ``∂u/∂z`` every 5 timesteps.

output_writer = OutputWriter(problem, "barotropic_zeta_refraction.h5"; overwrite=true)
add_output_variables!(
    output_writer;
    u=OutputVariables.u(problem),
    v=OutputVariables.v(problem),
    w=OutputVariables.w(problem),
    b=OutputVariables.b(problem),
    ∂u∂z=OutputVariables.∂u∂z(problem),
)
write_attributes!(output_writer; f=1, N²=N², Ro=Ro, L=L, H=H)
write!(output_writer)

# Run the simulation for 10 inertial periods.

for _ in 1:10
    advance!(ts, 5)
    write!(output_writer)
end

# ## Visualise the solution
# We'll plot the across-filament vertical shear and the wave energy density
# ```math
# \mathcal{E} = \frac{1}{2}\left(u^2 + v^2 + w^2 + b^2/N^2\right)
# ```
# in the top 25% of the domain.

z_idx = (3 * NZ ÷ 4):NZ

output = h5open("barotropic_zeta_refraction.h5", "r") do h5
    u = h5["u"][:, z_idx, :]
    v = h5["v"][:, z_idx, :]
    w = h5["w"][:, z_idx, :]
    b = h5["b"][:, z_idx, :]
    E = (u .^ 2 .+ v .^ 2 .+ w .^ 2 .+ b .^ 2 / N²) / 2
    (
        E=E,
        ∂u∂z=h5["∂u∂z"][:, z_idx, :],
        time=read_dataset(h5, "time"),
        x=read_dataset(h5, "x"),
        z=h5["z"][z_idx],
    )
end;

#-

n = Observable(1)
∂u∂zₙ = @lift output[:∂u∂z][:, :, $n]
Eₙ = @lift output[:E][:, :, $n]
title = @lift @sprintf "t = %.2f inertial periods" output[:time][$n] / 2π

fig = Figure(; size=(1200, 400))
Label(fig[1, 1:2], title; tellwidth=false)
ax_∂u∂z = Axis(fig[2, 1]; ylabel="z")
ax_E = Axis(fig[2, 2]; xlabel="x", ylabel="z")
linkaxes!(ax_∂u∂z, ax_E)

cf_∂u∂z = heatmap!(
    ax_∂u∂z, output[:x], output[:z], ∂u∂zₙ; colormap=:balance, colorrange=(-1000, 1000)
)
Colorbar(
    fig[3, 1], cf_∂u∂z; vertical=false, label=L"\partial u/\partial z", labelpadding=10
)

cf_E = heatmap!(
    ax_E, output[:x], output[:z], Eₙ; colormap=Reverse(:grays), colorrange=(0, 0.5)
)
Colorbar(fig[3, 2], cf_E; vertical=false, label=L"\mathcal{E}", labelpadding=10)

record(fig, "barotropic_zeta_refraction.mp4", 1:length(output[:time]); framerate=10) do i
    n[] = i
end
nothing #hide

![](barotropic_zeta_refraction.mp4)

# ## Performance with and without preconditioning
# Let's time how long it takes to run 10 inertial periods with and without preconditioning.
# First create a clean problem.

function setup_problem()
    grid = Grid(NX, NZ, (-LX / 2, LX / 2), LZ)
    domain = Domain(grid; dealias_x=dealias)
    xgrid, zgrid = gridpoints(grid)
    background_flow = BackgroundFlow(grid)
    background_flow.Vx .= Vx.(xgrid);
    background_flow.Bz .= N²;
    problem = Problem(domain, background_flow)
    set_ζ!(problem; u=u₀)
    return domain, problem
end;

# ### Without preconditioning
# Create the timestep and advance 1 step before timing.

domain, problem = setup_problem()
ts = Timestepper(problem, 2π / 50)
advance!(ts)

#-

@time advance!(ts, 500)

# ### With preconditioning

domain, problem = setup_problem()
preconditioner = DiagonalQuadraticPreconditioner(domain, ω₀², ω₁²)
ts = Timestepper(problem, 2π / 50, preconditioner)
advance!(ts)

#-

@time advance!(ts, 500)
