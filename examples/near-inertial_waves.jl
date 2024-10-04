# # [Near-inertial waves](@id near-inertial_waves_example)
#
# This simple example demonstrates how to set up and run a simulation with `SawyerEliassenSolver`.
#
# ## Setting up a problem
# First load the SawyerEliassenSolver and the other packages we need.

using CairoMakie
using HDF5
using Printf
using Statistics

using SawyerEliassenSolver

# ### Grid
# The first step is to create a grid. The grids are periodic in ``x`` and bounded in ``z``.
# Let's simulate some near-inertial waves in a domain of depth 1 km and width 250 km.

NX, NZ = 128, 64
const LX, LZ = 250, 1
grid = Grid(NX, NZ, (-LX / 2, LX / 2), LZ)

# Some remarks on the grid:
#  * Here we use the default floating point precision of `Float64`.
#  * Alternatively, we could use `Float32` by creating the grid with `grid = Grid(Float32, NX, NZ, (-LX/2,LX/2), LZ)`.
#  * The choice of precision here propagates through the entire simulation.
#  * We specified the ``x`` extent with a tuple of bounds but the ``z`` extent with a single depth value. Passing a the width / depth is equivalent to passing a tuple of `(0, LX)` / `(-LZ, 0)` for ``x`` and ``z`` respectively.
#
# ### Background Flow
# With a grid we can now create a background flow.
# When creating the background flow let's set the Coriolis frequency to a typical midlatitude value of ``10^{-4}\mathrm{s}^{-1}``.

const f = 1e-4
background_flow = BackgroundFlow(grid, f)

# The background flow is initialised to zeros. For our simple example, let's use a constant background stratification of ``10^{-4}\mathrm{s}^{-2}``

const N² = 1e-4
background_flow.Bz .= N²;

# ### Domain
# Now we create a `Domain`. The `Domain` object contains not just the grid but also details of the spectral domain and the transforms between the different representations.
# When creating a `Domain` we can specify the number of high wavenumbers to zero-out inorder to dealias the products with the background flow. In this example the background
# flow is constant and so we do not need to worry about dealiasing.

domain = Domain(grid)

# ### Problem
# The next step is to define a `Problem`. A problem consists of a domain, background flow and forcing. We leave forcing for a later example and thus create the problem with

problem = Problem(domain, background_flow)

# The problem object also contains the current state of the simulation.

# ## Initial conditions
#
# Here we simulate a pair of waves propagating horizontally in different directions.
# If the streamfunction of a wave is ``\psi = a \cos(k x - \omega t + \phi) \sin m z`` then the polarisation relations give
# ```math
# \begin{align*}
#   u & = -m a \cos(k x - \omega t + \phi) \cos m z \\
#   v & = -(f / \omega) m a \sin(k x - \omega t + \phi) \cos m z \\
#   w & = -k a \sin(k x - \omega t + \phi) \sin m z \\
#   b & = (N^2 / \omega) k a \cos(k x - \omega t + \phi) \sin m z
# \end{align*}
# ```

ω(k, m) = √(f^2 * m^2 + N² * k^2) / √(k^2 + m^2)
function u(x, z, t, W)
    return -W.m * W.a * cos(W.k * x - ω(W.k, W.m) * t + W.φ) * cos(W.m * z)
end
function v(x, z, t, W)
    return -(f / ω(W.k, W.m)) *
           W.m *
           W.a *
           sin(W.k * x - ω(W.k, W.m) * t + W.φ) *
           cos(W.m * z)
end
function w(x, z, t, W)
    return -W.k * W.a * sin(W.k * x - ω(W.k, W.m) * t + W.φ) * sin(W.m * z)
end
function b(x, z, t, W)
    return N² / ω(W.k, W.m) *
           W.k *
           W.a *
           cos(W.k * x - ω(W.k, W.m) * t + W.φ) *
           sin(W.m * z)
end

W₁ = (; k=2 * π / LX, m=π / LZ, a=0.1, φ=0)
W₂ = (; k=-12 * π / LX, m=4 * π / LZ, a=0.01, φ=π / 2);
ω(W₁.k, W₁.m) / f, ω(W₂.k, W₂.m) / f

# A couple of utility functions allow us to set initial conditions from the velocity and buoyancy.
# First we specify the initial conditions as functions of `x` and `z`.

u₀(x, z) = u(x, z, 0, W₁) + u(x, z, 0, W₂)
v₀(x, z) = v(x, z, 0, W₁) + v(x, z, 0, W₂)
w₀(x, z) = w(x, z, 0, W₁) + w(x, z, 0, W₂)
b₀(x, z) = b(x, z, 0, W₁) + b(x, z, 0, W₂);

# Then we use [`set_ζ!`](@ref) and [`set_vb!`](@ref) to set the initial conditions.

set_ζ!(problem; u=u₀, w=w₀)
set_vb!(problem; v=v₀, b=b₀)

# Note that `set_vb!` also sets the state variable `ζₜ` using ``\zeta_t = b_x - fv_z``.
# The functions `set_v!` and `set_b!` can be used to set `v` and `b` separately without computing `ζₜ`.
# Furthermore, with `v` and `b` specified the function [`compute_ζₜ!`](@ref) can be used to compute and set `ζₜ`.

# ## Timestepper
# The next step is to create a timestepper. We use a timestep of `Δt = 2πf⁻¹/100`.

Δt = 2π / f / 100
ts = Timestepper(problem, Δt)

# ## Output
# We can write output to an HDF5 file using an `OutputWriter`.

output_writer = OutputWriter(problem, "near_inertial_waves.h5"; overwrite=true)

# Here we'll just save the `u` component of the velocity.
# Common output variables are defined in the submodule [`OutputVariables`](@ref)

u_output = OutputVariables.u(problem)
add_output_variables!(output_writer; u=u_output)

# We can also write some attributes to the output file.

write_attributes!(output_writer; f=f, N²=N²)

# To save the current state of the simulation to the output file we use the [`write!`](@ref) function.
# i.e. to save the initial conditions

write!(output_writer)

# ## Running the simulation
# Finally, it is time to run the simulation. We use [`advance!`](@ref) to timestep the simulation.
# Let's run the simulation for 1000 timesteps saving output every 10 timesteps.

for _ in 1:100
    advance!(ts, 10)
    write!(output_writer)
end

# ## Visualise the solution
# Read in the output we saved.

output = h5open("near_inertial_waves.h5", "r") do h5
    (
        u=read_dataset(h5, "u"),
        time=read_dataset(h5, "time"),
        x=read_dataset(h5, "x"),
        z=read_dataset(h5, "z"),
    )
end;

# Make a video of `u`.

n = Observable(1)
uₙ = @lift output[:u][:, :, $n]
title = @lift @sprintf "t = %.2f inertial periods" output[:time][$n] * f / 2π

fig = Figure(; size=(1200, 400))
Label(fig[1, 1:2], title; tellwidth=false)
ax = Axis(fig[2, 1]; xlabel="x [km]", ylabel="z [km]")

cf = heatmap!(ax, output[:x], output[:z], uₙ; colormap=:balance, colorrange=(-0.5, 0.5))
Colorbar(fig[2, 2], cf; label="u [ms⁻¹]", labelpadding=10)

record(fig, "near-inertial_waves.mp4", 1:length(output[:time]); framerate=10) do i
    n[] = i
end
nothing #hide

# ![](near-inertial_waves.mp4)
# We have a mode 1 wave propagating to the right and a mode 4 wave propagating to the left.

# ## Error analysis
# We can compute the root mean square error between the theoretical solution and the simulated
# solution after 1 inertial period

x_grid, z_grid = gridpoints(domain)
u_theory = [
    u(x, z, 2π / f, W₁) + u(x, z, 2π / f, W₂)
    for x in xgridpoints(domain), z in zgridpoints(domain)
];

function run_one_inertial_period(problem, nsteps::Int)
    set_ζ!(problem; u=u₀, w=w₀)
    set_vb!(problem; v=v₀, b=b₀)
    Δt = 2π / f / nsteps
    ts = Timestepper(problem, Δt)
    for _ in 1:nsteps
        advance!(ts)
    end
    output = OutputVariables.u(problem)
    compute!(output)
    return output
end

function rms_error(simulated, theory)
    squared_error = (simulated .- theory) .^ 2
    return sqrt(mean(squared_error))
end;

# Make a plot of the root mean square error as a function of the timestep.
# Our time stepping scheme is third-order accurate and we see that the error scales as `Δt³`.
# The lower bound for the error is determined by the tolerance of the iterative solver.

fig = Figure(; size=(1200, 400))
ax = Axis(fig[1, 1]; xlabel="Δt [Inertial periods]", ylabel="RMS Error", xscale=log10, yscale=log10)
number_of_steps = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]

for nsteps in number_of_steps
    problem = Problem(domain, background_flow)
    u_sim = run_one_inertial_period(problem, nsteps)
    rms = rms_error(u_sim, u_theory)
    scatter!(ax, 1 / nsteps, rms, color=:red, marker=:cross)
end

lines!(ax, 1 ./ number_of_steps, 1e3 ./ number_of_steps .^ 3; color=:black, label="Δt³", linewidth=2)
axislegend(position=:lt)

current_figure() #hide
fig
