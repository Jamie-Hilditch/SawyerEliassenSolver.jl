# # [Near-inertial waves](@id near-inertial_waves_example)
#
# This simple example demonstrates how to set up and run a simulation with `SawyerEliassenSolver`.
#
# ## Setting up a problem
# First load the SawyerEliassenSolver objects and functions

using SawyerEliassenSolver

# ### Grid
# The first step is to create a grid. The grids are periodic in ``x`` and bounded in ``z``.
# Let's simulate some near-inertial waves in a domain of depth 1km and width 100km.

NX, NZ = 32, 16
const LX, LZ = 100, 1
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
#   v & = -(f / \omega) a m \sin(k x - \omega t + \phi) \cos m z \\
#   w & = -k a \sin(k x - \omega t + \phi) \sin m z \\
#   b & = -N^2 / \omega k a \cos(k x - \omega t + \phi) \sin m z
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
    return -N² / ω(W.k, W.m) *
           W.k *
           W.a *
           cos(W.k * x - ω(W.k, W.m) * t + W.φ) *
           sin(W.m * z)
end

W₁ = (; k=2 * π / LX, m=π / LZ, a=1, φ=0)
W₂ = (; k=-12 * π / LX, m=5 * π / LZ, a=0.1, φ=π / 2);

# A couple of utility functions allow us to set initial conditions from the velocity and buoyancy.
# First we specify the initial conditions as functions of `x` and `z`.

u₀(x, z) = u(x, z, 0, W₁) + u(x, z, 0, W₂)
v₀(x, z) = v(x, z, 0, W₁) + v(x, z, 0, W₂)
w₀(x, z) = w(x, z, 0, W₁) + w(x, z, 0, W₂)
b₀(x, z) = b(x, z, 0, W₁) + b(x, z, 0, W₂);

# Then we use `set_ζ!` and `set_vb!` to set the initial conditions.

set_ζ!(problem; u=u₀, w=w₀)
set_vb!(problem; v=v₀, b=b₀)

# Note that `set_vb!` also sets the state variable `ζₜ` using ``\zeta_t = b_x - fv_z``.
# The functions `set_v!` and `set_b!` can be used to set `v` and `b` separately without computing `ζₜ`.
# Furthermore, with `v` and `b` specified the function `compute_ζₜ!` can be used to compute and set `ζₜ`.

# ## Timestepper
# The next step is to create a timestepper. We use a timestep of `Δt = 2πf⁻¹/100`.

Δt = 2π / f / 100
ts = Timestepper(problem,Δt)

# ## Running the simulation
# We can now run the simulation for 1000 timesteps.

advance!(ts, 1000)
