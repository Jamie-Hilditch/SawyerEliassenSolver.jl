# # [Forced secondary circulation](@id forced_secondary_circulation_example)
#
# This example demonstrates how to apply RHS forcing and how to create output variables.
#
# With the flow initially at rest we will apply forcing to the RHS of the vorticity equation and generate a secondary circulation.
# THe forcing function will take the form
# ```math
# \mathcal{F}(x,z,t) = \mathcal{S}(x,z)\frac{1}{\tau}\sigma\left(\frac{t}{\tau}\right)
# ```
# where the spatial structure is given by a Gaussian
# ```math
# \mathcal{S}(x,z) = \exp\left(-\frac{x^2}{2L^2} + \frac{z^2}{2H^2}\right)
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
# `zeta` forcing can be applied in four different ways:
# - [`PointwisePhysicalForcing`](@ref)
# - [`GlobalPhysicalForcing`](@ref)
# - [`PointwiseSpectralForcing`](@ref)
# - [`GlobalSpectralForcing`](@ref)
# Here we implement Pointwise Physical Forcing where we specify the forcing as a function of `x`, `z` and `t`.

const L::Float64 = 10
const H::Float64 = 0.05

function forcing(x, z, t)
    return exp(-x^2 / (2 * L^2) + z^2 / (2 * H^2)) * t * exp(-t^2)
end
