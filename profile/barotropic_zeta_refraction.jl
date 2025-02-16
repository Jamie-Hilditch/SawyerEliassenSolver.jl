# More or less the same as the barotropic zeta refraction example

using Printf
using Profile
using PProf
using SawyerEliassenSolver

using SpecialFunctions: erf

const ALLOCS_FILEPATH = joinpath(@__DIR__, "barotropic_zeta_refraction.allocs.pb.gz")

NX, NZ = 128, 512
const LX, LZ = 2, 0.1
grid = Grid(NX, NZ, (-LX / 2, LX / 2), LZ)

const N²::Float64 = 200^2
const Ro = 0.5
const L = 0.2
Vx(x) = Ro * (-1 + 2 / L * exp(-π * x^2 / L^2));

xgrid, zgrid = gridpoints(grid)
background_flow = BackgroundFlow(grid)
background_flow.Vx .= Vx.(xgrid);
background_flow.Bz .= N²;

const dealias = 10
domain = Domain(grid; dealias_x=dealias)

const H = 2e-3
u₀(x, z) = 0.5 * (1 + erf(3 + 2 * z / H));

const ω₀² = 1.0
const ω₁² = N²

function setup_problem()
    grid = Grid(NX, NZ, (-LX / 2, LX / 2), LZ)
    domain = Domain(grid; dealias_x=dealias)
    xgrid, zgrid = gridpoints(grid)
    background_flow = BackgroundFlow(grid)
    background_flow.Vx .= Vx.(xgrid)
    background_flow.Bz .= N²
    problem = Problem(domain, background_flow)
    set_ζ!(problem; u=u₀)
    return domain, problem
end;

function setup_timestepper(use_preconditioner::Bool=false)
    domain, problem = setup_problem()
    if use_preconditioner
        preconditioner = DiagonalQuadraticPreconditioner(domain, ω₀², ω₁²)
        ts = Timestepper(problem, 2π / 50, preconditioner)
    else
        ts = Timestepper(problem, 2π / 50)
    end
    return ts
end;

function run_simulation(ts)
    advance!(ts, 10)
    return nothing
end;

# ### With preconditioning
# Create the timestep and advance 1 step before timing.

ts = setup_timestepper(true)
advance!(ts)

Profile.clear()
Profile.clear_malloc_data()
Profile.Allocs.@profile run_simulation(ts)
PProf.Allocs.pprof(; web=false, out=ALLOCS_FILEPATH)
