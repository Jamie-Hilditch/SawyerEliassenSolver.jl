using Printf: @sprintf
using Profile
using PProf

using SawyerEliassenSolver

const NX = 256
const NZ = 128
const LX = 1000.0
const LZ = 1.0
const N = 100.0

# output
const PROFILING_FILEPATH = joinpath(@__DIR__, "barotropic_niws.profile.pb.gz")
const ALLOCS_FILEPATH = joinpath(@__DIR__, "barotropic_niws.allocs.pb.gz")

function setup_simulation()
    # create domain
    grid = Grid(NX, NZ, LX, LZ)
    domain = Domain(grid)

    # set background flow
    bg = BackgroundFlow(grid)
    bg.Bz .= N^2

    # setup the simulation
    problem = Problem(domain, bg)

    function initial_conditions(x, z)
        sin(2 * π * x / LX) * cos(5 * π * z / LZ) +
        1 / 2 * sin(5 * 2 * π * x / LX) * cos(π * z / LZ)
    end

    set_ζ!(problem, u=initial_conditions)

    # timestepper
    preconditioner = DiagonalQuadraticPreconditioner(domain, 1.0, N^2)
    timestepper = Timestepper(problem, 2π / 20, preconditioner)

    return timestepper
end


# profile backtraces
ts = setup_simulation()
advance!(ts)

Profile.init(; n=10^8)
Profile.clear()
Profile.clear_malloc_data()
@profile advance!(ts, 100)
pprof(; web=false, out=PROFILING_FILEPATH)

# profile memory allocations
ts = setup_simulation()
advance!(ts)

Profile.clear()
Profile.clear_malloc_data()
Profile.Allocs.@profile advance!(ts, 100)
PProf.Allocs.pprof(; web=false, out=ALLOCS_FILEPATH)
