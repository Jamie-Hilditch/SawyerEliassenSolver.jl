using Printf: @sprintf
using Profile
using PProf

using SawyerEliassenSolver

const NX = 384
const NZ = 192
const LX = 1000.
const LZ = 1.
const Π = 100.

# timestepping
const NPERIODS = 5

# output 
const OUTPUT_FILENAME = Base.Filesystem.joinpath(@__DIR__,"barotropic_niws.h5")
const PROFILING_FILEPATH = Base.Filesystem.joinpath(@__DIR__,"barotropic_niws.profile.pb.gz")
const ALLOCS_FILEPATH = Base.Filesystem.joinpath(@__DIR__,"barotropic_niws.allocs.pb.gz")

initial_conditions(x,z) = @. sin(2*π*x/LX)*cos(5*π*z'/LZ) + 1/2*sin(5*2*π*x/LX)*cos(π*z'/LZ)

group_name(Δt) = @sprintf "%.3f" Δt/π

function run_simulation(Δt :: Float64, OUTPUT_CADENCE :: Int)
    grid = Grid(NX,NZ,LX,LZ)

    # set background flow
    bg = BackgroundFlow(grid)
    @. bg.Bz = Π^2
    
    # setup the simulation
    problem = Problem(grid,bg)
    timestepper, state = setup_simulation(problem,Δt)

    # set initial conditions
    Tools.ics_from_u!(state,problem,initial_conditions(grid.x,grid.z))

    # construct a file handler and save initial conditions
    output_writer = Tools.OutputWriter(OUTPUT_FILENAME,problem,state,
        :u,:w,:v,group=group_name(Δt),overwrite=true)
    output_writer()

    # run simulation
    @info @sprintf "Beginning Simulation with Δt = %.3fπ" Δt/π
    while state.clock.t < 2*π*NPERIODS
        timestepper(OUTPUT_CADENCE)
        output_writer()
        @info "Completed iteration $(state.clock.iteration)"
    end 

    # close the output file
    finalize(output_writer)
end

# const Δt = [π/10,π/20,π/40]
# const OUTPUT_CADENCE = [5,10,20]

# run first simulation. Now everything is compiled
run_simulation(π/10,5)

run_simulation(π/14,7)

# run with higher time resolution and profile backtraces
Profile.init(n = 10^8)
Profile.clear()
Profile.clear_malloc_data()
@profile run_simulation(π/20,10)
pprof(web=false,out=PROFILING_FILEPATH)

run_simulation(π/30,15)

# run with even higher resolution and profile memory allocations
Profile.clear()
Profile.clear_malloc_data()
Profile.Allocs.@profile run_simulation(π/40,20)
PProf.Allocs.pprof(web=false,out=ALLOCS_FILEPATH)

# run final simutional with yet higher resolution providing .mem files
Profile.clear()
Profile.clear_malloc_data()
run_simulation(π/80,40)