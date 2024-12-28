module SuperVFM
using CUDA, Adapt
using Dates
using StaticArrays
using LinearAlgebra

export Run

include("VF_boundary.jl")
include("VF_timestep.jl")
include("VF_initial_condition.jl")
include("VF_cdata.jl")
include("VF_derivatives.jl")
include("VF_general.jl")
include("VF_misc.jl")




function Run(SimParams::SimulationParams)
    print_banner()
    print_boundary_info(
        SimParams.boundary_x,
        SimParams.boundary_y,
        SimParams.boundary_z)
    print_GPU_info()

    pcount = getInitPcount(SimParams.initf,SimParams.δ)
    println("-: pcount is now at $pcount")

    #!!![NOTE] Function here to determine the number of threads and blocks to be used
    nthreads = min(pcount, 1024)
    nblocks = cld(pcount,nthreads)
    #!!!

    f = CUDA.zeros(Float32, 39, pcount)
    fint = CUDA.zeros(Int32, 3, pcount)

    initialiseVortex!(f,fint,pcount,SimParams.initf; threads=nthreads, blocks=nblocks)

    @time ghostp!(f,fint,pcount, SimParams; nthreads, nblocks)
    @time s_dot = get_deriv_1(f, pcount; nthreads, nblocks)
    t = 0; #Simulation time
    for it ∈ 1:SimParams.nsteps

        t += SimParams.dt
    end
    return f
end


end # module VortexFilament
