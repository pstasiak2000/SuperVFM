module SuperVFM
using CUDA, Adapt
using Dates
using StaticArrays
using LinearAlgebra
using Plots

export Run

#Define the max number of threads that can be run per block
const max_threads_per_block = CUDA.attribute(CUDA.device(), CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)

include("VF_cdata.jl")
include("VF_boundary.jl")
include("VF_timestep.jl")
include("VF_initial_condition.jl")
include("VF_derivatives.jl")
include("VF_general.jl")
include("VF_misc.jl")




function Run(SimParams::SimulationParams)
    # print_banner()
    # print_GPU_info()
    # print_boundary_info(
    #     SimParams.boundary_x,
    #     SimParams.boundary_y,
    #     SimParams.boundary_z)
    

    #!!![NOTE] Function here to determine the number of threads and blocks to be used
    nthreads = min(pcount, 1024)
    nblocks = cld(pcount, nthreads)
    #!!!
    #Check the timestep here
    @assert check_timestep(SimParams) "Timestep is too large dt=$(SimParams.dt)"

    #Initialise the vortex arrays [VF_initial_condition.jl]
    f, fint, pcount, nthreads, nblocks = (SimParams.initf)()
    
    @info "Computing the ghost points" #[VF_boundary.jl]
    @time ghostp!(f, fint, pcount, SimParams; nthreads, nblocks) 

    t = 0 #Simulation time
 
    x_pos = zeros(2,SimParams.nsteps)
    for it ∈ 1:SimParams.nsteps

        calc_fil_motion!(f, fint, pcount, SimParams::SimulationParams; nthreads=1, nblocks=1)
        fCPU = Array(f) 
        t += SimParams.dt

        x_pos[1,it] = t
        x_pos[2,it] = fCPU[1,1][1]
    end
    return f, x_pos
end


end # module VortexFilament
