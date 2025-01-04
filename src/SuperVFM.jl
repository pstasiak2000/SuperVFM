module SuperVFM
using CUDA, Adapt
using Dates
using StaticArrays
using LinearAlgebra
using Plots
import Printf: @sprintf
export Run

#Define the max number of threads that can be run per block
const max_threads_per_block = CUDA.attribute(CUDA.device(), CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)

include("VF_cdata.jl")
include("VF_boundary.jl")
include("VF_timestep.jl")
include("VF_initial_condition.jl")
include("VF_derivatives.jl")
include("VF_general.jl")
include("VF_output.jl")
include("VF_misc.jl")




function Run(SimParams::SimulationParams)
    # print_banner()
    # print_GPU_info()
    # print_boundary_info(
    #     SimParams.boundary_x,
    #     SimParams.boundary_y,
    #     SimParams.boundary_z)
    

    #Check the timestep here
    @assert check_timestep(SimParams) "Timestep is too large dt=$(SimParams.dt)"

    #Initialise the vortex arrays [VF_initial_condition.jl]
    f, fint, pcount, nthreads, nblocks = (SimParams.initf)(SimParams.δ)

    u = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)
    u1 = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)
    u2 = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)
    
    t = 0.0 #Simulation time
    for it ∈ 1:SimParams.nsteps

        #Find the right number of threads and blocks
        nthreads, nblocks = redefineThreads_Blocks(pcount)

   
        #Calculate the velocity and timestep
        calc_fil_motion!(f, u, u1, u2, fint, pcount, SimParams::SimulationParams; nthreads=1, nblocks=1)
        # fCPU = Array(f)

        t += SimParams.dt
        print_info(f, SimParams, pcount, it)
    end
    return f
end


end # module VortexFilament
