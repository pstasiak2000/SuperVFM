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

const eps32 = eps(0.0f0) #Machine epsilon of 0.0 in 32 bits
const ZeroVector = SVector{3,Float32}(0.0f0,0.0f0,0.0f0) #Zero vector

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

    u_loc = CUDA.fill(SVector{3,Float32}(0,0,0),pcount) #Local superfluid velocity
    u_sup = CUDA.fill(SVector{3,Float32}(0,0,0),pcount) #Total superfluid velocity

    u = CUDA.fill(SVector{3,Float32}(0,0,0),pcount) #Vortex velocity
    u1 = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)#Previous vortex velocity (AB2)
    u2 = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)#Previous vortex velocity (AB2)
    
    #Determines the empty particles based on f_infront and stores them as a CUDA Boolean array. 
    Empty = vec(CUDA.reduce(+, fint, dims=1) .== 0)

    #The ghost particles infront and behind
    ghosti = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)
    ghostb = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)


    t = 0.0 #Simulation time
    print_info(f, SimParams, pcount, 0)

    # x_pos = zeros(Float32, 2, SimParams.nsteps)
    for it ∈ 1:SimParams.nsteps

        #Find the right number of threads and blocks if pcount changes
        nthreads, nblocks = redefineThreads_Blocks(pcount)

        ghostp!(ghosti, ghostb, f, fint, pcount, SimParams.box_size; nthreads, nblocks)
        
        #Compute the superfluid velocity v_s
        @. u_loc = (SimParams.velocity)(f, ghosti, ghostb, Empty, SimParams.κ, SimParams.corea)
        u_sup .= u_loc

        #Calculate the velocities of the filaments
        @. u = calc_fil_motion(f, ghosti, ghostb, u_sup, Empty)
        
        #Timestep the new positions
        @. f += SimParams.dt * timestep(u, u1, u2, Empty)

        t += SimParams.dt
        # x_pos[1,it] = t
        # x_pos[2,it] = Array(f)[1,1][1]
        print_info(f, SimParams, pcount, it)
    end
    return f#, x_pos
end


end # module VortexFilament
