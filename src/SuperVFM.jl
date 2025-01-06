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
    print_banner()
    print_GPU_info()
    print_boundary_info(
        SimParams.boundary_x,
        SimParams.boundary_y,
        SimParams.boundary_z)
    
    print_filamentmodel_info(SimParams.FilamentModel)
    #Check the timestep here
    @assert check_timestep(SimParams) "Timestep is too large dt=$(SimParams.dt)"

    #Initialise the vortex arrays [VF_initial_condition.jl]
    f, fint, pcount, nthreads, nblocks = (SimParams.initf)(SimParams.δ)
    f_curv = CUDA.zeros(Float32, pcount) #Vortex filament curvature

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

    normal_velocity = CUDA.fill(SimParams.normal_velocity, pcount)

    t = 0.0 #Simulation time
    print_info_header()

    f_out = []
    itCount = 0
    for it ∈ 1:SimParams.nsteps

        #Find the right number of threads and blocks if pcount changes
        nthreads, nblocks = redefineThreads_Blocks(pcount)

        ghostp!(ghosti, ghostb, f, fint, pcount, SimParams.box_size; nthreads, nblocks)
        
        #Compute the superfluid velocity v_s
        @. u_loc = (SimParams.velocity)(f, ghosti, ghostb, Empty, SimParams.κ, SimParams.corea)
        u_sup .= u_loc

        #Calculate the velocities of the filaments
        @. u = (SimParams.FilamentModel)(f, ghosti, ghostb, u_sup, normal_velocity, Empty)
        #Timestep the new positions
        @. f += SimParams.dt * timestep(u, u1, u2, Empty)
        u2 .= u1
        u1 .= u
        t += SimParams.dt

        fCPU = Array(f)
        # x_pos[1,it] = t
        # x_pos[2,it] = fCPU[1,1][1]
        # x_pos[3,it] = fCPU[1,1][2]
        # x_pos[4,it] = fCPU[1,1][3]
        print_info(f, ghosti, ghostb, u, Empty, SimParams, pcount, it)
        if mod(it, SimParams.shots) == 0
            itCount += 1
            push!(f_out,fCPU)
        end
    end
    println("")
    printstyled("Simulation finished!\n", bold=:true, color=:green)
    return f_out
end


end # module VortexFilament
