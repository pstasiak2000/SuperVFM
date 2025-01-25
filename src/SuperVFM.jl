module SuperVFM
using CUDA, Adapt
using Dates
using StaticArrays
using LinearAlgebra
using Interpolations
using Plots
using Unitful
using DelimitedFiles
import Printf: @sprintf

export Run

#Define the max number of threads that can be run per block
# const max_threads_per_block = CUDA.attribute(CUDA.device(), CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)

const eps32 = eps(0.0f0) #Machine epsilon of 0.0 in 32 bits
const ZeroVector = SVector{3,Float32}(0.0f0,0.0f0,0.0f0) #Zero vector

#Unit vectors to avoid constant definitions in the code
const e_x = cu(SVector{3,Float32}(1.0f0,0.0f0,0.0f0)) #unit vector in the x direction
const e_y = cu(SVector{3,Float32}(0.0f0,1.0f0,0.0f0)) #unit vector in the y direction
const e_z = cu(SVector{3,Float32}(0.0f0,0.0f0,1.0f0)) #unit vector in the z direction

struct cpu end; export cpu 
struct gpu end; export gpu


include("VF_cdata.jl")
include("VF_boundary.jl")
include("VF_timestep.jl")
include("VF_initial_condition.jl")
include("VF_derivatives.jl")
include("VF_general.jl")
include("VF_output.jl")
include("VF_misc.jl")




function Run(::cpu,SimParams::SimulationParams)
    print_banner()
    print_boundary_info(
        SimParams.boundary_x,
        SimParams.boundary_y,
        SimParams.boundary_z)
    print_filamentmodel_info(SimParams.FilamentModel)
    
     #Check the timestep here
     @assert check_timestep(SimParams) "Timestep is too large dt=$(SimParams.dt)"
     printstyled("Timestep check passed!\n", bold=:true, color=:green)   
    @warn "The cpu version of the code is still currently in development!"
    return nothing, nothing
end

function Run(::gpu,SimParams::SimulationParams)
    print_banner()
    print_GPU_info()
    print_boundary_info(
        SimParams.boundary_x,
        SimParams.boundary_y,
        SimParams.boundary_z)
    
    print_filamentmodel_info(SimParams.FilamentModel)
    
    #Check the timestep here
    @assert check_timestep(SimParams) "Timestep is too large dt=$(SimParams.dt)"
    printstyled("Timestep check passed!\n", bold=:true, color=:green)

    #Initialise the vortex arrays [VF_initial_condition.jl]
    f, fint, pcount, nthreads, nblocks = (SimParams.initf)(SimParams)
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

    t = 0.0f0 #Simulation time
    

    f_out = []
    tt = [t]
    push!(f_out,Array(f)) #Initial configuration

    itCount = 0
    print_info_header()
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
        
        enforce_boundary!(f,SimParams.boundary_x,SimParams.boundary_y,SimParams.boundary_z)
        
        
        if mod(it, SimParams.shots) == 0
            print_info(f, ghosti, ghostb, u, Empty, SimParams, pcount, it)
            itCount += 1
            push!(f_out,Array(f))
            push!(tt,t)
        end
    end
    println("")
    printstyled("Simulation finished!\n", bold=:true, color=:green)
    return f_out, tt
end

end # module VortexFilament
