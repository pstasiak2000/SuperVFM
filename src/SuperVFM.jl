module SuperVFM

using KernelAbstractions
using Adapt
using Dates
using StaticArrays
using LinearAlgebra
using Interpolations
using Unitful
using DelimitedFiles
using PrecompileTools
import Printf: @sprintf

export Run

export CPU
s = @doc(KernelAbstractions.CPU)
@doc s.text[1] KernelAbstractions.CPU

const eps32 = eps(0.0f0) #Machine epsilon of 0.0 in 32 bits
const ZeroVector = SVector{3,Float32}(0.0f0,0.0f0,0.0f0) #Zero vector

#Unit vectors to avoid constant definitions in the code
const e_x = SVector{3,Float32}(1.0f0,0.0f0,0.0f0) #unit vector in the x direction
const e_y = SVector{3,Float32}(0.0f0,1.0f0,0.0f0) #unit vector in the y direction
const e_z = SVector{3,Float32}(0.0f0,0.0f0,1.0f0) #unit vector in the z direction


include("VF_cdata.jl")
include("VF_boundary.jl")
include("VF_timestep.jl")
include("VF_initial_condition.jl")
include("VF_derivatives.jl")
include("VF_general.jl")
include("VF_output.jl")
include("VF_misc.jl")


function Run(SP::SimulationParams{S,T}) where {S,T}
    print_banner(SP)
    print_boundary_info(SP)
    print_filamentmodel_info(SP.IO, SP.FilamentModel)

    
    @assert check_timestep(SP) "Timestep is too large dt=$(SP.dt)"
    printstyled(SP.IO,"Timestep check passed!\n", bold=:true, color=:green)

    f, f_infront, f_behind, pcount = initialiseVortex(SP)

    u_loc = allocate(SP.backend, SVector{3,T}, pcount)
    u_sup = allocate(SP.backend, SVector{3,T}, pcount)

    u_mf = allocate(SP.backend, SVector{3,T}, pcount)
    u = allocate(SP.backend, SVector{3,T}, pcount)

    ### Zero the initial velocities -- important!
    u1 = KernelAbstractions.zeros(SP.backend, SVector{3,T}, pcount)
    u2 = KernelAbstractions.zeros(SP.backend, SVector{3,T}, pcount)

    
    ###  Initialise the time
    t = 0.0

    ### Initialise the time vector
    tt = []; push!(tt, t);
    itC = 0

    print_info_header(SP.IO)
    for it ∈ 1:SP.nsteps
        
        compute_filament_velocity!(u, u_mf, u_loc, u_sup,SP.FilamentModel, SP;
            f, f_infront, f_behind, pcount, SP.normal_velocity)

        timestep!(f, u, u1, u2, f_infront, pcount, SP)
        t += SP.dt

        enforce_boundary!(f, SP.boundary_x, SP.boundary_y, SP.boundary_z; f_infront, pcount, SP)

        ### NEED TO IMPLEMENT BOUNDARY FORCING
        if mod(it, SP.shots) == 0
            print_info(u, f, f_infront, f_behind, pcount, SP, it)
            itC += 1
            push!(tt,t)
        end
    end
    return f, tt
end

#     itCount = 0
#     print_info_header()
#     for it ∈ 1:SimParams.nsteps
        
#         #Compute the superfluid velocity v_s
#         @. u_loc = (SimParams.velocity)(f, ghosti, ghostb, Empty, SimParams.κ, SimParams.corea)
#         u_sup .= u_loc

#         #Calculate the velocities of the filaments
#         @. u = (SimParams.FilamentModel)(f, ghosti, ghostb, u_sup, normal_velocity, Empty)
#         #Timestep the new positions
#         @. f += SimParams.dt * timestep(u, u1, u2, Empty)
#         u2 .= u1
#         u1 .= u
#         t += SimParams.dt
        
#         enforce_boundary!(f,SimParams.boundary_x,SimParams.boundary_y,SimParams.boundary_z)
        
        
#         if mod(it, SimParams.shots) == 0
#             print_info(f, ghosti, ghostb, u, Empty, SimParams, pcount, it)
#             itCount += 1
#             push!(f_out,Array(f))
#             push!(tt,t)
#         end
#     end
#     println("")
#     printstyled("Simulation finished!\n", bold=:true, color=:green)
#     return f_out, tt
# end


include("precompile.jl")

end # module VortexFilament
