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

#Directory of the simulationfile
const base_dir = pwd()

include("VF_cdata.jl")
include("VF_boundary.jl")
include("VF_timestep.jl")
include("VF_initial_condition.jl")
include("VF_derivatives.jl")
include("VF_general.jl")
include("VF_output.jl")
include("VF_misc.jl")

### Other packages
include("VortexTools.jl")


function Run(SP::SimulationParams{S,T}) where {S,T}
    print_banner(SP)
    print_boundary_info(SP)
    print_filamentmodel_info(SP.IO, SP.FilamentModel)
    generate_initial_files(SP)

    @assert check_timestep(SP) "Timestep is too large dt=$(SP.dt)"
    printstyled(SP.IO,"Timestep check passed!\n", bold=:true, color=:green)

    f, f_infront, f_behind, pcount = initialiseVortex(SP)

    u_loc = KernelAbstractions.zeros(SP.backend, SVector{3,T}, pcount)
    u_sup = KernelAbstractions.zeros(SP.backend, SVector{3,T}, pcount)

    u = KernelAbstractions.zeros(SP.backend, SVector{3,T}, pcount)
    u_mf = KernelAbstractions.zeros(SP.backend, SVector{3,T}, pcount)

    curv = KernelAbstractions.zeros(SP.backend,T, pcount)

    ### Zero the initial velocities -- important!
    u1 = KernelAbstractions.zeros(SP.backend, SVector{3,T}, pcount)
    u2 = KernelAbstractions.zeros(SP.backend, SVector{3,T}, pcount)

    
    ###  Initialise the time
    t = 0.0
    itC = 0

    save_vortex(itC; pcount, t, f, f_infront, curv, u, u_mf, u_loc, u_sup)

    print_info_header(SP.IO)

    ### Start the loop here
    for it ∈ 1:SP.nsteps
        
        #Compute superfluid velocity and the velocity of filaments
        compute_filament_velocity!(u, u_mf, u_loc, u_sup,SP.FilamentModel, SP;
            f, f_infront, f_behind, pcount, SP.normal_velocity)

        #Timestep the filaments
        timestep!(f, u, u1, u2, f_infront, pcount, SP)
        t += SP.dt

        #Enforce the boundary conditions
        enforce_boundary!(f, SP.boundary_x, SP.boundary_y, SP.boundary_z; f_infront, pcount, SP)

        #Printing and writing to file
        if mod(it, SP.shots) == 0
            itC += 1
            curv = print_info(it, SP; pcount, f, f_infront, f_behind, curv, u)

            save_vortex(itC; pcount, t, f, f_infront, curv, u, u_mf, u_loc, u_sup)
            flush(SP.IO)
        end
        
        #Quit the loop if there are not enough vortex points
        check_active_pcount(f_infront)
    end
    return f, t
end

include("precompile.jl")

end # module VortexFilament
