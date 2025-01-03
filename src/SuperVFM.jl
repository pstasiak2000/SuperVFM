module SuperVFM
using CUDA, Adapt
using Dates
using StaticArrays
using LinearAlgebra
using Plots

export Run



include("VF_cdata.jl")
include("VF_boundary.jl")
include("VF_timestep.jl")
include("VF_initial_condition.jl")
include("VF_derivatives.jl")
include("VF_general.jl")
include("VF_misc.jl")




function Run(SimParams::SimulationParams)
    use_animation = true
    # print_banner()
    # print_GPU_info()
    # print_boundary_info(
    #     SimParams.boundary_x,
    #     SimParams.boundary_y,
    #     SimParams.boundary_z)
    pcount = getInitPcount(SimParams.initf, SimParams.δ)
    println("-: pcount is now at $pcount")

    #!!![NOTE] Function here to determine the number of threads and blocks to be used
    nthreads = min(pcount, 1024)
    nblocks = cld(pcount, nthreads)
    #!!!
    f = CUDA.fill(SVector{3,Float32}(0, 0, 0), 12, pcount)  #Contains the vector components of the vortex filaments
    fScal = CUDA.zeros(Float32, 3, pcount)             #Contains the scalar components of the vortex filaments
    fint = CUDA.zeros(Int32, 3, pcount)                 #Contains the scalar integer components of the vortex filaments

    #Check the timestep here
    @assert check_timestep(SimParams)

    initialiseVortex!(f, fint, pcount, SimParams.initf; threads=nthreads, blocks=nblocks)
    f_init = f;
    @info "Computing the ghost points"
    @time ghostp!(f, fint, pcount, SimParams; nthreads, nblocks)

    t = 0 #Simulation time
    for it ∈ 1:SimParams.nsteps
        calc_fil_motion!(f, fint, pcount, SimParams::SimulationParams; nthreads, nblocks)
        t += SimParams.dt
    end
    return f
end


end # module VortexFilament
