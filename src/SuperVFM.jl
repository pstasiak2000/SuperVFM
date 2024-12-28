module SuperVFM
using CUDA, Adapt
using Dates

export Run

include("VF_timestep.jl")
include("VF_initial_condition.jl")
include("VF_cdata.jl")
include("VF_general.jl")
include("VF_boundary.jl")
include("VF_misc.jl")




function Run(SimParams::SimulationParams)
    banner_print()
    print_GPU_info()

    pcount = getInitPcount(SimParams.initf,SimParams.δ)
    println("-: pcount is now at $pcount")

    f = CUDA.zeros(Float32, 39, pcount)
    fint = CUDA.zeros(Int32, 3, pcount)

    initialiseVortex!(f,fint,pcount,SimParams.initf)
    return f
end


end # module VortexFilament
