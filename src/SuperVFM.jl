module SuperVFM
using CUDA, Adapt
using Dates

export Run

include("VF_general.jl")
include("initial_condition.jl")



function Run(IC,δ)
    banner_print()
    print_GPU_info()

    pcount = getInitPcount(IC,δ)
    println("-: pcount is now at $pcount                          ")

    f = CUDA.zeros(Float32, 39, pcount)
    fint = CUDA.zeros(Int32, 3, pcount)

    initialiseVortex!(f,fint,pcount,IC)
    return f
end


end # module VortexFilament
