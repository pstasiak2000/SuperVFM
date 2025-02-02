### This script contains the initialisation of vortex filaments 

VORTEX_CONFIGS = readdir(joinpath(@__DIR__, "VF_configs"))


for config ∈ VORTEX_CONFIGS
    include(joinpath(@__DIR__,"VF_configs",config))
end

"""
    initialiseVortex(SP::SimulationParams{S,T}) where {S,T}

Uses `SP.backend` to initialise the vortex structure according to the initial condition `SP.initf`.
"""
function initialiseVortex(SP::SimulationParams{S,T}) where {S,T}
    @assert supertype(typeof(SP.initf)) == InitCond "Invalid initial condition"

    ### Get the initial vortex filament itCount
    pcount = getInitpcount(SP.initf, SP)

    ### Print vortex information
    printVortexBanner(SP.initf, SP)
    println(SP.IO,"-> pcount is now at $pcount")

    ### Initialise empty vectors for filaments
    f = allocate(SP.backend, SVector{3,T}, pcount)
    f_infront = allocate(SP.backend, S, pcount)
    f_behind = allocate(SP.backend, S, pcount) 

    ### Kernel call to initialise the vortices
    kernel! = initVortex_kernel!(SP.backend,SP.workergroupsize)
    kernel!(f, f_infront, f_behind, pcount, SP.initf, ndrange=pcount)
    println(SP.IO,"--------------------------------------------------------")
    return f, f_infront, f_behind, pcount
end
