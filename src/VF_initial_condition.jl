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

    pcount = getInitpcount(SP.initf, SP)
    printVortexBanner(SP.initf, SP)
    println(SP.IO,"-> pcount is now at $pcount")

    f = allocate(SP.backend, SVector{3,T}, pcount)
    fint = allocate(SP.backend, S, 3, pcount) 

    initVortex_kernel!(SP.backend,SP.workergroupsize)(f, fint, pcount, SP.initf, ndrange=pcount)
    println(SP.IO,"--------------------------------------------------------")
    return f, fint, pcount
end
