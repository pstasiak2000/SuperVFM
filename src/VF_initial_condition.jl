### This script contains the initialisation of vortex filaments 

VORTEX_CONFIGS = readdir(joinpath(@__DIR__, "VF_configs"))


printstyled("Loading vortex initial conditions\n", bold=:true, color=:yellow)
for config ∈ VORTEX_CONFIGS
	println("Loading in $config ...")
    include(joinpath(@__DIR__,"VF_configs",config))
end
println("Done!")

"""
    initialiseVortex(SP::SimulationParams{S,T}) where {S,T}

Uses `SP.backend` to initialise the vortex structure according to the initial condition `SP.initf`.
"""
function initialiseVortex(SP::SimulationParams{S,T}) where {S,T}
    @assert supertype(typeof(SP.initf)) == InitCond "Invalid initial condition"

    pcount = getInitpcount(SP.initf, SP)
    printVortexBanner(SP.initf, SP)
    println("-> pcount is now at $pcount")

    f = allocate(SP.backend, SVector{3,T}, pcount)
    fint = allocate(SP.backend, S, 3, pcount) 

    initVortex_kernel!(SP.backend,SP.workergroupsize)(f, fint, pcount, SP.initf, ndrange=pcount)
    println("--------------------------------------------------------")
    return f, fint, pcount
end

# function (initf::InitCond)(SimParams::SimulationParams)
#     pcount = getInitPcount(initf, SimParams)
#     println("-: pcount is now at $pcount")

#     f = CUDA.fill(SVector{3,Float32}(0, 0, 0), pcount)  #Contains the vector components of the vortex filaments
#     fint = CUDA.zeros(Int32, 3, pcount)                 #Contains the scalar integer components of the vortex filaments

#     nthreads, nblocks = redefineThreads_Blocks(pcount)
#     @info "Using threads=$nthreads and blocks=$nblocks "
#     # @assert supertype(typeof(initf)) == InitCond "Invalid initial condition"
#     # config = launch_configuration(kernel.fun)
#     # threads = min(pcount, config.threads)
#     # blocks = cld(pcount,threads)
#     CUDA.@sync begin
#         @cuda threads=nthreads blocks=nblocks initVortex!(f,fint,pcount, initf)
#     end
#     return f, fint, pcount, nthreads, nblocks
# end

# function redefineThreads_Blocks(pcount)
#     nthreads = min(pcount,1024) #max_threads_per_block
#     nblocks = cld(pcount,nthreads)
#     return nthreads, nblocks
# end 
