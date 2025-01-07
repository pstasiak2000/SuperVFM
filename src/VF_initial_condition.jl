


VORTEX_CONFIGS = readdir("./src/VF_configs")

for config ∈ VORTEX_CONFIGS  
    include("./VF_configs/" * config)
end


function (initf::InitCond)(δ)
    pcount = getInitPcount(initf, δ)
    println("-: pcount is now at $pcount")

    f = CUDA.fill(SVector{3,Float32}(0, 0, 0), pcount)  #Contains the vector components of the vortex filaments
    fint = CUDA.zeros(Int32, 3, pcount)                 #Contains the scalar integer components of the vortex filaments

    nthreads, nblocks = redefineThreads_Blocks(pcount)
    @info "Using threads=$nthreads and blocks=$nblocks "
    # @assert supertype(typeof(initf)) == InitCond "Invalid initial condition"
    # config = launch_configuration(kernel.fun)
    # threads = min(pcount, config.threads)
    # blocks = cld(pcount,threads)
    CUDA.@sync begin
        @cuda threads=nthreads blocks=nblocks initVortex!(f,fint,pcount, initf)
    end
    return f, fint, pcount, nthreads, nblocks
end

function redefineThreads_Blocks(pcount)
    nthreads = min(pcount,1024) #max_threads_per_block
    nblocks = cld(pcount,nthreads)
    return nthreads, nblocks
end 