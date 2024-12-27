abstract type InitCond end #Super type of all initial conditions

include("./VFConfig/setup_single_line.jl")
include("./VFConfig/setup_single_ring.jl")


function initialiseVortex!(f,fint,pcount,IC::InitCond)
    @assert supertype(typeof(IC)) == InitCond
    kernel = @cuda launch=false init!(f,fint,pcount,IC)
    config = launch_configuration(kernel.fun)
    threads = min(pcount, config.threads)
    blocks = cld(pcount,threads)
    CUDA.@sync begin
        kernel(f,fint,pcount, IC; threads, blocks)
    end
end