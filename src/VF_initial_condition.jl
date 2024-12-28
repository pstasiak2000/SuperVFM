abstract type InitCond end #Super type of all initial conditions

include("./VFConfig/setup_single_line.jl")
include("./VFConfig/setup_single_ring.jl")


function initialiseVortex!(f,fint,pcount,initf::InitCond)
    @assert supertype(typeof(initf)) == InitCond "Invalid initial condition"
    kernel = @cuda launch=false init!(f,fint,pcount,initf)
    config = launch_configuration(kernel.fun)
    threads = min(pcount, config.threads)
    blocks = cld(pcount,threads)
    CUDA.@sync begin
        kernel(f,fint,pcount, initf; threads, blocks)
    end
end