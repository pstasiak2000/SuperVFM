abstract type InitCond end #Super type of all initial conditions


VORTEX_CONFIGS = readdir("./src/VF_configs")

for config ∈ VORTEX_CONFIGS  
    include("./VF_configs/" * config)
end


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