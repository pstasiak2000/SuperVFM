include("VF_velocity.jl") #Include the routines for computing the superfluid velocity



function calc_fil_motion!(f, fint, pcount, SimParams::SimulationParams; nthreads=1, nblocks=1)

    #Calculate the superfluid velocity v_s
    u_sup = calc_velocity(f, fint, pcount, SimParams.κ, SimParams.corea, SimParams.velocity; nthreads, nblocks)

    u = u_sup   

    CUDA.@sync begin
        @cuda threads=nthreads blocks=nblocks copy_to_f!(f, fint, u, 2, pcount)
        @cuda threads=nthreads blocks=nblocks timestep!(f, fint, pcount, SimParams.dt)
    end
end

function timestep!(f, fint, pcount, dt)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x

    for idx ∈ index:stride:pcount
        if fint[1,idx] == 0; continue; end;

        if(maximum(f[3,idx] == 0.f0))
            f[1,idx] += dt * f[2, idx]
        elseif(maximum(f[4,idx] == 0.f0))
            f[1,idx] += (3.f0/2.0f)*dt*f[2,idx] - (0.5f0)*dt*f[3,idx]
        else
            f[1,idx] += (23.f0/12.f0)*dt*f[2,idx] - (4.f0/3.f0)*dt*f[3,idx] + (5.0f0/12.0f0)*dt*f[4,idx]
        end
        f[4, idx] *= @. 0.f0; f[4,idx] += f[3,idx]
        f[3, idx] *= @. 0.f0; f[3,idx] += f[2,idx]
    end
    return nothing
end