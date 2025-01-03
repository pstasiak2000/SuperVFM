### Routines to calculate the special derivatives
function get_deriv_1(f, fint, pcount; nthreads=1, nblocks=1)
    s_dot = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)
    CUDA.@sync begin
        @cuda threads=nthreads blocks=nblocks deriv_1_kernel(s_dot, pcount, f, fint)
    end
    return s_dot
end

function deriv_1_kernel(s_dot, pcount, f, fint)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    disti = 0.0; distb = 0.0;
    for idx ∈ index:stride:pcount
        if fint[1,idx] == 0; continue; end;
        # fil = @SVector [f[1,idx],f[2,idx],f[3,idx]]
        # fil_i = @SVector [f[33,idx],f[34,idx],f[35,idx]]
        # fil_b = @SVector [f[36,idx],f[37,idx],f[38,idx]]

        disti = norm(f[1,idx]-f[11,idx])
        distb = norm(f[1,idx]-f[12,idx])

        s_dot[idx] += @. (distb*f[11,idx])
        s_dot[idx] += @. (disti-distb)*f[1,idx]
        s_dot[idx] -= @. disti * f[12,idx]
        s_dot[idx] /= @. 2*disti*distb
    end
    return nothing
end

function get_deriv_2(f, fint, pcount; nthreads=1, nblocks=1)
    s_ddot = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)
    CUDA.@sync begin
        @cuda threads=nthreads blocks=nblocks deriv_2_kernel(s_ddot, pcount, f, fint)
    end
    return s_ddot
end

function deriv_2_kernel(s_ddot, pcount, f, fint)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    disti = 0.0; distb = 0.0;
    tmp = @SVector [0,0,0]
    for idx ∈ index:stride:pcount
        if fint[1,idx] == 0; continue; end;
        # fil = @SVector [f[1,idx],f[2,idx],f[3,idx]]
        # fil_i = @SVector [f[33,idx],f[34,idx],f[35,idx]]
        # fil_b = @SVector [f[36,idx],f[37,idx],f[38,idx]]

        disti = norm(f[1,idx]-f[11,idx])
        distb = norm(f[1,idx]-f[12,idx])

        s_ddot[idx] += @. f[11,idx] / (disti*(disti+distb))
        s_ddot[idx] += @. f[12,idx] / (distb*(disti+distb))
        s_ddot[idx] -= @. f[1,idx] / (disti*distb)
        s_ddot[idx] *= @. 2.f0

    end
    return nothing    
end


# function dist_gen(a,b)
#     dist = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)
#     return dist
# end