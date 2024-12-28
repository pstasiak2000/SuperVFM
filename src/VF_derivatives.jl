### Routines to calculate the special derivatives
"""
    get_deriv_1(f,pcount; nthreads=1, nblocks=1)

Computes the first derivative of the vortex filaments using a second order adaptive mesh finite difference scheme.

```math
\\frac{d\\mathbf{s}_i}{d\\xi} = \\frac{\\ell_{i-1}\\mathbf{s}_{i+1} + (\\ell_{i+1} - \\ell_{i-1})\\mathbf{s}_i + \\ell_{i+1}\\mathbf{s}_{i-1}}{2\\ell_{i+1}\\ell_{i-1}} + {\\cal O}(\\ell^2) 
```
"""
function get_deriv_1(f,pcount; nthreads=1, nblocks=1)
    s_dot = CUDA.zeros(3,pcount)
    CUDA.@sync begin
        @cuda threads=nthreads blocks=nblocks deriv_1_kernel(s_dot, pcount, f)
    end
    return s_dot
end

function deriv_1_kernel(s_dot, pcount,f)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    disti = 0.0; distb = 0.0;
    tmp = @SVector [0,0,0]
    for idx ∈ index:stride:pcount
        fil = @SVector [f[1,idx],f[2,idx],f[3,idx]]
        fil_i = @SVector [f[33,idx],f[34,idx],f[35,idx]]
        fil_b = @SVector [f[36,idx],f[37,idx],f[38,idx]]

        disti = dist_gen(fil,fil_i)
        distb = dist_gen(fil,fil_b)

        tmp = ((distb.*fil_i) + (disti-distb).*fil - (disti.*fil_b))./(2*disti*distb)

        s_dot[1,idx] = tmp[1]
        s_dot[2,idx] = tmp[2]
        s_dot[3,idx] = tmp[3]
    end
    return nothing
end

function get_deriv_2(f,pcount; nthreads=1, nblocks=1)
    s_ddot = CUDA.zeros(3,pcount)
    CUDA.@sync begin
        @cuda threads=nthreads blocks=nblocks deriv_2_kernel(s_ddot, pcount, f)
    end
    return s_ddot
end

function deriv_2_kernel(s_ddot, pcount, f)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    disti = 0.0; distb = 0.0;
    tmp = @SVector [0,0,0]
    for idx ∈ index:stride:pcount
        fil = @SVector [f[1,idx],f[2,idx],f[3,idx]]
        fil_i = @SVector [f[33,idx],f[34,idx],f[35,idx]]
        fil_b = @SVector [f[36,idx],f[37,idx],f[38,idx]]

        disti = dist_gen(fil,fil_i)
        distb = dist_gen(fil,fil_b)

        tmp = 2.0.*(fil_i./(disti*(disti+distb)) - fil./(disti*distb) + fil_b./(distb*(disti+distb)))

        s_ddot[1,idx] = tmp[1]
        s_ddot[2,idx] = tmp[2]
        s_ddot[3,idx] = tmp[3]
    end
    return nothing    
end


function dist_gen(a,b)
    dist = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)
    return dist
end