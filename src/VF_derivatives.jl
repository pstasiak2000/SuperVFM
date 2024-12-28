### Routes to calculate the special derivatives

function get_deriv_1(f,pcount; nthreads=1, nblocks=1)
    s_dot = CUDA.zeros(pcount)
    CUDA.@sync begin
        @cuda threads=nthreads blocks=nblocks deriv_1_kernel(s_dot, pcount, f)
    end
    return s_dot
end

function deriv_1_kernel(s_dot, pcount,f)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for idx ∈ index:stride:pcount
        fil = @SVector [f[1,idx],f[2,idx],f[3,idx]]
        fil_i = @SVector [f[33,idx],f[34,idx],f[35,idx]]
        fil_b = @SVector [f[36,idx],f[37,idx],f[38,idx]]

        # disti = norm(fil - fil_i); distb = norm(fil - fil_b)
    end
    return nothing
end