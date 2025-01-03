
struct LIA <: VelocityMode
end
Adapt.@adapt_structure LIA
export LIA

#Don't compute anything if type is set to nothing
function calc_velocity!(f, pcount, ::Nothing; nthreads=1, nblocks=1)
    return nothing   
end

"""
    calc_velocity(f, pcount, ::LIA; nthreads=1, nblocks=1)

Computes the superfluid velocity using the Local Induction Approximation (LIA)

```math
\\mathbf{v}_i = \\beta \\mathbf{s}'_i \\times \\mathbf{s}''_i  
```
"""
function calc_velocity(f, fint, pcount, κ, corea, ::LIA; nthreads=1, nblocks=1)

    f_dot = get_deriv_1(f, fint, pcount; nthreads, nblocks)
    f_ddot =  get_deriv_2(f, fint, pcount; nthreads, nblocks)
    
    # println(norm(Array(f_ddot)[1]))

    curv = sqrt.(CUDA.map(dot, f_ddot, f_ddot)) 
    @cuda threads=nthreads blocks=nblocks check_zero_curvature!(curv,pcount)

    beta = (κ/(4*π)) .* log10.(4.6f0 * curv / corea)
    
    # println(Array(beta)[1])
    # println(Array(curv)[1])

    u = CUDA.map(.*, beta, CUDA.map(cross, f_dot, f_ddot)) #This is the local velocity
    @cuda threads=nthreads blocks=nblocks copy_to_f!(f, fint, u, 7, pcount)
    
    # println(Array(CUDA.map(cross, f_dot, f_ddot))[1])
    @show Array(CUDA.map(cross, f_dot, f_ddot))[10]
    return u  
end

function copy_to_f!(f, fint, vec, row, pcount)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for idx ∈ index:stride:pcount
        if fint[1,idx] == 0; continue; end;
        f[row,idx] *= @. 0.f0
        f[row,idx] += vec[idx]
    end
    return nothing
end

function check_zero_curvature!(curv, pcount)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for idx ∈ index:stride:pcount
        if(curv[idx] < eps(0.0f0))
            curv[idx] = eps(0.0f0)
        else
            curv[idx] ^= -1
        end
    end
    return nothing
end