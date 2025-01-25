export LIA
struct LIA <: VelocityMode end
Adapt.@adapt_structure LIA


#Don't compute anything if type is set to nothing
function (Velocity::Nothing)(f, ghosti, ghostb, Empty::Bool, κ, corea)
    return ZeroVector
end

"""
    (Velocity::LIA)(f, ghosti, ghostb, Empty, κ, corea)

Computes the superfluid velocity using the Local Induction Approximation (LIA)

```math
\\mathbf{v}_i = \\beta \\mathbf{s}'_i \\times \\mathbf{s}''_i  
```
"""
function (Velocity::LIA)(f, ghosti, ghostb, Empty::Bool, κ, corea)   
    if Empty
        return ZeroVector
    else
        f_dot = get_deriv_1(f, ghosti, ghostb, false)
        f_ddot = get_deriv_2(f, ghosti, ghostb, false)
        
        curv = sqrt(dot(f_ddot,f_ddot))
        if(curv < eps32)
            curv = eps32
        else
            curv ^= -1
        end 
        beta = (κ/Float32(4*π)) * log10(4.6f0 * curv / corea)
        u_loc = beta * cross(f_dot,f_ddot)
        return u_loc
    end 
end




### NO LONGER USED:

# """
#     calc_velocity(f, pcount, ::LIA; nthreads=1, nblocks=1)

# Computes the superfluid velocity using the Local Induction Approximation (LIA)

# ```math
# \\mathbf{v}_i = \\beta \\mathbf{s}'_i \\times \\mathbf{s}''_i  
# ```
# """
# function calc_velocity(f, fint, pcount, κ, corea, box_size,::LIA; nthreads=1, nblocks=1)
#     ghosti = similar(f)
#     ghostb = similar(f)

#     #Computing the ghost points #[VF_boundary.jl]
#     ghostp!(ghosti, ghostb, f, fint, pcount, box_size) 
  
#     f_dot = get_deriv_1(f, ghosti, ghostb, pcount)
#     f_ddot =  get_deriv_2(f, ghosti, ghostb, pcount)

#     curv = sqrt.(CUDA.map(dot, f_ddot, f_ddot)) 
#     @cuda threads=nthreads blocks=nblocks check_zero_curvature!(curv,pcount)

#     beta = (κ/(4*π)) .* log10.(4.6f0 * curv / corea)
    
#     # println(Array(beta)[1])
#     # println(Array(curv)[1])

#     u = CUDA.map(.*, beta, CUDA.map(cross, f_dot, f_ddot)) #This is the local velocity
    
#     # println(Array(CUDA.map(cross, f_dot, f_ddot))[1])
#     # @show Array(CUDA.map(cross, f_dot, f_ddot))[10]
#     return u  
# end

# function check_zero_curvature!(curv, pcount)
#     index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#     stride = gridDim().x * blockDim().x
#     for idx ∈ index:stride:pcount
#         if(curv[idx] < eps(0.0f0))
#             curv[idx] = eps(0.0f0)
#         else
#             curv[idx] ^= -1
#         end
#     end
#     return nothing
# end