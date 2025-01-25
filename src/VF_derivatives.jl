### Routines to calculate the special derivatives
"""
    get_deriv_1(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray, Empty::Bool)

Computes the first order derivative using a 2nd order finite difference adaptive scheme using ghost points, skipping filaments that are labelled as empty.
"""
function get_deriv_1(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray, Empty::Bool)
    # s_dot = CUDA.fill(SVector{3,Float32}(0,0,0),length(f))
    if Empty #Skip the empty particles
        return ZeroVector
    else
        disti = norm(f - ghosti)
        distb = norm(f - ghostb)
        s_dot = ((distb*ghosti) + (disti-distb)*f - disti * ghostb) / (2*disti*distb)
        return s_dot
    end
end

"""
    get_deriv_2(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray, Empty::Bool)

Computes the first order derivative using a 2nd order finite difference adaptive scheme, skipping filaments that are labelled as empty.
"""
function get_deriv_2(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray, Empty::Bool)
    # s_dot = CUDA.fill(SVector{3,Float32}(0,0,0),length(f))
    if Empty #Skip the empty particles
        return ZeroVector
    else
        disti = norm(f - ghosti)
        distb = norm(f - ghostb)
        s_ddot = 2.0f0 * (ghosti/(disti*(disti+distb)) + ghostb/(distb*(disti+distb)) - f/(disti*distb))
        return s_ddot
    end
end

# function get_deriv_2(f, ghosti, ghostb, Empty::Bool)
#     if Empty #Skip the empty particles
#         return ZeroVector
#     else
#     # s_ddot = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)
#         disti = norm(f-ghosti)
#         distb = norm(f-ghostb)

#         @. s_ddot = 2.0f0* (ghosti/(disti*(disti+distb)) + ghostb/(distb*(disti+distb)) - f/(disti*distb))
#         return s_ddot
#     end
# end
# function deriv_1_kernel(s_dot, pcount, f, fint)
#     index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#     stride = gridDim().x * blockDim().x
#     disti = 0.0; distb = 0.0;
#     for idx ∈ index:stride:pcount
#         if fint[1,idx] == 0; continue; end;
#         # fil = @SVector [f[1,idx],f[2,idx],f[3,idx]]
#         # fil_i = @SVector [f[33,idx],f[34,idx],f[35,idx]]
#         # fil_b = @SVector [f[36,idx],f[37,idx],f[38,idx]]

#         disti = norm(f[1,idx]-f[11,idx])
#         distb = norm(f[1,idx]-f[12,idx])

#         s_dot[idx] += @. (distb*f[11,idx])
#         s_dot[idx] += @. (disti-distb)*f[1,idx]
#         s_dot[idx] -= @. disti * f[12,idx]
#         s_dot[idx] /= @. 2*disti*distb
#     end
#     return nothing
# end

# function deriv_2_kernel(s_ddot, pcount, f, fint)
#     index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#     stride = gridDim().x * blockDim().x
#     disti = 0.0; distb = 0.0;
#     tmp = @SVector [0,0,0]
#     for idx ∈ index:stride:pcount
#         if fint[1,idx] == 0; continue; end;
#         # fil = @SVector [f[1,idx],f[2,idx],f[3,idx]]
#         # fil_i = @SVector [f[33,idx],f[34,idx],f[35,idx]]
#         # fil_b = @SVector [f[36,idx],f[37,idx],f[38,idx]]

#         disti = norm(f[1,idx]-f[11,idx])
#         distb = norm(f[1,idx]-f[12,idx])

#         s_ddot[idx] += @. f[11,idx] / (disti*(disti+distb))
#         s_ddot[idx] += @. f[12,idx] / (distb*(disti+distb))
#         s_ddot[idx] -= @. f[1,idx] / (disti*distb)
#         s_ddot[idx] *= @. 2.f0

#     end
#     return nothing    
# end


# function dist_gen(a,b)
#     dist = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)
#     return dist
# end