export LIA
struct LIA <: VelocityMode end
Adapt.@adapt_structure LIA


"""
    compute_velocity!(u_loc, u_sup, ::LIA, SP::SimulationParams{S,T},; kwargs...) where {S,T}

Computes the superfluid velocity using the Local Induction Approximation (LIA)

```math
\\mathbf{v}_i = \\beta \\mathbf{s}'_i \\times \\mathbf{s}''_i  
```
"""
function compute_velocity!(u_loc, u_sup, ::LIA, SP::SimulationParams{S,T}; kwargs...) where {S,T}
    (; f, f_infront, f_behind, pcount) = (; kwargs...)
    (; ghosti, ghostb) = (; kwargs...)

    kernel! = compute_LIA_kernel!(SP.backend,SP.workergroupsize)
    kernel!(u_loc, f, f_infront, ghosti, ghostb, SP.κ, SP.corea, ndrange=pcount)
    u_sup .= u_loc
    return nothing
end
"""
    compute_LIA_kernel!(u_loc, f, ghosti, ghostb, κ, corea)

Kernel to compute the superfluid velocity.
"""
@kernel function compute_LIA_kernel!(u_loc, f, f_infront, ghosti, ghostb, κ, corea)
    Idx = @index(Global, Linear)
    if f_infront[Idx] != 0
        f_dot = get_deriv_1(f[Idx], ghosti[Idx], ghostb[Idx])
        f_ddot = get_deriv_2(f[Idx], ghosti[Idx], ghostb[Idx])

        curv = sqrt(dot(f_ddot,f_ddot))
        if(curv < eps32)
            curv = eps32
        else
            curv ^= -1
        end
        beta = (κ/Float32(4*π)) * log10(4.6f0 * curv / corea)
        u_loc[Idx] = beta * cross(f_dot,f_ddot)
    end
end
# """
#     (Velocity::LIA)(f, ghosti, ghostb, Empty, κ, corea)

# Computes the superfluid velocity using the Local Induction Approximation (LIA)

# ```math
# \\mathbf{v}_i = \\beta \\mathbf{s}'_i \\times \\mathbf{s}''_i  
# ```
# """
# function (Velocity::LIA)(f, ghosti, ghostb, Empty::Bool, κ, corea)   
#     if Empty
#         return ZeroVector
#     else
#         f_dot = get_deriv_1(f, ghosti, ghostb, false)
#         f_ddot = get_deriv_2(f, ghosti, ghostb, false)
        
#         curv = sqrt(dot(f_ddot,f_ddot))
#         if(curv < eps32)
#             curv = eps32
#         else
#             curv ^= -1
#         end 
#         beta = (κ/Float32(4*π)) * log10(4.6f0 * curv / corea)
#         u_loc = beta * cross(f_dot,f_ddot)
#         return u_loc
#     end 
# end
