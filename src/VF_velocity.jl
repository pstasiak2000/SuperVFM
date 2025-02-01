export LIA
struct LIA <: VelocityMode end
# Adapt.@adapt_structure LIA


# #Don't compute anything if type is set to nothing
# function (Velocity::Nothing)(f, ghosti, ghostb, Empty::Bool, κ, corea)
#     return ZeroVector
# end

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
