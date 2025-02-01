include("VF_velocity.jl") #Include the routines for computing the superfluid velocity

export ZeroTemperature
export SchwarzModel

### VORTEX FILAMENT MODELS

#! Zero temperature model - no dissipation mechanism
struct ZeroTemperature <: FilamentModel end
Adapt.@adapt_structure ZeroTemperature

print_filamentmodel_info(io::IO,::ZeroTemperature) = print(io,"Using the "), printstyled(io,"zero temperature model\n\n", bold=:true, color=:yellow)



# function (FM::ZeroTemperature)(f, ghosti, ghostb, u_sup, normal_velocity, Empty::Bool)
#     if Empty
#         return ZeroVector
#     else
#         # f_dot = get_deriv_1(f, ghosti, ghostb, Empty)
#         return u_sup
#     end
# end
### ---------------------------------------------------------------------


#! Schwarz model - one way coupled model of finite temperature
struct SchwarzModel{A} <: FilamentModel
    α1::A
    α2::A
end
Adapt.@adapt_structure SchwarzModel
SchwarzModel(α1::AbstractFloat,α2::AbstractFloat) = SchwarzModel{Float32}(α1,α2)

print_filamentmodel_info(io::IO, FM::SchwarzModel) = print(io, "Using the "), printstyled(io, "Schwarz model with α=$(FM.α1) and α'=$(FM.α2)\n\n", bold=:true, color=:yellow)


# function (FM::SchwarzModel)(f, ghosti, ghostb, u_sup, normal_velocity, Empty::Bool)
#     if Empty
#         return ZeroVector
#     else
#         f_dot = get_deriv_1(f, ghosti, ghostb, false)
#         u_mf = FM.α1 * cross(f_dot,normal_velocity-u_sup) - FM.α2*cross(f_dot,cross(f_dot,normal_velocity - u_sup))
#         return u_sup + u_mf
#     end
# end

# #Employing the second order AB2 Scheme
# function timestep(u, u1, u2, Empty::Bool)
#     if Empty
#         return ZeroVector
#     else
#         if(maximum(norm(u1) == 0.f0))
#             return u
#         elseif(maximum(norm(u2) == 0.f0))
#             return (3.f0/2.0f0)*u - (0.5f0)*u1
#         else
#             return (23.f0/12.f0)*u - (4.f0/3.f0)*u1 + (5.0f0/12.0f0)*u2
#         end
#     end
# end
