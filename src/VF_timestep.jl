include("VF_velocity.jl") #Include the routines for computing the superfluid velocity

### VORTEX FILAMENT MODELS

#! Zero temperature model - no dissipation mechanism
struct ZeroTemperature <: FilamentModel
end
Adapt.@adapt_structure ZeroTemperature
export ZeroTemperature
print_filamentmodel_info(::ZeroTemperature) = print("Using the "), printstyled("zero temperature model\n\n", bold=:true, color=:yellow)

function (FM::ZeroTemperature)(f, ghosti, ghostb, u_sup, normal_velocity, Empty::Bool)
    if Empty
        return ZeroVector
    else
        # f_dot = get_deriv_1(f, ghosti, ghostb, Empty)
        return u_sup
    end
end
### ---------------------------------------------------------------------


#! Schwarz model - one way coupled model of finite temperature
struct SchwarzModel{A} <: FilamentModel
    α1::A
    α2::A
end
Adapt.@adapt_structure SchwarzModel
SchwarzModel(α1::AbstractFloat,α2::AbstractFloat) = SchwarzModel{Float32}(α1,α2)

export SchwarzModel

print_filamentmodel_info(FM::SchwarzModel) = print("Using the "), printstyled("Schwarz model with α=$(FM.α1) and α'=$(FM.α2)\n\n", bold=:true, color=:yellow)


function (FM::SchwarzModel)(f, ghosti, ghostb, u_sup, normal_velocity, Empty::Bool)
    if Empty
        return ZeroVector
    else
        f_dot = get_deriv_1(f, ghosti, ghostb, false)
        u_mf = FM.α1 * cross(f_dot,normal_velocity-u_sup) - FM.α2*cross(f_dot,cross(f_dot,normal_velocity - u_sup))
        return u_sup + u_mf
    end
end

#Employing the second order AB2 Scheme
function timestep(u, u1, u2, Empty::Bool)
    if Empty
        return ZeroVector
    else
        if(maximum(norm(u1) == 0.f0))
            return u
        elseif(maximum(norm(u2) == 0.f0))
            return (3.f0/2.0f0)*u - (0.5f0)*u1
        else
            return (23.f0/12.f0)*u - (4.f0/3.f0)*u1 + (5.0f0/12.0f0)*u2
        end
    end
end



### NO LONGER USED:

# function calc_fil_motion!(f, u, u1, u2, fint, pcount, SimParams::SimulationParams; nthreads=1, nblocks=1)

#     #Calculate the superfluid velocity v_s
#     u_sup = calc_velocity(f, fint, pcount, SimParams.κ, SimParams.corea, SimParams.box_size, SimParams.velocity; nthreads, nblocks)

#     u = u_sup   

#     timestep!(f,u, u1, u2, SimParams.dt)
# end

# function timestep!(f, u, u1, u2, dt)
#     if(maximum(norm.(u1) == 0.f0))
#         @. f += dt * u
#     elseif(maximum(norm.(u2) == 0.f0))
#         @. f += (3.f0/2.0f)*dt*u - (0.5f0)*dt*u1
#     else
#         @. f += (23.f0/12.f0)*dt*u - (4.f0/3.f0)*dt*u1 + (5.0f0/12.0f0)*dt*u2
#     end
#     u2 .= u1
#     u1 .= u
# end

# function timestep!(f, fint, pcount, dt)
#     index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#     stride = gridDim().x * blockDim().x

#     for idx ∈ index:stride:pcount
#         if fint[1,idx] == 0; continue; end;

#         if(maximum(f[3,idx] == 0.f0))
#             f[1,idx] += dt * f[2, idx]
#         elseif(maximum(f[4,idx] == 0.f0))
#             f[1,idx] += (3.f0/2.0f)*dt*f[2,idx] - (0.5f0)*dt*f[3,idx]
#         else
#             f[1,idx] += (23.f0/12.f0)*dt*f[2,idx] - (4.f0/3.f0)*dt*f[3,idx] + (5.0f0/12.0f0)*dt*f[4,idx]
#         end
#         f[4, idx] *= @. 0.f0; f[4,idx] += f[3,idx]
#         f[3, idx] *= @. 0.f0; f[3,idx] += f[2,idx]
#     end
#     return nothing
# end