include("VF_velocity.jl") #Include the routines for computing the superfluid velocity

export ZeroTemperature
export SchwarzModel


######## VORTEX FILAMENT MODELS

###! Zero temperature model - no dissipation mechanism
struct ZeroTemperature <: FilamentModel end
Adapt.@adapt_structure ZeroTemperature

print_filamentmodel_info(io::IO,::ZeroTemperature) = print(io,"Using the "), printstyled(io,"zero temperature model\n\n", bold=:true, color=:yellow)


###! Schwarz model - one way coupled model of finite temperature
struct SchwarzModel{A} <: FilamentModel
    α1::A
    α2::A
end
Adapt.@adapt_structure SchwarzModel
SchwarzModel(α1::AbstractFloat,α2::AbstractFloat) = SchwarzModel{Float32}(α1,α2)

print_filamentmodel_info(io::IO, FM::SchwarzModel) = print(io, "Using the "), printstyled(io, "Schwarz model with α=$(FM.α1) and α'=$(FM.α2)\n\n", bold=:true, color=:yellow)


################################################################################
"""
    compute_filament_velocity!(u, u_loc, u_sup, ::ZeroTemperature, SP::SimulationParams{S,T}; kwargs...) where {S,T}

Compute vortex filament velocities in the zero temperature limit. Filaments are evolved by
```math
    \\frac{d\\mathbf{s}}{dt} = \\mathbf{v}_s
```
"""
function compute_filament_velocity!(u, u_mf, u_loc, u_sup, ::ZeroTemperature, SP::SimulationParams{S,T}; kwargs...) where {S,T}
    f, f_infront, f_behind, pcount = (; kwargs...) 

    ### Compute the ghost points
    ghosti, ghostb = ghostp(f, f_infront, f_behind, pcount, SP)

    ### Compute the superfluid velocity
    compute_velocity!(u_loc, u_sup, SP.velocity, SP; f, f_infront, f_behind, pcount, ghosti, ghostb)

    u_mf .*= 0.0
    u .= u_sup
    return ghosti, ghostb
end

"""
    compute_filament_velocity!(u, u_loc, u_sup, FM::SchwarzModel, SP::SimulationParams{S,T}; kwargs...) where {S,T}

Compute vortex filament velocities using the Schwarz model.

```math
    \\frac{d\\mathbf{s}}{dt} = \\mathbf{v}_s + \\alpha\\left[\\mathbf{s}'\\times\\left(\\mathbf{v}_n - \\mathbf{v}_s\\right)\\right] - \\alpha'\\left(\\mathbf{s}'\\times\\left[\\mathbf{s}'\\times\\left(\\mathbf{v}_n - \\mathbf{v}_s\\right)\\right] \\right)
```
"""
function compute_filament_velocity!(u, u_mf, u_loc, u_sup, FM::SchwarzModel, SP::SimulationParams{S,T}; kwargs...) where {S,T}
    (; f, f_infront, f_behind, pcount) = (; kwargs...) 
    (; normal_velocity) = (; kwargs...)

    ### Compute the ghost points
    ghosti, ghostb = ghostp(f, f_infront, f_behind, pcount, SP)

    ### Compute the superfluid velocity
    compute_velocity!(u_loc, u_sup, SP.velocity, SP; f, f_infront, f_behind, pcount, ghosti, ghostb)

    kernel! = SchwarzModel_kernel!(SP.backend,SP.workergroupsize)
    kernel!(u_mf, u_sup, f, f_infront, ghosti, ghostb, FM, normal_velocity, ndrange=pcount)

    @. u = u_sup + u_mf
    return ghosti, ghostb
end

"""
    SchwarzModel_kernel!(u_mf, u_sup, f, f_infront, ghosti, ghostb, FM::SchwarzModel, normal_velocity)

Kernel launch for the Schwarz model.
"""
@kernel function SchwarzModel_kernel!(u_mf, u_sup, f, f_infront, ghosti, ghostb, FM::SchwarzModel, normal_velocity)
    Idx = @index(Global,Linear)
    if f_infront[Idx] != 0
        f_dot = get_deriv_1(f[Idx], ghosti[Idx], ghostb[Idx])
        u_mf[Idx] = FM.α1 * cross(f_dot, normal_velocity - u_sup[Idx]) - FM.α2*cross(f_dot,cross(f_dot,normal_velocity - u_sup[Idx]))
    end
end


"""
    timestep!(f, u, u1, u2, f_infront, pcount, SP::SimulationParams)

Perform a single timestep using the second order Adams-Bashforth method.
"""
function timestep!(f, u, u1, u2, f_infront, pcount, SP::SimulationParams)
    kernel! = timestep_kernel!(SP.backend, SP.workergroupsize)
    kernel!(f, u, u1, u2, f_infront, SP.dt, ndrange=pcount)
    return nothing
end


"""
    timestep_kernel!(f, u, u1, u2, f_infront, dt)

Kernel launch for timestep.
"""
@kernel function timestep_kernel!(f, u, u1, u2, f_infront, dt)
    Idx = @index(Global, Linear)
    if f_infront[Idx] != 0
        if maximum(norm(u1[Idx])) == 0
            f[Idx] += dt * u[Idx]          
        elseif maximum(norm(u2[Idx])) == 0
            f[Idx] += dt*((3.0f0/2.0f0)*u[Idx] - 0.5f0*u1[Idx])
        else
            f[Idx] += dt*((23.f0/12.f0)*u[Idx] - (4.f0/3.f00)*u1[Idx] + (5.f0/12.f0)*u2[Idx])
        end
        u2[Idx] = u1[Idx]
        u1[Idx] = u[Idx]
    end
end
