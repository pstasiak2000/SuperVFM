#Contains all of the routines used to enforce the boundary conditions within this file
#Support for PERIODIC boundary conditions only


export PeriodicBoundary
export OpenBoundary

#---------------------------------------
#--- Periodic boundary identifier
#---------------------------------------
struct PeriodicBoundary{A} <: BoundaryType
    dim::Int
end
Adapt.@adapt_structure PeriodicBoundary

"""
    PeriodicBoundary(dim::Int)

Initialises a periodic boundary in the chosen direction selected by `dim`. Vortex loops that exceed the box size (typically ``2π``) are looped back periodically to the other side of the box.Set `dim=1` for the ``x`` direction, `dim=2` for the ``y`` direction and `dim=3` fir the ``z`` direction. 

Example usage:
```julia
    boundary_x = PeriodicBoundary(1)
    boundary_y = PeriodicBoundary(2)
    boundary_z = PeriodicBoundary(3)
```

!!! warning
    As of v1.0.2, only the periodic boundary condition is fully implemented and working correctly, open boundary conditions and solid walls will be implemented in a future release.
"""
PeriodicBoundary(dim::Int) = PeriodicBoundary{Int32}(dim)




#--------------------------------------
#--- Open boundary identifier
#---------------------------------------
struct OpenBoundary{A} <: BoundaryType
    dim::Int
end
Adapt.@adapt_structure OpenBoundary

"""
    OpenBoundary(dim::Int)

Initialises an open boundary in the chosen direction. Vortex loops that exceed the box size (typically ``2π``) are not restricted. Set `dim=1` for the ``x`` direction, `dim=2` for the ``y`` direction and `dim=3` fir the ``z`` direction.

Example usage:
```julia
    boundary_x = OpenBoundary(1)
    boundary_y = OpenBoundary(2)
    boundary_z = OpenBoundary(3)
```
"""
OpenBoundary(dim::Int) = OpenBoundary{Int32}(dim)

function print_boundary(io::IO, BC::OpenBoundary)
    if BC.dim == 1
        ax = "x"
    end
    if BC.dim == 2
        ax = "y"
    end
    if BC.dim == 3
        ax = "z"
    end
    print(io, "boundary_$ax: "); 
    printstyled(io, "open\n", color=:blue);
    return nothing
end

function print_boundary(io::IO, BC::PeriodicBoundary)
    if BC.dim == 1
        ax = "x"
    end
    if BC.dim == 2
        ax = "y"
    end
    if BC.dim == 3
        ax = "z"
    end
    print(io, "boundary_$ax: ")
    printstyled(io, "periodic\n", color=:blue)
    return nothing
end


#######################################################################

""" 
    ghostp(f, f_infront, f_behind, pcount, SP::SimulationParams{S,T}) where {S,T}

Compute ghost points.
"""
function ghostp(f, f_infront, f_behind, pcount, SP::SimulationParams{S,T}) where {S,T}
    ghosti = allocate(SP.backend, SVector{3,T}, pcount)
    ghostb = allocate(SP.backend, SVector{3,T}, pcount)

    kernel = ghostp_kernel!(SP.backend, SP.workergroupsize)
    kernel(ghosti, ghostb, f, f_infront, f_behind, SP.box_size, ndrange=pcount)
    return ghosti, ghostb
end

"""
     ghostp!(ghosti, ghostb; kwargs...)

In place variant of `ghostp`.
"""
function ghostp!(ghosti, ghostb; kwargs...)
    f, f_infront, f_behind, pcount, SP = (;kwargs...)
    kernel = ghostp_kernel!(SP.backend, SP.workergroupsize)
    kernel(ghosti, ghostb, f, f_infront, f_behind, SP.box_size, ndrange=pcount)
end

"""
    ghostp_kernel!(ghosti, ghostb, f, fint, box_size)

Kernel for computation of ghost points.
"""
@kernel function ghostp_kernel!(ghosti, ghostb, f, f_infront, f_behind, box_size)
    Idx = @index(Global, Linear)
    if f_infront[Idx] != 0 #Ignore empty particles
        ghosti[Idx] = f[f_infront[Idx]]
        ghostb[Idx] = f[f_behind[Idx]]


        e_x = SVector{3,eltype(box_size)}(1,0,0)
        e_y = SVector{3,eltype(box_size)}(0,1,0)
        e_z = SVector{3,eltype(box_size)}(0,0,1)
        
        ### x direction
        c  = 1
        if sum((f[Idx] - ghosti[Idx]) .* e_x) > box_size[c] / 2
            ghosti[Idx] += box_size[c] * e_x
        elseif sum((f[Idx] - ghosti[Idx]) .* e_x) < -box_size[c] / 2
            ghosti[Idx] -= box_size[c] * e_x
        end

        if sum((f[Idx] - ghostb[Idx]) .* e_x) > box_size[c] / 2
            ghostb[Idx] += box_size[c] * e_x
        elseif sum((f[Idx] - ghostb[Idx]) .* e_x) < -box_size[c] / 2
            ghostb[Idx] -= box_size[c] * e_x
        end

        ### y direction
        c = 2
        if sum((f[Idx] - ghosti[Idx]) .* e_y) > box_size[c] / 2
            ghosti[Idx] += box_size[c] * e_y
        elseif sum((f[Idx] - ghosti[Idx]) .* e_y) < -box_size[c] / 2
            ghosti[Idx] -= box_size[c] * e_y
        end

        if sum((f[Idx] - ghostb[Idx]) .* e_y) > box_size[c] / 2
            ghostb[Idx] += box_size[c] * e_y
        elseif sum((f[Idx] - ghostb[Idx]) .* e_y) < -box_size[c] / 2
            ghostb[Idx] -= box_size[c] * e_y
        end

        ### z direction
        c = 3
        if sum((f[Idx] - ghosti[Idx]) .* e_z) > box_size[c] / 2
            ghosti[Idx] += box_size[c] * e_z
        elseif sum((f[Idx] - ghosti[Idx]) .* e_z) < -box_size[c] / 2
            ghosti[Idx] -= box_size[c] * e_z
        end

        if sum((f[Idx] - ghostb[Idx]) .* e_z) > box_size[c] / 2
            ghostb[Idx] += box_size[c] * e_z
        elseif sum((f[Idx] - ghostb[Idx]) .* e_z) < -box_size[c] / 2
            ghostb[Idx] -= box_size[c] * e_z
        end
    end
end


# #######################################################################


# function boundary(f, ::Val{1})
#     s = SVector{3,Float32}(2π,0,0)
#     if norm(f.*e_x) > π
#         f -= s
#     elseif norm(f.*e_x) < -π
#         f += s
#     end
#     return f
# end

# function boundary(f, ::Val{2})
#     s = SVector{3,Float32}(0,2π,0)
#     if norm(f.*e_y) > π
#         f -= s
#     elseif norm(f.*e_y) < -π
#         f += s
#     end
#     return f
# end

# function boundary(f, ::Val{3})
#     s = SVector{3,Float32}(0,0,2π)
#     if norm(f.*e_z) > π
#         f -= s
#     elseif norm(f.*e_z) < -π
#         f += s
#     end
#     return f
# end

# function enforce_boundary!(f,boundary_x,boundary_y,boundary_z)
#     f .= boundary.(f, Val(boundary_x.dim))
#     f .= boundary.(f, Val(boundary_y.dim))
#     f .= boundary.(f, Val(boundary_z.dim))
#    return nothing
# end
