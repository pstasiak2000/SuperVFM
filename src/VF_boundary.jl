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

function print_boundary(BC::Boundary) where {Boundary<:BoundaryType}
    if BC.dim == 1
        ax = "x"
    end
    if BC.dim == 2
        ax = "y"
    end
    if BC.dim == 3
        ax = "z"
    end
    print_BC(ax, BC)
    return nothing
end

print_BC(ax::AbstractString, ::OpenBoundary) = print("boundary_$ax: ");
printstyled("open\n", color=:blue);

function print_boundary(BC::PeriodicBoundary)
    if BC.dim == 1
        ax = "x"
    end
    if BC.dim == 2
        ax = "y"
    end
    if BC.dim == 3
        ax = "z"
    end
    print("boundary_$ax: ")
    printstyled("periodic\n", color=:blue)
    return nothing
end


#######################################################################
function ghostp(f, fint, pcount, SP::SimulationParams{S,T}) where {S,T}
    ghosti = allocate(SP.backend, SVector{3,T}, pcount)
    ghostb = allocate(SP.backend, SVector{3,T}, pcount)

    kernel = ghostp_kernel!(SP.backend, SP.workergroupsize)
    kernel(ghosti, ghostb, f, fint, SP.box_size, ndrange=pcount)
    return ghosti, ghostb
end

@kernel function ghostp_kernel!(ghosti, ghostb, f, fint, box_size)
    Idx = @index(Global, Linear)
    if fint[1, Idx] == 0
        ghosti[Idx] = ZeroVector
        ghostb[Idx] = ZeroVector
    else
        ghosti[Idx] = f[fint[1, Idx]]
        ghostb[Idx] = f[fint[2, Idx]]
    end


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




# function ghostp!(ghosti,ghostb,f,fint,pcount,box_size; nthreads=1, nblocks=1)
#     Id = CuArray([#Defines the identity matrix as three static arrays
#         SVector{3,Float32}(1,0,0),
#         SVector{3,Float32}(0,1,0),
#         SVector{3,Float32}(0,0,1)
#     ])
#     CUDA.@sync begin
#         @cuda threads=nthreads blocks=nblocks ghostp_Kernel!(ghosti,ghostb,f,fint,pcount,box_size,Id)
#     end
#     return nothing
# end

# # #Computes the ghost points infront and behind
# function ghostp_Kernel!(ghosti,ghostb,f,fint,pcount,box_size,Id)
#     index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
#     stride = gridDim().x * blockDim().x
#     for Idx ∈ index:stride:pcount 
#         # #Ghost point infront - 11 is infront
#         ghosti[Idx] = f[fint[1,Idx]]

#         # #Ghost point behind - 12 is behind
#         ghostb[Idx] = f[fint[2,Idx]]

#         #Periodic fixing for static arrays
#         for c ∈ 1:3
#             #Wrapping the points infront
#             if sum((f[Idx] - ghosti[Idx]) .* Id[c])  > box_size[c]/2
#                 ghosti[Idx] += box_size[c] * Id[c]
#             elseif sum((f[Idx] - ghosti[Idx]).* Id[c])  < -box_size[c]/2
#                 ghosti[Idx] -= box_size[c] * Id[c]
#             end

#             #Wrapping the points behind
#             if sum((f[Idx] - ghostb[Idx]) .* Id[c])  > box_size[c]/2
#                 ghostb[Idx] += box_size[c] * Id[c]
#             elseif sum((f[Idx] - ghostb[Idx]).* Id[c])  < -box_size[c]/2
#                 ghostb[Idx] -= box_size[c] * Id[c]
#             end
#         end
#     end
#     return nothing
# end


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
