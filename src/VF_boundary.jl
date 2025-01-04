#Contains all of the routines used to enforce the boundary conditions within this file
#Support for PERIODIC boundary conditions only

#---------------------------------------
#--- Periodic boundary identifier
#---------------------------------------
struct PeriodicBoundary{A} <: BoundaryType
    name::A
end
Adapt.@adapt_structure PeriodicBoundary

"""
    PeriodicBoundary()

Initialises a periodic boundary in the chosen direction. Vortex loops that exceed the box size (typically ``2π``) are looped back periodically to the other side of the box.

Example usage:
```julia
    boundary_x = PeriodicBoundary()
```
"""
function PeriodicBoundary()
    return PeriodicBoundary{String}("periodic")
end
export PeriodicBoundary
#--------------------------------------
#--- Open boundary identifier
#---------------------------------------
struct OpenBoundary{A} <: BoundaryType
    name::A
end
Adapt.@adapt_structure OpenBoundary

"""
    OpenBoundary()

Initialises an open boundary in the chosen direction. Vortex loops that exceed the box size (typically ``2π``) are not restricted. 

Example usage:
```julia
    boundary_x = OpenBoundary()
```
"""
function OpenBoundary()
    return OpenBoundary{String}("open")
end
export OpenBoundary

#######################################################################

function ghostp!(ghosti,ghostb,f,fint,pcount,box_size; nthreads=1, nblocks=1)
    Id = CuArray([#Defines the identity matrix as three static arrays
        SVector{3,Float32}(1,0,0),
        SVector{3,Float32}(0,1,0),
        SVector{3,Float32}(0,0,1)
    ])
    CUDA.@sync begin
        @cuda threads=nthreads blocks=nblocks ghostp_Kernel!(ghosti,ghostb,f,fint,pcount,box_size,Id)
    end
    return nothing
end

# #Computes the ghost points infront and behind
function ghostp_Kernel!(ghosti,ghostb,f,fint,pcount,box_size,Id)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for idx ∈ index:stride:pcount 
        # #Ghost point infront - 11 is infront
        ghosti[idx] = f[fint[1,idx]]

        # #Ghost point behind - 12 is behind
        ghostb[idx] = f[fint[2,idx]]


        #Periodic fixing for static arrays
        for c ∈ 1:3
            #Wrapping the points infront
            if sum((f[idx] - ghosti[idx]) .* Id[c])  > box_size[c]/2
                ghosti[idx] += box_size[c] * Id[c]
            elseif sum((f[idx] - ghosti[idx]).* Id[c])  < -box_size[c]/2
                ghosti[idx] -= box_size[c] * Id[c]
            end
            
            #Wrapping the points behind
            if sum((f[idx] - ghostb[idx]) .* Id[c])  > box_size[c]/2
                ghostb[idx] += box_size[c] * Id[c]
            elseif sum((f[idx] - ghostb[idx]).* Id[c])  < -box_size[c]/2
                ghostb[idx] -= box_size[c] * Id[c]
            end
        end
    end
    return nothing
end