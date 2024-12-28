#Contains all of the routines used to enforce the boundary conditions within this file
#Support for PERIODIC boundary conditions only

abstract type BoundaryType end #Super type of all boundary conditions

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

function ghostp!(f,fint,pcount,SimParams; nthreads=1, nblocks=1)
    CUDA.@sync begin
        @cuda threads=nthreads blocks=nblocks ghostp_Kernel!(f,fint,pcount,SimParams.box_size)
    end
end

#Computes the ghost points infront and behind
function ghostp_Kernel!(f,fint,pcount,box_size)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for idx ∈ index:stride:pcount 
        #Ghost point infront
        f[33,idx] = f[1,fint[1,idx]]
        f[34,idx] = f[2,fint[1,idx]]
        f[35,idx] = f[3,fint[1,idx]]

        #Ghost point behind
        
    end
    return nothing
end