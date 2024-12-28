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

function print_boundary_info(boundary_x,boundary_y,boundary_z)
    println("--------------------------------------------------------")
    println("                Boundary Informataion                   ")
    println("--------------------------------------------------------")
    println("boundary_x: $(boundary_x.name)")
    println("boundary_y: $(boundary_y.name)")
    println("boundary_z: $(boundary_z.name)")
end


#Computes the ghost points infront and behind
function ghostp!(f,fint,pcount,boundary_x,boundary_y,boundary_z)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for idx ∈ index:stride:pcount 
        f[33,idx] = f[1,fint[1,idx]]
        f[34,idx] = f[2,fint[2,idx]]
        f[35,idx] = f[3,fint[3,idx]]
    end
end