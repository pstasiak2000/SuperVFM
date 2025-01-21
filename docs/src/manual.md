# Manual

Here everything will be listed, the routines for computing 

## Vortex filament structure on GPU
In the standard implementations of the Vortex Filament Method (VFM), filaments are constructed by a Fortran 'type' array. Unfortunately, though structure types are native to Julia, structures defined and adapted using the `Adapt` package are required to be immutable. The implementation of the VFM **requires** structures to mutate according to the number of vortex points, the value will fluctuate to maintain discretisation between filaments. 

!!! tip
    The GPU array can be converted to a filament structure by calling `convert(VFArr::CUDA.CuMatrix{Float32},VFArrInt::CUDA.CuMatrix{Int32},VFStruct::VortexFilaments)`. [WIP!!!]



## Boundary Conditions
As of v1.0.2, only periodic boundary conditions and open boundary conditions are supported.
```julia
    #Initialises all open boundary conditions 
    boundary_x = OpenBoundary(1)#x direction
    boundary_y = OpenBoundary(2)#y direction 
    boundary_z = OpenBoundary(3)#z direction

    #Initialises all periodic boundary conditions 
    boundary_x = PeriodicBoundary(1)#x direction
    boundary_y = PeriodicBoundary(2)#y direction
    boundary_z = PeriodicBoundary(3)#z direction
```
Periodic and open boundary conditions can be mixed. For example, if we define a helical or straight line vortex, we require periodic boundary conditions in the ``z`` direction, however we may also set open boundary conditions in the ``x`` and ``y`` directions.
```julia
    #Set a straight line vortex - this condition  requires  periodic boundary in z
    IC = SingleLine() 

    #Sert open conditions in x,y and periodic in z 
    boundary_x = OpenBoundary(1)#x direction
    boundary_y = OpenBoundary(2)#y direction
    boundary_z = PeriodicBoundary(3)#z direction 
```
!!! warn
    Vortex initial conditions will throw an error if not supplied with the correct boundary condition!

## Vortex Configurations

Here is a list of all of the vortex configurations that can be used for the vortex filament method. Specific parameters relating to the initialisation of the vortex line are stored in the corresponding structure blocks.

##### List of implemented initial conditions
```@contents
Pages = Main.SUBSECTION_PAGES
Depth=1
```




## Derivatives 

The Vortex Filament requires the computation of spatial derivatives. The following methods compute these derivatives using finite difference adaptive techniques.

First order derivative using a second order adaptive-mesh finite difference scheme:
```math
\frac{d\mathbf{s}_i}{d\xi} = \frac{\ell_{i-1}\mathbf{s}_{i+1} + (\ell_{i+1} - \ell_{i-1})\mathbf{s}_i + \ell_{i+1}\mathbf{s}_{i-1}}{2\ell_{i+1}\ell_{i-1}} + {\cal O}(\ell^2) 
```

Second order derivative using a second order adaptive-mesh finite difference scheme:
```math
\frac{d^2 \mathbf{s}_i}{d \xi^2}=\frac{2\mathbf{s}_{i+1}}{\ell_{i+1}(\ell_{i+1}+\ell_{i-  1})}-\frac{2\mathbf{s}_i}{\ell_{i+1}\ell_{i-1}}+\frac{2\mathbf{s}_{i-1}}{\ell_{i-1}(\ell_{i+  1}+\ell_{i-1})}+{\cal O}(\ell^2)
```
