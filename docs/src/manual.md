# Manual

Here everything will be listed, the routines for computing 

## Vortex filament structure on GPU
In the standard implementations of the Vortex Filament Method (VFM), filaments are constructed by a Fortran 'type' array. Unfortunately, though structure types are native to Julia, structures defined and adapted using the `Adapt` package are required to be immutable. The implementation of the VFM **requires** structures to mutate according to the number of vortex points, the value will fluctuate to maintain discretisation between filaments. 

!!! tip
    The GPU array can be converted to a filament structure by calling `convert(VFArr::CUDA.CuMatrix{Float32},VFArrInt::CUDA.CuMatrix{Int32},VFStruct::VortexFilaments)`. [WIP!!!]



## Boundary Conditions
The following methods are the currently implemented boundary conditions in the vortex filament method solver.

```@docs
OpenBoundary(dims::Int)
```

```@docs
PeriodicBoundary(dims::Int)
```
!!! warning
    As of v1.0.2, only the periodic boundary condition is fully implemented and working correctly, open boundary conditions and solid walls will be implemented in a future release.

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
