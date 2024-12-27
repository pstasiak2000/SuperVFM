# Manual

Here everything will be listed, the routines for computing 

## Vortex filament structure on GPU
In the standard implementations of the Vortex Filament Method (VFM), filaments are constructed by a Fortran 'type' array. Unfortunately, though structure types are native to Julia, structures defined and adapted using the `Adapt` package are required to be immutable. The implementation of the VFM **requires** structures to mutate according to the number of vortex points, the value will fluctuate to maintain discretisation between filaments. 

!!! tip
<<<<<<< HEAD
    The GPU array can be converted to a filament structure by calling `convert(VFArr::CUDA.CuMatrix{Float32},VFArrInt::CUDA.CuMatrix{Int32},VFStruct::VortexFilaments)`. [WIP!!!]
=======
    The GPU array can be converted to a filament structure by calling `convert(VFArr::CUDA.CuMatrix{Float32},VFArrInt::CUDA.CuMatrix{Int32},VFStruct::VortexFilaments)`.
>>>>>>> b28cc71f9ab0bd7658202bc3e7a22387ef04881e

#### Float32 Array
| Rows |            Variable Description                            |
|:------:|:-------------------------------------------------------- |
| 1-3  | `xyz` position of the vortex point            |
| 4-6  | `u` velocity of the vortex point | 
| 7-9  |  `u1` stored velocity for Adams-Bashforth| 
| 10-12 |  `u1` stored velocity for Adams-Bashforth|
| 13-15 | `f_u_mf` velocity due to mutual friction |
| 16-18 | `f_u_sup` velocity due to background flow|
| 19-21 | `f_u_loc` velocity due to local contribution|
| 22-24 | `f_u_pll` velocity due to flow across tangent|
| 25-27 | `f_f_mf` mutual friction of the vortex point|
| 28-30 | `f_u_n` normal fluid velocity |
| 31 | `f_curv` curvature of vortex point |
| 32 | `f_stretch` helicity density ``\mathbf{s}'\cdot\mathbf{v}_n(\mathbf{s})``|
| 33-35 | `f_ghosti` ghost point infront |
| 36-38 | `f_ghostb` ghost point behind  |
| 39    | `f_closestd` distance to closest filament|

#### Int32 Array
| Rows |            Variable Description                            |
|:------:|:-------------------------------------------------------- |
| 1  | `f_infront` vortex point infront            |
| 2  | `f_behind` vortex point behind   |       
| 3  | `f_closest` vortex point closest |



## Variables


## Vortex Configurations

Here is a list of all of the vortex configurations that can be used for the vortex filament method. Specific parameters relating to the initialisation of the vortex line are stored in the corresponding structure blocks.

##### List of implemented initial conditions
```@contents
Pages = Main.SUBSECTION_PAGES
Depth=1
```
```math
\int x dx
```

```@example 1
a = 1
b = 2
a + b
```