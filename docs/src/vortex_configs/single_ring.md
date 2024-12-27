# Single vortex ring

*Initialises a single, vortex ring of size ``R``.*
<<<<<<< HEAD

![ Single vortex ring](../assets/single_ring.svg)

## Usage
```julia
struct SingleRing{A} <: InitCond
    Radius::A
end
```
Vortex ring in the ``yz``-plane, propagating in the ``x`` direction. Use the field `Radius` to set the radius of the vortex ring.

```julia
#Set an initial condition of ring size 1.0
IC = VortexFilament.SingleRing(1.0f0)
```
=======
>>>>>>> b28cc71f9ab0bd7658202bc3e7a22387ef04881e
