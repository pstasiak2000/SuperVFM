# Library

## Public API
```@docs
OpenBoundary(dims::Int)
```

```@docs
PeriodicBoundary(dims::Int)
```
!!! warning
    As of v1.0.2, only the periodic boundary condition is fully implemented and working correctly, open boundary conditions and solid walls will be implemented in a future release.

```@docs
Base.show(io::IO, SimParams::SimulationParams)
```

```@docs
GetSchwarzTempCoeffs(Temp::Real)
```

## Internal API