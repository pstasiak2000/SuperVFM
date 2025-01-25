# Library

## Public API

### Dimensional and dimensionless variables

```@docs
GetSchwarzTempCoeffs(Temp::Real)
```

```@docs
get_density(T::AbstractFloat)
```

```@docs
get_density(T::Unitful.Temperature)
```

```@docs
print_density_data(io::IO=stdout)
```

```@docs
get_dynamic_viscosity(T::AbstractFloat)
```

```@docs
get_dynamic_viscosity(T::Unitful.Temperature)
```

```@docs
print_dynamic_visosity_data(io::IO=stdout)
```

```@docs
check_timestep(SimParams::SimulationParams)
```

### Boundary Methods

```@docs
OpenBoundary(dims::Int)
```

```@docs
PeriodicBoundary(dims::Int)
```
!!! warning
    As of v1.0.2, only the periodic boundary condition is fully implemented and working correctly, open boundary conditions and solid walls will be implemented in a future release.

###

```@docs
Base.show(io::IO, SimParams::SimulationParams)
```



## Internal API