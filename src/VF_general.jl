export check_timestep

"""
    KA.zeros(BE::Backend,::Type{SVector{S,T}},dims::Tuple) where {S,T}

Initialise an array of static vectors of size `S` and type `T`
"""
KernelAbstractions.zeros(BE::Backend,::Type{SVector{S,T}},dims::Tuple) where {S,T} = T(0.0) * allocate(BE,SVector{S,T},dims)

"""
    KA.zeros(BE::Backend,::Type{SVector{S,T}},N::Int) where {S,T}

Initialise an array of static vectors of size `S` and type `T`
"""
KernelAbstractions.zeros(BE::Backend,::Type{SVector{S,T}},N::Int) where {S,T} = T(0.0) * allocate(BE,SVector{S,T},N)



"""
    check_timestep(SP::SimulationParams)

Checks if the current timestep in `SP` is small enough to resolve the smallest Kelvin waves. Returns true if the timestep ``\\Delta t < \\Delta t_{max}``.

The maximum timestep is given by
```math
    \\Delta t_{max} = \\frac{(\\delta/2)^2}{\\kappa\\log(\\delta/(2\\pi a_0))}
```
"""
function check_timestep(SP::SimulationParams)
    dt = SP.dt
    println("Performing  timestep check...")
    dt_max = ((SP.δ/2.0f0)^2)/(SP.κ * log10(SP.δ/(Float32(2π) * SP.corea)))
    @info "Max timestep is $dt_max"
    return dt < dt_max
end