



"""
    convert_to_struct(VFArr::CUDA.CuMatrix{Float32},
            VFArrInt::CUDA.CuMatrix{Int32})

Converts GPU vortex filament array into the CPU structure for saving to file.
"""
function convert_to_struct(VFArr::CUDA.CuMatrix{Float32},
                 VFArrInt::CUDA.CuMatrix{Int32})
    return nothing
end

"""
    check_timestep(SimParams::SimulationParams)

Checks if the current timestep in `SimParams` is small enough to resolve the smallest Kelvin waves. Returns true if it is true.

The maximum timestep is given by
```math
    \\Delta t_{max} = \\frac{(\\delta/2)^2}{\\kappa *\\log(\\delta/(2\\pi a_0))}
```
"""
function check_timestep(SimParams::SimulationParams)
    δ = SimParams.δ
    corea = SimParams.corea
    dt = SimParams.dt
    println("Performing  timestep check...")
    dt_max = ((δ/2.0f0)^2)/(SimParams.κ * log10(δ/(Float32(2π) * corea)))
    @info "Max timestep is $dt_max"
    return dt < dt_max
end