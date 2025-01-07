



"""
    convert_to_struct(VFArr::CUDA.CuMatrix{Float32},
            VFArrInt::CUDA.CuMatrix{Int32})

Converts GPU vortex filament array into the CPU structure for saving to file.
"""
function convert_to_struct(VFArr::CUDA.CuMatrix{Float32},
                 VFArrInt::CUDA.CuMatrix{Int32})
    return nothing
end


function check_timestep(SimParams::SimulationParams)
    δ = SimParams.δ
    κ = SimParams.κ
    corea = SimParams.corea
    dt = SimParams.dt
    println("Performing  timestep check...")
    dt_max = ((δ/2.0f0)^2)/(κ * log10(δ/(Float32(2π) * corea)))
    @info "Max timestep is $dt_max"
    return dt < dt_max
end