struct SimulationParams{A,B}
    shots::A
    nsteps::A
    δ::B
    box_size::Tuple{B,B,B}
    velocity::Union{VelocityMode,Nothing} #For testing only
    initf::InitCond
    boundary::Tuple{String,String,String}
    corea::B
    dt::B
end
Adapt.@adapt_structure SimulationParams

function SimulationParams(;
    shots=1,
    nsteps=1,
    δ=0.1,
    box_size=(2π,2π,2π),
    velocity=nothing,
    initf=SingleLine(),
    boundary=("periodic","periodic","periodic"),
    corea=Float32(6.29e-7),
    dt=0.1)
    return SimulationParams{Int32,Float32}(shots,nsteps,δ,box_size,velocity,initf,boundary,corea,dt)
end

export SimulationParams 