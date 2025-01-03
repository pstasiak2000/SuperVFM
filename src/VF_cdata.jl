abstract type InitCond end #Super type of all initial conditions
abstract type VelocityMode end #Super type for computing the velocity
abstract type BoundaryType end #Super type of all boundary conditions

struct SimulationParams{A,B}
    shots::A
    nsteps::A
    δ::B
    box_size::Tuple{B,B,B}
    velocity::Union{VelocityMode,Nothing} #For testing only
    initf::InitCond
    boundary_x::BoundaryType
    boundary_y::BoundaryType
    boundary_z::BoundaryType
    corea::B
    ν_0::B
    Γ::B
    κ::B
    dt::B
end
Adapt.@adapt_structure SimulationParams

function SimulationParams(;
    shots=1,
    nsteps=1,
    δ=0.1f0,
    box_size=(2π,2π,2π),
    velocity=nothing,
    initf=SingleLine(),
    boundary_x=nothing, #open boundary conditions unless specified
    boundary_y=nothing, #open boundary conditions unless specified
    boundary_z=nothing, #open boundary conditions unless specified
    corea=Float32(6.29e-7),
    ν_0=0.04f0,
    Γ=1.0f0,
    dt=0.1f0
    )
    κ = Γ * ν_0
    return SimulationParams{Int32,Float32}(shots,nsteps,δ,box_size,velocity,initf,boundary_x,boundary_y,boundary_z,corea,ν_0,Γ,κ,dt)
end

export SimulationParams 