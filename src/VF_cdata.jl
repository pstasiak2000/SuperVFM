abstract type InitCond end #Super type of all initial conditions
abstract type VelocityMode end #Super type for computing the velocity
abstract type FilamentModel end #Super type for the filament models
abstract type BoundaryType end #Super type of all boundary conditions

### Dimensional constants
const a₀ = 1e-10u"m"        #Helium-4 vortex core radius
const m₄ = 6.65e-24u"g"    #Mass of helium-4 atom
const κ = 9.97e-8u"m^2/s"   #Quantum of circulation 1

### Dimensional parameters of the simulation
mutable struct DimensionalParams
    T::Unitful.Temperature
    D::Unitful.Length
    λ::Unitful.Length
    τ::Unitful.Time
    ρ_n::typeof(1.0u"g/cm^3")
    ρ_s::typeof(1.0u"g/cm^3")
    ν_n::typeof(1.0u"m^2/s")
end

function DimensionalParams(;
        T=0.0u"K",
        D=1.0u"mm")

    λ = D/2π;
    τ = 0.0u"s" #Sets initially - re-calculated based on ν_0 when setting up the simulation parameters.
    
    # Read the normal and superfluid densities and dynamic viscosity
    ρ_n, ρ_s = get_density(T)
    η = get_dynamic_viscosity(T)

    #Compute the normal fluid kinematic viscosity
    ν_n = uconvert(u"m^2/s", η / ρ_n)
    return DimensionalParams(T,D,λ,τ,ρ_n,ρ_s,ν_n)
end

struct SimulationParams{A,B}
    shots::A
    nsteps::A
    δ::B
    box_size::Tuple{B,B,B}
    velocity::Union{VelocityMode,Nothing} #For testing only
    FilamentModel::FilamentModel
    initf::InitCond
    boundary_x::BoundaryType
    boundary_y::BoundaryType
    boundary_z::BoundaryType
    normal_velocity::SVector{3,B}
    corea::B
    ν_0::B
    Γ::B
    κ::B
    dt::B
end
Adapt.@adapt_structure SimulationParams

function SimulationParams(DimParams::DimensionalParams;
    shots=1,
    nsteps=1,
    δ=0.1f0,
    box_size=(2π,2π,2π),
    velocity=nothing,
    FilamentModel=ZeroTemperature(),
    initf=SingleLine(),
    boundary_x=nothing, #open boundary conditions unless specified
    boundary_y=nothing, #open boundary conditions unless specified
    boundary_z=nothing, #open boundary conditions unless specified
    normal_velocity=[0.0f0,0.0f0,0.0f0],
    ν_0=0.04f0,
    dt=0.1f0
    )
    #Compute the non-dimensional paramater Γ
    Γ = SuperVFM.κ / DimParams.ν_n

    #Compute the non-dimensional quantum of circulation κ
    κ = Γ * ν_0

    #Compute the non-dimensional vortex core radius
    corea = a₀/DimParams.λ

    #Compute the characteristic time scale
    DimParams.τ = uconvert(u"s",((DimParams.λ)^2)*(ν_0/DimParams.ν_n))

    norm_vel = @SVector [normal_velocity[1],normal_velocity[2],normal_velocity[3]]


    return SimulationParams{Int32,Float32}(shots,nsteps,δ,box_size,velocity,FilamentModel,initf,boundary_x,boundary_y,boundary_z,norm_vel,corea,ν_0,Γ,κ,dt)
end

export SimulationParams 

"""
    get_density(T::AbstractFloat)

Returns the normal fluid density ``ρ_n`` and superfluid density ``ρ_s`` for a given temperature ``T`` in arbitrary units of temperature.

Requires ``0.0\\leq T \\leq T_{\\lambda} = 2.178``
"""
function get_density(T::AbstractFloat)
    @assert T ≥ 0.0 "Invalid temperature: T cannot be negative"
    @assert T ≤ 2.17 "Invalid temperature: T cannot be above the transition temperature"
    data = readdlm(joinpath(@__DIR__,"data","density.txt"),skipstart=3)
    row = findall(view(data,:,1) .== T)
    if isempty(row)
        @error "Interpolation of values not yet implemented, please choose another value. The full list can be viewed by calling SuperVFM.print_density_data() "
        ρ_n = 0.0u"g/cm^3"
        ρ_s = 0.0u"g/cm^3"
    else
        ρ_n = reshape(data[row,3])u"g/cm^3"
        ρ_s = reshape(data[row,2])u"g/cm^3"
    end
    return ρ_n, ρ_s
end

"""
    get_density(T::Unitful.Temperature)

Returns the normal fluid density ``ρ_n`` and superfluid density ``ρ_s`` for a given temperature ``T`` in arbitrary units of temperature.

Requires ``0.0\\leq T \\leq T_{\\lambda} = 2.178``
"""
get_density(T::Unitful.Temperature) = get_density(ustrip(uconvert(u"K",T)))

"""
    get_dynamic_viscosity(T::AbstractFloat)

Returns the dynamic viscosity ``η`` at temperature ``T``.

Requires ``0.8\\leq T \\leq T_{\\lambda} = 2.178``
"""
function get_dynamic_viscosity(T::AbstractFloat)
    @assert T > 0.8 "Invalid temperature: Required T>0.8K"
    @assert T < 2.17 "Invalid temperature: T cannot be above the transition temperature"
    data = readdlm(joinpath(@__DIR__,"data","dynamic_viscosity.txt"),skipstart=0)
    row = findall(view(data,:,1) .== T)
    if isempty(row)
        @error "Interpolation of values not yet implemented, please choose another value. The full list can be viewed by calling SuperVFM.print_dynamic_visosity_data()"
        η = 0.0u"Pa*s"
    else
        η = reshape(data[row,2])u"Pa*s"
    end
    return η
end

"""
    get_dynamic_viscosity(T::Unitful.Temperature)

Returns the dynamic viscosity ``η`` at temperature ``T``.

Requires ``0.8\\leq T \\leq T_{\\lambda} = 2.178``
"""
get_dynamic_viscosity(T::Unitful.Temperature) = get_dynamic_viscosity(ustrip(uconvert(u"K",T)))


"""
    print_density_data(io::IO=stdout)

Prints the normal fluid ``ρ_n`` and superfluid ``ρ_s``densities to the IO buffer. Defaults to `stdout`.
"""
function print_density_data(io::IO=stdout)
    open(joinpath(@__DIR__,"data","density.txt")) do file
        while ! eof(file)
            println(io,readline(file))
        end
    end 
end

"""
    print_density_data(io::IO=stdout)

Prints the dynamic viscosity η to the IO buffer. Defaults to `stdout`.
"""
function print_dynamic_visosity_data(io::IO=stdout)
    open(joinpath(@__DIR__,"data","dynamic_viscosity.txt")) do file
        while ! eof(file)
            println(io,readline(file))
        end
    end 
end