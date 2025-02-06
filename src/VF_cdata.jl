abstract type InitCond end #Super type of all initial conditions
abstract type VelocityMode end #Super type for computing the velocity
abstract type FilamentModel end #Super type for the filament models
abstract type BoundaryType end #Super type of all boundary conditions

export SimulationParams 
export DimensionalParams
export get_density
export get_dynamic_viscosity
export print_density_data
export print_dynamic_visosity_data

### Dimensional constants
const a₀ = 1e-10u"m"        #Helium-4 vortex core radius
const m₄ = 6.65e-24u"g"     #Mass of helium-4 atom
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

    #Default to T=1.9K if the temperature is too low for viscosity measurements
    #Typically this is only used for pure superfluid simulations
    #ρ_temp and η simply set the characteristic time scale for the simulation!
    if T < 0.8u"K"
        η = get_dynamic_viscosity(1.9u"K") #Default to T=1.9K for characteristics
        ρ_temp = get_density(1.9u"K")
        ν_n = uconvert(u"m^2/s", η / ρ_temp[1])
    else
        η = get_dynamic_viscosity(T)
        ν_n = uconvert(u"m^2/s", η / ρ_n)
    end

    


    return DimensionalParams(T,D,λ,τ,ρ_n,ρ_s,ν_n)
end

struct SimulationParams{A,B}
    backend::Any
    workergroupsize::Int64
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
    IO::IO
end
Adapt.@adapt_structure SimulationParams

function SimulationParams{S,T}(DimParams::DimensionalParams;
    backend=CPU(),
    workergroupsize=64,
    shots=1,
    nsteps=1,
    δ=0.1f0,
    box_size=(2π,2π,2π),
    velocity=LIA(),
    FilamentModel=ZeroTemperature(),
    initf=SingleLine(),
    boundary_x=PeriodicBoundary{S}(1), #Periodic boundary conditions unless specified
    boundary_y=PeriodicBoundary{S}(2), #Periodic boundary conditions unless specified
    boundary_z=PeriodicBoundary{S}(3), #Periodic boundary conditions unless specified
    normal_velocity=[0.0f0,0.0f0,0.0f0],
    ν_0=0.04,
    dt=0.001,
    IO=stdout
    ) where {S,T}
    #Compute the non-dimensional paramater Γ
    Γ = SuperVFM.κ / DimParams.ν_n

    #Compute the non-dimensional quantum of circulation κ
    κ = Γ * ν_0

    #Compute the non-dimensional vortex core radius
    corea = a₀/DimParams.λ

    #Compute the characteristic time scale
    DimParams.τ = uconvert(u"s",((DimParams.λ)^2)*(ν_0/DimParams.ν_n))

    norm_vel = @SVector [normal_velocity[1],normal_velocity[2],normal_velocity[3]]


    return SimulationParams{S,T}(backend,workergroupsize,shots,nsteps,δ,box_size,velocity,FilamentModel,initf,boundary_x,boundary_y,boundary_z,norm_vel,corea,ν_0,Γ,κ,dt,IO)
end

### Obtain the precision of simulation parameters
# Base.precision(::SimulationParams{S,T},::Type{Int}) where {S,T} = S
# Base.precision(::SimulationParams{S,T},::Type{AbstractFloat}) where {S,T} = T

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
        interp_ρs = linear_interpolation(data[:,1],data[:,2]);
        interp_ρn = linear_interpolation(data[:,1],data[:,2]);
        ρ_n = interp_ρn(T)u"g/cm^3"
        ρ_s = interp_ρs(T)u"g/cm^3"
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

Requires ``0.8\\leq T < T_{\\lambda} = 2.178``
"""
function get_dynamic_viscosity(T::AbstractFloat)
    @assert T ≥ 0.8 "Invalid temperature: T cannot be negative"
    @assert T < 2.17 "Invalid temperature: T cannot be above the transition temperature"
    data = readdlm(joinpath(@__DIR__,"data","dynamic_viscosity.txt"),skipstart=3)
    row = findall(view(data,:,1) .== T)
    if isempty(row)
        interp_η = linear_interpolation(data[:,1],data[:,2]);
        η = interp_η(T)u"Pa*s"
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
    print_dynamic_visosity_data(io::IO=stdout)

Prints the dynamic viscosity η to the IO buffer. Defaults to `stdout`.
"""
function print_dynamic_visosity_data(io::IO=stdout)
    open(joinpath(@__DIR__,"data","dynamic_viscosity.txt")) do file
        while ! eof(file)
            println(io,readline(file))
        end
    end 
end

"""
    GetSchwarzTempCoeffs(Temp::Real)

Obtain the Schwarz α and α' parameters from observational data using the temperature. For values that lie within the observation data, compute a cubic spline interpolation to extract parameters for the desired temperature.
"""
function GetSchwarzTempCoeffs(Temp::Real)
    @assert Temp >= 0 "Temperature cannot be negative"
    @assert Temp <= 2.17 "Temperature cannot exceed T_c"
    ObsData = [
        0.00 0.000 0.000e-00;
        1.30 0.034 1.383e-02;
        1.35 0.042 1.543E-02;
        1.40 0.051 1.668E-02;
        1.45 0.061 1.746E-02;
        1.50 0.072 1.766E-02;
        1.55 0.084 1.721E-02;
        1.60 0.097 1.608E-02;
        1.65 0.111 1.437E-02;
        1.70 0.126 1.225E-02;
        1.75 0.142 1.003E-02;
        1.80 0.160 8.211E-03;
        1.85 0.181 7.438E-03;
        1.90 0.206 8.340E-03;
        1.95 0.236 1.079E-02;
        2.00 0.279 1.198E-02;
        2.02 0.302 1.097E-02;
        2.04 0.330 8.318E-03;
        2.06 0.366 3.018E-03;
        2.08 0.414 -6.690E-03;
        2.10 0.481 -2.412E-02]
    ObsIndex = findall(i -> (i == Temp), ObsData)
    
    if length(ObsIndex) == 0 && Temp < 2.00
        # @info "Cubic Spline interpolating α α'"
        interp_cubic_α1 = cubic_spline_interpolation((1.30:0.05:2.00),ObsData[2:16,2])
        interp_cubic_α2 = cubic_spline_interpolation((1.30:0.05:2.00),ObsData[2:16,3])
        
        α = (interp_cubic_α1(Temp),interp_cubic_α2(Temp))
    elseif length(ObsIndex) == 0 && Temp >= 2.00
        # @info "Cubic Spline interpolating α α'"
        interp_cubic_α1 = cubic_spline_interpolation((2.00:0.02:2.10),ObsData[16:21,2])
        interp_cubic_α2 = cubic_spline_interpolation((2.00:0.05:2.10),ObsData[16:21,3])    

        α = (interp_cubic_α1(Temp),interp_cubic_α2(Temp))
    else   
    ObsIndex = ObsIndex[1][1]
    α = (ObsData[ObsIndex, 2], ObsData[ObsIndex, 3])       
    end

    return α
end