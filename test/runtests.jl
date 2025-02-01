push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Unitful
using Printf
using Test
using Plots

make_animation = false
make_plot = false

#Vortex initial condition
IC = SingleHelix(0.2, 0.2, 2π)
# IC = SingleRing(0.25)



### Set the dimensional properties
DimParams = SuperVFM.DimensionalParams(;
    T=0.0u"K",
    D=0.1u"cm")


α = GetSchwarzTempCoeffs(ustrip(DimParams.T))

### Set the simulation parameters
PARAMS = SimulationParams(DimParams;
    backend=CPU(),
    shots=1,
    nsteps=100,
    δ=0.05f0,
    box_size=(2π, 2π, 2π),
    velocity=LIA(),
    # FilamentModel=SchwarzModel(α[1], α[2]),
    FilamentModel=ZeroTemperature(),
    initf=IC,
    boundary_x=PeriodicBoundary(1),
    boundary_y=PeriodicBoundary(2),
    boundary_z=PeriodicBoundary(3),
    normal_velocity=[0.0, 0.0, 0.0],
    ν_0=0.04,
    dt=1e-4
)

### Save parameters to file
open("parameterVF.txt","w") do io
    show(io,PARAMS)
end

@time f, tt = Run(PARAMS);