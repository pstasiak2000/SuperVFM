push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Unitful
using Printf
using Test
using Plots

make_animation = false
make_plot = false

### Set the device
dev = CPU();
# using CUDA; dev = CUDABackend()

### Set the precision (for GPU it is highly recommended to use single precision)
IntPrec = Int32;
FloatPrec = Float32;

### Vortex initial condition
# IC = SingleLine()
IC = SingleRing(1.0)
# IC = SingleHelix(0.2, 0.2, 2π)

### Set the dimensional properties
DimParams = SuperVFM.DimensionalParams(;
    T=1.9u"K",
    D=0.1u"cm")


α = GetSchwarzTempCoeffs(ustrip(DimParams.T))

### Set the simulation parameters
PARAMS = SimulationParams{IntPrec,FloatPrec}(DimParams;
    backend=dev,
    shots=100,
    nsteps=10000,
    δ=0.005f0,
    box_size=(2π, 2π, 2π),
    velocity=LIA(),
    FilamentModel=SchwarzModel(α[1], α[2]),
    # FilamentModel=ZeroTemperature(),
    initf=IC,
    boundary_x=PeriodicBoundary(1),
    boundary_y=PeriodicBoundary(2),
    boundary_z=PeriodicBoundary(3),
    normal_velocity=[0.0, 0.0, 0.0],
    ν_0=0.04,
    dt=1e-6,
    IO=IOBuffer()
)

### Save parameters to file
# open("parameterVF.txt","w") do io
#     show(io,PARAMS)
# end

@time Run(PARAMS);

# if make_plot    
#     let it = 1
#         plot_title = @sprintf "t = %4.2f" tt[it]
#         scatter(Tuple.(f[it]),
#             xlim=(-π, π), xlabel="x",
#             ylim=(-π, π), ylabel="y",
#             zlim=(-π, π), zlabel="z",
#             markerstrokewidth=0,
#             markersize=2,
#             linewidth=3,
#             title=plot_title,
#             label=nothing)
#     end
# end

# if make_animation
#     anim = @animate for it in eachindex(f)
#         plot_title = @sprintf "t = %4.2f" tt[it]
#         plt = scatter(Tuple.(f[it]),
#             xlim=(-π, π), xlabel="x",
#             ylim=(-π, π), ylabel="y",
#             zlim=(-π, π), zlabel="z",
#             markerstrokewidth=0.1,
#             markersize=2,
#             # linewidth=3,
#             label=nothing,
#             title=plot_title,
#             # camera=(0,90)
#         )
#         display(plt)
#     end
#     gif(anim, "animation.gif"; fps=30)
# end
