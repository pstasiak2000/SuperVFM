push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Unitful
using Printf
using Test
using Plots

make_animation = true
make_plot = false

### Set the device
# dev = CPU();
using CUDA; dev = CUDABackend()

### Set the precision (for GPU it is highly recommended to use single precision)
IntPrec = Int32;
FloatPrec = Float32;

### Vortex initial condition
# IC = SingleLine()
# IC = SingleRing(0.25)
# IC = SingleHelix(0.2, 0.2, 2π)
# IC = SimpleTrefoil{FloatPrec}(0.5)
IC = TorusKnot(17,25,0.5,2.0);

### Set the dimensional properties
DimParams = SuperVFM.DimensionalParams(;
    T=1.9u"K",
    D=0.1u"cm")


α = GetSchwarzTempCoeffs(ustrip(DimParams.T))

### Set the simulation parameters
PARAMS = SimulationParams{IntPrec,FloatPrec}(DimParams;
    backend=dev,
    shots=5000,
    nsteps=1500000,
    δ=0.05f0,
    box_size=(2π, 2π, 2π),
    velocity=LIA(),
    # FilamentModel=SchwarzModel(α[1], α[2]),
    FilamentModel=ZeroTemperature(),
    initf=IC,
    boundary_x=PeriodicBoundary{IntPrec}(1),
    boundary_y=PeriodicBoundary{IntPrec}(2),
    boundary_z=PeriodicBoundary{IntPrec}(3),
    normal_velocity=[0.0, 0.0, 0.0],
    ν_0=0.04,
    dt=1e-4,
)

### Save parameters to file
# open("parameterVF.txt","w") do io
#     show(io,PARAMS)
# end

@time f, tt = Run(PARAMS);



data = load_VF_file("OUTPUTS/VFdata/var.000000.log")
if make_plot    
    let it = 1
        plot_title = @sprintf "t = %4.2f" tt[it]
        scatter(data.xyz,
            xlim=(-π, π), xlabel="x",
            ylim=(-π, π), ylabel="y",
            zlim=(-π, π), zlabel="z",
            markerstrokewidth=0,
            markersize=0.5,
            linewidth=1,
            title=plot_title,
            label=nothing)
    end
end

DIR = readdir("OUTPUTS/VFdata")
DIR = DIR[2:end]

VF = Vector{Any}(undef,length(DIR))
for it in eachindex(DIR)
    VF[it] = load_VF_file(joinpath("OUTPUTS","VFdata",DIR[it]))
end


if make_animation
    anim = @animate for it in eachindex(VF)
        plot_title = @sprintf "t = %4.2f" VF[it].time
        plt = scatter(VF[it].xyz,
            xlim=(-π, π), xlabel="x",
            ylim=(-π, π), ylabel="y",
            zlim=(-π, π), zlabel="z",
            markerstrokewidth=0.1,
            markersize=0.5,
            # linewidth=3,
            label=nothing,
            title=plot_title,
            # camera=(0,90)
        )
        display(plt)
    end
    gif(anim, "animation.gif"; fps=30)
end
