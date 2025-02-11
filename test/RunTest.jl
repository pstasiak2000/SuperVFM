push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Unitful
using Printf
using Test
using Plots

make_animation = true
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
# IC = SimpleTrefoil{FloatPrec}(0.5)
# IC = TorusKnot(3,-2,0.5,2.0);

### Set the dimensional properties
DimParams = SuperVFM.DimensionalParams(;
    T=2.1u"K",
    D=0.1u"cm")


α = GetSchwarzTempCoeffs(ustrip(DimParams.T))

### Set the simulation parameters
PARAMS = SimulationParams{IntPrec,FloatPrec}(DimParams;
    backend=dev,
    shots=1000,
    nsteps=100000,
    δ=0.01f0,
    box_size=(2π, 2π, 2π),
    velocity=LIA(),
    FilamentModel=SchwarzModel(α[1], α[2]),
    # FilamentModel=ZeroTemperature(),
    initf=IC,
    boundary_x=PeriodicBoundary{IntPrec}(1),
    boundary_y=PeriodicBoundary{IntPrec}(2),
    boundary_z=PeriodicBoundary{IntPrec}(3),
    normal_velocity=[-2.0, 0.0, 0.0],
    ν_0=0.04,
    dt=1e-5,
)

### Save parameters to file
# open("parameterVF.txt","w") do io
#     show(io,PARAMS)
# end

@time f, tt = Run(PARAMS);
s = scatter(Tuple.(f),xlim=(-π,π),ylim=(-π,π),zlim=(-π,π))
display(s)


if make_plot
    let it = 70
        itstr = @sprintf "%06i" it
        data = load_VF_file("OUTPUTS/VFdata/var.$itstr.log")
        plot_title = @sprintf "t = %4.2f" data.time
        s = scatter(data.xyz,
            xlim=(-π, π), xlabel="x",
            ylim=(-π, π), ylabel="y",
            zlim=(-π, π), zlabel="z",
            markerstrokewidth=0,
            markersize=0.5,
            linewidth=1,
            title=plot_title,
            label=nothing)
            display(s)
    end
end




if make_animation
    DIR = readdir("OUTPUTS/VFdata")
    DIR = DIR[2:end]

    VF = Vector{Any}(undef, length(DIR))
    for it in eachindex(DIR)
        VF[it] = load_VF_file(joinpath("OUTPUTS", "VFdata", DIR[it]))
    end


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
