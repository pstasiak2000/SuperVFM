push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Printf
using Test
using Plots

#Vortex initial condition
IC = SingleHelix(0.2, 0.2, 2π)
# IC = SingleRing(0.25)

#Set the simulation parameters
PARAMS = SuperVFM.SimulationParams(;
    shots=100,
    nsteps=50000,
    δ=0.05f0,
    box_size=(2π, 2π, 2π),
    velocity=LIA(),
    FilamentModel=SchwarzModel(0.206, 8.34e-3),
    # FilamentModel=ZeroTemperature(),
    initf=IC,
    boundary_x=PeriodicBoundary(1),
    boundary_y=PeriodicBoundary(2),
    boundary_z=PeriodicBoundary(3),
    normal_velocity=[0.0, 0.0, 0.0],
    corea=6.29e-7,
    ν_0=0.04,
    Γ=4.8,
    dt=1e-4
)


@time f, tt = Run(cpu(),PARAMS);

let it = 1
    plot_title= @sprintf "t = %4.2f" tt[it]
    scatter(Tuple.(f[it]),
        xlim=(-π, π), xlabel="x",
        ylim=(-π, π), ylabel="y",
        zlim=(-π, π), zlabel="z",
        markerstrokewidth=0,
        markersize=2,
        linewidth=3,
        title=plot_title,
        label=nothing)
end

begin
    anim = @animate for it in eachindex(f)
        plot_title = @sprintf "t = %4.2f" tt[it]
        plt = scatter(Tuple.(f[it]),
            xlim=(-π, π), xlabel="x",
            ylim=(-π, π), ylabel="y",
            zlim=(-π, π), zlabel="z",
            markerstrokewidth=0.1,
            markersize=2,
            # linewidth=3,
            label=nothing,
            title=plot_title,
            # camera=(0,90)
        )
        display(plt)
    end
    gif(anim, "animation.gif"; fps=30)
end


