push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Test
using Plots

BOX_SIZE = (2π, 2π, 2π)

#Vortex initial condition
# IC = SingleHelix(0.2, 0.1, BOX_SIZE)
IC = SingleRing(0.25)

#Set the simulation parameters
PARAMS = SuperVFM.SimulationParams(;
        shots=1000,
        nsteps=200000,
        δ=0.05f0,
        box_size=BOX_SIZE,
        velocity=LIA(),
        # FilamentModel=SchwarzModel(0.206,8.34e-3),
        FilamentModel=ZeroTemperature(),
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


@time f = Run(PARAMS);

scatter(Tuple.(f[1]),
        xlim=(-π, π), xlabel="x",
        ylim=(-π, π), ylabel="y",
        zlim=(-π, π), zlabel="z",
        markerstrokewidth=0,
        markersize=2,
        linewidth=3,
        label=nothing)

anim = @animate for Vortex ∈ f[1:165]
        plt = scatter(Tuple.(Vortex),
        xlim=(-π,π), xlabel="x",
        ylim=(-π,π), ylabel="y",
        zlim=(-π,π), zlabel="z",
        markerstrokewidth=0.1,
        markersize=2,
        # linewidth=3,
        label=nothing,
        # camera=(0,90)
        )
        display(plt)
end
gif(anim, "animation.gif"; fps=30)



