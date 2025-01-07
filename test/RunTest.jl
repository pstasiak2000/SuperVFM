push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Test
using Plots

BOX_SIZE = (2π, 2π, 2π)

#Vortex initial condition
IC = SingleHelix(0.2, 0.1, BOX_SIZE)

#Set the simulation parameters
PARAMS = SuperVFM.SimulationParams(;
        shots=1000,
        nsteps=60000,
        δ=0.05f0,
        box_size=BOX_SIZE,
        velocity=LIA(),
        # FilamentModel=SchwarzModel(0.206,8.34e-3),
        FilamentModel=ZeroTemperature(),
        initf=IC,
        boundary_x=PeriodicBoundary(),
        boundary_y=PeriodicBoundary(),
        boundary_z=PeriodicBoundary(),
        normal_velocity=[0.0, 0.0, 1.0],
        corea=Float32(6.29e-7),
        ν_0=0.04f0,
        Γ=4.8f0,
        dt=1e-4 |> Float32
)


@time f = Run(PARAMS);

plot(Tuple.(f[1]),
        xlim=(-π, π), xlabel="x",
        ylim=(-π, π), ylabel="y",
        zlim=(-π, π), zlabel="z",
        linewidth=3,
        label=nothing)

anim = @animate for Vortex ∈ f
        plt = plot(Tuple.(Vortex),
        xlim=(-π,π), xlabel="x",
        ylim=(-π,π), ylabel="y",
        zlim=(-π,π), zlabel="z",
        linewidth=3,
        label=nothing,
        # camera=(0,90)
        )
        display(plt)
end
gif(anim, "animation.gif"; fps=30)



