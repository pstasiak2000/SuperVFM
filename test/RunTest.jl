push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Test
using Plots

#Vortex initial condition
IC = SingleRing(0.5f0)

#Set the simulatiprintln(Array(curv)[1])on parameters
PARAMS = SuperVFM.SimulationParams(;
        shots=2500,
        nsteps=150000,
        δ=0.05f0,
        box_size=(2π, 2π, 2π),
        velocity=LIA(),
        # FilamentModel=SchwarzModel(0.206,8.34e-3),
        FilamentModel=ZeroTemperature(),
        initf=IC,
        boundary_x=PeriodicBoundary(),
        boundary_y=PeriodicBoundary(),
        boundary_z=PeriodicBoundary(),
        normal_velocity=[0.0,0.0,0.0],
        corea=Float32(6.29e-7),
        ν_0=0.04f0,
        Γ=4.8f0,
        dt=1e-4 |> Float32
)


@time f = Run(PARAMS);

# plt = plot(xlim=(-π,π),ylim=(-π,π),zlim=(-π,π))
# plt = plot()
anim = @animate for Vortex ∈ f
        plt = plot(Tuple.(Vortex),
        xlim=(-π,π), xlabel="x",
        ylim=(-π,π), ylabel="y",
        zlim=(-π,π), zlabel="z",
        linewidth=3,
        label="R=$(IC.Radius)")
        display(plt)
end
gif(anim, "animation.gif"; fps=20)
# fCPU = Array(f)[1,1]
# fCPU_init = Array(f_init)[1,1]
# plot(Tuple.(fCPU[:,1]),
#         linewidth=5,
#         xlimits=(-π,π), xlabel="x",
#         ylimits=(-π,π), ylabel="y",
#         zlimits=(-π,π), zlabel="z",
#         label=false)
#  plot!(Tuple.(fCPU_init[:,1]),
#         linewidth=5,
#         xlimits=(-π,π), xlabel="x",
#         ylimits=(-π,π), ylabel="y",
#         zlimits=(-π,π), zlabel="z",
#         label=false)       


# println("Done")



