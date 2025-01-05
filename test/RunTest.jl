push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using Test

#Vortex initial condition
IC = SingleRing(1.0f0)

#Set the simulatiprintln(Array(curv)[1])on parameters
PARAMS = SuperVFM.SimulationParams(;
        shots=25,
        nsteps=1000,
        δ=0.05f0,
        box_size=(2π, 2π, 2π),
        velocity=LIA(),
        initf=IC,
        boundary_x=PeriodicBoundary(),
        boundary_y=PeriodicBoundary(),
        boundary_z=PeriodicBoundary(),
        corea=Float32(6.29e-7),
        ν_0=0.04f0,
        Γ=4.8f0,
        dt=1e-6 |> Float32
)

@time f, x_pos = Run(PARAMS);



plot(x_pos[1,:],x_pos[2,:])

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



