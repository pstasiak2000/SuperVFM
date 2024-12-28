push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM

#Vortex initial condition
IC = SingleRing(1.2)

#Set the simulation parameters
PARAMS = SuperVFM.SimulationParams(;
        shots=1,
        nsteps=1,
        δ=0.01,
        box_size=(2π,2π,2π),
        velocity=nothing,
        initf=IC,
        boundary=("periodic","periodic","periodic"),
        corea=Float32(6.29e-7),
        dt=0.01
)

f = Run(PARAMS)

fCPU = Array(f')

# plot(fCPU[:,1],fCPU[:,2],fCPU[:,3],
#         linewidth=5,
#         xlimits=(-π,π), xlabel="x",
#         ylimits=(-π,π), ylabel="y",
#         zlimits=(-π,π), zlabel="z",
#         label=false)

println("Done")


 
