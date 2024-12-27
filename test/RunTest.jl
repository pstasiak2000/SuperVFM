push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM

δ=0.001
IC = SingleRing(1.2)


f = Run(IC,δ)

fCPU = Array(f')

# plot(fCPU[:,1],fCPU[:,2],fCPU[:,3],
#         linewidth=5,
#         xlimits=(-π,π),
#         ylimits=(-π,π),
#         zlimits=(-π,π),
#         label=false)

# println("Done")


 
