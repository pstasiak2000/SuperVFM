push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using CUDA
using StaticArrays
using LinearAlgebra
#Add here all the tests

IC = SingleRing(0.5f0)

SimParams = SuperVFM.SimulationParams(;
shots=1,
nsteps=1,
δ=0.01f0,
box_size=(2π,2π,2π),
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
