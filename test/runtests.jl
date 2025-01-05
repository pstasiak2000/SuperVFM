push!(LOAD_PATH, "../") #Load the source path 
using SuperVFM
using CUDA
using StaticArrays
using LinearAlgebra
#Add here all the tests

IC = SingleRing(0.5f0)

SimParams = SuperVFM.SimulationParams(;
shots=1,
nsteps=10000,
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

f, fint, pcount, nthreads, nblocks = (IC)(SimParams.δ)
fScal = CUDA.zeros(Float32, 3, pcount) #Contains the scalar components of the vortex filaments

ghosti = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)
ghostb = CUDA.fill(SVector{3,Float32}(0,0,0),pcount)

Empty = (CUDA.reduce(+, fint, dims=1) .== 0)'

@info "Computing the ghost points" #[VF_boundary.jl]
SuperVFM.ghostp!(ghosti, ghostb, f, fint, pcount, SimParams.box_size; nthreads, nblocks)


begin
    f, fint, pcount, nthreads, nblocks = (IC)(SimParams.δ)
    fScal = CUDA.zeros(Float32, 3, pcount) #Contains the scalar components of the vortex filaments

    @info "Computing the ghost points" #[VF_boundary.jl]
    @time SuperVFM.ghostp!(f, fint, pcount, SimParams; nthreads, nblocks)

    t = 0 #Simulation time
    
    x_pos = zeros(2,SimParams.nsteps)

    f_dot = SuperVFM.get_deriv_1(f,fint, pcount; nthreads, nblocks)

    t = 0 #Simulation time 
    nsteps = 1000
    x_pos = zeros(2,nsteps)
    for it ∈ 1:nsteps#SimParams.nsteps

        SuperVFM.calc_fil_motion!(f, fint, pcount, SimParams::SimulationParams; nthreads, nblocks)
        fCPU = Array(f) 
        t += SimParams.dt

        f_dot = SuperVFM.get_deriv_1(f, fint, pcount; nthreads, nblocks)   |> Array
        f_ddot = SuperVFM.get_deriv_2(f, fint, pcount; nthreads, nblocks)   |> Array
        
        u = cross.(f_dot,f_ddot)
        
        x_pos[1,it] = t
        x_pos[2,it] = u[1][2]
    end
    plot(x_pos[1,:],x_pos[2,:])
end