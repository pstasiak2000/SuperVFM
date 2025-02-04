var documenterSearchIndex = {"docs":
[{"location":"api/#Library","page":"Library","title":"Library","text":"","category":"section"},{"location":"api/#Public-API","page":"Library","title":"Public API","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [SuperVFM]\nPrivate = false\nOrder = [:function, :type]","category":"page"},{"location":"api/#SuperVFM.GetSchwarzTempCoeffs-Tuple{Real}","page":"Library","title":"SuperVFM.GetSchwarzTempCoeffs","text":"GetSchwarzTempCoeffs(Temp::Real)\n\nObtain the Schwarz α and α' parameters from observational data using the temperature. For values that lie within the observation data, compute a cubic spline interpolation to extract parameters for the desired temperature.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.check_timestep-Tuple{SimulationParams}","page":"Library","title":"SuperVFM.check_timestep","text":"check_timestep(SP::SimulationParams)\n\nChecks if the current timestep in SP is small enough to resolve the smallest Kelvin waves. Returns true if the timestep Delta t  Delta t_max.\n\nThe maximum timestep is given by\n\n    Delta t_max = frac(delta2)^2kappalog(delta(2pi a_0))\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_density-Tuple{AbstractFloat}","page":"Library","title":"SuperVFM.get_density","text":"get_density(T::AbstractFloat)\n\nReturns the normal fluid density ρ_n and superfluid density ρ_s for a given temperature T in arbitrary units of temperature.\n\nRequires 00leq T leq T_lambda = 2178\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_density-Tuple{Union{Unitful.Quantity{T, 𝚯, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝚯, U}} where {L, S}} where {T, U}}","page":"Library","title":"SuperVFM.get_density","text":"get_density(T::Unitful.Temperature)\n\nReturns the normal fluid density ρ_n and superfluid density ρ_s for a given temperature T in arbitrary units of temperature.\n\nRequires 00leq T leq T_lambda = 2178\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_dynamic_viscosity-Tuple{AbstractFloat}","page":"Library","title":"SuperVFM.get_dynamic_viscosity","text":"get_dynamic_viscosity(T::AbstractFloat)\n\nReturns the dynamic viscosity η at temperature T.\n\nRequires 08leq T leq T_lambda = 2178\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_dynamic_viscosity-Tuple{Union{Unitful.Quantity{T, 𝚯, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝚯, U}} where {L, S}} where {T, U}}","page":"Library","title":"SuperVFM.get_dynamic_viscosity","text":"get_dynamic_viscosity(T::Unitful.Temperature)\n\nReturns the dynamic viscosity η at temperature T.\n\nRequires 08leq T leq T_lambda = 2178\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.print_characteristics-Tuple{IO, SimulationParams}","page":"Library","title":"SuperVFM.print_characteristics","text":"print_characteristics(io::IO,SimParams::SimulationParams)\n\nPrints the characteristic time and length scales of the simulation in dimensional units\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.print_density_data","page":"Library","title":"SuperVFM.print_density_data","text":"print_density_data(io::IO=stdout)\n\nPrints the normal fluid ρ_n and superfluid ρ_sdensities to the IO buffer. Defaults to stdout.\n\n\n\n\n\n","category":"function"},{"location":"api/#SuperVFM.print_dynamic_visosity_data","page":"Library","title":"SuperVFM.print_dynamic_visosity_data","text":"print_dynamic_visosity_data(io::IO=stdout)\n\nPrints the dynamic viscosity η to the IO buffer. Defaults to stdout.\n\n\n\n\n\n","category":"function"},{"location":"api/#KernelAbstractions.CPU","page":"Library","title":"KernelAbstractions.CPU","text":"CPU(; static=false)\n\nInstantiate a CPU (multi-threaded) backend.\n\nOptions:\n\nstatic: Uses a static thread assignment, this can be beneficial for NUMA aware code. Defaults to false.\n\n\n\n\n\n","category":"type"},{"location":"api/#SuperVFM.OpenBoundary","page":"Library","title":"SuperVFM.OpenBoundary","text":"OpenBoundary(dim::Int)\n\nInitialises an open boundary in the chosen direction. Vortex loops that exceed the box size (typically 2π) are not restricted. Set dim=1 for the x direction, dim=2 for the y direction and dim=3 fir the z direction.\n\nExample usage:\n\n    boundary_x = OpenBoundary(1)\n    boundary_y = OpenBoundary(2)\n    boundary_z = OpenBoundary(3)\n\n\n\n\n\n","category":"type"},{"location":"api/#SuperVFM.PeriodicBoundary","page":"Library","title":"SuperVFM.PeriodicBoundary","text":"PeriodicBoundary(dim::Int)\n\nInitialises a periodic boundary in the chosen direction selected by dim. Vortex loops that exceed the box size (typically 2π) are looped back periodically to the other side of the box.Set dim=1 for the x direction, dim=2 for the y direction and dim=3 fir the z direction. \n\nExample usage:\n\n    boundary_x = PeriodicBoundary(1)\n    boundary_y = PeriodicBoundary(2)\n    boundary_z = PeriodicBoundary(3)\n\nwarning: Warning\nAs of v1.0.2, only the periodic boundary condition is fully implemented and working correctly, open boundary conditions and solid walls will be implemented in a future release.\n\n\n\n\n\n","category":"type"},{"location":"api/#Internal-API","page":"Library","title":"Internal API","text":"","category":"section"},{"location":"api/","page":"Library","title":"Library","text":"Modules = [SuperVFM]\nPublic = false\nOrder = [:function, :type]","category":"page"},{"location":"api/#Base.show-Tuple{IO, SimulationParams}","page":"Library","title":"Base.show","text":"Base.show(io::IO,SimParams::SimulationParams)\n\nLists the simulation parameters stored in SimParams in a stylistic way with a key.\n\n\n\n\n\n","category":"method"},{"location":"api/#KernelAbstractions.zeros-Union{Tuple{T}, Tuple{S}, Tuple{KernelAbstractions.Backend, Type{StaticArraysCore.SVector{S, T}}, Int64}} where {S, T}","page":"Library","title":"KernelAbstractions.zeros","text":"KA.zeros(BE::Backend,::Type{SVector{S,T}},N::Int) where {S,T}\n\nInitialise an array of static vectors of size S and type T\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.SchwarzModel_kernel!-Tuple{Any}","page":"Library","title":"SuperVFM.SchwarzModel_kernel!","text":"SchwarzModel_kernel!(u_mf, u_sup, f, f_infront, ghosti, ghostb, FM::SchwarzModel, normal_velocity)\n\nKernel launch for the Schwarz model.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.all_periodic_enforce_kernel!-Tuple{Any}","page":"Library","title":"SuperVFM.all_periodic_enforce_kernel!","text":"all_periodic_enforce_kernel!(f, f_infront)\n\nLaunch kernel to enforce periodic boundary conditions in all dimensions.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.check_output_folder_structure-Tuple{IO}","page":"Library","title":"SuperVFM.check_output_folder_structure","text":"check_output_folder_structure(io::IO)\n\nCreates the output folder structure if it did not exist already. Prints result to buffer.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.compute_LIA_kernel!-Tuple{Any}","page":"Library","title":"SuperVFM.compute_LIA_kernel!","text":"compute_LIA_kernel!(u_loc, f, ghosti, ghostb, κ, corea)\n\nKernel to compute the superfluid velocity.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.compute_filament_velocity!-Union{Tuple{T}, Tuple{S}, Tuple{Any, Any, Any, Any, SchwarzModel, SimulationParams{S, T}}} where {S, T}","page":"Library","title":"SuperVFM.compute_filament_velocity!","text":"compute_filament_velocity!(u, u_loc, u_sup, FM::SchwarzModel, SP::SimulationParams{S,T}; kwargs...) where {S,T}\n\nCompute vortex filament velocities using the Schwarz model.\n\n    fracdmathbfsdt = \n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.compute_filament_velocity!-Union{Tuple{T}, Tuple{S}, Tuple{Any, Any, Any, Any, ZeroTemperature, SimulationParams{S, T}}} where {S, T}","page":"Library","title":"SuperVFM.compute_filament_velocity!","text":"compute_filament_velocity!(u, u_loc, u_sup, ::ZeroTemperature, SP::SimulationParams{S,T}; kwargs...) where {S,T}\n\nCompute vortex filament velocities in the zero temperature limit.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.compute_velocity!-Union{Tuple{T}, Tuple{S}, Tuple{Any, Any, LIA, SimulationParams{S, T}}} where {S, T}","page":"Library","title":"SuperVFM.compute_velocity!","text":"compute_velocity!(u_loc, u_sup, ::LIA, SP::SimulationParams{S,T},; kwargs...) where {S,T}\n\nComputes the superfluid velocity using the Local Induction Approximation (LIA)\n\nmathbfv_i = beta mathbfs_i times mathbfs_i  \n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.create_info_file-Union{Tuple{SimulationParams{S, T}}, Tuple{T}, Tuple{S}} where {S, T}","page":"Library","title":"SuperVFM.create_info_file","text":"create_info_file(::SimulationParams{S,T}) where {S,T}\n\nCreate file to print the simulation precision of floating point and integers. To be used by vortex reading methods for correct parsing.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.enforce_boundary!-Tuple{Any, OpenBoundary, OpenBoundary, OpenBoundary}","page":"Library","title":"SuperVFM.enforce_boundary!","text":"enforce_boundary!(f, boundary_x::OpenBoundary, boundary_y::OpenBoundary, boundary_z::OpenBoundary; kwargs...)\n\nEnforces open boundary conditions across all 3 dimensions.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.enforce_boundary!-Tuple{Any, PeriodicBoundary, PeriodicBoundary, PeriodicBoundary}","page":"Library","title":"SuperVFM.enforce_boundary!","text":"enforce_boundary!(f, boundary_x::PeriodicBoundary, boundary_y::PeriodicBoundary, boundary_z::PeriodicBoundary; kwargs...)\n\nEnforces periodic boundary conditions across all 3 dimensions.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_curvature!-Tuple{Any}","page":"Library","title":"SuperVFM.get_curvature!","text":"get_curvature!(curv; kwargs...)\n\nCalculate the vortex line curvature zeta.\n\n    zeta = mathbfs = leftfracdmathbfsdxiright\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_curvature_kernel!-Tuple{Any}","page":"Library","title":"SuperVFM.get_curvature_kernel!","text":"get_curvature_kernel!(curv, f, f_infront, ghosti, ghostb)\n\nKernel launch for computing vortex line curvature.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_deriv_1-Tuple{AbstractArray, AbstractArray, AbstractArray}","page":"Library","title":"SuperVFM.get_deriv_1","text":"get_deriv_1(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray)\n\nComputes the first order derivative using a 2nd order finite difference adaptive scheme using ghost points.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_deriv_2-Tuple{AbstractArray, AbstractArray, AbstractArray}","page":"Library","title":"SuperVFM.get_deriv_2","text":"get_deriv_2(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray)\n\nComputes the first order derivative using a 2nd order finite difference adaptive scheme using ghost points.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_Δξ-Union{Tuple{T}, Tuple{S}, Tuple{Any, Any, Any, Any, SimulationParams{S, T}}} where {S, T}","page":"Library","title":"SuperVFM.get_Δξ","text":"get_Δξ(f, ghosti, f_infront, pcount, SP::SimulationParams{S,T}) where {S,T}\n\nCompute the seperation distance Deltaxi between each filament.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.get_Δξ_kernel!-Tuple{Any}","page":"Library","title":"SuperVFM.get_Δξ_kernel!","text":"get_Δξ_kernel!(Δξ, f, ghosti, f_infront)\n\nLaunch kernel to compute Deltaxi.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.ghostp!-Tuple{Any, Any}","page":"Library","title":"SuperVFM.ghostp!","text":" ghostp!(ghosti, ghostb; kwargs...)\n\nIn place variant of ghostp.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.ghostp-Union{Tuple{T}, Tuple{S}, Tuple{Any, Any, Any, Any, SimulationParams{S, T}}} where {S, T}","page":"Library","title":"SuperVFM.ghostp","text":"ghostp(f, f_infront, f_behind, pcount, SP::SimulationParams{S,T}) where {S,T}\n\nCompute ghost points.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.ghostp_kernel!-Tuple{Any}","page":"Library","title":"SuperVFM.ghostp_kernel!","text":"ghostp_kernel!(ghosti, ghostb, f, fint, box_size)\n\nKernel for computation of ghost points.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.initVortex_kernel!-Tuple{Any}","page":"Library","title":"SuperVFM.initVortex_kernel!","text":"initVortex_kernel!(f, f_infront, f_behind, pcount, initf::SimpleTrefoil)\n\nLaunch kernel to initialise a simple trefoil.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.initialiseVortex-Union{Tuple{SimulationParams{S, T}}, Tuple{T}, Tuple{S}} where {S, T}","page":"Library","title":"SuperVFM.initialiseVortex","text":"initialiseVortex(SP::SimulationParams{S,T}) where {S,T}\n\nUses SP.backend to initialise the vortex structure according to the initial condition SP.initf.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.save_vortex-Tuple{Any}","page":"Library","title":"SuperVFM.save_vortex","text":"save_vortex(it;kwargs...)\n\nSave the vortex and related information to file.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.static2field-Union{Tuple{T}, Tuple{S}, Tuple{AbstractArray{StaticArraysCore.SVector{S, T}}, Int64}} where {S, T<:AbstractFloat}","page":"Library","title":"SuperVFM.static2field","text":"static2field(field::AbstractArray{SVector{S,T}}, N::Int) where {S<:Int,T<:AbstractFloat}\n\nConvert array of static vectors of size S and type T of length N, to a matrix of with dimensions (S,N). Return type is an child type of AbstractMatrix.\n\njulia> N = 128;\n\njulia> f = rand(SVector{3,Float32},N);\n\njulia> typeof(static2field(f,N)) <: AbstractMatrix\ntrue\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.timestep!-Tuple{Any, Any, Any, Any, Any, Any, SimulationParams}","page":"Library","title":"SuperVFM.timestep!","text":"timestep!(f, u, u1, u2, f_infront, pcount, SP::SimulationParams)\n\nPerform a single timestep using the second order Adams-Bashforth method.\n\n\n\n\n\n","category":"method"},{"location":"api/#SuperVFM.timestep_kernel!-Tuple{Any}","page":"Library","title":"SuperVFM.timestep_kernel!","text":"timestep_kernel!(f, u, u1, u2, f_infront, dt)\n\nKernel launch for timestep.\n\n\n\n\n\n","category":"method"},{"location":"misc/#Miscelleneus-functions","page":"Miscelleneus functions","title":"Miscelleneus functions","text":"","category":"section"},{"location":"misc/","page":"Miscelleneus functions","title":"Miscelleneus functions","text":"<!– @docs convert(VFArr::CUDA.CuMatrix{Float32},VFArrInt::CUDA.CuMatrix{Int32},VFStruct::VortexFilaments) –>","category":"page"},{"location":"vortex_configs/single_helix/#Single-helix","page":"Single helix","title":"Single helix","text":"","category":"section"},{"location":"vortex_configs/single_helix/","page":"Single helix","title":"Single helix","text":"Initialises a single, helix of amplitude A and wavenumber k.","category":"page"},{"location":"vortex_configs/single_helix/","page":"Single helix","title":"Single helix","text":"(Image:  Single vortex ring)","category":"page"},{"location":"vortex_configs/single_helix/#Usage","page":"Single helix","title":"Usage","text":"","category":"section"},{"location":"vortex_configs/single_helix/","page":"Single helix","title":"Single helix","text":"struct SingleHelix{A} <: InitCond\n    A_KW::A #Amplitude of the Kelvin wave\n    b_KW::A #Wavenumber of the Kelvin wave\n    box_length_z::A #Length in vertical z-direction\nend","category":"page"},{"location":"vortex_configs/single_helix/","page":"Single helix","title":"Single helix","text":"Helical vortex oriented in the z direction, of size box_length_z, taken to be the size of the box in z.","category":"page"},{"location":"vortex_configs/single_helix/","page":"Single helix","title":"Single helix","text":"Parameters are as follows:","category":"page"},{"location":"vortex_configs/single_helix/","page":"Single helix","title":"Single helix","text":"A_KW is the amplitude scaled by 2π, A_KW = A2π","category":"page"},{"location":"vortex_configs/single_helix/","page":"Single helix","title":"Single helix","text":"b_KW is the wavelength λ scaled by 2π, b_kw = λ2π = 1k","category":"page"},{"location":"vortex_configs/single_helix/","page":"Single helix","title":"Single helix","text":"Example:","category":"page"},{"location":"vortex_configs/single_helix/","page":"Single helix","title":"Single helix","text":"#Set a helix with amplitude A=0.2 and k=10\nIC = SingleHelix(0.2, 0.1, 2π)","category":"page"},{"location":"vortex_configs/single_ring/#Single-vortex-ring","page":"Single vortex ring","title":"Single vortex ring","text":"","category":"section"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"Initialises a single, vortex ring of size R.","category":"page"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"(Image:  Single vortex ring)","category":"page"},{"location":"vortex_configs/single_ring/#Usage","page":"Single vortex ring","title":"Usage","text":"","category":"section"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"struct SingleRing{A} <: InitCond\n    Radius::A\nend","category":"page"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"Vortex ring in the yz-plane, propagating in the x direction. Use the field Radius to set the radius of the vortex ring.","category":"page"},{"location":"vortex_configs/single_ring/","page":"Single vortex ring","title":"Single vortex ring","text":"#Set an initial condition of ring size 1.0\nIC = VortexFilament.SingleRing(1.0f0)","category":"page"},{"location":"manual/#Manual","page":"Manual","title":"Manual","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"Here everything will be listed, the routines for computing ","category":"page"},{"location":"manual/#Vortex-filament-structure-on-GPU","page":"Manual","title":"Vortex filament structure on GPU","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"In the standard implementations of the Vortex Filament Method (VFM), filaments are constructed by a Fortran 'type' array. Unfortunately, though structure types are native to Julia, structures defined and adapted using the Adapt package are required to be immutable. The implementation of the VFM requires structures to mutate according to the number of vortex points, the value will fluctuate to maintain discretisation between filaments. ","category":"page"},{"location":"manual/#Boundary-Conditions","page":"Manual","title":"Boundary Conditions","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"As of v1.0.2, only periodic boundary conditions and open boundary conditions are supported.","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"#Initialises all open boundary conditions \nboundary_x = OpenBoundary(1) #x direction\nboundary_y = OpenBoundary(2) #y direction \nboundary_z = OpenBoundary(3) #z direction\n\n#Initialises all periodic boundary conditions \nboundary_x = PeriodicBoundary(1) #x direction\nboundary_y = PeriodicBoundary(2) #y direction\nboundary_z = PeriodicBoundary(3) #z direction","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"Periodic and open boundary conditions can be mixed. For example, if we define a helical or straight line vortex, we require periodic boundary conditions in the z direction, however we may also set open boundary conditions in the x and y directions.","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"#Set a straight line vortex - this condition  requires  periodic boundary in z\nIC = SingleLine() \n\n#Set open conditions in x,y and periodic in z \nboundary_x = OpenBoundary(1)#x direction\nboundary_y = OpenBoundary(2)#y direction\nboundary_z = PeriodicBoundary(3)#z direction ","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"warning: Warning\nVortex initial conditions will throw an error if not supplied with the correct boundary condition!","category":"page"},{"location":"manual/#Vortex-Configurations","page":"Manual","title":"Vortex Configurations","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"Here is a list of all of the vortex configurations that can be used for the vortex filament method. Specific parameters relating to the initialisation of the vortex line are stored in the corresponding structure blocks.","category":"page"},{"location":"manual/#List-of-implemented-initial-conditions","page":"Manual","title":"List of implemented initial conditions","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"Pages = Main.SUBSECTION_PAGES\nDepth=1","category":"page"},{"location":"manual/#Derivatives","page":"Manual","title":"Derivatives","text":"","category":"section"},{"location":"manual/","page":"Manual","title":"Manual","text":"The Vortex Filament requires the computation of spatial derivatives. The following methods compute these derivatives using finite difference adaptive techniques.","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"First order derivative using a second order adaptive-mesh finite difference scheme:","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"fracdmathbfs_idxi = fracell_i-1mathbfs_i+1 + (ell_i+1 - ell_i-1)mathbfs_i + ell_i+1mathbfs_i-12ell_i+1ell_i-1 + cal O(ell^2) ","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"Second order derivative using a second order adaptive-mesh finite difference scheme:","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"fracd^2 mathbfs_id xi^2=frac2mathbfs_i+1ell_i+1(ell_i+1+ell_i-  1)-frac2mathbfs_iell_i+1ell_i-1+frac2mathbfs_i-1ell_i-1(ell_i+  1+ell_i-1)+cal O(ell^2)","category":"page"},{"location":"#SuperVFM.jl","page":"Home","title":"SuperVFM.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A little bit of text here to see the changes.","category":"page"},{"location":"#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"vortex_configs/straight_line/#Straight-line","page":"Straight line","title":"Straight line","text":"","category":"section"},{"location":"vortex_configs/straight_line/","page":"Straight line","title":"Straight line","text":"Initialises a single, straight line vortex at the centre of the computational box.","category":"page"},{"location":"vortex_configs/straight_line/","page":"Straight line","title":"Straight line","text":"(Image:  Straight line vortex)","category":"page"},{"location":"vortex_configs/straight_line/#Usage","page":"Straight line","title":"Usage","text":"","category":"section"},{"location":"vortex_configs/straight_line/","page":"Straight line","title":"Straight line","text":"struct SingleLine <: InitCond end","category":"page"},{"location":"vortex_configs/straight_line/","page":"Straight line","title":"Straight line","text":"To set up a single vortex line, simply pass SingleLine().","category":"page"},{"location":"vortex_configs/straight_line/","page":"Straight line","title":"Straight line","text":"IC = SingleLine()","category":"page"}]
}
