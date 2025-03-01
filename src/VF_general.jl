### List of useful functions to be used within the package.

export check_timestep

"""
    static2field(field::AbstractArray{SVector{S,T}}, N::Int) where {S<:Int,T<:AbstractFloat} 

Convert array of static vectors of size `S` and type `T` of length `N`, to a matrix of with dimensions `(S,N)`. Return type is an child type of `AbstractMatrix`.

```
julia> N = 128;

julia> f = rand(SVector{3,Float32},N);

julia> typeof(static2field(f,N)) <: AbstractMatrix
true
```
"""
static2field(field::AbstractArray{SVector{S,T}}, N::Int) where {S,T<:AbstractFloat} = reshape(reinterpret(T, reshape(field,N)),S, :)



"""
    check_output_folder_structure(io::IO)

Creates the output folder structure if it did not exist already. Prints result to buffer.
"""
function check_output_folder_structure(io::IO)
    if isdir(joinpath(base_dir,"OUTPUTS")) && isdir(joinpath(base_dir,"OUTPUTS","VFdata"))
        println(io,"OUTPUTS directory already exists!")
    else
        println(io,"OUTPUTS directory not found. Creating output directory...")
        mkpath(joinpath(base_dir,"OUTPUTS","VFdata"))
        #Add more diectories here if needed
    end
end


"""
    check_active_pcount(f_infront)

Exits the run if there are not enough vortex points, if there are less than 5 points (needed to properly resolve the derivatives).
"""
function check_active_pcount(f_infront)
    @assert sum(f_infront .> 0) >= 5 "pcount too small: quitting loop..."
end

"""
    KA.zeros(BE::Backend,::Type{SVector{S,T}},N::Int) where {S,T}

Initialise an array of static vectors of size `S` and type `T`
"""
function KernelAbstractions.zeros(BE::Backend,::Type{SVector{S,T}},N::Int64) where {S,T}
    @kernel function zero_kernel!(arr)
        Idx = @index(Global, Linear)
        arr[Idx] = ZeroVector
    end

    arr = allocate(BE,SVector{S,T},N)
    zero_kernel!(BE,64)(arr,ndrange=N)
    return arr
end

"""
    get_Δξ(f, ghosti, f_infront, pcount, SP::SimulationParams{S,T}) where {S,T}

Compute the seperation distance ``\\Delta\\xi`` between each filament.
"""
function get_Δξ(f, ghosti, f_infront, pcount, SP::SimulationParams{S,T}) where {S,T}
    Δξ = KernelAbstractions.zeros(SP.backend, T, pcount) 
    kernel! = get_Δξ_kernel!(SP.backend, SP.workergroupsize)
    kernel!(Δξ, f, ghosti, f_infront, ndrange=pcount)
    return Δξ
end

"""
    get_Δξ_kernel!(Δξ, f, ghosti, f_infront)

Launch kernel to compute ``\\Delta\\xi``.
"""
@kernel function get_Δξ_kernel!(Δξ, f, ghosti, f_infront)
    I = @index(Global, Linear)
    if f_infront[I] != 0
        Δξ[I] = norm(f[I]-ghosti[I])
    end
end

"""
    get_curvature!(curv; kwargs...)

Calculate the vortex line curvature ``\\zeta``.

```math
    \\zeta = |\\mathbf{s}''| = \\left|\\frac{d\\mathbf{s}}{d\\xi}\\right|
```
"""
function get_curvature!(curv; kwargs...)
    (; f, f_infront, ghosti, ghostb, pcount, SP) = (; kwargs...)

    kernel! = get_curvature_kernel!(SP.backend,SP.workergroupsize)
    kernel!(curv, f, f_infront, ghosti, ghostb, ndrange=pcount)
    return nothing
end

"""
    get_curvature_kernel!(curv, f, f_infront, ghosti, ghostb)

Kernel launch for computing vortex line curvature.
"""
@kernel function get_curvature_kernel!(curv, f, f_infront, ghosti, ghostb)
    Idx = @index(Global, Linear)
    if f_infront[Idx] != 0 
        curv[Idx] = norm(get_deriv_2(f[Idx],ghosti[Idx], ghostb[Idx]))
    end
end

"""
    check_timestep(SP::SimulationParams)

Checks if the current timestep in `SP` is small enough to resolve the smallest Kelvin waves. Returns true if the timestep ``\\Delta t < \\Delta t_{max}``.

The maximum timestep is given by
```math
    \\Delta t_{max} = \\frac{(\\delta/2)^2}{\\kappa\\log(\\delta/(2\\pi a_0))}
```
"""
function check_timestep(SP::SimulationParams)
    dt = SP.dt
    println(SP.IO,"Performing  timestep check...")
    dt_max = ((SP.δ/2.0f0)^2)/(SP.κ * log10(SP.δ/(Float32(2π) * SP.corea)))
    printstyled(SP.IO,"Max timestep is $dt_max \n", color=:blue)
    return dt < dt_max
end