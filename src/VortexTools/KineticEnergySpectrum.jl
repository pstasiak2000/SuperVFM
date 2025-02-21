### Computes the kinetic energy spectrum from vortex data
export computeEnergySpectrum
export generateMapping
export Anisotropic
export Isotropic

struct Anisotropic end
struct Isotropic end



"""
    computeEnergySpectrum(VortexData, kMap, ::Isotropic)

Computes the superfluid kinetic energy spectrum ``E(k)`` using the isotropic approximation formula:

```math
E(k) = \\frac{\\rho_s \\kappa^2}{4\\pi^2}\\int_{\\mathcal{L}_1}\\int_{\\mathcal{L}_2}d\\xi_1 d\\xi_2 \\mathbf{s}'(\\xi_1)\\cdot\\mathbf{s}'(\\xi_2)\\frac{\\sin\\left(k|\\mathbf{s}(\\xi_2) -\\mathbf{s}(\\xi_1)| \\right)}{k|\\mathbf{s}(\\xi_2) -\\mathbf{s}(\\xi_1)| \\right)},
```
such that the total superfluid kinetic energy ``E`` is 

```math
E = \\int_{0}^{\\infty}E(k)dk.
```
"""
function computeEnergySpectrum(VortexData, kMap, ::Isotropic)
    Nk = length(kMap)

    dev = get_backend(VortexData.xyz)
    workergroupsize = 1024
    Ekin = KernelAbstractions.zeros(CPU(), Float32, Nk)

    ### Compute f_behind for ghost filaments
    T = eltype(eltype(VortexData.xyz))
    L = length(VortexData.xyz)

    f = allocate(dev,SVector{3,T},L)
    f = convert(typeof(f),VortexData.xyz)
    f_infront = VortexData.front
    f_behind = allocate(dev, Int32, L)

    kernel! = get_f_behind_from_infront_kernel!(dev, workergroupsize)
    kernel!(f_behind, f_infront, ndrange=L)

    ghosti = allocate(dev, SVector{3,T}, L)
    ghostb = allocate(dev, SVector{3,T}, L)

    ghostp!(ghosti, ghostb; f, f_infront, f_behind, dev, workergroupsize)

    ∇f = get_deriv_1.(f, ghosti, ghostb)
    Δξ = norm.(f - ghosti)

    PREF = KernelAbstractions.zeros(dev, T, L, L)
    dRdR = dot.(Δξ,Δξ')
    ∇2f = dot.(∇f,∇f')
    Δf = norm.(reshape(f, :, 1) .- reshape(f, 1, :)) .+ eps(Float32)

    @. PREF = ∇2f * dRdR
    dEk = allocate(dev,Float32,L,L)

    for ik ∈ eachindex(kMap)
        @. dEk = PREF * sin(ik*Δf)/(ik*Δf)
        Ekin[ik] = sum(dEk) / (4π^2) 
    end
    return Ekin
end


"""
    computeEnergySpectrum(VortexData, kMap, ::Anisotropic)

Full anisotropic computation of the superfluid kinetic energy spectrum ``E(k)``. The method uses the vortex data read from file and uses kMap to loop through integer `k`-shells and the assosciated `k` vectors to each shell and saves the final array as a 1D k vector.
!!! warning
    The full anistropic form is very slow, and should only be used sparingly when the isotropic version does not suffice. Only use this for the time steps necessary! Consider running on GPU for vastly improved performance.

The full anisotropic formula for the 3D energy spectrum ``E(\\mathbf{k})``is given by:

```math
E(\\mathbf{k}) = \\frac{\\rho_s \\kappa^2}{16\\pi^3}\\frac{1}{k^2}\\int_{\\mathcal{L}_1}\\int_{\\mathcal{L}_2}d\\xi_1 d\\xi_2 \\mathbf{s}'(\\xi_1)\\cdot\\mathbf{s}'(\\xi_2) e^{i\\mathbf{k}\\cdot\\left[\\mathbf{s}(\\xi_2) - \\mathbf{s}(\\xi_1) \\right]}.
```
The 1D spectrum is computed by integration over concentric shells

```math
E(k) = \\sum\\limits_{k-1<\\left\\lVert \\mathbf{k}\\right\\rVert}\\leq k E(\\mathbf{k}),
```

such that the total superfluid kinetic energy ``E`` is

```math
E = \\int_{\\hat{\\Omega}}E(\\mathbf{k})d^3\\mathbf{k} = \\int_{0}^{\\infty}E(k)dk.
```
"""
function computeEnergySpectrum(VortexData, kMap, ::Anisotropic)
    Nk = length(kMap)

    dev = get_backend(VortexData.xyz)
    workergroupsize = 1024
    Ekin = KernelAbstractions.zeros(CPU(), Float32, Nk)

    ### Compute f_behind for ghost filaments
    T = eltype(eltype(VortexData.xyz))
    L = length(VortexData.xyz)

    f = allocate(dev,SVector{3,T},L)
    f = convert(typeof(f),VortexData.xyz)
    f_infront = VortexData.front
    f_behind = allocate(dev, Int32, L)

    kernel! = get_f_behind_from_infront_kernel!(dev, workergroupsize)
    kernel!(f_behind, f_infront, ndrange=L)

    ghosti = allocate(dev, SVector{3,T}, L)
    ghostb = allocate(dev, SVector{3,T}, L)

    ghostp!(ghosti, ghostb; f, f_infront, f_behind, dev, workergroupsize)

    ∇f = get_deriv_1.(f, ghosti, ghostb)
    Δξ = norm.(f - ghosti)

    PREF = KernelAbstractions.zeros(dev, T, L, L)
    
    dRdR = dot.(Δξ,Δξ')
    ∇2f = dot.(∇f,∇f')
    Δf = reshape(f, :, 1) .- reshape(f, 1, :)

    @. PREF = ∇2f * dRdR
    dEk = allocate(dev,Float32,L,L)
    kernel! = anisotropic_spectrumKernel!(dev,workergroupsize)

    for ik ∈ eachindex(kMap)
        println("Computing: $ik")
        @time for iMap=1:kMap[ik].n_map       
            kernel!(dEk,kMap[ik].kk[iMap], PREF, Δf,ndrange=size(dEk))
            Ekin[ik] += sum(dEk) ./ (ik^2)
        end
    end
    
    return Ekin./16π^3
end

@kernel function anisotropic_spectrumKernel!(dEk, kVec, PREF, Δf)
    Idx, Idy = @index(Global,NTuple)
    @inbounds dEk[Idx,Idy] = @fastmath real(exp(1im * dot(kVec,Δf[Idx,Idy])))
    @inbounds dEk[Idx,Idy] *= PREF[Idx,Idy]
end

struct ComboMapping{A,B}
    ik::A
    n_map::A
    kk::B
end
Adapt.@adapt_structure ComboMapping

"""
    generateMapping(N::Int)

Generate a vector of combination k-maps. A combination map with `id=n` contains all of the 3D k-vectors ``(k_x,k_y,k_z)`` which reside within two consective integer spherical shells of size `n`.

The mapping partititions the vector as a set ``\\Omega = \\lbrace \\Omega_k : 1\\leq k \\leq M = N/2\\rbrace`` where ``N`` is the resolution and 
```math
\\Omega_k = \\lbrace \\mathbf{k} = (k_x,k_y,k_z) : k-1 < ||\\mathbf{k}||_{2} \\leq k \\rbrace
```
"""
function generateMapping(N::Int)
    kMap = [];
    if (N == 128)
        max_kmap = 50586
    elseif (N == 256)
        max_kmap = 205502
    elseif (N = 512)
        max_kmap = 827402
    end

    n_maps = zeros(Int32, div(N, 2))
    kk = zeros(Int32, 3, max_kmap, div(N, 2))

    kSpan = -div(N, 2):div(N, 2)
    for ikz ∈ kSpan
        @inbounds for iky ∈ kSpan, ikx ∈ kSpan
            ik = Int32(ceil(sqrt(ikx^2 + iky^2 + ikz^2)))
            if (ik > 0 && ik <= div(N, 2))
                n_maps[ik] += 1
                kk[1, n_maps[ik], ik] = ikx
                kk[2, n_maps[ik], ik] = iky
                kk[3, n_maps[ik], ik] = ikz
            end
        end
    end

    for ik ∈ 1:div(N,2)
        slice = reinterpret(SVector{3,Int32},kk[:,1:n_maps[ik],ik])'
        push!(kMap,ComboMapping(Int32(ik),n_maps[ik],slice))
    end

    return kMap
end


@kernel function get_f_behind_from_infront_kernel!(f_behind, f_infront)
    Idx = @index(Global, Linear)
    if f_infront[Idx] != 0
        f_behind[f_infront[Idx]] = Idx
    end
end