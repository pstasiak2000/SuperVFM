### Computes the kinetic energy spectrum from vortex data
export computeEnergySpectrum
export generateMapping

struct Anisotropic end
struct Isotropic end

export Anisotropic
export Isotropic

"""
    computeEnergySpectrum(VortexData, kMap, ::Anisotropic)

Anisotropic computation of the superfluid kinetic energy spectrum ``E(k)``. The method uses the vortex data read from file and uses kMap to loop through integer `k`-shells and the assosciated `k` vectors to each shell.
!!! warn
    The full anistropic form is very slow, and should only be used sparingly when the isotropic version does not suffice. Only use this for the time steps necessary! Consider running on GPU for vastly improved performance.

The full anistropic formula is given by:
```math
    List the formula here
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
    
    # dRdR = Kronecker.kronecker(Δξ,Δξ')
    dRdR = dot.(Δξ,Δξ')
    ∇2f = dot.(∇f,∇f')
    Δf = reshape(f, :, 1) .- reshape(f, 1, :)

    @. PREF = ∇2f * dRdR
    dEk = allocate(dev,Float32,L,L)
    kernel! = spectrumKernel!(dev,workergroupsize)

    for ik ∈ 1:eachindex(kMap)
        println("Computing: $ik")
        @time for iMap=1:kMap[ik].n_map       
            kernel!(dEk,kMap[ik].kk[iMap], PREF, Δf,ndrange=size(dEk))
            Ekin[ik] += sum(dEk) ./ (ik^2)
        end
    end
    
    return Ekin./16π^3
end

@kernel function spectrumKernel!(dEk, kVec, PREF, Δf)
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

function generateMapping(N)

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