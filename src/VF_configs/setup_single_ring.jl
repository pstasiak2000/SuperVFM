export SingleRing

struct SingleRing{A} <: InitCond
    Radius::A
end
Adapt.@adapt_structure SingleRing


function getInitpcount(initf::SingleRing, SP::SimulationParams{S,T}) where {S,T}
    return Int64(ceil(2π*initf.Radius/(0.75*SP.δ)))
end

function printVortexBanner(initf::SingleRing,SP::SimulationParams)
    println("--------------------------------------------------------")
    println("------------- Initialising vortex ring  --------0-------")
    println("--------------------------------------------------------")
    println("Changing size of pcount to fit with box_length and δ    ")
    println("-> Radius of ring: R=$(initf.Radius)")
    println("-> δ=$(SP.δ)                                            ")
end

"""
    initVortex_kernel!(f, fint, pcount, initf::SingleRing)

Launch kernel to initialise a single ring vortex.
"""
@kernel function initVortex_kernel!(f, fint, pcount, initf::SingleRing)
    Idx = @index(Global, Linear)
    f[Idx] = @SVector [
        0.0,
        initf.Radius * cos(π * Float32((2*Idx - 1))/Float32(pcount)),
        initf.Radius * sin(π * Float32((2*Idx - 1))/Float32(pcount))
    ]
    if Idx == 1
        fint[2,Idx] = pcount
        fint[1,Idx] = Idx+1
    elseif Idx == pcount
        fint[2,Idx] = Idx - 1
        fint[1,Idx] = 1
    else
        fint[2,Idx] = Idx - 1
        fint[1,Idx] = Idx + 1
    end
end