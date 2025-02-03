export SingleRing

struct SingleRing{A} <: InitCond
    Radius::A
end
Adapt.@adapt_structure SingleRing


function getInitpcount(initf::SingleRing, SP::SimulationParams{S,T}) where {S,T}
    return Int64(ceil(2π*initf.Radius/(0.75*SP.δ)))
end

function printVortexBanner(initf::SingleRing,SP::SimulationParams)
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"------------- Initialising vortex ring  ----------------")
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"Changing size of pcount to fit with box_length and δ    ")
    println(SP.IO,"-> Radius of ring: R=$(initf.Radius)")
    println(SP.IO,"-> δ=$(SP.δ)                                            ")
end

"""
    initVortex_kernel!(f, f_infront, f_behind, pcount, initf::SingleRing)

Launch kernel to initialise a single ring vortex.
"""
@kernel function initVortex_kernel!(f, f_infront, f_behind, pcount, initf::SingleRing)
    Idx = @index(Global, Linear)
    f[Idx] = @SVector [
        0.0,
        initf.Radius * cos(π * Float32((2*Idx - 1))/Float32(pcount)),
        initf.Radius * sin(π * Float32((2*Idx - 1))/Float32(pcount))
    ]
    if Idx == 1
        f_behind[Idx] = pcount
        f_infront[Idx] = Idx+1
    elseif Idx == pcount
        f_behind[Idx] = Idx - 1
        f_infront[Idx] = 1
    else
        f_behind[Idx] = Idx - 1
        f_infront[Idx] = Idx + 1
    end
end