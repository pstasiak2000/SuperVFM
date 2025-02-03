export SimpleTrefoil

struct SimpleTrefoil{A} <: InitCond
    scale::A
end
Adapt.@adapt_structure SimpleTrefoil


function getInitpcount(initf::SimpleTrefoil, SP::SimulationParams{S,T}) where {S,T}
    trefoil_length = 28.8621 #Numerically computed for scale=1
    line_length = initf.scale * trefoil_length
    return Int64(ceil(line_length/(0.75*SP.δ)))
end

function printVortexBanner(initf::SimpleTrefoil,SP::SimulationParams)
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"------------- Initialising vortex ring  ----------------")
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"Changing size of pcount to fit with box_length and δ    ")
    println(SP.IO,"Drawing a simple trefoil scaled from blueprint")
    println(SP.IO,"-> scale=$(initf.scale)")
    println(SP.IO,"-> δ=$(SP.δ)                                            ")
end

"""
    initVortex_kernel!(f, f_infront, f_behind, pcount, initf::SimpleTrefoil)

Launch kernel to initialise a simple trefoil.
"""
@kernel function initVortex_kernel!(f, f_infront, f_behind, pcount, initf::SimpleTrefoil)
    Idx = @index(Global, Linear)

    t = π * Float32(2*Idx - 1)/Float32(pcount)
    f[Idx] = @SVector [
        initf.scale * (cos(t) + 2.0*cos(2.0*t)),
        initf.scale * (sin(t) - 2.0*sin(2.0*t)),
        initf.scale * (-1.0*sin(3.0*t))
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