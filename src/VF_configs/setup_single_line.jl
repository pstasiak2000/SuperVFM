### Singleton for straight line

export SingleLine

struct SingleLine <: InitCond end
Adapt.@adapt_structure SingleLine


function getInitpcount(::SingleLine, SP::SimulationParams{S,T}) where {S,T}
    return Int64(round((SP.box_size[3])/(0.75*SP.δ)))
end

function printVortexBanner(::SingleLine,SP::SimulationParams)
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"----------- Initialising straight line vortex ----------")
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"   Changing size of pcount to fit with box_length and δ ")
    println(SP.IO,"-> δ=$(SP.δ)                                            ")
end

"""
    initVortex_kernel!(f, f_infront, f_behind, pcount, ::SingleLine)

Kernel call to initialise a single vortex line.
"""
@kernel function initVortex_kernel!(f, f_infront, f_behind, pcount, ::SingleLine)
    Idx = @index(Global, Linear)
    f[Idx] = @SVector [0.0,0.0, -π + (Float32(Idx)-0.5)*2π/Float32(pcount)]

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