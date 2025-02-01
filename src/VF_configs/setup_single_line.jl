### Singleton for straight line

export SingleLine

struct SingleLine <: InitCond end
Adapt.@adapt_structure SingleLine


function getInitpcount(::SingleLine, SP::SimulationParams{S,T}) where {S,T}
    return Int64(round((SP.box_size[3])/(0.75*SP.δ)))
end

function printVortexBanner(::SingleLine,SP::SimulationParams)
    println("--------------------------------------------------------")
    println("----------- Initialising straight line vortex ----------")
    println("--------------------------------------------------------")
    println("   Changing size of pcount to fit with box_length and δ ")
    println("-> δ=$(SP.δ)                                            ")
end


@kernel function initVortex_kernel!(f, fint, pcount, ::SingleLine)
    Idx = @index(Global, Linear)
    f[Idx] = @SVector [0.0,0.0, -π + (Float32(Idx)-0.5)*2π/Float32(pcount)]

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