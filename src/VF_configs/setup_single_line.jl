### Singleton for straight line

export SingleLine
struct SingleLine <: InitCond end
Adapt.@adapt_structure SingleLine


function getInitPcount(::SingleLine,SimParams::SimulationParams)
    println("--------------------------------------------------------")
    println("----------- Initialising straight line vortex ----------")
    println("--------------------------------------------------------")
    println("Changing size of pcount to fit with box_length and δ ")
    println("-: δ=$(SimParams.δ)                                     ")
    return Int32(round((2π)/(0.75*SimParams.δ)))
end

#Generate the structure
function initVortex!(f,fint,pcount,::SingleLine)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for idx ∈ index:stride:pcount
        f[idx] += @SVector [0.0,0.0,-π + (Float32(idx)-0.5)*2π/pcount]
        # f[1,idx][1] = 0.0
        # f[1,idx][2] = 0.0
        # f[1,idx][3] = -π + (Float32(idx)-0.5)*2π/pcount

        if idx == 1 #The first element
            fint[2,idx] = pcount
            fint[1,idx] = idx+1
        elseif idx == pcount #The last element
            fint[2,idx] = idx - 1
            fint[1,idx] = 1
        else 
            fint[2,idx] = idx - 1
            fint[1,idx] = idx + 1
        end
    end
    return nothing
end