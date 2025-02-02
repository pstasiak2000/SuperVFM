export SingleRing

struct SingleRing{A} <: InitCond
    Radius::A
end
Adapt.@adapt_structure SingleRing

#Obtains initial number of vortex points
function getInitPcount(initf::SingleRing,SimParams::SimulationParams)
    println("--------------------------------------------------------")
    println("----------- Initialising vortex ring  ------------------")
    println("--------------------------------------------------------")
    println("Changing size of pcount to fit with box_length and δ ")
    println("Radius of ring: R=$(initf.Radius)")
    println("-: δ=$(SimParams.δ)                                     ")
    tmp = ceil(2π*initf.Radius/(0.75*SimParams.δ))
    return Int32(tmp)
end

#Initialises the vortex configuration
function initVortex!(f,fint,pcount,initf::SingleRing)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for idx ∈ index:stride:pcount
        f[idx] += @SVector [
            0.0,
            initf.Radius * cos(π * Float32((2*idx - 1))/Float32(pcount)),
            initf.Radius * sin(π * Float32((2*idx - 1))/Float32(pcount))
        ]
        # f[1,idx] = 0.0
        # f[2,idx] = initf.Radius * sin(π * Float32((2*idx - 1))/Float32(pcount))
        # f[3,idx] = initf.Radius * cos(π * Float32((2*idx - 1))/Float32(pcount))

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