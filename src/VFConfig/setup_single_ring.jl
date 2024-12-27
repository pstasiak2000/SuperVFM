export SingleRing

struct SingleRing{A} <: InitCond
    Radius::A
end
Adapt.@adapt_structure SingleRing

#Obtains initial number of vortex points
function getInitPcount(IC::SingleRing,δ)
    println("-----------------------------------------------------")
    println("-------- Initialising straight line vortex ----------")
    println("-----------------------------------------------------")
    println("Changing size of pcount to fit with box_length and δ ")
    println("-: δ=$δ                                              ")
    tmp = ceil(2π*IC.Radius/(0.75*δ))
    return Int32(tmp)
end

#Initialises the vortex configuration
function init!(f,fint,pcount,IC::SingleRing)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    for idx ∈ index:stride:pcount
        f[1,idx] = 0.0
        f[2,idx] = IC.Radius * sin(π * Float32((2*idx - 1))/Float32(pcount))
        f[3,idx] = IC.Radius * cos(π * Float32((2*idx - 1))/Float32(pcount))

        if index == 1 #The first element
            fint[2,idx] = pcount
            fint[1,idx] = idx+1
        elseif index == pcount #The last element
            fint[2,idx] = idx - 1
            fint[1,idx] = 1
        else 
            fint[2,idx] = idx - 1
            fint[1,idx] = idx + 1
        end
    end
    return nothing
end