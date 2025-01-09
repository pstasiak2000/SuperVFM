export SingleHelix

struct SingleHelix{A} <: InitCond
    A_KW::A #Amplitude of the Kelvin wave
    b_KW::A #Wavenumber of the Kelvin wave
    box_length_z::A #Length in vertical z-direction
end
Adapt.@adapt_structure SingleHelix

function SingleHelix(A_KW,b_KW,box_size)
    return SingleHelix{Float32}(A_KW,b_KW,box_size[3])
end

#Obtains initial number of vortex points
function getInitPcount(initf::SingleHelix,SimParams::SimulationParams)
    println("--------------------------------------------------------")
    println("----------- Initialising single helix  -----------------")
    println("--------------------------------------------------------")
    println("Changing size of pcount to fit with box_length and δ ")
    println("Amplitude of wave: A/2π=$(initf.A_KW)")
    println("Wavenumber $(initf.b_KW)")
    println("-: δ=$(SimParams.δ)                                     ")
    @assert SimParams.boundary_z == PeriodicBoundary(3) "Periodic boundary conditions in z required"
    tmp = ceil(initf.box_length_z*(sqrt(((initf.A_KW/initf.b_KW)^2)+1.0f0))/(0.75*SimParams.δ))
    return Int32(tmp)
end

#Initialises the vortex configuration
function initVortex!(f,fint,pcount,initf::SingleHelix)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x

    step_KW = (initf.box_length_z)*(sqrt(((initf.A_KW/initf.b_KW)^2)+1.0f0))/pcount

    for idx ∈ index:stride:pcount
        f[idx] += @SVector [
        initf.A_KW*cos((idx-1)*step_KW/(sqrt((initf.A_KW)^2 + (initf.b_KW)^2))),
        initf.A_KW*sin((idx-1)*step_KW/(sqrt((initf.A_KW)^2 + (initf.b_KW)^2))),
        initf.b_KW*(idx-1)*step_KW/(sqrt((initf.A_KW)^2 + (initf.b_KW)^2)) - initf.box_length_z/2
        ]

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