export SingleHelix

struct SingleHelix{A} <: InitCond
    A_KW::A #Amplitude of the Kelvin wave
    b_KW::A #Wavenumber of the Kelvin wave
    box_length_z::A #Length in vertical z-direction
end
Adapt.@adapt_structure SingleHelix


function getInitpcount(initf::SingleHelix, SP::SimulationParams{S,T}) where {S,T}
    @assert SP.boundary_z == PeriodicBoundary(3) "Periodic boundary conditions in z required"
    return Int64(ceil(initf.box_length_z*(sqrt(((initf.A_KW/initf.b_KW)^2)+1.0))/(0.75*SP.δ)))
end


function printVortexBanner(initf::SingleHelix,SP::SimulationParams)
    local k = Int32(1/initf.b_KW)
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"----------- Initialising single helix  -----------------")
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"Changing size of pcount to fit with box_length and δ "   )
    println(SP.IO,"-> Amplitude of wave: A/2π=$(initf.A_KW)")
    println(SP.IO,"-> Wavenumber k=$(k)")
    println(SP.IO,"-> δ=$(SP.δ)         ")
end



"""
    initVortex_kernel!(f, fint, pcount, initf::SingleHelix)

Launch kernel to initialise a single ring vortex.
"""
@kernel function initVortex_kernel!(f, fint, pcount, initf::SingleHelix)
    Idx = @index(Global, Linear)

    step_KW = (initf.box_length_z)*(sqrt(((initf.A_KW/initf.b_KW)^2)+1.0f0))/Float32(pcount) 

    f[Idx] = @SVector [
        initf.A_KW*cos((Idx-1)*step_KW/(sqrt((initf.A_KW)^2 + (initf.b_KW)^2))),
        initf.A_KW*sin((Idx-1)*step_KW/(sqrt((initf.A_KW)^2 + (initf.b_KW)^2))),
        initf.b_KW*(Idx-1)*step_KW/(sqrt((initf.A_KW)^2 + (initf.b_KW)^2)) - initf.box_length_z/2
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