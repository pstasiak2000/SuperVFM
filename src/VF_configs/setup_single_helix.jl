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
    initVortex_kernel!(f, f_infront, f_behind, pcount, initf::SingleHelix)

Launch kernel to initialise a single he_behi_infront
"""
@kernel function initVortex_kernel!(f, f_infront, f_behind, pcount, initf::SingleHelix)
    Idx = @index(Global, Linear)

    step_KW = (initf.box_length_z)*(sqrt(((initf.A_KW/initf.b_KW)^2)+1.0f0))/Float32(pcount) 

    f[Idx] = @SVector [
        initf.A_KW*cos((Idx-1)*step_KW/(sqrt((initf.A_KW)^2 + (initf.b_KW)^2))),
        initf.A_KW*sin((Idx-1)*step_KW/(sqrt((initf.A_KW)^2 + (initf.b_KW)^2))),
        initf.b_KW*(Idx-1)*step_KW/(sqrt((initf.A_KW)^2 + (initf.b_KW)^2)) - initf.box_length_z/2
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