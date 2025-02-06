export TorusKnot

struct TorusKnot{A,B} <: InitCond
    p::A
    q::A
    a_knot::B
    b_knot::B
end
Adapt.@adapt_structure TorusKnot

function getInitpcount(initf::TorusKnot, SP::SimulationParams{S,T}) where {S,T}
    ϵ = 0.01; #Accuracy of vortex length computation
    N_knot = 50; #Initial number of points

    ### Computing the vortex length to a given accuracy
    vortex_length = 0.0 
    vortex_length_new = Inf
    while abs(vortex_length_new - vortex_length) > ϵ
        vortex_length = vortex_length_new 
        N_knot *= 2; #Double the number of points each time

        t = LinRange(0,2π,N_knot)
        r1 = @. initf.a_knot * cos(initf.q* t) + initf.b_knot

        X = allocate(CPU(),Float64,N_knot,3)

        X[:,1] = @. r1 * cos(initf.p * t)
        X[:,2] = @. r1 * sin(initf.p * t)
        X[:,3] = @. initf.a_knot*sin(initf.q * t)

        ΔX = X[1:end-1,:] - circshift(X[1:end-1,:],(1,0))
        Δξ = sqrt.(sum(abs2,ΔX,dims=2))
        vortex_length_new = sum(Δξ)
    end

    step_knot = 0.75*SP.δ
    return Int64(ceil(vortex_length/step_knot))
end

function printVortexBanner(initf::TorusKnot,SP::SimulationParams)
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"------------ Initialising torus knot  ----------0-------")
    println(SP.IO,"--------------------------------------------------------")
    println(SP.IO,"Changing size of pcount to fit with box_length and δ    ")
    println(SP.IO,"Initialising a torus knot with parameters (p,q)         ")
    println(SP.IO,"-> p=$(initf.p)")
    println(SP.IO,"-> q=$(initf.q)")
    println(SP.IO,"-> a_knot=$(initf.a_knot)")
    println(SP.IO,"-> b_knot=$(initf.b_knot)")
    println(SP.IO,"-> δ=$(SP.δ)")
end

"""
    initVortex_kernel!(f, f_infront, f_behind, pcount, initf::TorusKnot)

Launch kernel to initialise a torus knot.
"""
@kernel function initVortex_kernel!(f, f_infront, f_behind, pcount, initf::TorusKnot)
    Idx = @index(Global, Linear)

    t = π * Float32(2*Idx - 1)/Float32(pcount)
    r1 = initf.a_knot*cos(initf.q*t) + initf.b_knot
    f[Idx] = @SVector [
        r1*cos(initf.p * t),
        r1*sin(initf.p * t),
        initf.a_knot*sin(initf.q*t)
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