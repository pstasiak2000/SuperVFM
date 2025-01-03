using LinearAlgebra
using BenchmarkTools
using Plots

struct SingleRing{A<:AbstractFloat}
    Radius::A
end

function (initf::SingleRing)(δ::T, PR) where T
    pcount = PR.int(ceil((2π * initf.Radius)/(0.75*δ)))
    f = fill(zeros(T,3), pcount)
    fint = zeros(PR.int,2,pcount)

    for id ∈ eachindex(f)
        f[id] = [
        0.0,
        initf.Radius * cos(π * Float64(2*id - 1) / pcount),
        initf.Radius * sin(π * Float64(2*id - 1) / pcount)
        ]
        if(id==1)
            fint[1,id] = id+1;
            fint[2,id] = pcount;
        elseif(id==pcount)
            fint[1,id] = 1;
            fint[2,id] = id-1;
        else
            fint[1,id] = id + 1;
            fint[2,id] = id - 1;
        end
    end
    return f, fint, pcount
end

function ghostp(f,fint)
    ghosti = similar(f)
    ghostb = similar(f)

    for id ∈ eachindex(f)
        ghosti[id] = f[fint[1,id]] 
        ghostb[id] = f[fint[2,id]]
    end
    return ghosti, ghostb
end


function get_deriv_1(f, ghosti, ghostb)
    f_dot = similar(f)

    disti = @. norm(f - ghosti)
    distb = @. norm(f - ghostb)

    @. f_dot = ((distb * ghosti) + ((disti-distb)*f) - (disti * ghostb)) / (2*disti*distb)
    return f_dot
end

function get_deriv_2(f, ghosti, ghostb)
    f_ddot = similar(f)

    disti = @. norm(f - ghosti)
    distb = @. norm(f - ghostb)

    @. f_ddot = 2 * (ghosti/(disti*(disti+distb)) - f/(disti*distb) + ghostb/(distb*(disti+distb)))
    return f_ddot
end

function timestep!(f, u, u1, u2, dt, it)
    if it==1
        @. f += dt * u
    elseif it==2
        @. f += @. 1.5f0*dt*u - 0.5f0*dt*u1
    else
        @. f += (23f0/12)*dt*u - (16f0/12f0)*dt*u1 + (5f0/12f0)*dt*u2
    end
    u2 .= u1
    u1 .= u
    return nothing
end
begin
    PRECISIONS = []
    push!(PRECISIONS,(int=Int32,fl=Float32)) #Single precision
    push!(PRECISIONS,(int=Int64,fl=Float64)) #Double precision

    plt = plot(ylabel= "x position of ring", xlabel="t")
    for PR in PRECISIONS
        δ=0.01              |> PR.fl
        κ = 1.0             |> PR.fl
        a_0 = 6.28e-7       |> PR.fl
        nsteps = 5000       |> PR.int
        dt = 1e-6       |> PR.fl

        initf = SingleRing{PR.fl}(1.0)

        dt_max = PR.fl(((δ/2.0)^2)/(κ * log10(δ/(2π * a_0))))
        @assert dt < dt_max "Timestep is too large dt=$dt: \n dt needs to be smaller than $dt_max"
        @info "Condition satisfied"  
        typeof(dt_max)

        f, fint, pcount = initf(δ, PR)
        
        f_init = f;
        u = similar(f)
        u1 = similar(f)
        u2 = similar(f)

        #Time and x position of the ring
        x_pos = zeros(PR.fl,2,nsteps)
        t = 0 |> PR.fl

        for it ∈ 1:nsteps

            ghosti, ghostb = ghostp(f,fint)

            f_dot = get_deriv_1(f, ghosti, ghostb)
            f_ddot = get_deriv_2(f, ghosti, ghostb)

            curv = @. 1/norm(f_ddot)
            beta = @. PR.fl.((κ/(4π)*log10(4.6*curv/a_0)))
            @. u = beta * cross(f_dot, f_ddot)

            timestep!(f, u, u1, u2, dt, it)
            x_pos[1,it] = it*dt
            x_pos[2,it] = u[1][1]
            
            # println(typeof(f))
        end
        plot!(plt, x_pos[1,:],x_pos[2,:],label="$(PR.fl)",linewidth=1)
        display(plt)
    end
end
"""
begin
    δ=0.01          
    κ = 1.0
    a_0 = 6.28e-7
    nsteps = 1
    dt = 0.000005

    initf = SingleRing(0.1)

    dt_max = ((δ/2.0)^2)/(κ * log10(δ/(2π * a_0)))
    @assert dt < dt_max "Timestep is too large dt=$dt: \n dt needs to be smaller than $dt_max"
    @info "Condition satisfied"

    f, fint, pcount = initf(δ)

    f_init = f;
    u = similar(f)
    u1 = similar(f)
    u2 = similar(f)

    x_pos = zeros(Float64,2,nsteps)

    t = 0;
    anim = @animate for it ∈ 1:nsteps

        ghosti, ghostb = ghostp(f,fint)

        f_dot = get_deriv_1(f, ghosti, ghostb)
        f_ddot = get_deriv_2(f, ghosti, ghostb)

        curv = @. 1/norm(f_ddot)
        beta = @. κ/(4π)*log10(4.6*curv/a_0)

        u .= @. beta * cross(f_dot, f_ddot)

        timestep!(f, u, u1, u2, dt, it)
        x_pos[1,it] = it*dt
        x_pos[2,it] = f[1][1]
        
        println(it)

        plot(Tuple.(f), xlim=(-π,π), xlabel=("x"),
        ylim=(-π,π), ylabel=("y"),
        zlim=(-π,π), zlabel=("z"))
    end every 1000
    gif(anim, "anim_fps30.gif", fps=60)
# plot(x_pos[1,:],x_pos[2,:])
end

begin
    plot(Tuple.(f), xlim=(-π,π), xlabel=("x"),
    ylim=(-π,π), ylabel=("y"),
    zlim=(-π,π), zlabel=("z"))

    plot!(Tuple.(f_init), xlim=(-π,π), xlabel=("x"),
    ylim=(-π,π), ylabel=("y"),
    zlim=(-π,π), zlabel=("z"))
  
end
"""