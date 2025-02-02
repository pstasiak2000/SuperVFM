function print_info_header(io::IO)
    printstyled(io,"--var--------t--------pcount--------recon-----wall_recon---avg_d-----length--------maxu---------maxdu-------curv------removed\n", bold=:true)
    return nothing
end

function print_info(u, f, f_infront, f_behind, pcount, SP::SimulationParams{S,T}, it) where {S,T}
    ghosti, ghostb = ghostp(f, f_infront, f_behind, pcount, SP) 

    pcountx = sum(f_infront .> 0)

    itstr = @sprintf "%06d" div(it, SP.shots)
    t = @sprintf "%10.7f" it * SP.dt
    pcount_str = @sprintf "%5i" pcountx
    recon_str = @sprintf "%5i" 0
    wall_recon_str = @sprintf "%5i" 0

    L = sum(get_Δξ(f, ghosti, f_infront, pcount, SP))
    avg_d_str = @sprintf "%1.4f" L / (pcountx * SP.δ)
    length_str = @sprintf "%9.6f" L

    u_max_str = @sprintf "%1.5f" maximum(norm.(u))
    du_max_str = @sprintf "%1.5f" 0.0f0

    curv = KernelAbstractions.zeros(SP.backend, T, pcount)
    get_curvature!(curv; f, f_infront, ghosti, ghostb, pcount, SP)
    curv_str = @sprintf "%3.2f" sum(curv) / pcountx

    removed_str = @sprintf "%8i" 0
    println(SP.IO,itstr * "   " *
            t * "   " *
            pcount_str * "       " *
            recon_str * "       " *
            wall_recon_str * "       " *
            avg_d_str * "    " *
            length_str * "    " *
            u_max_str * "       " *
            du_max_str * "      " *
            curv_str * "     " *
            removed_str
    )
end
