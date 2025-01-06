function print_info_header()
    printstyled("--var--------t--------pcount--------recon-----wall_recon---avg_d-----length--------maxu---------maxdu-------curv------removed\n", bold=:true)
    return nothing
end

function print_info(f, ghosti, ghostb, u, Empty, SimParams, pcount, it)
    if it == 0
        println()
    elseif mod(it, SimParams.shots) == 0
        itstr = @sprintf "%06d" div(it, SimParams.shots)
        t = @sprintf "%10.7f" it * SimParams.dt
        pcount_str = @sprintf "%5i" pcount
        recon_str = @sprintf "%5i" 0
        wall_recon_str = @sprintf "%5i" 0

        d = sum(norm.(f-ghosti))
        avg_d_str = @sprintf "%1.4f" d/(pcount*SimParams.δ)
        length_str = @sprintf "%4.6f" d

        u_max_str = @sprintf "%1.5f" maximum(norm.(u))
        du_max_str = @sprintf "%1.5f" 0.0f0

        curv = norm.(get_deriv_2.(f, ghosti, ghostb, Empty))
        curv_str = @sprintf "%3.2f" sum(curv)/pcount

        removed_str = @sprintf "%8i" 0
        println(itstr          * "   "*
                t              * "   "* 
                pcount_str     * "       "*
                recon_str      * "       "*
                wall_recon_str * "       "*
                avg_d_str      * "    "*
                length_str     * "    "*
                u_max_str      * "       "*
                du_max_str     * "      "*
                curv_str       * "     "*
                removed_str   
                )

    end
end