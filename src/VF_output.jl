function print_info(f, SimParams, pcount, it)
    if it == 1
        println("--var--------t--------pcount--------recon-----wall_recon---avg_d-----length--------maxu---------maxdu-------curv------removed")
    end

    if mod(it, SimParams.shots) == 0
        itstr = @sprintf "%06d" div(it, SimParams.shots)
        t = @sprintf "%3.7f" it * SimParams.dt
        pcount_str = @sprintf "%5i" pcount
        recon_str = @sprintf "%5i" 0
        wall_recon_str = @sprintf "%5i" 0
        println(itstr          * "    "*
                t              * "   "* 
                pcount_str     * "       "*
                recon_str      * "       "*
                wall_recon_str * "    "    
                )

    end
end