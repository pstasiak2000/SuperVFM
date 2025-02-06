"""
    save_vortex(it;kwargs...)

Save the vortex and related information to file.
"""
function save_vortex(it;kwargs...)
    (; pcount, t) = (; kwargs...)
    (; f, f_infront, curv) = (; kwargs...)
    (; u, u_mf, u_loc, u_sup) = (; kwargs...)

    itstr = @sprintf "%06d" it
    open(joinpath(base_dir,"OUTPUTS","VFdata","var." * itstr * ".log"),"w") do io
        write(io, t)
        write(io, pcount)
        write(io, Array(f))
        write(io, Array(f_infront))
        write(io, Array(u))
        write(io, Array(u_mf))
        write(io, Array(u_loc))
        write(io, Array(u_sup))
        write(io, Array(curv))
    end   
end


function print_info_header(io::IO)
    header_string = "--var--------t--------pcount--------recon-----wall_recon---avg_d-----length--------maxu---------maxdu-------curv------removed\n"

    ### Print to file
    open(joinpath(base_dir,"OUTPUTS","VFdata","ts.log"),"w") do file
        printstyled(file,header_string, bold=:true)
    end

    ### Print to simulation buffer (by default stdout)
    printstyled(io,header_string, bold=:true)
    return nothing
end

function print_info(it, SP::SimulationParams{S,T};kwargs...) where {S,T}
    (; pcount) = (; kwargs...)
    (; f, f_infront, f_behind, curv) = (; kwargs...)
    (; u) = (; kwargs...)

    ghosti, ghostb = ghostp(f, f_infront, f_behind, pcount, SP) 

    pcountx = sum(f_infront .> 0)

    itstr = @sprintf "%06d" div(it, SP.shots)
    t = @sprintf "%11.7f" it * SP.dt
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

    output_string = itstr * "   " * t * "   " * pcount_str * "       " * recon_str * "       " * wall_recon_str * "       " * avg_d_str * "    " * length_str * "    " * u_max_str * "       " * du_max_str * "      " * curv_str * "     " * removed_str

    ### Print to file
    open(joinpath(base_dir,"OUTPUTS","VFdata","ts.log"),"a") do file
        println(file,output_string)
    end

    ### Print to simulation buffer (by default stdout)
    println(SP.IO,output_string)
    return curv
end

"""
    generate_initial_files(SP::SimulationParams)

Generates the output folder structure and the initialisation files.
"""
function generate_initial_files(SP::SimulationParams)
    check_output_folder_structure(SP.IO)
    create_info_file(SP)
    
    open(joinpath(base_dir,"OUTPUTS","parameterVF.txt"),"w") do io
        show(io,SP)
    end
end 

"""
    create_info_file(::SimulationParams{S,T}) where {S,T}

Create file to print the simulation precision of floating point and integers. To be used by vortex reading methods for correct parsing.
"""
function create_info_file(::SimulationParams{S,T}) where {S,T}
    filename = "precision.info"
    open(joinpath(base_dir,"OUTPUTS",filename),"w") do io
        println(io,T) #Write integer precision
        println(io,S) #Write floating point precision
    end
end