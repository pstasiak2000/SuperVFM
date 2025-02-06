export print_characteristics
export GetSchwarzTempCoeffs

function print_banner(SP::SimulationParams)
    ### Gets the system independent username to print
    function get_username()
        varnames = ["LOGNAME", "USER", "LNAME",  "USERNAME"]
        for varname in varnames
            haskey(ENV, varname) && return ENV[varname]
        end
        return nothing
    end
    
    hostname = gethostname()  #Get the hostname
    username = get_username() #Get the username
    date = now()              #Get the current date and time
    docsite = raw"https://pstasiak2000.github.io/SuperVFM/stable/"

    printstyled(SP.IO,"                                             \n", bold=:true, color=:blue)
    printstyled(SP.IO,"        _   _     __     __         _        \n", bold=:true, color=:blue)
    printstyled(SP.IO,"       | \\ | |_   \\ \\   / /__  _ __| |_   \n", bold=:true, color=:blue)
    printstyled(SP.IO,"       |  \\| | | | \\ \\ / / _ \\| '__| __| \n", bold=:true, color=:blue)
    printstyled(SP.IO,"       | |\\  | |_| |\\ V / (_) | |  | |_    \n", bold=:true, color=:blue)
    printstyled(SP.IO,"       |_| \\_|\\__,_| \\_/ \\___/|_|   \\__|\n", bold=:true, color=:blue)
    printstyled(SP.IO,"                                             \n", bold=:true, color=:blue)
    printstyled(SP.IO,"                                             \n", bold=:true, color=:blue)

    println(SP.IO,"user info:  $(username)@$(hostname)")
    println(SP.IO,"Launched on $(Dates.day(date))/$(Dates.month(date))/$(Dates.year(date)) @ $(Dates.Time(date))")
    print(SP.IO,"Documentation available at: ")
    printstyled(SP.IO,docsite * "\n", underline=:true)
    println(SP.IO,"--------------------------------------------------------")
end

# function print_GPU_info()
#     device = CUDA.device()  # Get the current GPU device

#     # Access properties
#     device_name = CUDA.name(device)
#     num_multiprocessors = CUDA.attribute(device, CUDA.DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT)
#     max_threads_per_mp = CUDA.attribute(device, CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR)
#     max_threads_per_block = CUDA.attribute(device, CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)
#     warp_size = CUDA.attribute(device, CUDA.DEVICE_ATTRIBUTE_WARP_SIZE)
#     println("========================================================")
#     printstyled("                   GPU Device info                      \n", bold=:true)
#     println("========================================================")
#     print("Device Name: ")
#     printstyled("$device_name\n", color=:blue, bold=:true)
#     print("Total Multiprocessors: ")
#     printstyled("$num_multiprocessors\n", color=:blue, bold=:true)
#     print("Maximum Threads per Multiprocessor: ")
#     printstyled("$max_threads_per_mp\n", color=:blue, bold=:true)
#     print("Maximum Threads per Block: ")
#     printstyled("$max_threads_per_block\n", color=:blue, bold=:true)
#     print("Warp Size: ")
#     printstyled("$warp_size\n", color=:blue, bold=:true)
#     println("========================================================")
# end

function print_boundary_info(SP::SimulationParams)
    println(SP.IO,"========================================================")
    printstyled(SP.IO,"                 Boundary Information                   \n", bold=:true)
    println(SP.IO,"========================================================")
    print_boundary(SP.IO, SP.boundary_x)
    print_boundary(SP.IO, SP.boundary_y)
    print_boundary(SP.IO, SP.boundary_z)
    println(SP.IO,"========================================================")
end


"""
    Base.show(io::IO,SimParams::SimulationParams)

Lists the simulation parameters stored in `SimParams` in a stylistic way with a key.
"""
Base.show(io::IO, SimParams::SimulationParams) = list_parameters(io,SimParams)

function list_parameters(io::IO,SimParams::SimulationParams{S,T}) where {S,T}
    printstyled(io,"                           SIMULATION PARAMETERS                                   \n", bold=:true)
    println(io,"------------------------------------------------------------------------------------")

    total_args = length(fieldnames(typeof(SimParams)))
    for arg in fieldnames(typeof(SimParams))

        
        argname = @sprintf "%20s = " arg
        argval = @sprintf "%2s" getfield(SimParams, arg)

        #Choose the colouring based on type
        if eltype(getfield(SimParams, arg)) == T
            col = :blue
            printstyled(io,"|□        ", color=col)
        elseif eltype(getfield(SimParams, arg)) == S
            col = :green
            printstyled(io,"|○        ", color=col)
        elseif typeof(getfield(SimParams, arg)) <: BoundaryType
            col = :red
            printstyled(io,"|△        ", color=col)
        elseif typeof(getfield(SimParams, arg)) <: InitCond
            col = :yellow
            printstyled(io,"|★        ", color=col)
        elseif typeof(getfield(SimParams, arg)) <: Union{VelocityMode,FilamentModel}
            col = :magenta
            printstyled(io,"|⋄        ", color=col)
        else
            col = :white
            printstyled(io,"|x        ", color=col)
        end


        printstyled(io,argname)
        printstyled(io,join([argval, "\n"]), bold=:true, color=col)

    end
    # printstyled("                                                                                     \n", underline=:true)
    println(io,"------------------------------------------------------------------------------------")
    printstyled(io,"Key:\n", italic=:true)
    printstyled(io,"   ○ $T   / Vector{$T}   / Tuple{$T}\n", italic=:true, color=:green)
    printstyled(io,"   □ $S / Vector{$S} / Tuple{$S} \n", italic=:true, color=:blue)
    printstyled(io,"   ★ Initial vortex configuration \n", italic=:true, color=:yellow)
    printstyled(io,"   △ Boundary conditions \n", italic=:true, color=:red)
    printstyled(io,"   ⋄ Velocity and filament models \n", italic=:true, color=:magenta)
    printstyled(io,"   x Other  \n", italic=:true, color=:white)
    return nothing
end



"""
    print_characteristics(io::IO,SimParams::SimulationParams)

Prints the characteristic time and length scales of the simulation in dimensional units
"""
function print_characteristics(io::IO,SimParams::SimulationParams)
    return nothing
end