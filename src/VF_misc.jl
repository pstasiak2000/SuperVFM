function print_banner()
    hostname = gethostname() #Get the hostname
    username = ENV["USER"]   #Get the username
    date = now()             #Get the current date and time
    docsite = raw"https://pstasiak2000.github.io/SuperVFM/stable/"

    printstyled("                                             \n", bold=:true, color=:blue)
    printstyled("        _   _     __     __         _        \n", bold=:true, color=:blue)
    printstyled("       | \\ | |_   \\ \\   / /__  _ __| |_   \n", bold=:true, color=:blue)
    printstyled("       |  \\| | | | \\ \\ / / _ \\| '__| __| \n", bold=:true, color=:blue)
    printstyled("       | |\\  | |_| |\\ V / (_) | |  | |_    \n", bold=:true, color=:blue)
    printstyled("       |_| \\_|\\__,_| \\_/ \\___/|_|   \\__|\n", bold=:true, color=:blue)
    printstyled("                                             \n", bold=:true, color=:blue)
    printstyled("                                             \n", bold=:true, color=:blue)

    println("user info:  $(username)@$(hostname)")
    println("Launched on $(Dates.day(date))/$(Dates.month(date))/$(Dates.year(date)) @ $(Dates.Time(date))")
    print("Documentation available at: ")
    printstyled(docsite * "\n", underline=:true)
    println("--------------------------------------------------------")
end

function print_GPU_info()
    device = CUDA.device()  # Get the current GPU device

    # Access properties
    device_name = CUDA.name(device)
    num_multiprocessors = CUDA.attribute(device, CUDA.DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT)
    max_threads_per_mp = CUDA.attribute(device, CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR)
    max_threads_per_block = CUDA.attribute(device, CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)
    warp_size = CUDA.attribute(device, CUDA.DEVICE_ATTRIBUTE_WARP_SIZE)
    println("========================================================")
    printstyled("                   GPU Device info                      \n", bold=:true)
    println("========================================================")
    print("Device Name: ")
    printstyled("$device_name\n", color=:blue, bold=:true)
    print("Total Multiprocessors: ")
    printstyled("$num_multiprocessors\n", color=:blue, bold=:true)
    print("Maximum Threads per Multiprocessor: ")
    printstyled("$max_threads_per_mp\n", color=:blue, bold=:true)
    print("Maximum Threads per Block: ")
    printstyled("$max_threads_per_block\n", color=:blue, bold=:true)
    print("Warp Size: ")
    printstyled("$warp_size\n", color=:blue, bold=:true)
    println("========================================================")
end

function print_boundary_info(boundary_x, boundary_y, boundary_z)
    println("========================================================")
    printstyled("                 Boundary Information                   \n", bold=:true)
    println("========================================================")
    print_boundary(boundary_x)
    print_boundary(boundary_y)
    print_boundary(boundary_z)
    println("========================================================")
end


"""
    Base.show(io::IO,SimParams::SimulationParams)

Lists the simulation parameters stored in `SimParams` in a stylistic way with a key.
"""
Base.show(io::IO, SimParams::SimulationParams) = list_parameters(io,SimParams)


function list_parameters(io::IO,SimParams::SimulationParams)
    printstyled(io,"                           SIMULATION PARAMETERS                                   \n", bold=:true, underline=:true)
    # println("------------------------------------------------------------------------------------")

    total_args = length(fieldnames(typeof(SimParams)))
    for arg in fieldnames(typeof(SimParams))

        print(io,"□        ")
        argname = @sprintf "%20s = " arg
        argval = @sprintf "%2s" getfield(SimParams, arg)

        #Choose the colouring based on type
        if eltype(getfield(SimParams, arg)) == Float32
            col = :blue

        elseif eltype(getfield(SimParams, arg)) == Int32
            col = :green
        elseif typeof(getfield(SimParams, arg)) <: BoundaryType
            col = :red
        elseif typeof(getfield(SimParams, arg)) <: InitCond
            col = :yellow
        elseif typeof(getfield(SimParams, arg)) <: Union{VelocityMode,FilamentModel}
            col = :magenta
        else
            col = :white
        end


        printstyled(io,argname)
        printstyled(io,join([argval, "\n"]), bold=:true, color=col)

    end
    # printstyled("                                                                                     \n", underline=:true)
    println(io,"------------------------------------------------------------------------------------")
    printstyled(io,"Key:\n", italic=:true)
    printstyled(io,"   • Int32   / Vector{Int32}   / Tuple{Int32}\n", italic=:true, color=:green)
    printstyled(io,"   • Float32 / Vector{Float32} / Tuple{Float32} \n", italic=:true, color=:blue)
    printstyled(io,"   • Boundary conditions \n", italic=:true, color=:red)
    printstyled(io,"   • Velocity and filament models \n", italic=:true, color=:magenta)
    return nothing
end
