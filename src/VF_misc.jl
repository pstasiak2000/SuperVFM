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

export GetTempCoeffs

"""
    GetTempCoeffs(Temp::Real)

Obtain the Schwarz α and α' parameters from observational data using the temperature. For values that lie within the observation data, compute a cubic spline interpolation to extract parameters for the desired temperature.
"""
function GetTempCoeffs(Temp::Real)
    @assert Temp >= 0 "Temperature cannot be negative"
    @assert Temp <= 2.17 "Temperature cannot exceed T_c"
    ObsData = [
        0.00 0.000 0.000e-00;
        1.30 0.034 1.383e-02;
        1.35 0.042 1.543E-02;
        1.40 0.051 1.668E-02;
        1.45 0.061 1.746E-02;
        1.50 0.072 1.766E-02;
        1.55 0.084 1.721E-02;
        1.60 0.097 1.608E-02;
        1.65 0.111 1.437E-02;
        1.70 0.126 1.225E-02;
        1.75 0.142 1.003E-02;
        1.80 0.160 8.211E-03;
        1.85 0.181 7.438E-03;
        1.90 0.206 8.340E-03;
        1.95 0.236 1.079E-02;
        2.00 0.279 1.198E-02;
        2.02 0.302 1.097E-02;
        2.04 0.330 8.318E-03;
        2.06 0.366 3.018E-03;
        2.08 0.414 -6.690E-03;
        2.10 0.481 -2.412E-02]
    ObsIndex = findall(i -> (i == Temp), ObsData)
    
    if length(ObsIndex) == 0 && Temp < 2.00
        @info "Cubic Spline interpolating α α'"
        interp_cubic_α1 = cubic_spline_interpolation((1.30:0.05:2.00),ObsData[2:16,2])
        interp_cubic_α2 = cubic_spline_interpolation((1.30:0.05:2.00),ObsData[2:16,3])
        
        α = (interp_cubic_α1(Temp),interp_cubic_α2(Temp))
    elseif length(ObsIndex) == 0 && Temp >= 2.00
        @info "Cubic Spline interpolating α α'"
        interp_cubic_α1 = cubic_spline_interpolation((2.00:0.02:2.10),ObsData[16:21,2])
        interp_cubic_α2 = cubic_spline_interpolation((2.00:0.05:2.10),ObsData[16:21,3])    

        α = (interp_cubic_α1(Temp),interp_cubic_α2(Temp))
    else   
    ObsIndex = ObsIndex[1][1]
    α = (ObsData[ObsIndex, 2], ObsData[ObsIndex, 3])       
    end

    return α
end
