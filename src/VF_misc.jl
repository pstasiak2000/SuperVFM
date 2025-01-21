function print_banner()
    hostname = gethostname() #Get the hostname
    username = ENV["USER"]   #Get the username
    date = now()             #Get the current date and time
    docsite = raw"https://pstasiak2000.github.io/SuperVFM/stable/"   

    printstyled("                                             \n",bold=:true,color=:blue)
    printstyled("        _   _     __     __         _        \n",bold=:true,color=:blue)
    printstyled("       | \\ | |_   \\ \\   / /__  _ __| |_   \n",bold=:true,color=:blue)
    printstyled("       |  \\| | | | \\ \\ / / _ \\| '__| __| \n",bold=:true,color=:blue)
    printstyled("       | |\\  | |_| |\\ V / (_) | |  | |_    \n",bold=:true,color=:blue)
    printstyled("       |_| \\_|\\__,_| \\_/ \\___/|_|   \\__|\n",bold=:true,color=:blue)
    printstyled("                                             \n",bold=:true,color=:blue)
    printstyled("                                             \n",bold=:true,color=:blue)
    
    println("user info:  $(username)@$(hostname)")
    println("Launched on $(Dates.day(date))/$(Dates.month(date))/$(Dates.year(date)) @ $(Dates.Time(date))")
    print("Documentation available at: "); printstyled(docsite * "\n", underline=:true)
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
    printstyled("                   GPU Device info                      \n",bold=:true)
    println("========================================================")
    print("Device Name: ");
            printstyled("$device_name\n", color=:blue, bold=:true)
    print("Total Multiprocessors: ");
            printstyled("$num_multiprocessors\n", color=:blue, bold=:true)
    print("Maximum Threads per Multiprocessor: ");
            printstyled("$max_threads_per_mp\n", color=:blue, bold=:true)
    print("Maximum Threads per Block: ");
            printstyled("$max_threads_per_block\n", color=:blue, bold=:true)
    print("Warp Size: ");
            printstyled("$warp_size\n", color=:blue, bold=:true)
    println("========================================================")
end

function print_boundary_info(boundary_x,boundary_y,boundary_z)
    println("========================================================")
    printstyled("                 Boundary Information                   \n", bold=:true)
    println("========================================================")
    print_boundary(boundary_x)
    print_boundary(boundary_y)
    print_boundary(boundary_z)
    println("========================================================")
end
