function banner_print()
    hostname = gethostname() #Get the hostname
    username = ENV["USER"]   #Get the username
    date = now()             #Get the current date and time

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
    println("---------------------------------------------")
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
    println("                   GPU Device info                      ")
    println("========================================================")
    println("Device Name: $device_name")
    println("Total Multiprocessors: $num_multiprocessors")
    println("Maximum Threads per Multiprocessor: $max_threads_per_mp")
    println("Maximum Threads per Block: $max_threads_per_block")
    println("Warp Size: $warp_size")
    println("========================================================")
end
