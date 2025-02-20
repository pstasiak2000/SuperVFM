export load_VF_file

const VecN{N,T} = NTuple{N,T} where {N,T}

struct VFData{A,B,C}
    time::Float64
    num_particles::Int32  # = Np
    num_loops::Int32      # = Nl

    # Connects every particle to a neighbouring particle in the same vortex
    # loop. A "removed" particle has a value of 0.
    front::A   # [Np]

    # Indices of reference particles for each vortex loop.
    # Each vortex has a different reference particle (which can be any particle
    # inside the vortex).
    # Knowing the index of a vortex particle and the connectivities in the
    # `front` vector allows to reconstruct a whole vortex loop.
    loops::A   # [Nl]

    # Number of particles in each loop.
    num_particles_per_loop::A  # [Nl]

    # ID of each vortex.
    # Useful when a single vortex is broken into separate objects due to
    # periodicity.
    vortex_id::A  # [Nl]

    xyz::C # [3, Np]
    u::C  # [3, Np]
    vs_loc::C #[3, Np]
    vs_sup::C #[3, Np]
    # u_pll :: Matrix{Float64} #[3, Np]
    # u_n  :: Matrix{Float64}  # [3, Np]
    # f_mf :: Matrix{Float64}  # [3, Np]
    u_mf::C #[3, Np]
    v_curv::B
    # v_stretch :: Vector{Float64}
end
Adapt.@adapt_structure VFData

function read_vector(io, N, T::DataType)
    v = Array{T}(undef, N)
    read!(io, v)
end


"""
Identify a single vortex loop.

The input `n_start` indicates from which particle index to start looking for
non-visited particles.

Returns the index of a reference particle in the loop and marks all particles in
the loop as visited. Also returns the number of particles in the loop.

If no loop is found (because all particles have already been visited), returns
zero.
"""
function identify_loop!(visited, front::Vector{T},
    n_start=one(T))::Tuple{T,T} where {T}
    Np = length(front)
    @assert n_start >= one(T)
    num_particles = zero(T)  # number of particles in the loop

    # Determine reference particle.
    # Skip particles that have already been visited, as well as removed
    # particles (such that front = 0).
    nref = n_start
    while nref <= Np && (visited[nref] || front[nref] == 0)
        nref += one(T)
    end

    # No unvisited loop could be found.
    if nref > Np
        return zero(T), num_particles
    end

    # Visit all particles in the loop containing the reference particle.
    n = nref
    while !visited[n]
        visited[n] = true
        num_particles += one(T)
        n = front[n]
    end

    (nref, num_particles)::Tuple{T,T}
end


"""
Identify reference particles for different vortex loops.
"""
function identify_loops(front::Vector{T}) where {T}
    Np = length(front)
    visited = falses(Np)  # a given particle has been "visited"?
    loops = similar(front, 0)  # one particle index per vortex loop
    particles_per_loop = similar(loops)  # number of particles in each loop

    nref = zero(T)
    while true
        # Identify a single loop.
        nref, num_particles = identify_loop!(visited, front, nref + one(T))
        nref == zero(T) && break
        push!(loops, nref)
        push!(particles_per_loop, num_particles)
    end

    loops, particles_per_loop
end


function load_VF_file(io::IO)
    time = read(io, Float64)
    time = round(time, digits=3)

    num_particles = read(io, Int64)

    buf = Array{VecN{3,Float32}}(undef, num_particles)
    
    # buf = Array{SVector{3,Float32}}(undef,num_particles)
    read!(io, buf)
    xyz = buf

    # xyz += fill(VecN{3,Float32}(π,π,π),num_particles)
    front = read_vector(io, num_particles, Int32)

    loops, particles_per_loop = identify_loops(front)
    num_loops = length(loops)
    total_particles_in_loops = sum(particles_per_loop)
    # @info "Found $num_loops vortex loops ($total_particles_in_loops particles)"

    vortex_id = collect(1:num_loops)

    buf = Array{VecN{3,Float32}}(undef, num_particles)
    read!(io, buf)
    u = buf

    buf = Array{VecN{3,Float32}}(undef, num_particles)
    read!(io, buf)
    u_mf = buf

    buf = Array{VecN{3,Float32}}(undef, num_particles)
    read!(io, buf)
    vs_loc = buf

    buf = Array{VecN{3,Float32}}(undef, num_particles)
    read!(io, buf)
    vs_sup = buf
    
    
    v_curv = read_vector(io, num_particles, Float32)

    @assert eof(io)
    return VFData{Vector{Int32},Vector{Float32},Vector{VecN{3,Float32}}}(time, num_particles, num_loops, front, loops, particles_per_loop,
        vortex_id, xyz, u, vs_loc, vs_sup, u_mf, v_curv)
end


load_VF_file(filename::AbstractString) = open(load_VF_file, filename, "r")


# ### Read in the vortex tool scripts
# VortexToolScripts = readdir(joinpath(@__DIR__, "VortexTools"))
# for scripts ∈ VortexToolScripts
#     include(joinpath(@__DIR__,"VortexTools",scripts))
# end




