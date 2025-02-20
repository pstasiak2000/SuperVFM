### All routines that alter the geometry of the vortex filaments/flux tubes should be contained within this file
### The main routines here insert and remove particles to maintain a roughly constant resolution along filaments

# function pinsert!(; kwargs...)
#     (; pcount) = (; kwargs...)
#     (; f, f_infront, f_behind, ghosti, ghostb) = (; kwargs...)
#     (; u, u1, u2) = (; kwargs...)
#     (; SP) = (; kwargs...)
#     Δξ = get_Δξ(f, ghosti, f_infront, pcount)
    

# end

# @kernel function pinsert_kernel!(f, f_infront, f_behind, ghosti, ghostb, u, u1, u2, δ)

# end

"""
    premove()
"""
function premove!(; kwargs...)
    (; pcount) = (; kwargs...)
    (; f, f_infront, f_behind, ghosti, ghostb) = (; kwargs...)
    (; u, u1, u2) = (; kwargs...)
    (; SP) = (; kwargs...)

    removed = 0
    pre_remove = sum(f_infront .== 0)
    
    remove_vector = KernelAbstractions.zeros(SP.backend,Int32,pcount) .+ Int32(1)
   
    pcount_active = sum(f_infront .> 0)

    kernel! = premove_kernel!(SP.backend,SP.workergroupsize)
    kernel!(remove_vector, f, f_infront, f_behind, u, u1, u2, SP.δ, ndrange=pcount)
    KernelAbstractions.synchronize(SP.backend)

    post_remove = sum(f_infront .== 0)
    removed = post_remove - pre_remove
    return removed
end

@kernel function premove_kernel!(remove_vector, f, f_infront, f_behind, u, u1, u2, δ)
    Idx = @index(Global, Linear)
    if f_infront[Idx] != 0
        infront = f_infront[Idx]
        distii = norm(f[Idx] - f[f_infront[infront]])
        if (distii < 0.99 * δ)
            tinfront = f_infront[f_infront[Idx]]
            f_behind[tinfront] = Idx
            f_infront[Idx] = tinfront
            remove_vector[infront] = 0
        end
    end

    @synchronize

    if remove_vector[Idx] == 0
        clear_particle(f, f_infront, f_behind, u, u1, u2, Idx)
    end
end


function clear_particle(f, f_infront, f_behind, u, u1, u2, i)
    f[i] *= 0.0
    f_infront[i] = 0
    f_behind[i] = 0
    u[i] *= 0.0
    u1[i] *= 0.0
    u2[i] *= 0.0
end