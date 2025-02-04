### Routines to calculate the special derivatives

"""
    get_deriv_1(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray)

Computes the first order derivative using a 2nd order finite difference adaptive scheme using ghost points.
"""
function get_deriv_1(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray)
    disti = norm(f - ghosti)
    distb = norm(f - ghostb)
    s_dot = ((distb * ghosti) + (disti - distb) * f - disti * ghostb) / (2 * disti * distb)
    return s_dot
end

"""
    get_deriv_2(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray)

Computes the first order derivative using a 2nd order finite difference adaptive scheme using ghost points.
"""
function get_deriv_2(f::AbstractArray, ghosti::AbstractArray, ghostb::AbstractArray)
        disti = norm(f - ghosti)
        distb = norm(f - ghostb)
        s_ddot = 2.0f0 * (ghosti/(disti*(disti+distb)) + ghostb/(distb*(disti+distb)) - f/(disti*distb))
        return s_ddot
end