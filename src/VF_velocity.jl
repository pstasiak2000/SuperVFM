abstract type VelocityMode end #Super type for computing the velocity


struct LIA <: VelocityMode
end
Adapt.@adapt_structure LIA
export LIA

#Don't compute anything if type is set to nothing
function calc_velocity!(f, pcount, ::Nothing; nthreads=1, nblocks=1)
    return nothing   
end

"""
    calc_velocity(f, pcount, ::LIA; nthreads=1, nblocks=1)

Computes the superfluid velocity using the Local Induction Approximation (LIA)

```math
\\mathbf{v}_i = \\beta \\mathbf{s}'_i \\times \\mathbf{s}''_i  
```
"""
function calc_velocity!(f, pcount, ::LIA; nthreads=1, nblocks=1)
    f_dot = get_deriv_1(f, pcount; nthreads, nblocks)
    f_ddot =  get_deriv_2(f, pcount; nthreads, nblocks)
    
    curv = CUDA.reduce(+, f_ddot.^2)
    return nothing   
end