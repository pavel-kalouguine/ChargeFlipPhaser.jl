
"""
    F0WaasKirf

Callable object representing the non-dispersive part ``f_0`` of the atomic
scattering factor for an element or ion. The representation follows 
the formulas from this article:
New Analytical Scattering Factor Functions for Free Atoms and Ions for Free 
Atoms and Ions, D. Waasmaier & A. Kirfel, Acta Cryst. (1995). A51, 416-413                      

The non-dispersive part ``f_0`` of the atomic scattering factor is a	    
function of the selected element and of ``\\kappa=\\sin(\\theta)/\\lambda``, 
where ``\\lambda`` is the photon wavelengh and ``\\theta`` is incident angle.		    
This function can be approximated by a function:			    

```math
f_0(\\kappa) = c + \\sum_{i=1}^5 a_i \\exp(-b_i*(\\kappa^2)) 				    
```

# Fields
- `sym`: The chemical symbol of the element or ion (e.g., "Fe", "Fe2+", "O2-")
- `z`: The atomic number of the element or ion (e.g., 26 for Fe)
- `a`: The vector of coefficients for the exponential terms in the scattering factor
- `b`: The vector of coefficients for the exponential decay in the scattering factor
- `c`: The constant term in the scattering factor

The object should be called with a single real argument ``κ = sin(θ) / λ``.
"""
struct F0WaasKirf
    sym::String
    z::Int
    a::SVector{5, Float64}
    b::SVector{5, Float64}
    c::Float64
end


(f::F0WaasKirf)(κ::Float64)=f.c+sum(f.a .* exp.(-κ^2*f.b))

"""
    WeightedF0

A callable object representing the weighted non-dispersive part of the atomic 
scattering factor for a material of a given composition.

# Fields
- `v`: A vector of tuples, where each tuple contains the `F0WaasKirf` object for 
an atom or ion and the fraction of this type of atoms in the material.

# Constructors
- `WeightedF0(c::Vector{Tuple{String, Float64}})`: Creates a `WeightedF0` object 
from a vector of tuples, where each tuple contains the chemical symbol of the 
atom or ion and its fraction in the material.

# Example
The object should be called with a single real argument ``κ = sin(θ) / λ``.
```julia-repl
julia> w=WeightedF0([("As", 1.0), ("Ga", 1.0)]); [w(κ) for κ in 0.0:0.2:1.0]
6-element Vector{Float64}:
 32.00163144733964
 25.58541327156792
 19.004924495418177
 13.806263061220022
 10.184215575459087
  8.03425166926575
```
"""
struct WeightedF0
    v::Vector{Tuple{F0WaasKirf,Float64}}
    function WeightedF0(c::Vector{Tuple{String, Float64}})
        n = sum(t[2] for t in c)
        new([(wktab[t[1]], t[2]/n) for t in c])
    end
end

(w::WeightedF0)(κ::Float64)::Float64 = sqrt(sum((t[1](κ))^2*t[2] for t in w.v))

waaskirf_filepath() = joinpath(@__DIR__, "..", "data", "f0_WaasKirf.dat")

function load_waaskirf()::Dict{String, F0WaasKirf}
    t=readlines(waaskirf_filepath())
    tab=Dict{String, F0WaasKirf}()
    for i=70:4:910
        h=split(t[i])
        z=parse(Int, h[2])
        sym=h[3]
        t[i+1]=="#N 11" || throw(ArgumentError("wrong format"))
        t[i+2]=="#L a1  a2  a3  a4  a5  c  b1  b2  b3  b4  b5" || throw(ArgumentError("wrong format"))
        v=[parse(Float64, x) for x in split(t[i+3])]
        a=v[1:5]
        c=v[6]
        b=v[7:11]
        tab[sym]=F0WaasKirf(sym, z, a, b, c)
    end
    tab
end

const wktab=load_waaskirf()
