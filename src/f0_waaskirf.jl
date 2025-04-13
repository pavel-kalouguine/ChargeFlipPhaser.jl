struct F0WaasKirf
    sym::String
    z::Int
    a::SVector{5, Float64}
    b::SVector{5, Float64}
    c
end

(f::F0WaasKirf)(k::Float64)=f.c+sum(f.a .* exp.(-k^2*f.b))

struct WeightedF0
    v::Vector{Tuple{F0WaasKirf,Float64}}
    function WeightedF0(c::Vector{Tuple{String, Float64}})
        n = sum(t[2] for t in c)
        new([(wktab[t[1]], t[2]/n) for t in c])
    end
end

(w::WeightedF0)(k::Float64)= sqrt(sum((t[1](k))^2*t[2] for t in w.v))

waaskirf_filepath() = joinpath(@__DIR__, "..", "data", "f0_WaasKirf.dat")

function load_waaskirf()
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
