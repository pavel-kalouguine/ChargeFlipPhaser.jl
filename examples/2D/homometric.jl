module Homometric
using StaticArrays, ChargeFlipPhaser, SpaceGroups

# Form-factor of three delta peaks at the vertices of an equilateral triangle of edge sqrt(3)*a
trio(a::Float64, k::SVector{2,Int}) = map(α -> exp.(2π * im * a * (cos(α) * k[1] + sin(α) * k[2])), (0:2:4) * π / 3) |> sum
# Formfactor of a well-known example of a homometric structure (9 delta peaks)
# The pattern is a convolution of two groups of three δ-peaks at the vertices of an equilateral triangle
# with the edge length of 0.2618033 and -0.1, respectively.
# The ratio of the edge lengths is 2.6180339887, which is the square of the golden ratio,
# since this number is the most poorly apporximated by rationals.
homometric(k::SVector{2,Int}) = trio(0.2618033, k) * trio(-0.1, k)

function generate_difdata(r::Real)
    G = SpaceGroupQuotient{2,Int}(Vector{SpaceGroupElement{2,Int}}()) # trivial space group
    md = SMatrix{2,2,Float64}([1.0 0.0; 0.0 1.0]) # The lattice unit cell happen to be square
    dd = DiffractionData(G, md)
    n = Int(floor(r))
    for kx = 0:n, ky = -n:n
        if (kx != 0 || ky > 0) && (kx^2 + ky^2 < r^2)
            k = SVector{2,Int}([kx, ky])
            f = homometric(k)
            i = add_peak!(dd, k, abs(f)^2)
        end
    end
    dd
end


end