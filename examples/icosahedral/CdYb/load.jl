module CdYb
using DataFrames, CSV, StaticArrays, LinearAlgebra
using ChargeFlipPhaser
include(joinpath(@__DIR__, "..", "icosahedron.jl"))
using .Icosahedron

const datafilename = "cdyb_integrated_bck"
const datafilepath = joinpath(@__DIR__, datafilename)
const a = 5.689 # Icosahedral lattice parameter, Angstrems
const composition = [("Cd", 5.7), ("Yb", 1.0)]
const B_factor = 1.8 # Debye-Waller exponent, Angstrem^2
G = PI # The symmetry group

# Functions used for the sanity check
function kpar(k::SVector{6, Int})
    base = Icosahedron.epar * (π / a)
    base * k
end
function kper(k::SVector{6, Int})
    base = Icosahedron.eper * (π / a)
    base * k
end

dd = DiffractionData(G, Icosahedron.epar * (π / a))
# Check the consistency of the metric data
inconsistency = metric_data_inconsistency(dd)
println("Metric data inconsistency: $inconsistency")

t = CSV.File(datafilepath; header=7, skipto=9, delim=' ', ignorerepeated=true)
dt = DataFrame(t)
ratio_q = Vector{Float64}()
for r in eachrow(dt)
    k = SVector{6,Int}(Int(r.n1), Int(r.n2), Int(r.n3), Int(r.n4), Int(r.n5), Int(r.n6))
    I = r.I
    dI = r[Symbol("Sigma(I)")]
    if I > 0.0 && I > dI
        n=add_peak!(dd, k, I)
        if n==0
            println("Warning: Peak with the wavevector $k already exists in the diffraction data.")
        end
    end
    # Sanity check:
    qpar = r.Qpar * (sqrt(2) * π / a)
    push!(ratio_q, qpar / norm(kpar(k)))
end
println("Sanity check: $(minimum(ratio_q)) <= q_par/k_par <= $(maximum(ratio_q))")

# Compute the formfactors to apply to the diffraction data
atomic_formfactor = WeightedF0(composition)
formfactors=ones(Float64, length(dd.bps))
# First, divide by the average atomic formfactor and the Debye-Waller factor to sharpen the atomic surfaces
for (i, bp) in enumerate(dd.bps)
    q = physicalnorm(bp.o.aps[1].k, dd)
    κ = q / (4 * π) # sin(θ)/λ
    formfactors[i] = 1.0/(atomic_formfactor(κ) * exp(-B_factor * κ^2))
end
# Then, apply the windoning function in the parallel and in the perp space
windowing_function=ball_autocorr
lim_kpar = maximum(norm(kpar(bp.o.aps[1].k)) for bp in dd.bps)*(1.0 + 1.0 / length(dd.bps))
lim_kper = maximum(norm(kper(bp.o.aps[1].k)) for bp in dd.bps)*(1.0 + 1.0 / length(dd.bps))
for (i, bp) in enumerate(dd.bps)
    formfactors[i] *= windowing_function(norm(kpar(bp.o.aps[1].k)) / lim_kpar)
    formfactors[i] *= windowing_function(norm(kper(bp.o.aps[1].k)) / lim_kper)
end
end