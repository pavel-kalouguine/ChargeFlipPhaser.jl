mutable struct SweepDown <: AbstractPhasingAlgorithm
    aux::Vector{Float64} # Auxiliairy storage for sorting the density
    fraction_flipped::Float64 # The fraction of voxels flipped
    decrement::Float64 # The decrement factor for `fraction_flipped`
end

SweepDown(; fraction_flipped=0.8, decrement=0.995) =
    SweepDown(Vector{Float64}(), fraction_flipped, decrement)

function flip_charge!(rho::Vector{Float64}, algorithm::SweepDown)
    length(algorithm.aux) == length(rho) || resize!(algorithm.aux, length(rho))
    copyto!(algorithm.aux, rho)
    sort!(algorithm.aux)
    rho0 = algorithm.aux[Int(ceil(length(rho) * algorithm.fraction_flipped))]
    map!(x -> x < rho0 ? 2 * rho0 - x : x, rho, rho)
    algorithm.fraction_flipped *= algorithm.decrement
end

function flip_amplitudes!(wa::WorkingAmplitudes, ::SweepDown)
    wa.f_r .= wa.a_r .* sign.(real.(wa.f_r_back))
    wa.f_c .= wa.a_c .* (wa.f_c_back ./ abs.(wa.f_c_back))
end