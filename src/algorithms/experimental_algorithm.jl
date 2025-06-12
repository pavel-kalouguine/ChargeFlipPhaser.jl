mutable struct ExperimentalAlgorithm <: AbstractPhasingAlgorithm
end

function flip_charge!(ρ::Vector{Float64}, algorithm::ExperimentalAlgorithm)
    ρ0 = median(ρ)
    for i in 1:length(ρ)
        if ρ[i] < ρ0
            ρ[i] += 4*rand()*(ρ0-ρ[i])
        end
    end
end

function flip_amplitudes!(wa::WorkingAmplitudes, algorithm::ExperimentalAlgorithm)
    for i in 1:length(wa.f_r)
        wa.f_r[i] = wa.a_r[i] * sign(real(wa.f_r_back[i]))
    end
    for i in 1:length(wa.f_c)
        wa.f_c[i] = wa.a_c[i] * (wa.f_c_back[i]/abs(wa.f_c_back[i]))
    end
end