mutable struct RandomizedMedianFlip <: AbstractPhasingAlgorithm
    threshold::Float64 # Coefficient roughly proportional to the probability of making a phase flip
    decrement::Float64 # Decrement factor for the threshold after each iteration
end

function RandomizedMedianFlip(; threshold::Float64=1.0, decrement::Float64=0.998)
    return RandomizedMedianFlip(threshold, decrement)
end


function flip_charge!(ρ::Vector{Float64}, ::RandomizedMedianFlip)
    ρ0 = median(ρ)
    map!(x -> x < ρ0 ? -x : x, ρ, ρ) # Flip the density below the median
end


function flip_amplitudes!(wa::WorkingAmplitudes, algorithm::RandomizedMedianFlip)
    # Assign the new phases
    for i in 1:length(wa.f_r)
        f_in=wa.f_r[i]
        f_out=wa.f_r_back[i]
        df=f_out-f_in
        if abs(f_in) < 2*rand()*algorithm.threshold*abs(df)
            ϕ = real(f_out/f_in)>0 ? sign(wa.f_r[i]) : -sign(wa.f_r[i])
            wa.f_r[i]=ϕ * wa.a_r[i] 
        end
    end
    for i in 1:length(wa.f_c)
        f_in=wa.f_c[i]
        f_out=wa.f_c_back[i]
        df=f_out-f_in
        if abs(f_in) < 2*rand()*algorithm.threshold*abs(df)
            ϕ = f_out/abs(f_out)
            wa.f_c[i]=ϕ * wa.a_c[i]
        end
    end
    algorithm.threshold*= algorithm.decrement # Decrease the threshold for the next iteration
end