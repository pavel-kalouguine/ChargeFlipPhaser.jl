
function find_injective_projector(dd::DiffractionData{N,D,T}, sparseness::Real=8.0)::Tuple{SVector{N,T},T} where {N,D,T<:Integer}
    # Finds the vector v such as the scalar products of v with all wave vectors in the 
    # diffraction data are distinct
    all_k = [k for k in keys(dd.k_to_bp)]
    # Compute the covariance matrix of the wave vectors in dd
    M = MMatrix{N,N,Float64}(zeros(Float64, N, N))
    for k in all_k
        M .+= k * k'
    end
    M = Symmetric(M) # Forces M to be symmetric
    Minv = inv(M)
    rmax = maximum(k' * Minv * k for k in all_k)
    M *= rmax
    S = M^-0.5 # The scaling matrix
    maxproj = sparseness * length(all_k) # Maximal possible absolute value of the scalar product k⋅v
    println("Finding an injective projector...")
    i=0
    while true
        i+=1
        r = SVector{N}(randn(N))
        r *= (maxproj / norm(r)) # make a random vector of length maxproj
        v = T.(round.(S * r)) # Round the scaled r to the nearest integer vector
        xx = [v ⋅ k for k in all_k] # Compute the scalar products
        if length(unique(xx)) == length(xx) # Check if all scalar products are distinct
            println("Projector $v found after $i iterations")
            return v, maximum(abs.(xx))
        end
    end    
end

# Returns the formfactor used for the synthetic data
function formfactor(dd::DiffractionData, windowing_function::Function)
    q_max = maximum(physicalnorm(dd, bp.o.aps[1].k) for bp in dd.bps) # The maximum norm of the physical wavevector
    q_max *= (1.0 + 1.0 / length(dd.bps)) # Scale the q_max up not to lose the data of the peak with the biggest q
    function ff(q::Real)
        if q > q_max
            return 0.0
        else
            return windowing_function(q / q_max)
        end
    end
    return ff
end

# Returns the formfactor used for the real diffraction data for a given composition and Debye Waller
# factor
function formfactor(dd::DiffractionData, windowing_function::Function, composition::Vector{Tuple{String,R}}, b_factor::Real) where {R<:Real}
    q_max = maximum(physicalnorm(bp.o.aps[1].k, dd) for bp in dd.bps) # The maximum norm of the physical wavevector
    q_max *= (1.0 + 1.0 / length(dd.bps)) # Scale the q_max up not to lose the data of the peak with the biggest q
    atomic_formfactor = WeightedF0(composition)
    function ff(q::Real)
        if q > q_max
            return 0.0
        else
            κ = q / (4 * π) # sin(θ)/λ
            return windowing_function(q / q_max) / (atomic_formfactor(κ) * exp(-b_factor * κ^2))
        end
    end
    return ff
end


# Utility structure for batch construction of a complex sparse matrix
struct SparseData
    irows::Vector{Int}
    icols::Vector{Int}
    vals::Vector{ComplexF64}
end

function SparseData()
    return SparseData(Vector{Int}(), Vector{Int}(), Vector{ComplexF64}())
end

# Utility function returning the per-column numbers of non-zero elements in a CSC sparse matrix
nnz_per_column(m::SparseMatrixCSC)::Vector{Int} = diff(m.colptr)

struct Phaser{N}
    dd::DiffractionData
    ff::Function # The formfactor
    real_orbits::Vector{Int} # The indices of the real orbits in the diffraction data
    complex_orbits::Vector{Int} # The indices of the complex orbits in the diffraction data
    v::SVector{N,Int} # The vector of the injective projector
    numamps::Int # The number of 1D amplitudes, not counting the antipodes and the zero wavevector
    ampl::Vector{Float64} # The vector of the amplitudes, in the order of orbits in dd
    p2f_r::SparseMatrixCSC{ComplexF64,Int} # The phase factors for real orbits
    p2f_c1::SparseMatrixCSC{ComplexF64,Int} # The phase factors for complex orbits
    p2f_c2::SparseMatrixCSC{ComplexF64,Int} # The phase factors for complex orbits, conjugated
    f::Vector{ComplexF64} # The vector of phased amplitudes
end

function Phaser(dd::DiffractionData{N}, ff::Function) where {N}
    v, maxproj = find_injective_projector(dd)
    n = Int(ceil(log2(maxproj)))
    println("Using 2^$(n+1) sampling points")
    numamps = 2^n
    real_orbits = Vector{Int}()
    complex_orbits = Vector{Int}()
    ampl = zeros(Float64, length(dd.bps))
    for i in 1:length(dd.bps)
        bp = dd.bps[i]
        ap = bp.o.aps[1] # Orbits assumed non-empty
        q = physicalnorm(ap.k, dd)
        ampl[i] = sqrt(bp.I) * ff(q) # The amplitude
        if bp.o isa RealOrbit
            push!(real_orbits, i)
        else
            push!(complex_orbits, i)
        end
    end
    rdata = SparseData()
    cdata1 = SparseData()
    cdata2 = SparseData()
    for i in 1:length(real_orbits)
        o = dd.bps[real_orbits[i]].o # The orbit
        for ap in o.aps # The affine phase
            p = v ⋅ ap.k # projected wavevector
            if p > 0 # No need for antipodes here
                push!(rdata.irows, p + 1) # The row index
                push!(rdata.icols, i) # The column index
                push!(rdata.vals, exp(2π * im * (ap.ϕ))) # The value
            end
        end
    end
    for i in 1:length(complex_orbits)
        o = dd.bps[complex_orbits[i]].o # The orbit
        for ap in o.aps # The affine phase
            p1 = v ⋅ ap.k # projected wavevector
            p2 = -p1 # The antipode
            if p1 > 0
                push!(cdata1.irows, p1 + 1) # The row index
                push!(cdata1.icols, i) # The column index
                push!(cdata1.vals, exp(2π * im * (ap.ϕ))) # The value
            else
                push!(cdata2.irows, p2 + 1) # The row index
                push!(cdata2.icols, i) # The column index
                push!(cdata2.vals, exp(-2π * im * (ap.ϕ))) # The value
            end
        end
    end
    p2f_r = sparse(rdata.irows, rdata.icols, rdata.vals, numamps + 1, length(real_orbits)) # The real coefficients
    p2f_c1 = sparse(cdata1.irows, cdata1.icols, cdata1.vals, numamps + 1, length(complex_orbits)) # The complex coefficients
    p2f_c2 = sparse(cdata2.irows, cdata2.icols, cdata2.vals, numamps + 1, length(complex_orbits)) # The conjugate complex coefficients
    # Prepare the FFT plans for the direct and inverse FFT
    f = zeros(ComplexF64, numamps + 1) # The vector of phased amplitudes
    ρ = zeros(Float64, 2^(n + 1)) # The vector of the density
    f2ρ = plan_irfft(f, 2 * numamps)
    ρ2f = plan_rfft(ρ)
    return Phaser(dd, ff, real_orbits, complex_orbits, v, numamps, ampl, p2f_r, p2f_c1, p2f_c2, f)
end

const default_callbacks = Dict{String,Function}("go" => () -> nothing, "show" => (args...) -> nothing, "done" => () -> false)

function do_phasing!(phaser::Phaser; action::Dict{String,Function}=default_callbacks)
    println("Starting phasing...")
    f = copy(phaser.f) # Working with the copy since FFTW.mul! modifies the input
    aux = similar(f) # Auxiliary vector for accumulation
    ρ = zeros(Float64, 2 * phaser.numamps) # The vector of the density
    # Prepare the FFT plans for the direct and inverse FFT
    f2ρ = plan_irfft(f, 2 * phaser.numamps)
    ρ2f = plan_rfft(ρ)

    mul_r = nnz_per_column(phaser.p2f_r)
    mul_c = nnz_per_column(phaser.p2f_c1) + nnz_per_column(phaser.p2f_c2)

    a_r = phaser.ampl[phaser.real_orbits] # The amplitudes for the real orbits
    a_c = phaser.ampl[phaser.complex_orbits] # The amplitudes for the complex orbits

    # Create the real and complex phases and set their inital values randomly
    ϕ_r = rand([-1.0, 1.0], length(phaser.real_orbits))
    ϕ_c = exp.(2π * im * rand(length(phaser.complex_orbits)))

    set_num_threads(4) # TODO: Set the number of threads automatically

    # Create adjoint sparse matrices for the backprojection
    f2p_r = phaser.p2f_r'
    f2p_c1 = phaser.p2f_c1'
    f2p_c2 = phaser.p2f_c2'


    # Enter the phasing loop
    coeff = 1.0
    for i = 1:200 # TODO: pass the number of iterations
        f_r = a_r .* ϕ_r # The amplitudes for the real orbits
        f_c = a_c .* ϕ_c # The amplitudes for the complex orbits
        mul!(f, phaser.p2f_r, f_r) # Apply the real phases
        mul!(aux, phaser.p2f_c1, f_c) # Apply the complex phases
        f += aux # Add the complex phases to the real phases
        mul!(aux, phaser.p2f_c2, conj(f_c)) # Apply the conjugate complex phases
        f += aux # Add the conjugate complex phases to the real and complex phases
        phaser.f.=f # Save the copy of the phased amplitudes for viewers

        # Compute the density
        mul!(ρ, f2ρ, f) # Apply the inverse FFT

        # Call the callbacks
        action["show"](phaser, Dict("iteration"=>i, "limits"=>extrema(ρ).*(phaser.numamps)))
        action["go"]() # A potentially blocking call to wait for the user input        
        if(action["done"]()) 
            break
        end

    
        ρ0 = median(ρ)

        # Flip the density below the median
        map!(x -> x<ρ0 ? -x : x, ρ, ρ)

        # Compute new amplitudes from the modified density
        mul!(f, ρ2f, ρ) # Apply the FFT

        # Compute the backprojected amplitudes
        f_r_back = f2p_r * f ./ mul_r
        f_c_back = (f2p_c1 * f + conj(f2p_c2 * f)) ./ mul_c

        # Assign the new phases
        for i in 1:length(phaser.real_orbits)
            f_in=f_r[i]
            f_out=f_r_back[i]
            df=f_out-f_in
            if abs(f_in) < 2*rand()*coeff*abs(df)
                ϕ_r[i] = real(f_out/f_in)>0 ? ϕ_r[i] : -ϕ_r[i]
            end
        end
        for i in 1:length(phaser.complex_orbits)
            f_in=f_c[i]
            f_out=f_c_back[i]
            df=f_out-f_in
            if abs(f_in) < 2*rand()*coeff*abs(df)
                ϕ_c[i] = f_out/abs(f_out)
            end
        end
        coeff*=0.99
    end

end