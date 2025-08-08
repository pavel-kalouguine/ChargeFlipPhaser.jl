
# Utility function returning the per-column numbers of non-zero elements in a CSC sparse matrix
nnz_per_column(m::SparseMatrixCSC)::Vector{Int} = diff(m.colptr)


struct Phaser{N}
    dd::DiffractionData
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

function Phaser(dd::DiffractionData{N}, formfactors::Vector{Float64}; small_primes=(2,3,5,7)) where {N}
    v, maxproj = find_injective_projector(dd)
    numamps = nextprod(small_primes, maxproj+1)
    @info "Using $(2*numamps) sampling points"
    real_orbits = Vector{Int}()
    complex_orbits = Vector{Int}()
    ampl = zeros(Float64, length(dd.bps))
    for i in 1:length(dd.bps)
        bp = dd.bps[i]
        ap = bp.o.aps[1] # Orbits assumed non-empty
        q = physicalnorm(ap.k, dd)
        ampl[i] = sqrt(bp.I) * formfactors[i] # The amplitude
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
    ρ = zeros(Float64, 2 * numamps) # The vector of the density
    f2ρ = plan_irfft(f, 2 * numamps)
    ρ2f = plan_rfft(ρ)
    return Phaser(dd, real_orbits, complex_orbits, v, numamps, ampl, p2f_r, p2f_c1, p2f_c2, f)
end

struct WorkingAmplitudes
    a_r::Vector{Float64} # The observed amplitudes for the real orbits 
    a_c::Vector{Float64} # The observed amplitudes for the complex orbits
    f_r::Vector{Float64} # The working amplitudes for the real orbits
    f_c::Vector{ComplexF64} # The working amplitudes for the complex orbits
    f_r_back::Vector{ComplexF64} # The backprojected amplitudes for the real orbits
    f_c_back::Vector{ComplexF64} # The backprojected amplitudes for the complex orbits
end

function WorkingAmplitudes(phaser::Phaser)
    a_r = phaser.ampl[phaser.real_orbits]
    a_c = phaser.ampl[phaser.complex_orbits]
    f_r = a_r .* rand([-1.0, 1.0], length(phaser.real_orbits)) # Set the initial signes randomly
    f_c = a_c .* exp.(2π * im * rand(length(phaser.complex_orbits))) # Set the initial phases randomly
    f_r_back = Vector{ComplexF64}(undef, length(f_r)) # Note: backprojection returns complex amplitudes
    f_c_back = similar(f_c)
    return WorkingAmplitudes(a_r, a_c, f_r, f_c, f_r_back, f_c_back)
end

function set_amplitudes!(phaser::Phaser, wa::WorkingAmplitudes)
    mul!(phaser.f, phaser.p2f_r, wa.f_r) # Apply the real phases
    mul!(phaser.f, phaser.p2f_c1, wa.f_c, 1.0, 1.0) # Apply the complex phases and accumulate
    mul!(phaser.f, phaser.p2f_c2, conj(wa.f_c), 1.0, 1.0) # Apply the conjugate complex phases and accumulate
end


# Default lifecycle behavior is to do nothing
on_go(::AbstractHooks) = nothing
on_show(::AbstractHooks, ::Phaser, ::Dict) = nothing
is_done(::AbstractHooks) = false
function on_save(::AbstractHooks, saver::AbstractSaver,
    phaser::Phaser, wa::WorkingAmplitudes)
    save_result(saver, phaser, wa)
end
save_result(::AbstractSaver, ::Phaser, ::WorkingAmplitudes) = nothing
struct DefaultHooks <: AbstractHooks end
struct DefaultSaver <: AbstractSaver end



function do_phasing!(phaser::Phaser; algorithm::AbstractPhasingAlgorithm,
    hooks::AbstractHooks=DefaultHooks(), saver::AbstractSaver=DefaultSaver(),
    max_iterations::Int=1000)
    @info "Starting phasing..."
    f = similar(phaser.f) # Working space for amplitudes, needed since FFTW.mul! modifies the input
    ρ = zeros(Float64, 2 * phaser.numamps) # The vector of the density
    # Prepare the FFT plans for the direct and inverse FFT
    f2ρ = plan_irfft(f, 2 * phaser.numamps)
    ρ2f = plan_rfft(ρ)

    mul_r = nnz_per_column(phaser.p2f_r)
    mul_c = nnz_per_column(phaser.p2f_c1) + nnz_per_column(phaser.p2f_c2)

    wa = WorkingAmplitudes(phaser) # The working amplitudes

    set_num_threads(4) # TODO: Set the number of threads automatically

    # Create adjoint sparse matrices for the backprojection
    f2p_r = phaser.p2f_r'
    f2p_c1 = phaser.p2f_c1'
    f2p_c2 = phaser.p2f_c2'


    # Enter the phasing loop
    for i = 1:max_iterations
        set_amplitudes!(phaser, wa)
        f .= phaser.f # Copy the amplitudes to the working space

        # Compute the density
        mul!(ρ, f2ρ, f) # Apply the inverse FFT

        # Call the callbacks
        on_show(hooks, phaser, Dict("iteration" => i, "limits" => extrema(ρ) .* (phaser.numamps)))
        on_go(hooks) # A potentially blocking call to wait for the user input        
        if (is_done(hooks))
            break
        end

        flip_charge!(ρ, algorithm) # Flip the charge applying the given algorithm

        # Compute new amplitudes from the modified density
        mul!(f, ρ2f, ρ) # Apply the FFT

        # Compute the backprojected amplitudes
        wa.f_r_back .= f2p_r * f ./ mul_r
        wa.f_c_back .= (f2p_c1 * f + conj(f2p_c2 * f)) ./ mul_c

        flip_amplitudes!(wa, algorithm) # Flip the amplitudes applying the given algorithm
    end
    on_save(hooks, saver, phaser, wa)

end