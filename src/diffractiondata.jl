"""
   BraggPeaksOrbit{N,T<:Integer}
  
   A structure representing an orbit of Bragg peaks in a diffraction pattern.
   Type parameters:
    - `N`: The number of dimensions (e.g., 3 for 3D diffraction patterns).
    - `T`: The type of the integer used for indexing (e.g., `Int`).
   The `BraggPeaksOrbit` structure is used to store information about a specific
   Bragg peak in a diffraction pattern. It contains the orbit of wave vectors
   containing the Bragg peak, as well as the intensity of the peak.
   The structure contains the following fields:
   - `o`: A `PhysicalOrbit` object representing the orbit of wave vectors
        associated with the Bragg peaks.
   - `I`: The intensity of the Bragg peak.
   
"""
struct BraggPeaksOrbit{N,T<:Integer}
    o::PhysicalOrbit{N,T}
    I::AbstractFloat
end

"""
   DiffractionData{N,D,T<:Integer}
  
   A structure representing a diffraction pattern in N-dimensional space.
   Type parameters:
    - `N`: The number of dimensions of the reciprocal lattice.
    - `D`: The number of dimensions of the real space. For the periodic crystals, 
        `D` is equal to `N`.
    - `T`: The type of the integer used for indexing (e.g., `Int`).
   The `DiffractionData` structure is used to store information about a
   diffraction pattern, including the Bragg peaks and their intensities.
   The structure contains the following fields:
    - `G`: A `SpaceGroupQuotient` object representing the symmetry group of the
        structure.
    - `md`: The metric data of the crystal structure, represented by a DxN static 
        matrix. Given the integer `N`-dimensional vector `k`, the value of the corresponding
        diffraction wavevector is `md*k`. For the real physical data, the units of `md` are Å^-1
    - `bps`: A vector of `BraggPeaksOrbit` objects representing the orbits of Bragg peaks
        in the diffraction pattern.
    - `k_to_bp`: A dictionary mapping wave vectors to their corresponding
        `BraggPeaksOrbit` objects. This field is used to efficiently check if a
        wave vector is already present in the diffraction pattern and to retrieve
        the associated Bragg peak orbit when needed. It ensures that each wave
        vector is uniquely associated with a single Bragg peak orbit.
"""
struct DiffractionData{N,D,T<:Integer}
    G::SpaceGroupQuotient{N,T}
    md::SMatrix{D,N,Float64}
    bps::Vector{BraggPeaksOrbit{N,T}}
    k_to_bp::Dict{SVector{N,T},BraggPeaksOrbit{N,T}}
end

"""
   DiffractionData{N,T<:Integer}(G::SpaceGroupQuotient{N,T}, md::SMatrix{D,N,Float64})
  
   Constructor for the `DiffractionData` structure.
   Parameters:
    - `G`: A `SpaceGroupQuotient` object representing the symmetry group of the crystal
        structure.
    - `md`: A static matrix representing the metric data of the crystal structure. The matrix
        should be of size DxN, where D is the number of dimensions in real space and N is the
        order of the module of the Bragg peak vectors (the dimension of the space group). The 
        physical wavevector `q` corresponding to the integer vector of indices `k` is given by 
        `q = md * k`. 
   Returns a new `DiffractionData` object with an empty vector of Bragg
   peaks and the specified symmetry group.
"""
function DiffractionData(G::SpaceGroupQuotient{N,T}, md::SMatrix{D,N,Float64}) where {N,D,T<:Integer}
    DiffractionData{N,D,T}(G, md, BraggPeaksOrbit{N,T}[], Dict{SVector{N,T},BraggPeaksOrbit{N,T}}())
end

"""
    metric_data_inconsistency(dd::DiffractionData{N,D,T})::Float64 where {N,D,T<:Integer}
    Check for inconsistencies in the metric data of the `DiffractionData` object.
    Parameters:
    - `dd`: The `DiffractionData` object to check.
    For the metric data `md` to be consistent with the symmetry group, the matrix `r=md'*md`
    must be invariant with respect to the action of the symmetry group `G`. Namely, for every 
    element g∈G, the following condition must hold:
    `g.a'*r*g.a = r`
    This function computes the matrix `r` and returns the tolerance to which this condition is satisfied.
    Returns:
    - The ratio of the maximal absolute value of of the elements of the matrix `g.a'*r*g.a - r` to the 
    trace of `r`. This value should be close to zero for the metric data to be consistent with the symmetry group.
"""
function metric_data_inconsistency(dd::DiffractionData{N,D,T}) where {N,D,T<:Integer}
    r = dd.md' * dd.md
    max_inconsistency = 0.0
    trace_r = tr(r)

    for g in dd.G
        rg = g.a' * r * g.a
        inconsistency = maximum(abs.(rg .- r))
        max_inconsistency = max(max_inconsistency, inconsistency)
    end

    max_inconsistency / trace_r
end

"""
   physicalnorm(k::SVector{N, T}, dd::DiffractionData{N,D,T}) where {N,D,T<:Integer}
  
   Compute the norm of a physical wave vector corresponding to the integer vector `k`using 
   the metric data of the `DiffractionData` object. 
   Parameters:
    - `k`: An integer wave vector (a vector of Miller indices).
    - `dd`: The `DiffractionData` object containing the metric data.
   Returns:
    - The physical norm of the physical wave vector corresponding to `k`, calculated as `norm(dd.md*k)`.
"""
physicalnorm(k::SVector{N,T}, dd::DiffractionData{N,D,T}) where {N,D,T<:Integer} = norm(dd.md * k)

"""
   add_peak!(dd::DiffractionData{N,D,T}, k::SVector{N,T}, I::AbstractFloat)
  
   Add a new orbit of Bragg peaks to the `DiffractionData` object.
   Parameters:
    - `dd`: The `DiffractionData` object to which the Bragg peak will be added.
    - `k`: A vector representing the wave vector of the Bragg peak.
    - `I`: The intensity of the Bragg peak.
   This function creates a new `BraggPeaksOrbit` object 
    and adds it to the `DiffractionData` object.

   If the wave vector belongs to an extinct orbit, the function will
    throw an `ArgumentError`.

   Return value:
    - the number of Bragg peaks in the added orbit. If the wavevector is already
    taken into account, the function will return 0.     
"""
function add_peak!(dd::DiffractionData{N,D,T}, k::SVector{N,T}, I::AbstractFloat)::Int where {N,D,T<:Integer}
    # Create the orbit correponding to the wave vector k:
    o = make_orbit(k, dd.G)
    # Check if the orbit is extinct
    if o isa ExtinctOrbit
        throw(ArgumentError("The wavevector $k belongs to an extinct orbit."))
    end

    # Check if the wave vector is already present in the diffraction data
    if haskey(dd.k_to_bp, k)
        return 0
    end

    # Create a new Bragg peak orbit
    bp = BraggPeaksOrbit(o, I)

    # Add the wavevectors from the orbit to the dictionary
    for ap in o.aps
        dd.k_to_bp[ap.k] = bp
        # Add the antipode if the orbit is of complex type
        if o isa ComplexOrbit
            dd.k_to_bp[-ap.k] = bp
        end
    end
    # Add the Bragg peak to the DiffractionData object
    push!(dd.bps, bp)
    n=length(o.aps)
    o isa ComplexOrbit ? 2*n : n
end




function find_injective_projector(dd::DiffractionData{N,D,T}; num_candidates=10)::Tuple{SVector{N,T},T} where {N,D,T<:Integer}
    all_k = [k for k in keys(dd.k_to_bp)]
    # Compute the covariance matrix of the wave vectors in dd
    M = MMatrix{N,N,Float64}(zeros(Float64, N, N))
    for k in all_k
        #bpo = dd.k_to_bp[k] # the orbit of the Bragg peak
        #M .+= k * k' * bpo.I # peacks are weighted by their intensity
        M .+= k * k'
    end
    M = Symmetric(M) # Forces M to be symmetric
    Minv = inv(M)
    rmax = maximum(k' * Minv * k for k in all_k)
    M *= rmax
    scaling = M^0.5 # The scaling matrix

    B0 = basis_of_dense_packing(N) # The basis of the dense packing lattice

    @info "Finding an injective projector..."
    i = 0
    candidates=Vector{Tuple{SVector{N,T},T}}()
    while true
        i += 1
        R = svd(randn(N, N)).U # Random orthogonal matrix
        B = R * B0
        L = (Int.(round.(scaling * B)))' # The rows of L for the basis of a superlattice
        s = nothing
        try
            s = snf(L)
        catch e
            continue # SNF may fail because of integer overflow. Try again with a different random matrix
        end
        P=L*s.V 
        if !all(iszero,P[1:(end-1), end])
            continue # The projector may be non-injective due to rounding of L in the wrong way, try again
        end
        v = SVector{N, Int}(s.V[:, end]) # Take the last column of the right unimodular factor as the projector
        xx = [v ⋅ k for k in all_k] # Compute the scalar products
        if length(unique(xx)) == length(xx) # Check if all scalar products are distinct
            push!(candidates, (v, maximum(abs.(xx))))
            if length(candidates) >= num_candidates
                break # We have enough candidates, stop the search
            end
        end
    end
    # Find the candidate with the minimal maximal scalar product
    v, maxprod = argmin(x -> x[2], candidates) # Take the candidate with the minimal maximal scalar product
    @info "Projector $v found after $i iterations"
    @info "Maximal scalar product: $maxprod"
    return (v, maxprod)

end

function formfactors_synthetic(dd::DiffractionData, windowing_function::Function)
    q_max = maximum(physicalnorm(bp.o.aps[1].k, dd) for bp in dd.bps) # The maximum norm of the physical wavevector
    q_max *= (1.0 + 1.0 / length(dd.bps)) # Scale the q_max up not to lose the data of the peak with the biggest q
    ff = ones(Float64, length(dd.bps)) # The form factor is a constant function
    for i in 1:length(dd.bps)
        bp = dd.bps[i]
        ap = bp.o.aps[1] # Take the first vector of the orbit
        q = physicalnorm(ap.k, dd)
        ff[i] = windowing_function(q / q_max)
    end
    ff
end