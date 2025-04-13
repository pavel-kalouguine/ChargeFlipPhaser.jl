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
    o::PhysicalOrbit{N, T}
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
    k_to_bp::Dict{SVector{N,T}, BraggPeaksOrbit{N,T}}
end

"""
   DiffractionData{N,T<:Integer}(G::SpaceGroupQuotient{N,T}, md::SMatrix{D,N,Float64})
  
   Constructor for the `DiffractionData` structure.
   Parameters:
    - `G`: A `SpaceGroupQuotient` object representing the symmetry group of the crystal
        structure.
    - `md`: A static matrix representing the metric data of the crystal structure.
   Returns a new `DiffractionData` object with an empty vector of Bragg
   peaks and the specified symmetry group.
""" 
function DiffractionData(G::SpaceGroupQuotient{N,T}, md::SMatrix{D,N,Float64}) where {N,D,T<:Integer}
    DiffractionData{N,D,T}(G, md, BraggPeaksOrbit{N,T}[], Dict{SVector{N,T}, BraggPeaksOrbit{N,T}}())
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
physicalnorm(k::SVector{N, T}, dd::DiffractionData{N,D,T}) where {N,D,T<:Integer}=norm(dd.md*k) 

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
function add_peak!(dd::DiffractionData{N, D, T}, k::SVector{N,T}, I::AbstractFloat) where {N,D,T<:Integer}
    # Create the orbit correponding to the wave vector k:
    o=make_orbit(k, dd.G)
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
    length(dd.bps)
end