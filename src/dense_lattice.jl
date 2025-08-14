"""
    basis_of_dense_packing(n::Int)::Matrix{Float64}

Generate a basis for the lattice corresponding to the densest known packing of spheres in n-dimensional space.

# Arguments
- `n`: The dimension of the space (must be a positive integer).

# Returns
An `n×n` `Matrix{Float64}` where each column represents a basis vector of the lattice.

# Details
Produces the known densest packings in dimensions from 1 to 8. For higher dimensions, falls back to the Dₙ lattice.

The sphere packing uses spheres of radius 1 centered at the lattice points.

# Examples
```julia-repl
julia> basis_of_dense_packing(2)  # Returns basis for hexagonal packing in 2D
2×2 Matrix{Float64}:
 2.0  1.0
 0.0  1.73205
```
"""
function basis_of_dense_packing(n::Int)::Matrix{Float64}
    if n < 1
        throw(ArgumentError("Dimension n must be at least 1."))
    elseif n == 1
        return [2.0;;] # 1D case, just a single basis vector of length 2
    elseif n == 2
        return [2.0 1.0; 0.0 sqrt(3)] # 2D case, return basis vectors triangular lattice
    elseif n>= 6 && n <= 8
        return basis_e678(n) # For n = 6, 7, or
    else
        return basis_dn(n) # For all other dimensions, use the D_n lattice
    end
end

"""
    basis_dn(n::Int)::Matrix{Float64}

Generate a basis for the Dₙ lattice with nearest-neighbor distance 2.

# Arguments
- `n`: The dimension of the lattice (must be a positive integer).

# Returns
An `n×n` `Matrix{Float64}` where each column represents a basis vector of the Dₙ lattice.

# Details
The Dₙ lattice is defined such that the minimal distance between lattice points is 2. 

For `n ≥ 3`, this corresponds to the root lattice of the SO(2n) Lie algebra.

# Examples
```julia-repl
julia> ChargeFlipPhaser.basis_dn(3)  # Basis for the D₃ lattice (cubic FCC in 3D)
3×3 Matrix{Float64}:
 1.41421  1.41421  0.0
 0.0      1.41421  1.41421
 1.41421  0.0      1.41421
```
"""
function basis_dn(n::Int)::Matrix{Float64}
    M = zeros(Float64, n, n)
    for i in 1:n
        M[i, i] = 1              # Set diagonal to 1
        if i < n
            M[i, i+1] = 1        # Set superdiagonal to 1
        end
    end
    M[n, 1] = -(-1)^n                 # Set bottom-left corner to -1
    return M * sqrt(2) # Scale the basis vectors to have length 2
end


"""
    basis_e678(dim::Int)::Matrix{Float64}

Generate a basis for the E₆, E₇, or E₈ root lattice with nearest-neighbor distance 2.

# Arguments
- `dim`: The dimension of the lattice (must be 6, 7, or 8).

# Returns
A `dim×dim` `Matrix{Float64}` where each column represents a basis vector of the Eₙ lattice.

# Details
The Eₙ lattices are exceptional root lattices with minimal distance 2 between lattice points:
- **E₈**: The unique even unimodular lattice in 8D (densest sphere packing in ℝ⁸).
- **E₇**: A 7D lattice related to one of the exceptional Lie algebras.
- **E₆**: A 6D lattice related to one of the exceptional Lie algebras.
"""
function basis_e678(dim::Int)::Matrix{Float64}
    if dim < 6 || dim > 8
        error("Only dimension 6, 7 and 8 are allowed for the lattices of type E")
    end
    # The first 6, 7 and 8 columns of b8 form bases of E6, E7 and E8 respectively
    b8 = [
        -0.5 1.0 -1.0 0.0 0.0 0.0 0.0 0.0
        -0.5 1.0 1.0 -1.0 0.0 0.0 0.0 0.0
        -0.5 0.0 0.0 1.0 -1.0 0.0 0.0 0.0
        -0.5 0.0 0.0 0.0 1.0 -1.0 0.0 0.0
        -0.5 0.0 0.0 0.0 0.0 1.0 -1.0 0.0
        -0.5 0.0 0.0 0.0 0.0 0.0 1.0 -1.0
        -0.5 0.0 0.0 0.0 0.0 0.0 0.0 1.0
        -0.5 0.0 0.0 0.0 0.0 0.0 0.0 0.0
    ]
    _, R = qr(b8[:, 1:dim])
    sqrt(2) * R # Scale the basis vectors to have length 2
end