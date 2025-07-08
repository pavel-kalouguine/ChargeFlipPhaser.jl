"""
   basis_of_dense_packing(n::Int)::Matrix{Float64}
  
   Generate a basis for the lattice corresponding to dense packing of spheres of radius 1 in n-dimensional space.
   Produces the known densest packings in dimensions from 1 to 8, for higher dimensions falls back 
   to the D_n lattice.
   Parameters:
    - `n`: The dimension of the space.
   Returns:
    - An n×n Float64 matrix, where each column represents a basis vector of the lattice.
"""
function basis_of_dense_packing(n::Int)::Matrix{Float64}
    if n < 1
        throw(ArgumentError("Dimension n must be at least 1."))
    elseif n == 1
        return [2.0] # 1D case, just a single basis vector of length 2
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
  
   Generate a basis for the lattice D_n with the distance between the nearest neighbors equal to 2.
   Parameters:
    - `n`: The dimension of the lattice.
   Returns:
    - An n×n Float64 matrix, where each column represents a basis vector of the lattice D_n.
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
      
    Generate a basis for the E6, E7, or E8 lattice with the distance between the nearest neighbors equal to 2.
    Parameters:
     - `dim`: The dimension of the lattice (6, 7, or 8).
    Returns:
     - An n×n Float64 matrix, where each column represents a basis vector of the E_n lattice.
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