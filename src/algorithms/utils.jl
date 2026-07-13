# Utilities shared by the phasing algorithms

"""
    safe_phase(z::Complex{T}) where {T}

Return the unit-modulus complex number with the same phase as `z` (`z / abs(z)`),
or `one(z)` when `z` is zero, where the phase is undefined.
"""
@inline function safe_phase(z::Complex{T})::Complex{T} where {T}
    return iszero(z) ? one(z) : z / abs(z)
end
