# This file contains definitions of types used across multiples source files of the project,
# not necessarily exported in the main module.

# Abstract type for all phasing algorithms
abstract type AbstractPhasingAlgorithm end

# Utility structure for batch construction of a complex sparse matrix
struct SparseData
    irows::Vector{Int}
    icols::Vector{Int}
    vals::Vector{ComplexF64}
end

function SparseData()
    return SparseData(Vector{Int}(), Vector{Int}(), Vector{ComplexF64}())
end
