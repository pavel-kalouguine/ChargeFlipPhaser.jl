# This file contains definitions of types used across multiples source files of the project,
# not necessarily exported in the main module.

"""
    AbstractPhasingAlgorithm

Abstract supertype for all phasing algorithms in `ChargeFlipPhaser`.

To define a new algorithm, create a subtype:

    struct MyAlgorithm <: AbstractPhasingAlgorithm
        # fields
    end

and implement the following methods:

- `flip_charge!(rho::Vector{Float64}, algorithm::MyAlgorithm)`
- `flip_amplitudes!(wa::WorkingAmplitudes, algorithm::MyAlgorithm)`

These methods are part of the `ChargeFlipPhaser` API. They should be
either imported before definition:
```
    import ChargeFlipPhaser: flip_charge!, flip_amplitudes!
```
or defined with the module prefix:
```
    function ChargeFlipPhaser.flip_charge!(...)
```
"""
abstract type AbstractPhasingAlgorithm end

"""
    AbstractHooks

Abstract supertype for hooks used in the phasing process.
To define a new hooks type, create a subtype:

    struct MyHooks <: AbstractHooks
        # fields
    end

and implement the following methods:
- `on_go(hooks::MyHooks)`
- `on_show(hooks::MyHooks, phaser::Phaser,  ρ::Vector{Float64}, iteration::Int)`
- `is_done(hooks::MyHooks)`

These methods are part of the `ChargeFlipPhaser` API. They should be
either imported before definition:
```
    import ChargeFlipPhaser: on_go, on_show, is_done
```
or defined with the module prefix:
```
    function ChargeFlipPhaser.on_show(...)
```
"""
abstract type AbstractHooks end

"""
    AbstractSaver

Abstract supertype for saving results.
To define a new saver type, create a subtype:

    struct MySaver <: AbstractSaver
        # fields
    end

and implement the following method:
- `function save_result(saver::MySaver, phaser::Phaser, wa::WorkingAmplitudes)`

This method is part of the `ChargeFlipPhaser` API. It should be
either imported before definition:
```
    import ChargeFlipPhaser: save_result
```
or defined with the module prefix:
```
    function ChargeFlipPhaser.save_result(...)
```
"""
abstract type AbstractSaver end

# Utility structure for batch construction of a complex sparse matrix
struct SparseData
    irows::Vector{Int}
    icols::Vector{Int}
    vals::Vector{ComplexF64}
end

function SparseData()
    return SparseData(Vector{Int}(), Vector{Int}(), Vector{ComplexF64}())
end
