module ChargeFlipPhaser
import StaticArrays: SMatrix, SVector, MMatrix
import SparseArrays: sparse, nnz, SparseMatrixCSC
import LinearAlgebra: norm, transpose, inv, ⋅, Symmetric, tr, mul!, axpy!, det, svd
import SpaceGroups: RealOrbit, ComplexOrbit, PhysicalOrbit, ExtinctOrbit, SpaceGroupQuotient, make_orbit
import FFTW: plan_irfft, plan_rfft, mul!, set_num_threads, irfft
import Statistics: median
import NormalForms: snf
using Makie, GLMakie

export WeightedF0, DiffractionData, add_peak!, find_injective_projector,
    metric_data_inconsistency, physicalnorm, formfactor, PhasedData, 
    do_phasing!, ball_autocorr, Phaser,
    PhasingMonitor, Cut2D, add_panel!,
    AbstractHooks, DefaultHooks, MonitorHooks, display,
    AbstractPhasingAlgorithm, formfactors_synthetic

include("types.jl")
include("f0_waaskirf.jl")
include("diffractiondata.jl")
include("windowing.jl")
include("phaser.jl")
include("phasingmonitor.jl")
include("algorithms.jl")

end # module PhaserTmp
