using ChargeFlipPhaser, StaticArrays, LinearAlgebra

include("icosahedral/CdYb/load.jl")


ff=formfactor(CdYb.dd, ball_autocorr, CdYb.composition, CdYb.B_factor)
phaser=Phaser(CdYb.dd, ff)
callbacks = Dict{String,Function}("go" => () -> nothing, "show" => (phaser, extra) -> println(extra), "done" => () -> false)
do_phasing!(phaser, action=callbacks)
