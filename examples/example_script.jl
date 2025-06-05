using ChargeFlipPhaser, StaticArrays, LinearAlgebra

include("icosahedral/CdYb/load.jl")


phaser=Phaser(CdYb.dd, CdYb.formfactors)
callbacks = Dict{String,Function}("go" => () -> nothing, "show" => (phaser, extra) -> println(extra), "done" => () -> false)
do_phasing!(phaser, action=callbacks)
