using ChargeFlipPhaser, StaticArrays, LinearAlgebra, GLMakie

function example_monitor()
    include("icosahedral/CdYb/load.jl")

    ff = formfactor(CdYb.dd, ball_autocorr, CdYb.composition, CdYb.B_factor)
    phaser = Phaser(CdYb.dd, ff)
    pm=PhasingMonitor(phaser)

    # Add panels to the phasing monitor    
    # Five-fold axis 
    s5=Cut2D([0 1 1 1 1 -1; 1 0 0 0 0 0], [0.,0.,0.,0.,0.,0.], (843, 377))
    add_panel!(pm, (1, 1:3), s5, "Five fold axis", sqrt(5))

    # Two-fold axis passing through the origin
    s2_1=Cut2D([1 0 0 1 0 0; 0 0 1 0 0 -1], [0.,0.,0.,0.,0.,0.], (512, 512))
    add_panel!(pm, (2, 1), s2_1, "Two fold axis through the origin", 1)

    # Two-fold axis passing through the body center
    s2_2=Cut2D([1 0 0 1 0 0; 0 0 1 0 0 -1], [0.5,0.5,0.5,0.5,0.5,0.5], (512, 512))
    add_panel!(pm, (2, 2), s2_2, "Two fold axis through the body center", 1)

    # Three-fold axis
    s3=Cut2D([1 1 1 0 0 0; 0 0 0 1 1 1], [0.0,0.0,0.0,0.0,0.0,0.0], (512, 512))
    add_panel!(pm, (2, 3), s3, "Three fold axis", 1)

    Base.display(pm.fig)    
    do_phasing!(phaser, action=action(pm))
end

example_monitor()
