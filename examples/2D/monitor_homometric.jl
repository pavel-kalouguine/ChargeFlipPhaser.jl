using ChargeFlipPhaser, GLMakie

function example_monitor()
    include("homometric.jl")

    kmax=30
    dd = Homometric.dd

    # Create a phaser with the generated diffraction data and a form factor function
    formfactors= formfactors_synthetic(dd, ball_autocorr)

    phaser = Phaser(dd, formfactors)
    pm=PhasingMonitor(phaser)

    # Add a panel to the phasing monitor     
    s=Cut2D([2 0; 0 2], [0.,0.], (1024, 1024))
    add_panel!(pm, (1, 1), s, "Homometric structures example", 1.0)

    Base.display(pm.fig)
    do_phasing!(phaser, action=action(pm), algorithm=ExperimentalAlgorithm())
end

example_monitor()