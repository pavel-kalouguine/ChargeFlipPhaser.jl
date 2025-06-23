struct PhasingStatus
    phaser::Phaser
    extra::Dict
end

@enum MonitorState Created Running Paused

struct PhasingMonitor
    fig::Figure
    ops::Observable{PhasingStatus}
    c::Condition
    state::Observable{MonitorState}
end

function PhasingMonitor(phaser::Phaser)
    fig = Figure()
    ops = Observable(PhasingStatus(phaser, Dict("limits" => (0.0, 1.0))))
    c = Condition()
    state = Observable(Created)
    Label(fig[2, 1], lift(ps -> join(("$k = $v" for (k, v) in ps.extra), "\n"), ops), tellwidth=false)
    go_button = Button(fig[3, 1], label=lift(x -> if x == Created
        "START"
    elseif x == Running
        "PAUSE"
    else
        "RESUME"
    end, state), tellwidth=false)
    monitor = PhasingMonitor(fig, ops, c, state)
    on(go_button.clicks) do x
        if state[] == Created
            state[] = Running
        elseif state[] == Running
            state[] = Paused
        else
            state[] = Running
        end
        notify(c)
    end
    on(events(fig.scene).window_open) do x
        notify(c)
    end
    monitor
end

function display(pm::PhasingMonitor)
    Base.display(pm.fig)
end


struct MonitorHooks <: AbstractHooks
    pm::PhasingMonitor
end

function on_go(hooks::MonitorHooks)
    yield()
    hooks.pm.state[] == Running || wait(hooks.pm.c)
end

function on_show(hooks::MonitorHooks, phaser::Phaser, extra::Dict)
    hooks.pm.ops[] = PhasingStatus(phaser, extra)
end

function is_done(hooks::MonitorHooks)
    !events(hooks.pm.fig.scene).window_open[]
end

struct Cut2D{N}
    direction::SMatrix{2,N}
    origin::SVector{N,Float64}
    size::Tuple{Int,Int}
end

function Cut2D(dir::Matrix{Int}, orig::Vector{Float64}, sz::Tuple{Int,Int})
    M, N = size(dir)
    if M != 2
        throw(ErrorException("Cut2D: direction must be two-dimensional"))
    end
    if N != length(orig)
        throw(ErrorException("Cut2D: direction and origin must have the same dimension"))
    end
    Cut2D(SMatrix{2,N}(dir), SVector{N}(orig), sz)
end

struct Cutter
    phaser::Phaser
    cut::Cut2D
    fproj::Matrix{ComplexF64} # The projection of the amplitudes
    ρ::Matrix{Float64} # The cut of the density
    f2fproj1::SparseMatrixCSC{ComplexF64} # The sparse matrix to convert the 1D amplitudes to those of the cut
    f2fproj2::SparseMatrixCSC{ComplexF64} # Idem for the antipodes
    fproj2ρ # The irfft plan
end

alias(k::SVector{M,Int}, size::NTuple{M,Int}) where M =
    tuple((mod.(k, SVector{M,Int}(size...)) .+ 1)...)

function Cutter(phaser::Phaser, cut::Cut2D)
    m, n = cut.size
    mf = div(m, 2) + 1 # The first dimension in the frequencies domain
    fproj = zeros(ComplexF64, mf, n)
    ρ = zeros(Float64, m, n)
    cdata1 = SparseData()
    cdata2 = SparseData()
    # Loop over all wavevectors
    for k in keys(phaser.dd.k_to_bp)
        # Project the wavevector onto the cut
        i, j = alias(cut.direction * k, cut.size)
        if i > mf
            continue # Ignore the frequencies not used in irfft
        end
        # Flatten the indices
        ind = i + (j - 1) * mf
        p1 = phaser.v ⋅ k # The 1D projection of the wavevector
        p2 = -p1 # The antipode
        if p1 > 0
            push!(cdata1.irows, ind)
            push!(cdata1.icols, p1 + 1)
            push!(cdata1.vals, exp(2π * im * (k ⋅ cut.origin))) # The value 
        else
            push!(cdata2.irows, ind)
            push!(cdata2.icols, p2 + 1)
            push!(cdata2.vals, exp(-2π * im * (k ⋅ cut.origin))) # The value 
        end



    end
    f2fproj1 = sparse(cdata1.irows, cdata1.icols, cdata1.vals, n * mf, phaser.numamps + 1)
    f2fproj2 = sparse(cdata2.irows, cdata2.icols, cdata2.vals, n * mf, phaser.numamps + 1)
    fproj2ρ = plan_irfft(fproj, m)
    Cutter(phaser, cut, fproj, ρ, f2fproj1, f2fproj2, fproj2ρ)
end

function make_cut(cutter::Cutter)
    mul!(reshape(cutter.fproj, :), cutter.f2fproj1, cutter.phaser.f)
    mul!(reshape(cutter.fproj, :), cutter.f2fproj2, conj(cutter.phaser.f), 1.0, 1.0)
    mul!(cutter.ρ, cutter.fproj2ρ, cutter.fproj)
    cutter.ρ * length(cutter.fproj) # Normalize
end

function add_panel!(pm::PhasingMonitor, (i, j)::Tuple, cut::Cut2D, title::String, aspect::Real)
    phaser = pm.ops[].phaser # Get the phaser from the observable
    cutter = Cutter(phaser, cut)
    hmsig = lift(ps -> begin
        yield()
        make_cut(cutter)
    end, pm.ops)
    limsig = lift(ps -> ps.extra["limits"], pm.ops)
    panels = pm.fig[1, 1]
    ax = Axis(panels[i, j], title=title, aspect=aspect)
    hidedecorations!(ax)
    heatmap!(panels[i, j], hmsig, colormap=:jet, colorrange=limsig)
end