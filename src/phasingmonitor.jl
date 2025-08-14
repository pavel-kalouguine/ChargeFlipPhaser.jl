"""
    PhasingStatus

Container for the current status of the phasing process, used for GUI updates.

# Fields
- `phaser::Phaser`
  The phaser instance being monitored.

- `ρ::Vector{Float64}`
  The current charge density.

- `iteration::Int`
  The current iteration number.
"""
struct PhasingStatus
    phaser::Phaser
    ρ::Vector{Float64}
    iteration::Int
end

"""
    MonitorState

Enumeration representing the state of the phasing monitor GUI.

# Values
- `Created`: The monitor has been created but not started.
- `Running`: The monitor is actively running and updating.
- `Paused`: The monitor is paused and not updating.
"""
@enum MonitorState Created Running Paused


"""
    FormatLimits

Format the limits for display in the GUI.

# Arguments
- `limits::Tuple{Float64, Float64}`: The limits to format.

# Returns
- `String`: The formatted limits string.
"""
function format_limits(limits::Tuple{Float64, Float64})::String
    s=join((@sprintf("%.5g", a) for a in limits), ", ")
    "limits=($s)"
end

"""
    PhasingMonitor

GUI monitor for the phasing application.

# Fields
- `fig::Figure`  
  Makie `Figure` representing the monitor's GUI window.

- `ops::Observable{PhasingStatus}`  
  Observable triggering widget redraws whenever the phasing status changes.

- `olimits::Observable{Tuple{Float64, Float64}}`  
  Observable for the extrema of the charge density, used to update the corresponding label in the GUI.

- `c::Condition`  
  Condition object used to notify the main program about user events.

- `state::Observable{MonitorState}`  
  Current state of the monitor window (`Created`, `Running`, or `Paused`).
"""
struct PhasingMonitor
    fig::Figure
    ops::Observable{PhasingStatus}
    olimits::Observable{Tuple{Float64, Float64}}
    c::Condition
    state::Observable{MonitorState}
end

"""
    PhasingMonitor(phaser::Phaser) -> PhasingMonitor

Construct a monitor window for the phasing application associated with `phaser`.

The returned [`PhasingMonitor`](@ref) contains the basic GUI layout and controls,
but no panels for displaying 2D sections of the density.  
These panels should be added after construction using [`add_panel!`](@ref).

# Arguments
- `phaser::Phaser`: The phasing engine instance to be monitored.

# Returns
- `PhasingMonitor`: A new monitor window linked to the given `phaser`.
"""
function PhasingMonitor(phaser::Phaser)
    fig = Figure()
    ops = Observable(PhasingStatus(phaser, [0.0, 1.0], 0))
    olimits=lift(ps -> ps.phaser.numamps.*extrema(ps.ρ), ops)
    c = Condition()
    state = Observable(Created)
    Label(fig[2, 1], lift(ps -> "iteration=$(ps.iteration)", ops), tellwidth=false)
    Label(fig[3, 1], lift(limits -> format_limits(limits), olimits), tellwidth=false)
    go_button = Button(fig[4, 1], label=lift(x -> if x == Created
        "START"
    elseif x == Running
        "PAUSE"
    else
        "RESUME"
    end, state), tellwidth=false)
    monitor = PhasingMonitor(fig, ops, olimits, c, state)
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

"""
    display(pm::PhasingMonitor)

Display the GUI window of a phasing monitor.

This calls `Base.display` on the underlying Makie figure (`pm.fig`),
bringing up the monitor window for user interaction.

# Arguments
- `pm::PhasingMonitor`: The monitor whose GUI window should be displayed.
"""
function display(pm::PhasingMonitor)
    Base.display(pm.fig)
end

"""
    MonitorHooks <: AbstractHooks

Hooks for integrating the phasing workflow with a [`PhasingMonitor`](@ref) GUI.

`MonitorHooks` wraps a `PhasingMonitor` instance and can be passed to the
phasing engine so that GUI updates occur automatically at key points
in the workflow.

# Fields
- `pm::PhasingMonitor`  
  The GUI monitor window to be updated during phasing.
"""
struct MonitorHooks <: AbstractHooks
    pm::PhasingMonitor
end

function on_go(hooks::MonitorHooks)
    yield()
    hooks.pm.state[] == Running || wait(hooks.pm.c)
end

function on_show(hooks::MonitorHooks, phaser::Phaser,  ρ::Vector{Float64}, iteration::Int)
    hooks.pm.ops[] = PhasingStatus(phaser, ρ, iteration)
end

function is_done(hooks::MonitorHooks)
    !events(hooks.pm.fig.scene).window_open[]
end

"""
    Cut2D{N}

Representation of a 2D section (cut) of the charge density in an N-dimensional lattice.

# Type Parameters
- `N`: The dimension of the problem.

# Fields
- `direction::SMatrix{2,N}`  
  A 2×N matrix whose columns are the lattice vectors spanning the cutting plane.

- `origin::SVector{N,Float64}`  
  Coordinates of the origin of the cutting plane in N-dimensional space.

- `size::Tuple{Int,Int}`  
  The dimensions of the grid on which the density cut is computed.
"""
struct Cut2D{N}
    direction::SMatrix{2,N}
    origin::SVector{N,Float64}
    size::Tuple{Int,Int}
end

"""
    Cut2D(dir::Matrix{Int}, orig::Vector{Float64}, sz::Tuple{Int,Int}) -> Cut2D

Construct a `Cut2D` object from standard arrays.

# Arguments
- `dir::Matrix{Int}`: 2×N matrix whose columns span the cutting plane.
- `orig::Vector{Float64}`: Coordinates of the origin of the plane (length N).
- `sz::Tuple{Int,Int}`: Grid size on which the density cut will be computed.

# Returns
- `Cut2D{N}`: A 2D section of the charge density in N-dimensional space.

# Notes
- Throws an error if `dir` is not 2×N or if `length(orig) != N`.
"""
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

"""
    Cutter

Data structure for efficiently computing a 2D section (cut) of the charge density.

# Fields
- `phaser::Phaser`  
  The phasing engine providing the amplitudes and lattice information.

- `cut::Cut2D`  
  The 2D section specification.

- `fproj::Matrix{ComplexF64}`  
  Projection of the amplitudes onto the cutting plane.

- `ρ::Matrix{Float64}`  
  The computed density values on the 2D cut.

- `f2fproj1::SparseMatrixCSC{ComplexF64}`  
  Sparse matrix converting 1D amplitudes to the cut for the original directions.

- `f2fproj2::SparseMatrixCSC{ComplexF64}`  
  Sparse matrix for the antipodal directions.

- `fproj2ρ`  
  IRFFT plan for converting the projected amplitudes to real-space density.
"""
struct Cutter
    phaser::Phaser
    cut::Cut2D
    fproj::Matrix{ComplexF64} # The projection of the amplitudes
    ρ::Matrix{Float64} # The cut of the density
    f2fproj1::SparseMatrixCSC{ComplexF64} # The sparse matrix to convert the 1D amplitudes to those of the cut
    f2fproj2::SparseMatrixCSC{ComplexF64} # Idem for the antipodes
    fproj2ρ # The irfft plan
end

"""
    alias(k::SVector{M,Int}, size::NTuple{M,Int}) -> NTuple{M,Int}

Convert an integer vector `k` into 1-based indices suitable for FFT arrays.

# Arguments
- `k::SVector{M,Int}`: The input integer vector (e.g., a lattice index).
- `size::NTuple{M,Int}`: The dimensions of the FFT grid.

# Returns
- `NTuple{M,Int}`: 1-based indices corresponding to `k`, wrapped modulo the grid size.

# Notes
- This ensures that negative or out-of-bounds indices are correctly mapped
  into the valid FFT array range.
"""
alias(k::SVector{M,Int}, size::NTuple{M,Int}) where M =
    tuple((mod.(k, SVector{M,Int}(size...)) .+ 1)...)

"""
    Cutter(phaser::Phaser, cut::Cut2D) -> Cutter

Construct a `Cutter` object for efficiently computing a 2D section of the charge density.

# Arguments
- `phaser::Phaser`: The phasing engine providing the amplitudes and lattice information.
- `cut::Cut2D`: Specification of the 2D section (cut) to compute.

# Returns
- `Cutter`: A data structure containing precomputed projections and FFT plans
  for fast computation of the density cut.

"""    
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

"""
    add_panel!(pm::PhasingMonitor, (i, j)::Tuple, cut::Cut2D, title::String, aspect::Real)

Add a panel to the GUI monitor window to visualize a 2D section of the charge density.

The panel is placed at grid position `(i, j)` within the monitor layout and
is dynamically updated whenever the phasing status or charge-density limits change.

# Arguments
- `pm::PhasingMonitor`: The monitor window to which the panel will be added.
- `(i, j)::Tuple`: Grid coordinates (row, column) in the monitor's layout.
- `cut::Cut2D`: The specification of the 2D section to visualize.
- `title::String`: Title displayed above the panel.
- `aspect::Real`: Aspect ratio of the panel's axes.

# Notes
- The charge density cut is computed using a [`Cutter`](@ref) constructed from `pm`'s current `Phaser`.
- Panel content updates automatically via observables linked to `pm.ops` and `pm.olimits`.
- Decorations are hidden for a cleaner visual layout.
"""
function add_panel!(pm::PhasingMonitor, (i, j)::Tuple, cut::Cut2D, title::String, aspect::Real)
    phaser = pm.ops[].phaser # Get the phaser from the observable
    cutter = Cutter(phaser, cut)
    hmsig = lift(ps -> begin
        yield()
        make_cut(cutter)
    end, pm.ops)
    limsig = lift(limits -> limits, pm.olimits)
    panels = pm.fig[1, 1]
    ax = Axis(panels[i, j], title=title, aspect=aspect)
    hidedecorations!(ax)
    heatmap!(panels[i, j], hmsig, colormap=:jet, colorrange=limsig)
end