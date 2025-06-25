
struct CSVSaver <: AbstractSaver
    output_file_path::String
end

function save_result(saver::CSVSaver, phaser::Phaser, wa::WorkingAmplitudes)
    bps = phaser.dd.bps
    # Prepare the amplitudes
    ampl = zeros(ComplexF64, length(bps))
    for i in eachindex(wa.f_r)
        j = phaser.real_orbits[i]
        ϕ = bps[j].o.aps[1].ϕ
        ampl[j] = wa.f_r[i] * exp(2π * im * ϕ)
    end
    for i in eachindex(wa.f_c)
        j = phaser.complex_orbits[i]
        ϕ = bps[j].o.aps[1].ϕ
        ampl[j] = wa.f_c[i] * exp(2π * im * ϕ)
    end
    # Prepare the header
    dim = length(phaser.v)
    ks_header = join(("k$i" for i in 1:dim), ", ")
    header = "$ks_header, I, real(ampl), imag(ampl)\n"
    open(saver.output_file_path, "w") do io
        # Write header
        write(io, header)
        # Write data
        for j in eachindex(bps)
            I = bps[j].I
            a=ampl[j]
            k = bps[j].o.aps[1].k
            ks = join((x for x in k), ", ")
            values=join((@sprintf("%.5g", x) for x in (I, real(a), imag(a))), ", ")
            write(io, "$ks, $values\n")
        end
    end

    @info "Results saved in $(saver.output_file_path)"
end