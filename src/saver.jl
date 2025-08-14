
"""
    CSVSaver <: AbstractSaver

Type for saving results in CSV format.

# Fields
- `output_file_path`: The path to the output CSV file where results will be saved.
"""
struct CSVSaver <: AbstractSaver
    output_file_path::String
end

"""
    save_result(saver::CSVSaver, phaser::Phaser, wa::WorkingAmplitudes)

Save the results to a CSV file.

The output file has the following columns:
- `k1`, `k2`, ..., `kn`: The reflection indices.
- `I`: The peak intensity (the value is taken from the DiffractionData object).
- `real(ampl)`: The real part of the amplitude.
- `imag(ampl)`: The imaginary part of the amplitude.

Note that the absolute value of the amplitude is not necessarity equal to 
the square root of the intensity, because of the form-factors used in the phasing.

# Arguments
- `saver::CSVSaver`: The CSV saver instance.
- `phaser::Phaser`: The phaser instance containing the diffraction data.
- `wa::WorkingAmplitudes`: The current working amplitudes to save.
"""
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