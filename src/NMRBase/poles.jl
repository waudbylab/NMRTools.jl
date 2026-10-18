# Pole models of spectra ###########################################################################
#
# A frequency-domain spectrum is modelled as a sum of damped complex exponentials (Lorentzian lines)
# by band-wise harmonic inversion: the spectrum is split into overlapping bands, each band is
# inverse Fourier transformed to a short decimated time signal, and HSVD (`_hankelpoles`,
# `_amplitudes` in harmonic.jl) is run on it. Each pole gives a line's frequency, linewidth, complex
# amplitude (integral and phase) and height. Splitting into bands keeps the model order and SVD size
# small however many lines the spectrum contains.

"""
    Pole

One damped complex exponential (a Lorentzian line) from a harmonic-inversion model of a spectrum,
as returned by [`estimatepoles`](@ref).

# Fields
- `shift`: chemical shift (ppm)
- `linewidth`: full width at half height (Hz)
- `amplitude`: complex amplitude. `abs(amplitude)` is the line integral (spectrum intensity × points)
  and `angle(amplitude)` its phase, which is zero for a positive absorption line
- `height`: absorption peak height (spectrum intensity units)
- `snr`: `height` divided by the noise level
"""
struct Pole
    shift::Float64
    linewidth::Float64
    amplitude::ComplexF64
    height::Float64
    snr::Float64
end

"phase of a pole (radians); zero for a positive absorption line"
Base.angle(p::Pole) = angle(p.amplitude)

function Base.show(io::IO, p::Pole)
    return print(io, "Pole(", round(p.shift; digits=4), " ppm, ", round(p.linewidth; digits=2),
                 " Hz, phase ", round(rad2deg(angle(p)); digits=1), "°, SNR ",
                 round(p.snr; sigdigits=3), ")")
end

"""
    _referencetrace(spec, di) -> Vector{ComplexF64}

Complex trace along array dimension `di` of `spec`, in the imaginary unit given by that dimension's
`:mcindex`. Other imaginary units take their real component. Other dimensions are reduced so that
every signal contributes: a time dimension contributes its first point, any other dimension is
summed (for a frequency dimension the sum is proportional to its first time point).
"""
function _referencetrace(spec::NMRData, di::Integer)
    A = parent(spec)
    mcindex = metadata(dims(spec, di), :mcindex)
    N = MulticomplexNumbers.order(first(A))
    Ac = copy(A)
    w = MulticomplexNumbers.unsafe_ascomplex!(Ac, mcindex)   # N-1 component axes, then data axes
    v = w[ntuple(i -> i <= N - 1 ? 1 : Colon(), ndims(w))...]
    v = collect(ComplexF64, v)
    for (j, d) in enumerate(dims(spec))
        j == di && continue
        v = d isa TimeDimension ? selectdim(v, j, 1:1) : sum(v; dims=j)
    end
    return vec(collect(v))
end

"""
    _tracenoise(s; blocksize=64, q=0.2) -> Float64

rms noise (per real component) of a complex trace. The trace is split into blocks (at most
`blocksize` points, fewer for short traces, so that there are at least 16), a straight line
is removed from the real and imaginary parts of each block, and the noise is taken from a low
quantile `q` of the block standard deviations, so that blocks containing signal (including the
dispersive tails of an unphased spectrum) do not inflate it.
"""
function _tracenoise(s::AbstractVector{<:Complex}; blocksize::Integer=64, q::Real=0.2)
    n = length(s)
    blocksize = clamp(n ÷ 16, 8, blocksize)          # short traces (indirect dimensions): smaller blocks
    nblocks = n ÷ blocksize
    nblocks >= 4 || return std(real.(s))
    x = collect(1.0:blocksize) .- (blocksize + 1) / 2
    sxx = sum(abs2, x)
    σ = map(1:nblocks) do b
        y = view(s, ((b - 1) * blocksize + 1):(b * blocksize))
        r = y .- mean(y) .- x .* (sum(x .* y) / sxx)          # remove a complex straight line
        return sqrt(sum(abs2, r) / (2 * (blocksize - 2)))
    end
    # a low quantile of block standard deviations underestimates σ slightly; correct for Gaussian
    # noise with ~2·blocksize degrees of freedom per block
    z = sqrt(2) * erfinv(2q - 1)                               # standard normal quantile
    return quantile(σ, q) / (1 + z / sqrt(4 * blocksize))
end

"""
    estimatepoles(spec, [dim]; bandwidth=:auto, ncoef=16, snr=5, maxlinewidth=nothing, exclude=())
        -> Vector{Pole}

Model a frequency-domain spectrum along `dim` (default: the direct dimension) as a sum of Lorentzian
lines, by band-wise harmonic inversion (HSVD), and return one [`Pole`](@ref) per line, sorted by
decreasing chemical shift. The spectrum may be unphased: the phase of each line is part of the
result.

The spectrum is split into overlapping bands of `bandwidth` points. Each band is inverse Fourier
transformed to a short time signal (using only the part covering acquired, not zero-filled,
points) and fitted with `ncoef` damped exponentials. Poles are kept from the central half of each
band only, to avoid edge artefacts; bands wrap around the ends of the spectrum, which is periodic.
Bands with no signal above the threshold are skipped.

For multidimensional data, the trace along `dim` takes the first point of every other time
dimension and sums over the remaining dimensions.

# Keywords
- `bandwidth`: band size in points. `:auto` uses `n ÷ 32` (a power of two between 256 and 4096) for
  an `n`-point spectrum.
- `ncoef`: number of exponentials fitted per band.
- `snr`: minimum signal-to-noise ratio (peak height / noise) for a pole to be kept.
- `maxlinewidth`: largest linewidth (Hz) kept; broader poles model baseline rather than lines.
  Defaults to an eighth of the band width in Hz.
- `exclude`: chemical-shift intervals (ppm) to ignore, for example `[4.6 .. 4.9]` to drop water.

See also [`estimatephase`](@ref), [`autophase`](@ref).
"""
function estimatepoles(spec::NMRData, dim=1; bandwidth=:auto, ncoef::Integer=16,
                       snr::Real=5.0, maxlinewidth=nothing, exclude=())
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    d isa FrequencyDimension ||
        throw(ArgumentError("estimatepoles: dimension $di ($(typeof(d))) is not a " *
                            "FrequencyDimension; Fourier transform it first"))
    isnothing(metadata(d, :mcindex)) &&
        throw(NMRToolsError("estimatepoles: dimension $di has no :mcindex metadata"))

    s = _referencetrace(spec, di)
    n = length(s)
    n >= 16 || throw(ArgumentError("estimatepoles: too few points ($n)"))
    noise = _tracenoise(s)

    δ = collect(data(d))
    Δδ = δ[2] - δ[1]                                  # negative: high ppm first
    swhz = metadata(d, :swhz)
    Δf = isnothing(swhz) ? abs(Δδ) * metadata(d, :bf) * 1e-6 : swhz / n   # Hz per point (bf in Hz)
    td = metadata(d, :td)
    tdfrac = isnothing(td) ? 1.0 : clamp(td / n, 0.0, 1.0)

    bw = bandwidth === :auto ? clamp(nextpow(2, n ÷ 32), 256, 4096) : bandwidth
    nb = min(2 * (bw ÷ 2), 2 * (n ÷ 2))
    ncore = nb ÷ 2
    npad = (nb - ncore) ÷ 2
    fsdec = nb * Δf                                   # sampling rate of the decimated signal (Hz)
    lwmax = isnothing(maxlinewidth) ? fsdec / 8 : Float64(maxlinewidth)
    nuse = clamp(round(Int, nb * tdfrac), min(nb, 4 * ncoef + 2), nb)

    fwd = reverse(s)                                  # increasing frequency
    poles = Pole[]
    for c0 in 1:ncore:n
        core = c0:min(c0 + ncore - 1, n)
        maximum(abs, view(fwd, core)) < snr * noise && continue
        # the spectrum of a DFT is periodic, so bands wrap around at the edges and stay centred
        lo = c0 - npad
        band = mod1.(lo:(lo + nb - 1), n)

        y = FFTW.ifft(FFTW.ifftshift(fwd[band]))[1:nuse]   # band centre → zero frequency
        z = _hankelpoles(y, ncoef)
        c = _amplitudes(y, z)

        for k in eachindex(z)
            ρ = abs(z[k])
            0 < ρ < 1 || continue
            lw = -log(ρ) * fsdec / π
            lw <= lwmax || continue
            pos = lo + nb ÷ 2 + angle(z[k]) * nb / (2π)        # fractional index in `fwd`
            first(core) - 0.5 <= pos < last(core) + 0.5 || continue
            height = abs(c[k]) * (1 - ρ^nuse) / (1 - ρ)
            height >= snr * noise || continue
            shift = δ[1] + (n - pos) * Δδ                      # index n + 1 - pos in `δ`
            any(r -> shift in r, exclude) && continue
            push!(poles, Pole(shift, lw, nb * c[k], height, height / noise))
        end
    end
    return sort!(poles; by=p -> p.shift, rev=true)
end
