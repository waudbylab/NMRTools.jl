# Solvent suppression #############################################################################
#
# Removes an on- (or off-) resonance signal from the time-domain FID. Two methods:
#   :lowpass  — convolution difference (nmrPipe SOL): subtract a low-pass estimate of the FID. The
#               low-pass is a real-weighted moving average, so this runs natively on multicomplex
#               data (no reinterpretation). Fast; the standard choice for 2D+.
#   :hsvd     — fit the FID with damped exponentials (Hankel-SVD core) and subtract only those whose
#               frequency lies within a band of the target. Frequency-selective, minimal distortion
#               of neighbouring peaks; best for 1D where the cost is affordable.
#
# The target defaults to the carrier. A ppm target is moved on resonance with `frequencyshift`,
# filtered, then shifted back.

"relative weights of a symmetric low-pass window of half-width `fl` (`:boxcar`/`:sine`/`:sine²`)"
function _lowpasswindow(fl::Integer, shape::Symbol)
    m = -fl:fl
    if shape === :boxcar
        return ones(length(m))
    elseif shape === :sine
        return [cos(π * mi / (2 * (fl + 1))) for mi in m]
    elseif shape === :sine²
        return [cos(π * mi / (2 * (fl + 1)))^2 for mi in m]
    else
        throw(ArgumentError("solventfilter: shape must be :boxcar, :sine or :sine², got $shape"))
    end
end

"""
    _lowpassestimate(data, di, w, fl) -> Array

Low-pass estimate of `data` along array dimension `di` by a moving average with (relative) weights
`w` over a ±`fl` window, truncated and renormalised at the edges. Real-weighted, so it works directly
on multicomplex arrays.
"""
function _lowpassestimate(data::AbstractArray, di::Integer, w::AbstractVector, fl::Integer)
    D = ndims(data)
    perm = pushfirst!(deleteat!(collect(1:D), di), di)
    dp = permutedims(data, perm)
    L = size(dp, 1)
    mat = reshape(dp, L, :)
    out = similar(mat)
    for j in axes(mat, 2)
        col = view(mat, :, j)
        for n in 1:L
            lo = max(1, n - fl)
            hi = min(L, n + fl)
            acc = zero(eltype(mat))
            wsum = 0.0
            for m in lo:hi
                wt = w[m - n + fl + 1]
                acc += wt * col[m]
                wsum += wt
            end
            out[n, j] = acc / wsum
        end
    end
    return permutedims(reshape(out, size(dp)...), invperm(perm))
end

"resolve a `shift` spec (`:carrier` or a ppm value) to a frequency offset (Hz) from the carrier"
function _shifthz(shift, d)
    shift === :carrier && return 0.0
    shift isa Real ||
        throw(ArgumentError("solventfilter: shift must be :carrier or a chemical shift in ppm"))
    bf = metadata(d, :bf)
    offsetppm = metadata(d, :offsetppm)
    return (shift - offsetppm) * bf * 1e-6
end

"""
    solventfilter(spec, [dim]; shift=:carrier, method=:lowpass,
                  width=40.0, shape=:sine, band=nothing, ncoef=:auto)

Suppress a solvent (or any narrow on-/off-resonance signal) in the time domain. With no dimension
argument the directly detected dimension is used.

# Keywords
- `shift`: `:carrier` (default) or a chemical shift in ppm. A ppm target is moved on resonance via
  [`frequencyshift`](@ref), filtered, then shifted back.
- `method`: `:lowpass` (convolution difference, native, fast — good for 2D+) or `:hsvd`
  (frequency-selective subspace subtraction — best for 1D).
- `width`: approximate suppression notch width in Hz (larger ⇒ wider notch). Mapped to a low-pass
  filter length via the spectral width.
- `shape`: low-pass window, `:boxcar` | `:sine` | `:sine²` (`:lowpass` only).
- `band`: full frequency band (Hz) around the target to remove (`:hsvd` only; default `width`).
- `ncoef`: number of exponentials to fit (`:hsvd` only); `:auto` ⇒ `clamp(td ÷ 16, 16, 48)`. `:hsvd`
  runs an SVD per data vector, so it is intended for 1D (or a modest number of vectors); prefer
  `:lowpass` for 2D+.

Runs on time-domain data before `apodize`. Element type stays `Multicomplex`.
"""
function solventfilter(spec::NMRData, dim; shift=:carrier, method::Symbol=:lowpass,
                       width::Real=40.0, shape::Symbol=:sine, band=nothing, ncoef=:auto)
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    d isa TimeDimension ||
        throw(ArgumentError("solventfilter: dimension $di ($(typeof(d))) is not a TimeDimension"))
    mcindex = metadata(d, :mcindex)
    isnothing(mcindex) &&
        throw(ArgumentError("solventfilter: dimension $di has no :mcindex"))

    swhz = metadata(d, :swhz)
    Δ = _shifthz(shift, d)                       # Hz offset of target from carrier
    sp = Δ == 0 ? spec : frequencyshift(spec, di, Δ)   # move target on resonance

    if method === :lowpass
        fl = clamp(round(Int, swhz / width), 1, size(sp, di) ÷ 2)
        w = _lowpasswindow(fl, shape)
        lowfreq = _lowpassestimate(parent(sp), di, w, fl)
        out = rebuild(sp; data=parent(sp) .- lowfreq)
    elseif method === :hsvd
        dt = 1 / swhz
        bandhz = isnothing(band) ? Float64(width) : Float64(band)
        n = size(sp, di)
        nc = ncoef === :auto ? clamp(n ÷ 16, 16, 48) : ncoef
        kernel = _percolumn(x -> _hsvdsubtract(x; dt=dt, targethz=0.0, bandhz=bandhz, ncoef=nc))
        out = rebuild(sp; data=_maptransform(parent(sp), di, mcindex, kernel))
    else
        throw(ArgumentError("solventfilter: method must be :lowpass or :hsvd, got $method"))
    end

    return Δ == 0 ? out : frequencyshift(out, di, -Δ)
end

# default dimension = directly detected dimension (array dimension 1)
solventfilter(spec::NMRData; kw...) = solventfilter(spec, 1; kw...)

solventfilter(dim; kw...) = spec -> solventfilter(spec, dim; kw...)
solventfilter(; kw...) = spec -> solventfilter(spec; kw...)
function solventfilter(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => solventfilter(v, args...; kw...) for (k, v) in d)
end
