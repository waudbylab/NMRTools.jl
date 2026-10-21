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

"move the time origin of `A` along `di` forward by `g` points (a circular shift, by Fourier transform)"
function _shiftorigin!(A, mcindex, di, g)
    _mcfft!(A, mcindex, di)
    _removegroupdelay!(A, mcindex, di, g)
    return _mcfft!(A, mcindex, di; inverse=true)
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
    solventfilter(spec, [dim]; shift=:carrier, method=:lowpass, width=40.0, shape=:sine,
                  band=nothing, ncoef=:auto)

Suppress a solvent signal, or any narrow signal, in the time domain. With no dimension argument the
direct dimension is filtered.

# Keywords
- `shift`: `:carrier`, or the chemical shift (ppm) of the signal to remove. A signal away from the
  carrier is moved to it with [`frequencyshift`](@ref), filtered, and moved back.
- `method`: `:lowpass` subtracts a low-pass filtered copy of the FID (the convolution difference
  method of nmrPipe's `SOL`). It is fast and suits multidimensional data. `:hsvd` fits the FID with
  decaying oscillations and subtracts those within `band` of the target, which disturbs nearby
  peaks less but fits each FID separately, so it suits 1D data.
- `width`: approximate width (Hz) of the region suppressed by `:lowpass`. A larger value
  suppresses a wider region.
- `shape`: low-pass window for `:lowpass`: `:boxcar`, `:sine` or `:sine²`.
- `band`: width (Hz) of the region removed by `:hsvd`. Defaults to `width`.
- `ncoef`: number of oscillations fitted by `:hsvd`. `:auto` uses `clamp(td ÷ 16, 16, 48)`.

Use on time-domain data, before [`apodize`](@ref). The digital-filter group delay recorded by
[`loadfid`](@ref NMRTools.NMRIO.loadfid) is taken into account: the FID is filtered from t = 0,
with the build-up points before t = 0 treated as negative time, as in `ft`.

# Examples
```julia
fid |> solventfilter |> apodize |> zerofill |> ft
solventfilter(fid; shift=4.7, method=:hsvd, band=60)
```
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

    # Filter a FID whose time origin is its first point. The digital filter's build-up points
    # before t = 0 are not part of the decay: left in place, they bias the low-pass average and
    # the HSVD fit at the start of the array, and the error appears after the group-delay
    # correction in `ft` as a baseline ripple with period 1/(group delay). The delay is removed
    # here exactly as `ft` removes it, which moves the build-up points to the end of the array
    # (negative time), and restored after filtering.
    g = _groupdelaypoints(d)
    A = copy(parent(sp))
    g == 0 || _shiftorigin!(A, mcindex, di, g)
    if method === :lowpass
        fl = clamp(round(Int, swhz / width), 1, size(A, di) ÷ 2)
        w = _lowpasswindow(fl, shape)
        A .-= _lowpassestimate(A, di, w, fl)
    elseif method === :hsvd
        dt = 1 / swhz
        bandhz = isnothing(band) ? Float64(width) : Float64(band)
        nc = ncoef === :auto ? clamp(size(A, di) ÷ 16, 16, 48) : ncoef
        # fit the decay from t = 0 only: the build-up points now at the end of the array are not
        # part of it, and are left unchanged
        L = size(A, di) - ceil(Int, max(g, 0))
        hsvd(x) = _hsvdsubtract(x; dt=dt, targethz=0.0, bandhz=bandhz, ncoef=nc)
        kernel = _percolumn(x -> vcat(hsvd(x[1:L]), x[(L + 1):end]))
        A = _maptransform(A, di, mcindex, kernel)
    else
        throw(ArgumentError("solventfilter: method must be :lowpass or :hsvd, got $method"))
    end
    g == 0 || _shiftorigin!(A, mcindex, di, -g)
    out = rebuild(sp; data=A)

    return Δ == 0 ? out : frequencyshift(out, di, -Δ)
end

# default dimension = directly detected dimension (array dimension 1)
solventfilter(spec::NMRData; kw...) = solventfilter(spec, 1; kw...)

solventfilter(dim; kw...) = spec -> solventfilter(spec, dim; kw...)
solventfilter(; kw...) = spec -> solventfilter(spec; kw...)
function solventfilter(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => solventfilter(v, args...; kw...) for (k, v) in d)
end
