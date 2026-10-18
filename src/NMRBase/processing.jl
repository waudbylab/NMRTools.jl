# Processing pipeline ##############################################################################
#
# Composable, non-mutating processing functions for time-domain NMRData. Every function accepts an
# NMRData (and optionally a dimension identifier), returns a new NMRData, and provides a curried
# single-argument form for use with `|>` piping. Collection forms map over the values of a
# `Dict{Symbol,<:NMRData}`.
#
# The element type stays `Multicomplex` throughout; only `real`/`imag`/`abs`/`complex` convert.
# Each processing function reads `:mcindex` from the target dimension to determine which imaginary
# unit it should operate on.

# Mapping from dimension index to TimeDimension constructor (mirrors FDIM_CONSTRUCTORS)
const TDIM_CONSTRUCTORS = Dict(1 => T1Dim, 2 => T2Dim, 3 => T3Dim, 4 => T4Dim)

# Internal helpers #################################################################################

"""
    _fourierdims(spec) -> Vector{Int}

Return the array-dimension indices of every dimension that carries an `:mcindex` (i.e. the
`TimeDimension`s / `FrequencyDimension`s that participate in Fourier transformation). Non-Fourier
axes such as `TrelaxDim` or `X2Dim` have `:mcindex == nothing` and are skipped.
"""
function _fourierdims(spec::NMRData)
    return [i for (i, d) in enumerate(dims(spec)) if !isnothing(metadata(d, :mcindex))]
end

"reshape a vector `v` so it broadcasts along array dimension `di` of an `nd`-dimensional array"
_reshape_along(v, di, nd) = reshape(v, ntuple(i -> i == di ? length(v) : 1, nd))

"return a copy of a dimension tuple with the `di`-th element replaced by `newdim`"
function _replacedimtuple(dimstuple, di, newdim)
    # type the vector as the abstract supertype so a FrequencyDimension can replace a TimeDimension
    v = Vector{NMRDimension}([dimstuple...])
    v[di] = newdim
    return Tuple(v)
end

"""
    _rebuilddims(spec, newdata, newdims)

Rebuild `spec` with new data and dimensions, going through the `NMRData(A, dims; …)` constructor so
that DimensionalData's `format` runs on the new dimensions (a plain `rebuild` bypasses formatting and
leaves freshly-built dimensions as unformatted `AutoLookup`s).
"""
function _rebuilddims(spec, newdata, newdims)
    return NMRData(newdata, Tuple(newdims); refdims=spec.refdims, name=spec.name,
                   metadata=deepcopy(metadata(spec)), missingval=spec.missingval)
end

"""
    _setdimmetadata(spec, di, pairs...)

Return `spec` with the metadata of dimension `di` updated by `pairs` (`key => value`). The dimension
is rebuilt with a copy of its metadata, so the input spectrum is left unchanged.
"""
function _setdimmetadata(spec, di, pairs::Pair...)
    d = dims(spec, di)
    md = deepcopy(metadata(d))
    for (k, v) in pairs
        md[k] = v
    end
    newd = rebuild(d; val=rebuild(lookup(d); metadata=md))     # metadata lives on the lookup
    return rebuild(spec; dims=_replacedimtuple(dims(spec), di, newd))
end

"acquisition time (s) after t = 0 for a time dimension, used as the default window `tmax`"
function _acqtime(d)
    sw = metadata(d, :swhz)
    td = metadata(d, :td)
    if !isnothing(sw) && !isnothing(td)
        return (td - _groupdelaypoints(d)) / sw
    elseif !isnothing(sw)
        return length(d) / sw
    else
        return Inf
    end
end

"digital-filter group delay (points) still to be removed from a dimension by `ft`; 0 if none"
_groupdelaypoints(d) = d isa TimeDimension ? Float64(something(metadata(d, :grpdly), 0.0)) : 0.0

"""
    _removegroupdelay!(A, unit, di, g)

Remove a digital-filter group delay of `g` points from the spectrum `A` (immediately after the FFT
along array dimension `di`, before `fftshift`) by multiplying by the linear phase ramp
`exp(2πi·g·k/N)`, with `k` the signed frequency index. This is the Fourier shift theorem: it moves
the time origin forward by `g` points, so build-up points stored before t = 0 end up at negative
time, as in TopSpin and NMRPipe.
"""
function _removegroupdelay!(A::AbstractArray{<:Multicomplex}, unit::Integer, di, g::Real)
    N = size(A, di)
    k = FFTW.fftfreq(N, N)                        # 0, 1, …, -1: signed frequency index
    θ = 2π * g .* k ./ N
    imk = eltype(A)(imN(unit))
    pf = @. cos(θ) + sin(θ) * imk
    A .*= _reshape_along(pf, di, ndims(A))
    return A
end

# fill in a window's tmax when it was left unspecified (Inf), so sine/cosine windows have a length
_settmax(w::NullWindow, ::Real) = w
_settmax(w::UnknownWindow, ::Real) = w
_settmax(w::ExponentialWindow, t::Real) = isfinite(w.tmax) ? w : ExponentialWindow(w.lb, t)
_settmax(w::CosWindow, t::Real) = isfinite(w.tmax) ? w : CosWindow(t)
_settmax(w::Cos2Window, t::Real) = isfinite(w.tmax) ? w : Cos2Window(t)
function _settmax(w::GeneralSineWindow, t::Real)
    return isfinite(w.tmax) ? w : GeneralSineWindow(w.offset, w.endpoint, w.power, t)
end
_settmax(w::WindowFunction, ::Real) = w

"""
    _mcfft!(A, unit, di; inverse=false)

In-place Fourier transform of a multicomplex array `A` along array dimension `di`, treating the
imaginary unit `unit` (the dimension's `:mcindex`) as the complex `i`. Reinterprets the multicomplex
memory as `Complex` (via MulticomplexNumbers' unsafe reinterpretation, exactly as its FFTW extension
does), applies `FFTW.fft!`/`ifft!`, then reinterprets back. `ifft!` includes the `1/N` normalisation
so that `_mcfft!(inverse=true)` inverts `_mcfft!`.
"""
function _mcfft!(A::AbstractArray{<:Multicomplex{T,N,C}}, unit::Integer, di;
                 inverse=false) where {T,N,C}
    w = MulticomplexNumbers.unsafe_ascomplex!(A, unit)
    d = di .+ (N - 1)   # the complex reinterpretation prepends (N-1) collapsed dimensions
    try
        inverse ? FFTW.ifft!(w, d) : FFTW.fft!(w, d)
    finally
        MulticomplexNumbers.unsafe_fromcomplex!(A, unit)
    end
    return A
end

"""
    _buildfreqdim(d::TimeDimension, di, n) -> FrequencyDimension

Construct the `FrequencyDimension` (`F1Dim`…`F4Dim`, chosen by the array position `di`) produced by
Fourier transforming time dimension `d` to `n` points. The chemical-shift axis is computed from
`:bf`, `:swhz` and `:offsetppm`, decreasing (high ppm first) to match the convention used elsewhere
in NMRTools. All time-dimension metadata is carried over, including `:mcindex`.
"""
function _buildfreqdim(d, di, n)
    md = deepcopy(metadata(d).val)
    bf = md[:bf]
    swhz = md[:swhz]
    offsetppm = md[:offsetppm]
    swppm = 1e6 * swhz / bf
    md[:swppm] = swppm
    md[:tdzf] = n
    md[:npoints] = n
    md[:region] = missing
    md[:units] = "ppm"
    md[:grpdly] = 0.0                             # removed by the transform
    # decreasing ppm axis (a range → Regular lookup). `reverse ∘ fftshift` places the carrier (DC)
    # at data index n÷2, so the axis labels that point offsetppm; δ[1] is therefore one point below
    # the naive band edge offsetppm + swppm/2.
    Δ = swppm / n
    δ = range(offsetppm + (n ÷ 2 - 1) * Δ; step=-Δ, length=n)
    return FDIM_CONSTRUCTORS[di](δ; metadata=md)
end

"""
    _buildtimedim(d::FrequencyDimension, di, n) -> TimeDimension

Inverse of [`_buildfreqdim`](@ref): construct the `TimeDimension` (`T1Dim`…`T4Dim`, chosen by the
array position `di`) obtained by inverse Fourier transforming frequency dimension `d` to `n` points.
"""
function _buildtimedim(d, di, n)
    md = deepcopy(metadata(d).val)
    swhz = md[:swhz]
    dt = 1 / swhz
    md[:tdzf] = n
    md[:npoints] = n
    md[:units] = "s"
    md[:region] = missing
    t = range(0.0; step=dt, length=n)
    return TDIM_CONSTRUCTORS[di](t; metadata=md)
end

# apodize ##########################################################################################

"""
    apodize(spec, [dim], [window]; inverse=false)

Multiply the time-domain data by a window function. With no dimension argument every Fourier time
dimension is apodized; with no window, each dimension's `:window` metadata is used (no apodization if
it has none).

The window's length is set from the acquisition time when it is not given, so `Cos2Window()` or
`GeneralSineWindow(0.5, 1, 2)` can be used without a `tmax`. The window starts at t = 0;
digital-filter build-up points before t = 0 are left unwindowed. The window applied is recorded in
the dimension's `:window` metadata.

# Keywords
- `inverse`: divide by the window instead, undoing a previous apodization.

# Examples
```julia
loadfid("expt/1") |> apodize(ExponentialWindow(5.0)) |> zerofill |> ft
loadfid("hsqc/1") |> apodize(T1Dim, ExponentialWindow(5.0)) |> apodize(T2Dim, Cos2Window())
```

See also [`WindowFunction`](@ref), [`zerofill`](@ref).
"""
function apodize(spec::NMRData, dim, window::WindowFunction; inverse=false)
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    d isa TimeDimension ||
        throw(ArgumentError("apodize: dimension $di ($(typeof(d))) is not a TimeDimension"))

    w = _settmax(window, _acqtime(d))
    t = collect(data(d))
    wv = apod(max.(t, 0.0), w)
    wv[t .< 0] .= 1                               # digital-filter build-up before t = 0: unwindowed
    inverse && (wv = 1 ./ wv)

    newdata = parent(spec) .* _reshape_along(wv, di, ndims(spec))
    return _setdimmetadata(rebuild(spec; data=newdata), di, :window => w)
end

# specific dimension, window taken from the dimension's metadata
function apodize(spec::NMRData, dim; inverse=false)
    d = dims(spec, resolvedim(spec, dim))
    w = metadata(d, :window)
    isnothing(w) && (w = NullWindow())
    return apodize(spec, dim, w; inverse=inverse)
end

# all Fourier dimensions, explicit window
function apodize(spec::NMRData, window::WindowFunction; inverse=false)
    out = spec
    for di in _fourierdims(spec)
        out = apodize(out, di, window; inverse=inverse)
    end
    return out
end

# all Fourier dimensions, window from metadata
function apodize(spec::NMRData; inverse=false)
    out = spec
    for di in _fourierdims(spec)
        out = apodize(out, di; inverse=inverse)
    end
    return out
end

# curried forms for piping
apodize(; inverse=false) = spec -> apodize(spec; inverse=inverse)
apodize(window::WindowFunction; inverse=false) = spec -> apodize(spec, window; inverse=inverse)
apodize(dim; inverse=false) = spec -> apodize(spec, dim; inverse=inverse)
function apodize(dim, window::WindowFunction; inverse=false)
    return spec -> apodize(spec, dim, window; inverse=inverse)
end

# collection form
function apodize(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => apodize(v, args...; kw...) for (k, v) in d)
end

# zerofill #########################################################################################

"""
    zerofill(spec, [dim]; factor=2)

Append zeros to the end of the FID, multiplying the number of points by the integer `factor`. With
no dimension argument every Fourier time dimension is zero-filled. The `:tdzf` metadata is updated.

# Examples
```julia
fid |> zerofill |> ft                   # double every Fourier dimension
fid |> zerofill(T1Dim; factor=4) |> ft
```
"""
function zerofill(spec::NMRData, dim; factor=2)
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    d isa TimeDimension ||
        throw(ArgumentError("zerofill: dimension $di ($(typeof(d))) is not a TimeDimension"))
    mcindex = metadata(d, :mcindex)
    isnothing(mcindex) &&
        throw(ArgumentError("zerofill: dimension $di has no :mcindex and cannot be zero-filled"))

    n = size(spec, di)
    newn = factor * n
    nd = ndims(spec)
    M = eltype(parent(spec))
    newsz = ntuple(i -> i == di ? newn : size(spec, i), nd)
    newdata = zeros(M, newsz)
    frontidx = ntuple(i -> i == di ? (1:n) : Colon(), nd)
    newdata[frontidx...] = parent(spec)

    newdim = _extendtimedim(d, di, newn)
    newdims = _replacedimtuple(dims(spec), di, newdim)
    return _rebuilddims(spec, newdata, newdims)
end

function zerofill(spec::NMRData; kw...)
    out = spec
    for di in _fourierdims(spec)
        out = zerofill(out, di; kw...)
    end
    return out
end

zerofill(; kw...) = spec -> zerofill(spec; kw...)
zerofill(dim; kw...) = spec -> zerofill(spec, dim; kw...)
function zerofill(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => zerofill(v, args...; kw...) for (k, v) in d)
end

"extend a time dimension to `n` points, continuing the time axis at the same dwell"
function _extendtimedim(d, di, n; prepend::Integer=0)
    md = deepcopy(metadata(d).val)
    swhz = md[:swhz]
    dt = 1 / swhz
    md[:tdzf] = n
    # keeps negative times of any filter build-up; `prepend` points added before the first
    t = range(first(data(d)) - prepend * dt; step=dt, length=n)
    # points prepended to a FID that starts after t = 0 (negative group delay) move its start back
    prepend > 0 && haskey(md, :grpdly) && (md[:grpdly] = something(md[:grpdly], 0.0) + prepend)
    return TDIM_CONSTRUCTORS[di](t; metadata=md)
end

# Fourier transform ################################################################################

"""
    ft(spec, [dim]; inverse=false)

Fourier transform one dimension, or with no dimension argument every Fourier dimension in order
(direct dimension first). A time dimension becomes the corresponding frequency dimension, or the
reverse with `inverse=true`. The transform uses the dimension's imaginary unit (`:mcindex`), applies
no normalisation, as TopSpin does, and refreshes the noise estimate.

The spectrum is ordered from high to low chemical shift, as for data loaded with
[`loadnmr`](@ref NMRTools.NMRIO.loadnmr), and `ft(ft(spec); inverse=true)` recovers `spec`.

A digital-filter group delay recorded by [`loadfid`](@ref NMRTools.NMRIO.loadfid) (`:grpdly`) is
removed after the transform with a linear phase ramp, which moves the filter's build-up points to
negative time. The frequency dimension has no group delay left, so an inverse transform gives a FID
whose first point is at t = 0.

# Keywords
- `inverse`: inverse transform from the frequency domain to the time domain.

# Examples
```julia
spec = loadfid("expt/1") |> apodize(Cos2Window()) |> zerofill |> ft
spec2 = ft(spec, F1Dim; inverse=true)
```
"""
function ft(spec::NMRData, dim; inverse=false)
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    mcindex = metadata(d, :mcindex)
    isnothing(mcindex) && throw(NMRToolsError("cannot Fourier transform dimension $di " *
                                              "($(typeof(d))): no :mcindex metadata"))

    newdata = copy(parent(spec))
    if inverse
        newdata = reverse(newdata; dims=di)
        newdata = ifftshift(newdata, di)
        _mcfft!(newdata, mcindex, di; inverse=true)
        newdim = _buildtimedim(d, di, size(newdata, di))
    else
        _mcfft!(newdata, mcindex, di; inverse=false)
        g = _groupdelaypoints(d)
        g == 0 || _removegroupdelay!(newdata, mcindex, di, g)
        newdata = fftshift(newdata, di)
        newdata = reverse(newdata; dims=di)
        newdim = _buildfreqdim(d, di, size(newdata, di))
    end

    newdims = _replacedimtuple(dims(spec), di, newdim)
    # refresh the noise estimate: it only becomes physically meaningful in the frequency domain
    return estimatenoise(_rebuilddims(spec, newdata, newdims))
end

function ft(spec::NMRData; kw...)
    out = spec
    for di in _fourierdims(spec)
        out = ft(out, di; kw...)
    end
    return out
end

ft(; kw...) = spec -> ft(spec; kw...)
ft(dim; kw...) = spec -> ft(spec, dim; kw...)
function ft(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => ft(v, args...; kw...) for (k, v) in d)
end

# phase ############################################################################################

"""
    phase(spec, dim, p0, p1=0; pivot=nothing)
    phase(dim, p0, p1=0; pivot=nothing)

Apply a zero-order (`p0`) and first-order (`p1`) phase correction, in degrees, along dimension
`dim`. The phase applied to point `n` of `N` (counting from 0 at the first, highest-ppm point) is
`θ = p0 + p1·n/N`, multiplied in as `exp(-imₖ·θ)` with the dimension's imaginary unit. This is
TopSpin's convention, so `p0` and `p1` are its `PHC0` and `PHC1`. Works in the time or frequency
domain. The second form returns a function for use with `|>`.

For a frequency dimension the total correction applied so far is recorded in the dimension's
`:phase` metadata as `(p0=…, p1=…)`, with `p0` referred to the first point, and the noise estimate
is refreshed.

# Keywords
- `pivot`: position at which only `p0` acts, in ppm for a frequency dimension or seconds for a time
  dimension. Defaults to the first point.

# Examples
```julia
spec |> phase(F1Dim, 128.5) |> real
phase(spec, F1Dim, 40.0, -20.0; pivot=4.7)
```

See also [`autophase`](@ref), [`autophase0`](@ref).
"""
function phase(spec::NMRData, dim, p0::Real, p1::Real=0.0; pivot=nothing)
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    mcindex = metadata(d, :mcindex)
    isnothing(mcindex) &&
        throw(NMRToolsError("cannot phase dimension $di ($(typeof(d))): no :mcindex metadata"))

    M = eltype(parent(spec))
    imk = M(imN(mcindex))
    n = size(spec, di)

    # fractional point position; pivot defaults to the first point (index 0), as in TopSpin/nmrPipe
    ipiv = isnothing(pivot) ? 0.0 : _pivotindex(d, pivot)
    frac = ((0:(n - 1)) .- ipiv) ./ n
    θ = deg2rad(p0) .+ deg2rad(p1) .* frac
    pf = @. cos(θ) - sin(θ) * imk
    newdata = parent(spec) .* _reshape_along(pf, di, ndims(spec))
    newspec = rebuild(spec; data=newdata, metadata=deepcopy(metadata(spec)))
    d isa FrequencyDimension || return newspec

    # accumulate the total correction, with p0 referred to the first point
    prev = something(metadata(d, :phase), (p0=0.0, p1=0.0))
    total = (p0=_wrap180(prev.p0 + p0 - p1 * ipiv / n), p1=prev.p1 + p1)
    # the noise estimate depends on the real part, so refresh it
    return estimatenoise(_setdimmetadata(newspec, di, :phase => total))
end

"0-based fractional point index at which the dimension axis equals `pivot`"
function _pivotindex(d, pivot)
    ax = data(d)
    return (pivot - ax[1]) / (ax[2] - ax[1])
end

phase(dim, p0::Real, p1::Real=0.0; pivot=nothing) = spec -> phase(spec, dim, p0, p1; pivot)
# resolve the ambiguity between `phase(spec, dim, p0)` and the curried `phase(dim, p0, p1)`
phase(spec::NMRData, dim::Real, p0::Real; pivot=nothing) = phase(spec, dim, p0, 0.0; pivot)
function phase(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => phase(v, args...; kw...) for (k, v) in d)
end

# frequencyshift ###################################################################################

"""
    frequencyshift(spec, dim, Δ)
    frequencyshift(dim, Δ)

Shift time-domain data in frequency by `Δ` Hz, by multiplying the FID by `exp(-imₖ·2π·Δ·t)`: a
signal at `+Δ` Hz moves to the carrier. [`solventfilter`](@ref) uses this to suppress signals away
from the carrier; it is also useful for re-referencing in the time domain. The second form returns a
function for use with `|>`.

# Examples
```julia
fid |> frequencyshift(T1Dim, 250.0) |> ft
```
"""
function frequencyshift(spec::NMRData, dim, Δ::Real)
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    mcindex = metadata(d, :mcindex)
    isnothing(mcindex) &&
        throw(NMRToolsError("cannot frequency-shift dimension $di ($(typeof(d))): no :mcindex"))

    M = eltype(parent(spec))
    imk = M(imN(mcindex))
    t = collect(data(d))
    θ = -2π * Δ .* t
    pf = @. cos(θ) + sin(θ) * imk
    newdata = parent(spec) .* _reshape_along(pf, di, ndims(spec))
    return rebuild(spec; data=newdata)
end

frequencyshift(dim, Δ::Real) = spec -> frequencyshift(spec, dim, Δ)
function frequencyshift(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => frequencyshift(v, args...; kw...) for (k, v) in d)
end

# extract ##########################################################################################

"""
    extract(dim, range)
    extract(dim1 => range1, dim2 => range2, ...)

Return a pipeline-compatible closure that selects a region of a spectrum. `dim` is a dimension type
(for example `F1Dim`) and `range` an interval (for example `6.0 .. 10.0`). The noise estimate (`:noise`) of the
selected region is refreshed.

```julia
loadfid("1/ser") |> ft |> extract(F1Dim, 6.0 .. 10.0) |> phase(F1Dim, -50) |> real
```
"""
extract(dim, range) = spec -> _withnoise(spec[dim(range)])
function extract(pairs::Pair...)
    return spec -> _withnoise(spec[(p.first(p.second) for p in pairs)...])
end

"copy of `spec` with its own metadata and a fresh noise estimate"
_withnoise(spec) = estimatenoise(rebuild(spec; metadata=deepcopy(metadata(spec))))

# element-type conversions #########################################################################

"""
    real(spec::NMRData)

Extract the purely real component of a multicomplex spectrum (elementwise `realest`), returning an
`NMRData` with `Float64` element type. This is the final step of a processing pipeline.
"""
Base.real(spec::NMRData{<:Multicomplex}) = rebuild(spec; data=realest.(parent(spec)))

"""
    imag(spec::NMRData)

Return the imaginary part of a multicomplex spectrum: the coefficient of its highest-order imaginary
unit.
"""
Base.imag(spec::NMRData{<:Multicomplex}) = rebuild(spec; data=imag.(parent(spec)))

"""
    abs(spec::NMRData)

Return the magnitude of each point of a multicomplex spectrum, as `Float64` data.
"""
Base.abs(spec::NMRData{<:Multicomplex}) = rebuild(spec; data=abs.(parent(spec)))

"""
    complex(spec::NMRData)

Convert a first-order multicomplex spectrum (`Multicomplex{T,1,2}`) to Julia `Complex{T}`.
"""
function Base.complex(spec::NMRData{<:Multicomplex{T,1,2}}) where {T}
    z = [Complex(realest(m), imag(m)) for m in parent(spec)]
    return rebuild(spec; data=z)
end
