"""
    Peak{N}

A peak in an `N`-dimensional frequency-domain spectrum: where it is, how big it is, how
broad it is, and what it is called.

A peak records what was picked, not what a fit produced. Positions are in ppm and
linewidths are full widths at half maximum in Hz, both ordered like the dimensions of the
spectrum the peak came from, so for a 2D spectrum the first coordinate is the direct
dimension (`x`) and the second the indirect dimension (`y`). Linewidths are kept in Hz
because, unlike ppm, they mean the same thing on every nucleus and at every field; the
corresponding relaxation rate is `R2 = π * linewidth`.

Unknown quantities are `NaN`: a peak read from a Sparky list has a position and label but
no intensity or linewidth.

A peak list is simply a `Vector{Peak{N}}`, so `sort`, `filter`, `push!`, `deleteat!` and
indexing all work as usual. Peak lists are also Tables.jl tables, so `DataFrame(pks)`
works directly.

# Examples
```jldoctest
julia> p = Peak(8.21, 120.4; label="G10")      # a 2D peak with only a position and label
Peak{2} G10
  position (ppm): 8.21, 120.4
  linewidth (Hz): NA, NA
  intensity:      NA

julia> ppm(p)
(8.21, 120.4)

julia> label(Peak(p; label="G10N"))            # a copy with a new label
"G10N"
```

```julia
pks = detectpeaks(spec)                   # Vector{Peak{2}}
sort(pks; by=intensity, rev=true)
```

See also [`detectpeaks`](@ref), [`findpeak`](@ref), [`findnearest`](@ref), [`ppm`](@ref),
[`intensity`](@ref), [`linewidth`](@ref), [`loadpeaks`](@ref NMRTools.NMRIO.loadpeaks) and
[`savepeaks`](@ref NMRTools.NMRIO.savepeaks).
"""
struct Peak{N}
    position::NTuple{N,Float64}
    intensity::Float64
    linewidth::NTuple{N,Float64}
    label::String
end

function Peak(position::NTuple{N,Real}; intensity=NaN, linewidth=ntuple(_ -> NaN, N),
              label="") where {N}
    length(linewidth) == N ||
        throw(DimensionMismatch("a $N-dimensional peak needs $N linewidths, got $(length(linewidth))"))
    return Peak{N}(Float64.(position), Float64(intensity), Float64.(Tuple(linewidth)),
                   String(label))
end
Peak(position::Real...; kwargs...) = Peak(position; kwargs...)
Peak(position::AbstractVector{<:Real}; kwargs...) = Peak(Tuple(position); kwargs...)

"""
    Peak(p::Peak; position, intensity, linewidth, label)

A copy of `p` with any of its fields replaced, e.g. `Peak(p; label="G10")`.
"""
function Peak(p::Peak; position=p.position, intensity=p.intensity, linewidth=p.linewidth,
              label=p.label)
    return Peak(Tuple(position); intensity, linewidth, label)
end

Base.ndims(::Peak{N}) where {N} = N
Base.ndims(::Type{<:Peak{N}}) where {N} = N

"""
    ppm(peak)
    ppm(peak, n)

The position of a peak in ppm: a tuple ordered like the spectrum's dimensions, or the
`n`th coordinate on its own.
"""
ppm(p::Peak) = p.position
ppm(p::Peak, n::Integer) = p.position[n]

"""
    intensity(peak)

The height of a peak at its maximum, in the units of the spectrum it was picked from, or
`NaN` if unknown.
"""
intensity(p::Peak) = p.intensity

"""
    linewidth(peak)
    linewidth(peak, n)

The full width at half maximum of a peak in Hz: a tuple ordered like the spectrum's
dimensions, or the `n`th dimension's on its own. `NaN` where unknown. The corresponding
transverse relaxation rate is `π * linewidth(peak, n)`.
"""
linewidth(p::Peak) = p.linewidth
linewidth(p::Peak, n::Integer) = p.linewidth[n]

"""
    label(peak)

The label of a peak, e.g. an assignment such as `"G10"`.
"""
label(p::Peak) = p.label

function Base.show(io::IO, p::Peak)
    print(io, "Peak(", join(_fmt.(ppm(p)), ", "))
    isempty(label(p)) || print(io, "; label=", repr(label(p)))
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", p::Peak{N}) where {N}
    println(io, "Peak{$N}", isempty(label(p)) ? "" : " $(label(p))")
    println(io, "  position (ppm): ", join(_fmt.(ppm(p)), ", "))
    println(io, "  linewidth (Hz): ", join(_fmt.(linewidth(p)), ", "))
    return print(io, "  intensity:      ", _fmt(intensity(p)))
end

function Base.show(io::IO, ::MIME"text/plain", pks::AbstractVector{<:Peak{N}}) where {N}
    print(io, length(pks), "-element Vector{Peak{$N}}")
    isempty(pks) && return
    println(io, ":")
    names = _peakcolumnnames(N)
    rows = [[label(p); _fmt.(collect(ppm(p))); _fmt.(collect(linewidth(p)));
             _fmt(intensity(p))] for p in pks]
    header = ["label"; ["$n (ppm)" for n in names]; ["lw$n (Hz)" for n in names];
              "intensity"]
    widths = [maximum(length, [header[j]; [r[j] for r in rows]]) for j in eachindex(header)]
    println(io, "  ", join(lpad.(header, widths), "  "))
    lines = ["  " * join(lpad.(r, widths), "  ") for r in rows]
    nmax = displaysize(io)[1] - 6
    if length(lines) > nmax > 2
        half = nmax ÷ 2
        lines = [lines[1:half]; "  ⋮"; lines[(end - half + 1):end]]
    end
    return print(io, join(lines, "\n"))
end

_fmt(x::Real) = isnan(x) ? "NA" : string(round(x; sigdigits=6))

# Column names for each dimension, following NMRAnalysis.jl's `x` (direct), `y`, ...
_peakcolumnnames(N) = N ≤ 4 ? ("x", "y", "z", "a")[1:N] : ntuple(n -> "x$n", N)

# ---- Tables.jl ------------------------------------------------------------------------

Tables.istable(::Type{<:AbstractVector{<:Peak}}) = true
Tables.columnaccess(::Type{<:AbstractVector{<:Peak}}) = true

function Tables.columns(pks::AbstractVector{<:Peak{N}}) where {N}
    names = _peakcolumnnames(N)
    cols = Pair{Symbol,Vector}[:label => label.(pks)]
    for n in 1:N
        push!(cols, Symbol(names[n]) => [ppm(p, n) for p in pks])
    end
    for n in 1:N
        push!(cols, Symbol("lw", names[n]) => [linewidth(p, n) for p in pks])
    end
    push!(cols, :intensity => intensity.(pks))
    return (; cols...)
end

# ---- peak detection -------------------------------------------------------------------

"""
    detectpeaks(spec; snr=10, threshold=nothing, sign=:both) -> Vector{Peak{N}}
    detectpeaks(; kwargs...)

Find the peaks in a frequency-domain spectrum of any dimensionality.

The forms without the data argument return a function for use with `|>`.

A peak is a local maximum (or, for negative peaks, minimum) over its immediate neighbours
in every dimension whose height exceeds `threshold`, which by default is `snr` times the
noise level (estimated with [`estimatenoise`](@ref) if it has not been already). Positions
and heights are refined to better than the digital resolution by parabolic interpolation,
and linewidths are estimated from the half-height crossings in each dimension (`NaN` where
an overlapping peak or the edge of the spectrum gets in the way). Heights are measured from
the baseline, taken as the median of the spectrum, so a small offset does not produce
spurious peaks; the reported intensities are heights above that baseline.

`sign` may be `:both`, `:positive` or `:negative`. Peaks are returned strongest first and
are labelled `X1`, `X2`, … to match NMRAnalysis.jl.

Every dimension must be a `FrequencyDimension`. To pick peaks in one plane of a pseudo-2D
experiment, index it first, e.g. `detectpeaks(spec[:, 1])`.

# Examples
```julia
pks = detectpeaks(spec)                     # positive and negative peaks above 10σ
pks = detectpeaks(spec; snr=20, sign=:positive)
pks = detectpeaks(spec; threshold=1e6)      # an absolute threshold
```

See also [`findpeak`](@ref), [`Peak`](@ref).
"""
function detectpeaks(spec::NMRData{T,N}; snr=10, threshold=nothing,
                     sign=:both, snr_threshold=nothing) where {T,N}
    if !isnothing(snr_threshold)
        Base.depwarn("`detectpeaks(spec; snr_threshold)` is deprecated, use `snr` instead.",
                     :detectpeaks)
        snr = snr_threshold
    end
    _checkfrequencydims(spec, "detectpeaks")
    sign in (:both, :positive, :negative) ||
        throw(ArgumentError("sign must be :both, :positive or :negative, got $(repr(sign))"))
    y = _peakdata(spec)
    t = isnothing(threshold) ? snr * _noise(spec) : Float64(threshold)

    candidates = CartesianIndex{N}[]
    sign in (:both, :positive) && append!(candidates, _localmaxima(y, t))
    sign in (:both, :negative) && append!(candidates, _localmaxima(-y, t))

    pks = Peak{N}[_peakat(spec, y, I) for I in candidates]
    sort!(pks; by=p -> abs(intensity(p)), rev=true)
    return Peak{N}[Peak(p; label="X$i") for (i, p) in enumerate(pks)]
end

"""
    findpeak(spec, position; label="") -> Peak
    findpeak(position; label="")

The peak nearest to `position` (a tuple in ppm, ordered like the dimensions of `spec`),
found by climbing from the closest data point to the top of the hill it sits on, or to
the bottom of the valley if the spectrum is negative there. The position, height and
linewidths are refined just as in [`detectpeaks`](@ref).

The forms without the data argument return a function for use with `|>`.

This is what snapping a mouse click onto a peak means.

# Examples
```julia
p = findpeak(spec, (8.2, 120.0); label="G10")
p = findpeak(spec1d, 4.7)
```
"""
function findpeak(spec::NMRData{T,N}, position; label="") where {T,N}
    _checkfrequencydims(spec, "findpeak")
    δ = Tuple(position)
    length(δ) == N ||
        throw(DimensionMismatch("a position in a $N-dimensional spectrum needs $N coordinates"))
    y = _peakdata(spec)
    I = CartesianIndex(ntuple(n -> argmin(abs.(data(dims(spec, n)) .- δ[n])), N))
    s = y[I] < 0 ? -1 : 1
    I = _climb(s .* y, I)
    return Peak(_peakat(spec, y, I); label)
end

"""
    findnearest(peaks, position; scale=1) -> Union{Int,Nothing}

The index of the peak in `peaks` nearest to `position` (a tuple in ppm, or a number for
1D peaks), or `nothing` if the list is empty. This is what selecting a peak with a mouse
click means; the index can be passed straight to `deleteat!` or used to replace the peak.

`scale` weights the dimensions before distances are compared, as a tuple or a single
number, since a ppm of ¹⁵N is not a ppm of ¹H. A typical choice for an ¹H,¹⁵N spectrum is
`scale=(1, 0.2)`; a GUI can pass the inverse of the visible range in each dimension so that
distances are measured as they appear on screen.

# Examples
```julia
i = findnearest(pks, (8.2, 120.0); scale=(1, 0.2))
deleteat!(pks, i)
```
"""
function findnearest(pks::AbstractVector{<:Peak{N}}, position; scale=1) where {N}
    isempty(pks) && return nothing
    δ = Tuple(position)
    w = scale isa Real ? ntuple(_ -> scale, N) : Tuple(scale)
    length(δ) == N == length(w) ||
        throw(DimensionMismatch("$N-dimensional peaks need a position and scale with $N coordinates"))
    distance(p) = sum(abs2, (ppm(p) .- δ) .* w)
    return argmin(distance.(pks))
end

function _checkfrequencydims(spec, fname)
    all(d -> d isa FrequencyDimension, dims(spec)) ||
        throw(ArgumentError("$fname needs a spectrum whose dimensions are all frequency " *
                            "dimensions; for a pseudo-nD experiment, pick peaks in a " *
                            "single plane, e.g. spec[:, 1]"))
    return nothing
end

# Real data measured from the baseline, taken as the median of the spectrum, since most
# points of a spectrum are noise. Spectra commonly sit on a small offset, and measuring
# heights from zero would turn it into a forest of spurious peaks.
function _peakdata(spec)
    y = _realpart.(data(spec))
    return y .- median(y)
end
_realpart(x::Real) = Float64(x)
_realpart(x::Complex) = Float64(real(x))
_realpart(x::Multicomplex) = Float64(realest(x))

function _noise(spec)
    σ = get(metadata(spec), :noise, nothing)
    if isnothing(σ)
        σ = estimatenoise(spec)[:noise]
    end
    return σ
end

# All interior points that exceed `t` and are no lower than any of their neighbours. On a
# plateau, only the first point in storage order is kept.
function _localmaxima(y::AbstractArray{<:Real,N}, t) where {N}
    R = CartesianIndices(y)
    lo, hi = first(R) + oneunit(first(R)), last(R) - oneunit(first(R))
    box = CartesianIndices(ntuple(_ -> -1:1, N))
    centre = CartesianIndex(ntuple(_ -> 0, N))
    found = CartesianIndex{N}[]
    for I in lo:hi
        v = y[I]
        v > t || continue
        ismax = true
        for d in box
            d == centre && continue
            w = y[I + d]
            if w > v || (w == v && (I + d) < I)
                ismax = false
                break
            end
        end
        ismax && push!(found, I)
    end
    return found
end

# Steepest ascent from `I` to a local maximum of `y`.
function _climb(y::AbstractArray{<:Real,N}, I) where {N}
    R = CartesianIndices(y)
    box = CartesianIndices(ntuple(_ -> -1:1, N))
    while true
        best = I
        for d in box
            J = I + d
            J in R && y[J] > y[best] && (best = J)
        end
        best == I && return I
        I = best
    end
end

# A peak at grid point `I`, refined by fitting a parabola through it and its neighbours
# along each dimension, with linewidths from the half-height crossings.
function _peakat(spec::NMRData{T,N}, y, I::CartesianIndex{N}) where {T,N}
    v = y[I]
    s = v < 0 ? -1 : 1
    position = zeros(N)
    ratio = 1.0
    lw = fill(NaN, N)
    for n in 1:N
        ax = dims(spec, n)
        x = data(ax)
        i = I[n]
        e = CartesianIndex(ntuple(k -> k == n ? 1 : 0, N))
        δ = Float64(x[i])
        if 1 < i < length(x)
            a, c = s * y[I - e], s * y[I + e]
            b = s * v
            denom = a - 2b + c
            if denom < 0
                p = clamp(0.5 * (a - c) / denom, -0.5, 0.5)
                δ += p * (p > 0 ? x[i + 1] - x[i] : x[i] - x[i - 1])
                ratio *= (b - 0.25 * (a - c) * p) / b
            end
        end
        position[n] = δ
        lw[n] = _halfwidth(s .* view(y, ntuple(k -> k == n ? Colon() : I[k], N)...), x, i,
                           s * v * ratio) * ax[:bf] * 1e-6
    end
    return Peak(Tuple(position); intensity=v * ratio, linewidth=Tuple(lw))
end

# Full width (ppm) at half of `h` of the peak at index `i` of the 1D trace `y`, from
# linear interpolation of the half-height crossing on each side. A side that meets a
# valley or the edge of the trace before half height gives no answer, and the other side
# is doubled; if neither side gives an answer the width is NaN.
function _halfwidth(y, x, i, h)
    half = h / 2
    crossing(step) = begin
        j = i
        while true
            k = j + step
            checkbounds(Bool, y, k) || return NaN
            y[k] > y[j] && return NaN               # valley before half height
            if y[k] ≤ half
                f = (y[j] - half) / (y[j] - y[k])
                return Float64(x[j] + f * (x[k] - x[j]))
            end
            j = k
        end
    end
    lo, hi = crossing(-1), crossing(1)
    x0 = Float64(x[i])
    isnan(lo) && isnan(hi) && return NaN
    isnan(lo) && return 2abs(hi - x0)
    isnan(hi) && return 2abs(x0 - lo)
    return abs(hi - lo)
end

# ---- registration ---------------------------------------------------------------------

"""
    register!(peaks, spec; maxshift=nothing) -> peaks
    register(peaks, spec; maxshift=nothing) -> Vector{Peak{N}}
    register!(spec; maxshift=nothing)
    register(spec; maxshift=nothing)

Translate a peak list onto a spectrum: find the single offset, applied to every peak, that
best matches the list to the signals observed in `spec`, and move the peaks by it.
`register` returns a moved copy and leaves `peaks` unchanged. The offset found is reported
with `@info`.

The forms without the data argument return a function for use with `|>`.

This corrects for referencing differences, e.g. between a list predicted from BMRB shifts
(see [`loadbmrb`](@ref NMRTools.NMRIO.loadbmrb)) and a spectrum recorded under different
conditions. It does not move peaks individually.

The match is scored as ``\\sum_i \\log(1 + |S(\\delta_i + \\Delta)| / \\sigma)``: the absolute
intensity of the spectrum ``S`` at each shifted peak position, so that negative peaks
count as much as positive ones, in units of the noise ``\\sigma``. The logarithm gives
diminishing returns for very intense signals, so that many peaks landing on signal outweigh
one peak landing on a huge one. To tolerate small disagreements between individual peaks,
each point of the spectrum takes the largest absolute intensity among its immediate
neighbours.

Every whole-point offset within `maxshift` is tried, so the match is global rather than the
nearest local optimum, and the best is then refined to a fraction of a point. `maxshift` is
the largest offset considered in each dimension in ppm, as a tuple or a single number, and
defaults to 10% of the spectral width. Set it to `0` in a dimension to hold that dimension
fixed.

# Examples
```julia
pks = loadbmrb(15000, H1, N15)
register!(pks, spec)                        # shift the BMRB list onto the HSQC
register!(pks, spec; maxshift=(0.3, 2.0))   # search at most ±0.3 ppm ¹H, ±2 ppm ¹⁵N
register!(pks, spec; maxshift=(0.3, 0))     # shift ¹H only
```
"""
function register!(pks::AbstractVector{<:Peak{N}}, spec::NMRData{T,N};
                   maxshift=nothing) where {T,N}
    _checkfrequencydims(spec, "register!")
    isempty(pks) && return pks
    Δ = _registeroffset(pks, spec, maxshift)
    @info "Peak list shifted by $(join(_fmt.(Δ), ", ")) ppm"
    for i in eachindex(pks)
        pks[i] = Peak(pks[i]; position=ppm(pks[i]) .+ Δ)
    end
    return pks
end

register(pks::AbstractVector{<:Peak}, spec::NMRData; kwargs...) =
    register!(copy(pks), spec; kwargs...)

function _registeroffset(pks::AbstractVector{<:Peak{N}}, spec, maxshift) where {N}
    x = [data(dims(spec, n)) for n in 1:N]
    step = [Float64(x[n][2] - x[n][1]) for n in 1:N]
    width = [abs(Float64(last(x[n]) - first(x[n]))) for n in 1:N]
    limit = isnothing(maxshift) ? 0.1 .* width :
            maxshift isa Real ? fill(Float64(maxshift), N) : collect(Float64, maxshift)
    length(limit) == N ||
        throw(DimensionMismatch("maxshift needs $N values for a $N-dimensional spectrum"))
    K = round.(Int, limit ./ abs.(step))

    w = log1p.(_dilate(abs.(_peakdata(spec))) ./ _noise(spec))
    # each peak's nearest grid point
    base = [CartesianIndex(ntuple(n -> round(Int, (ppm(p, n) - x[n][1]) / step[n]) + 1, N))
            for p in pks]
    score(k) = sum(I -> checkbounds(Bool, w, I + k) ? w[I + k] : 0.0, base)

    best, kbest = -Inf, CartesianIndex(ntuple(_ -> 0, N))
    for k in CartesianIndices(ntuple(n -> (-K[n]):K[n], N))
        s = score(k)
        # on a tie, prefer the smaller offset
        if s > best || (s == best && sum(abs2, Tuple(k)) < sum(abs2, Tuple(kbest)))
            best, kbest = s, k
        end
    end

    # refine each dimension to a fraction of a point with a parabola through the neighbours
    return ntuple(N) do n
        e = CartesianIndex(ntuple(m -> m == n ? 1 : 0, N))
        a, b, c = score(kbest - e), best, score(kbest + e)
        denom = a - 2b + c
        f = K[n] > 0 && denom < 0 ? clamp(0.5 * (a - c) / denom, -0.5, 0.5) : 0.0
        return (kbest[n] + f) * step[n] + 0.0   # + 0.0 turns -0.0 into 0.0
    end
end

# The largest value among each point and its immediate neighbours.
function _dilate(a::AbstractArray{<:Real,N}) where {N}
    R = CartesianIndices(a)
    box = CartesianIndices(ntuple(_ -> -1:1, N))
    out = similar(a)
    for I in R
        m = a[I]
        for d in box
            J = I + d
            J in R && (m = max(m, a[J]))
        end
        out[I] = m
    end
    return out
end

# curried forms for piping
detectpeaks(; kwargs...) = spec -> detectpeaks(spec; kwargs...)
findpeak(position; kwargs...) = spec -> findpeak(spec, position; kwargs...)
register!(spec::NMRData; kwargs...) = pks -> register!(pks, spec; kwargs...)
register(spec::NMRData; kwargs...) = pks -> register(pks, spec; kwargs...)
