# Automatic phase correction #######################################################################
#
# Lines are found with a pole model of the spectrum (`estimatepoles`), which gives their positions,
# linewidths and signal-to-noise ratios. Poles are filtered and grouped first:
#
#   - Multiplets. J-evolution during echoes (for example excitation sculpting) gives the lines of a
#     multiplet antiphase dispersive components, which cancel only in the sum over the multiplet.
#     Poles closer than a coupling constant (plus their half-linewidths) are grouped.
#   - Broad poles that model baseline distortion or the tails of non-Lorentzian lines such as a
#     radiation-damped water signal. Poles much broader than the typical linewidth are dropped when
#     a narrower line lies within their linewidth; isolated broad poles are kept as genuine lines.
#
# Each group's phase is the argument of the spectrum's integral over a window around the group, less
# a linear baseline. The phases of the poles themselves are not used: a weak line's pole can absorb
# part of the dispersive tail of a strong neighbour, which shifts its phase by tens of degrees or more
# while barely changing its position.
#
# The phase correction θ(f) = p0 + p1·f (f the fractional position from the first point) makes
# every group absorptive: φᵣ − θ(fᵣ) ≡ 0 (mod π). The fit is a weighted circular regression on the
# doubled angle, so negative lines (φ = π) count as correctly phased rather than as phase errors:
#
#     maximise  |Σᵣ wᵣ exp(2i(φᵣ − p1·fᵣ))|   over p1,   then   p0 = ½ arg Σᵣ wᵣ exp(2i(φᵣ − p1·fᵣ))
#
# The two edges of the spectrum add two more observations. Far from the lines, a correctly phased
# spectrum has a flat real part, while the imaginary part carries the dispersion tails of the strong
# lines (largest for spectra with strong or broad signals, such as protein spectra with residual
# water). A small phase error leaks these tails into the real part as a tilt (zero order) or curl
# (first order) of the baseline at the edges. The phase that makes the real part vanish over each
# edge region is therefore a phase observation, weighted by its signal-to-noise ratio like the lines;
# the left and right edges constrain p0 and p0 + p1.
#
# Observations are weighted by inverse variance: a statistical part, 1/SNR², plus a systematic part.
# For line groups the systematic part is the scatter of their residuals about the fit, which
# reflects the accuracy of the pole model; for the edges it is a small fixed floor. This lets a few
# precise edge observations outweigh many individually biased line groups.
#
# p0 has a closed form for each p1, so p1 is found by a bounded one-dimensional scan with a penalty
# towards its expected value. The penalty is scaled by the observed phase scatter and divided by the
# effective number of observations, so it acts as a prior that matters only when p1 is poorly
# determined. The remaining 180° ambiguity is resolved by a majority vote of the line groups, each
# counting once, so that one dominant line such as residual water cannot decide it.

"""
    _phasegroups(poles, bf, coupling, maxwidthratio) -> Vector{Vector{Pole}}

Group poles (sorted by decreasing shift) for phasing: drop poles broader than `maxwidthratio` times
the SNR-weighted median linewidth that have a narrower pole within their linewidth, then merge
neighbours closer than `coupling` Hz plus their mean linewidth. `bf` is the spectrometer frequency
(Hz).
"""
function _phasegroups(poles::AbstractVector{Pole}, bf::Real, coupling::Real, maxwidthratio::Real)
    lw = [p.linewidth for p in poles]
    o = sortperm(lw)
    cw = cumsum([poles[i].snr for i in o])
    lwmedian = lw[o[findfirst(>=(cw[end] / 2), cw)]]

    # a broad pole is dropped only when a narrower pole lies within its linewidth: it then models the
    # tails of that line or the baseline under it. Isolated broad poles are genuine broad lines.
    hz(p, q) = abs(p.shift - q.shift) * bf * 1e-6
    underlying(p) = any(q -> q.linewidth <= maxwidthratio * lwmedian && hz(p, q) < p.linewidth,
                        poles)
    groups = Vector{Vector{Pole}}()
    for p in poles
        p.linewidth > maxwidthratio * lwmedian && underlying(p) && continue
        if !isempty(groups)
            q = last(groups[end])
            if (q.shift - p.shift) * bf * 1e-6 < coupling + 0.5 * (q.linewidth + p.linewidth)
                push!(groups[end], p)
                continue
            end
        end
        push!(groups, [p])
    end
    return groups
end

"""
    _edgeobservations(s, poles, δ, bf, edges, noise) -> Vector{NamedTuple}

Phase observations from the two edges of the complex trace `s`: for each edge, the region of
`edges` × the spectrum (at least 16 points) at that end, excluding points within four linewidths of a
pole. The observation is the phase at which the real part of the region's sum vanishes, with the
region's centre as a fraction from the first point and the signal-to-noise ratio of the sum.
"""
function _edgeobservations(s, poles, δ, bf, edges, noise)
    n = length(s)
    width = max(16, round(Int, edges * n))
    Δ = δ[2] - δ[1]
    nearline = falses(n)
    for p in poles
        c = 1 + (p.shift - δ[1]) / Δ
        h = 4 * p.linewidth / (abs(Δ) * bf * 1e-6)
        nearline[max(1, floor(Int, c - h)):min(n, ceil(Int, c + h))] .= true
    end
    obs = NamedTuple{(:edge, :phase, :position, :snr),Tuple{Symbol,Float64,Float64,Float64}}[]
    for (edge, r) in ((:left, 1:width), (:right, (n - width + 1):n))
        r = filter(i -> !nearline[i], collect(r))
        length(r) >= 8 || continue
        total = sum(view(s, r))
        push!(obs, (edge=edge, phase=angle(total) - π / 2, position=(mean(r) - 1) / n,
                    snr=abs(total) / (noise * sqrt(length(r)))))
    end
    return obs
end

"systematic phase uncertainty assumed for an edge observation (degrees)"
const _EDGEFLOOR = 2.0

"""
    _groupintegralphases(s, groups, δ, bf) -> Vector{Float64}

Phase of each line group measured from the spectrum: the integral of the complex trace `s` over a
window centred on the group, extending four times its largest linewidth beyond its outer lines but
no more than halfway to the neighbouring groups, less a straight complex baseline through the means
of the window's outer points (which removes the tails of distant lines).
"""
function _groupintegralphases(s, groups, δ, bf)
    n = length(s)
    Δ = δ[2] - δ[1]
    hzpt = abs(Δ) * bf * 1e-6
    idx(x) = 1 + (x - δ[1]) / Δ
    firsti = [idx(maximum(p.shift for p in g)) for g in groups]
    lasti = [idx(minimum(p.shift for p in g)) for g in groups]
    centre = [idx(_groupcentre(g)) for g in groups]
    reach = [4 * maximum(p.linewidth for p in g) / hzpt for g in groups]
    lo, hi = firsti .- reach, lasti .+ reach
    for i in 1:(length(groups) - 1)
        mid = (lasti[i] + firsti[i + 1]) / 2
        hi[i] = min(hi[i], mid)
        lo[i + 1] = max(lo[i + 1], mid)
    end
    return map(eachindex(groups)) do i
        h = max(min(centre[i] - lo[i], hi[i] - centre[i]), 4.0)
        r = max(1, round(Int, centre[i] - h)):min(n, round(Int, centre[i] + h))
        m = max(2, length(r) ÷ 20)
        a0, a1 = mean(view(s, r[1:m])), mean(view(s, r[(end - m + 1):end]))
        x = (0:(length(r) - 1)) .- (m - 1) / 2
        L = length(r) - m
        base = a0 .+ (a1 - a0) .* x ./ L
        return angle(sum(view(s, r) .- base))
    end
end

"amplitude-weighted chemical shift of a group of poles"
_groupcentre(g) = sum(abs(p.amplitude) * p.shift for p in g) / sum(abs(p.amplitude) for p in g)

"wrap an angle to (-π/2, π/2]: the residual of a phase taken modulo π"
_halfwrap(x) = angle(cis(2x)) / 2

"wrap an angle in degrees to (-180, 180]"
_wrap180(x) = rad2deg(angle(cis(deg2rad(x))))

"weighted doubled-angle resultant for first-order phase `p1r` (radians)"
_resultant(φ, f, w, p1r) = sum(w .* cis.(2 .* (φ .- p1r .* f)))

"optimal zero-order phase (radians, modulo π) for first-order phase `p1r`"
_bestp0(φ, f, w, p1r) = angle(_resultant(φ, f, w, p1r)) / 2

"""
    _scanp1(φ, f, w, p1range, p1expected, λ) -> Float64

First-order phase (radians) maximising the doubled-angle resultant within `p1range` (degrees), less a
penalty `λ/nₑ·((p1 − p1expected)/180)²`, where `nₑ = (Σw)²/Σw²` is the effective number of groups.
The caller scales `λ` by the observed phase scatter, so the penalty acts as a prior whose influence
falls as the fit becomes better determined.
"""
function _scanp1(φ, f, w, p1range, p1expected, λ)
    sw = sum(w)
    neff = sw^2 / sum(abs2, w)
    J(p1deg) = abs(_resultant(φ, f, w, deg2rad(p1deg))) / sw -
               λ / neff * ((p1deg - p1expected) / 180)^2
    lo, hi = p1range
    grid = range(lo, hi; length=max(2, ceil(Int, hi - lo) + 1))     # 1° steps
    best = grid[argmax(J.(grid))]
    res = optimize(x -> -J(x), max(lo, best - step(grid)), min(hi, best + step(grid)))
    return deg2rad(Optim.minimizer(res))
end

"""
    estimatephase(spec, [dim]; order=1, p1=0, maxdelay=20e-6, pivot=nothing, snr=10,
                  coupling=20, maxwidthratio=3, maxresidual=30, λ=20, edges=0.01, kw...)
        -> NamedTuple

Estimate the zero- and first-order phase correction (degrees) for dimension `dim` (default: the
direct dimension) of a frequency-domain spectrum. Returns
`(p0, p1, pivot, residual, nlines, rejected, edges)`, where `phase(spec, dim, p0, p1; pivot)` applies
the correction, `residual` is the weighted rms phase error (degrees) of the observations used,
`nlines` the number of line groups used, `rejected` lists the line groups left out of the fit as
`(shift, residual, relsnr)`, strongest first, where `relsnr` is the signal-to-noise ratio relative to
the strongest group used, and `edges` lists the edge observations as
`(edge, snr, residual, used)`. A large residual, a small number of groups, or a strong rejected group
means the estimate should not be trusted. [`autophase`](@ref) warns when strong groups are rejected.

Lines are found with [`estimatepoles`](@ref). Poles within a coupling constant of each other are
grouped, because the lines of a multiplet can carry antiphase dispersive components that cancel only
in the sum. Each group's phase is the argument of the integral of the spectrum over a window around
the group (up to four linewidths beyond its outer lines, and no more than halfway to the next group),
less a linear baseline through the ends of the window. Poles much broader than the typical linewidth are ignored when they overlap
a narrower line, since they then model its tails or the baseline. The correction is
the weighted circular regression of group phase against position, taken modulo 180° so that negative
peaks are not treated as phase errors. Groups with a residual above `maxresidual` degrees are
discarded and the fit repeated. The 180° ambiguity is resolved by a majority vote of the groups, so
that most lines are positive.

The edges of the spectrum add two observations. With the correct phase the real part is flat far
from the lines, while the imaginary part carries the dispersion tails of the strong lines; the phase
at which the real part of each edge region vanishes is used as a further observation, weighted by its
signal-to-noise ratio. These are most informative for spectra with strong or broad signals, such as
protein spectra with residual water. They are skipped when their signal-to-noise ratio is below
`snr`, and used only if they agree with the fit to the line groups within three times the combined
uncertainty of the groups' scatter and the edge's own 2° (and `maxresidual`), so that they refine
the line-based estimate rather than override it. Observations are weighted by inverse variance,
combining their signal-to-noise ratio with a systematic uncertainty: the scatter of the line groups
about the fit (a scaled median absolute residual, so that a few noise groups do not inflate it), or
2° for an edge.

# Keywords
- `order`: `1` fits `p0` and `p1`; `0` fits `p0` only, with the first-order phase fixed at `p1`.
- `p1`: expected first-order phase (degrees). Close to zero once the digital-filter group delay has
  been removed.
- `maxdelay`: largest timing error (seconds) to allow for. A first-order phase is a shift of the time
  origin: an error of τ gives a first-order phase of 360°·SW·τ across the spectral width SW, so the
  search covers `p1 ± 360·SW·maxdelay` degrees. The default of 20 µs covers the residual dead-time
  and receiver delays of most experiments; it gives ±100° for a 14 kHz ¹H spectrum and ±950° for a
  130 kHz ¹⁹F spectrum.
- `pivot`: chemical shift (ppm) at which `p0` is reported. Defaults to the first point, as for
  TopSpin's `PHC0`/`PHC1`.
- `snr`: minimum signal-to-noise ratio of lines used.
- `coupling`: poles closer than this (Hz, plus their mean linewidth) are grouped as one multiplet.
- `maxwidthratio`: poles broader than this multiple of the median linewidth are ignored when a
  narrower line lies within their linewidth.
- `maxresidual`: groups further than this (degrees) from the fitted phase are discarded.
- `λ`: strength of the prior pulling the first-order phase towards `p1`. It is scaled by the
  variance of the phase residuals and divided by the effective number of observations, so it matters
  when few observations determine the fit or their phases scatter, and not when many agree closely.
  Increase it when the lines cover only a small part of the spectrum, so `p1` is poorly determined.
- `edges`: fraction of the spectrum at each end used as an edge observation; `0` disables them.
  Points within four linewidths of a line are excluded.
- Other keywords (`exclude`, `bandwidth`, …) are passed to [`estimatepoles`](@ref). Use `exclude` to
  drop solvent or buffer lines.

See also [`autophase`](@ref), [`phase`](@ref).
"""
function estimatephase(spec::NMRData, dim=1; order::Integer=1, p1::Real=0.0,
                       maxdelay::Real=20e-6, pivot=nothing, snr::Real=10.0,
                       coupling::Real=20.0, maxwidthratio::Real=3.0, maxresidual::Real=30.0,
                       λ::Real=20.0, edges::Real=0.01, kw...)
    order in (0, 1) || throw(ArgumentError("estimatephase: order must be 0 or 1, got $order"))
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    poles = estimatepoles(spec, di; snr=snr, kw...)
    isempty(poles) &&
        throw(NMRToolsError("estimatephase: no lines found with signal-to-noise above $snr"))
    groups = _phasegroups(poles, metadata(d, :bf), coupling, maxwidthratio)
    # first-order search range from the timing uncertainty
    swhz = something(metadata(d, :swhz), abs(data(d)[2] - data(d)[1]) * length(d) *
                                         metadata(d, :bf) * 1e-6)
    p1max = 360 * swhz * maxdelay
    p1range = (p1 - p1max, p1 + p1max)

    δ = data(d)
    n = length(δ)
    a = [sum(p.amplitude for p in g) for g in groups]
    # group phases from integrals of the spectrum: pole phases are less reliable, since a weak
    # line's fitted pole can absorb the dispersive tail of a strong neighbour
    φ = _groupintegralphases(_referencetrace(spec, di), groups, δ, metadata(d, :bf))
    # position: amplitude-weighted centre of the group, as a fraction from the first point
    f = [(_groupcentre(g) - δ[1]) / (δ[2] - δ[1]) / n for g in groups]
    # phase variance ∝ 1/SNR²; capped so that a few strong lines cannot dominate
    groupsnr = [abs(a[i]) / sum(abs(p.amplitude) for p in g) * sum(p.snr for p in g)
                for (i, g) in enumerate(groups)]

    # edge observations (flat baseline), appended after the line groups
    nlinegroups = length(groups)
    edgeobs = []
    if edges > 0
        s = _referencetrace(spec, di)
        edgeobs = filter(e -> e.snr >= snr,
                         _edgeobservations(s, poles, δ, metadata(d, :bf), edges, _tracenoise(s)))
        append!(φ, [e.phase for e in edgeobs])
        append!(f, [e.position for e in edgeobs])
        append!(groupsnr, [e.snr for e in edgeobs])
    end
    isline = (1:length(φ)) .<= nlinegroups
    # Inverse-variance weights: statistical variance 1/SNR² plus a systematic term, the observed
    # scatter of the line groups (model error of the pole fit) or a fixed floor for the edges.
    weights(σlines) = [1 / (1 / groupsnr[i]^2 + (isline[i] ? σlines : deg2rad(_EDGEFLOOR))^2)
                       for i in eachindex(groupsnr)]
    # The prior on p1 is scaled by the variance of the phase residuals: weak when the lines agree
    # closely (well-determined fit), stronger when model errors make the phases scatter. Start from
    # a scatter `σ0` and update it from each fit. `candidates` are the observations allowed in.
    function fit(candidates, σ0)
        inl = copy(candidates)
        p0r, p1r = 0.0, deg2rad(p1)
        σ = σ0                           # scatter of the line groups about the fit
        w = weights(σ)
        for iter in 1:8
            p1r = order == 0 ? deg2rad(p1) :
                  _scanp1(φ[inl], f[inl], w[inl], p1range, p1, λ * σ^2)
            p0r = _bestp0(φ[inl], f[inl], w[inl], p1r)
            r = _halfwrap.(φ .- p0r .- p1r .* f)
            used = inl .& isline
            # robust scatter (scaled median absolute residual): a few noise groups with random
            # phases inside `maxresidual` would otherwise hold it high, keeping their weight up
            σnew = any(used) ? max(1.4826 * median(abs.(r[used])), deg2rad(1)) : σ
            newinl = candidates .& (abs.(r) .<= deg2rad(maxresidual))
            converged = (newinl == inl || !any(newinl)) && abs(σnew - σ) < deg2rad(0.5)
            σ = σnew
            w = weights(σ)
            any(newinl) && (inl = newinl)
            converged && iter > 1 && break
        end
        return p0r, p1r, inl, w, σ
    end

    # Fit the line groups first, from a scatter of 15°. An edge joins only if it agrees with that fit
    # within three times the combined uncertainty of the line groups and the edge floor (and
    # `maxresidual`), so that the edges refine the line-based estimate but cannot override it. The
    # refit starts from the scatter of the line groups: starting it from 15° again would give one
    # edge (2° floor) the weight of dozens of lines. Edges are wrong when the baseline has an offset (for example from a
    # wrong first-point factor) or when a strong line has a genuinely different phase, since its
    # dispersion tail reaches the edges too.
    p0r, p1r, inl, w, σlines = fit(isline, deg2rad(15))
    redge = _halfwrap.(φ .- p0r .- p1r .* f)
    gate = min(deg2rad(maxresidual), 3 * hypot(σlines, deg2rad(_EDGEFLOOR)))
    goodedge = .!isline .& (abs.(redge) .<= gate)
    any(goodedge) && ((p0r, p1r, inl, w, _) = fit(isline .| goodedge, σlines))

    # resolve the 180° ambiguity: most line groups positive (edges carry no sign information)
    lines = inl .& isline
    any(lines) || (lines = isline)
    sum(sign.(cos.(φ[lines] .- p0r .- p1r .* f[lines]))) < 0 && (p0r += π)

    rall = rad2deg.(_halfwrap.(φ .- p0r .- p1r .* f))
    residual = sqrt(sum(w[inl] .* rall[inl] .^ 2) / sum(w[inl]))
    # line groups left out of the fit, strongest first, with their phase disagreement (degrees)
    out = sort(findall(.!inl .& isline); by=i -> -groupsnr[i])
    maxused = maximum(groupsnr[inl .& isline]; init=maximum(groupsnr[inl]))
    rejected = [(shift=_groupcentre(groups[i]), residual=rall[i], relsnr=groupsnr[i] / maxused)
                for i in out]
    edgeinfo = [(edge=e.edge, snr=e.snr, residual=rall[nlinegroups + j],
                 used=inl[nlinegroups + j]) for (j, e) in enumerate(edgeobs)]

    p0deg, p1deg = rad2deg(p0r), rad2deg(p1r)
    isnothing(pivot) || (p0deg += p1deg * _pivotindex(d, pivot) / n)
    return (p0=_wrap180(p0deg), p1=p1deg, pivot=pivot, residual=residual,
            nlines=count(inl .& isline), rejected=rejected, edges=edgeinfo)
end

"""
    autophase(spec, [dim]; kw...)

Estimate the phase correction of dimension `dim` with [`estimatephase`](@ref) and apply it with
[`phase`](@ref). Without `dim`, every frequency dimension is phased in turn, direct dimension first;
a dimension in which no lines are found is left unphased, with a warning.
The total correction applied is recorded in each dimension's `:phase` metadata, and the estimate is
logged at debug level. Keywords are passed to `estimatephase`, for example `order=0` for zero-order
only, or `exclude=[4.6 .. 4.9]` to ignore water; without `dim` they apply to every dimension.

For multidimensional data the correction is estimated from a single trace along `dim` (see
[`estimatepoles`](@ref)) and applied to the whole dataset.

```julia
spec = loadfid("expt/1") |> apodize(ExponentialWindow(1.0)) |> zerofill |> ft |> autophase |> real
```
"""
function autophase(spec::NMRData, dim; kw...)
    est = estimatephase(spec, dim; kw...)
    ntotal = est.nlines + length(est.rejected)
    @debug "autophase: p0 = $(round(est.p0; digits=1))°, p1 = $(round(est.p1; digits=1))° " *
          "(residual $(round(est.residual; digits=1))° over $(est.nlines) of $ntotal line groups)"
    # warn about strong lines that could not be phased together with the rest
    strong = filter(r -> r.relsnr >= 0.02, est.rejected)
    # Only when the lines left out are a substantial part of the signal: in crowded spectra many
    # small overlapping groups scatter without affecting the result.
    if !isempty(strong) && length(strong) <= max(2, est.nlines ÷ 4)
        where = join(["$(round(r.shift; digits=2)) ppm ($(round(Int, r.residual))°)"
                      for r in first(strong, 3)], ", ")
        more = length(strong) > 3 ? " and $(length(strong) - 3) more" : ""
        @warn "autophase: strong lines left out of the fit, with their remaining phase error: " *
              "$where$more. They may have a genuinely different phase, or the first-order phase " *
              "may be far from the expected value; pass `p1` (and `maxdelay`) to set it."
    end
    return phase(spec, dim, est.p0, est.p1; pivot=est.pivot)
end

function autophase(spec::NMRData; kw...)
    fdims = _frequencydims(spec)
    isempty(fdims) &&
        throw(ArgumentError("autophase: no frequency dimensions to phase; Fourier transform first"))
    length(fdims) == 1 && return autophase(spec, only(fdims); kw...)
    # several dimensions: one without enough signal to phase (common for indirect dimensions) is
    # left as it is, with a warning, as long as at least one dimension can be phased
    out = spec
    failures = NMRToolsError[]
    for di in fdims
        try
            out = autophase(out, di; kw...)
        catch err
            err isa NMRToolsError || rethrow()
            push!(failures, err)
            @warn "autophase: dimension $di left unphased: $(err.message)"
        end
    end
    length(failures) == length(fdims) && throw(first(failures))
    return out
end

autophase(; kw...) = spec -> autophase(spec; kw...)
autophase(dim; kw...) = spec -> autophase(spec, dim; kw...)
function autophase(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => autophase(v, args...; kw...) for (k, v) in d)
end

"""
    autophase0(spec, [dim]; exclude=())

Correct the zero-order phase of dimension `dim` of a frequency-domain spectrum, or without `dim`
of every frequency dimension in turn. `p0` is chosen so that the integral of the dispersion mode over the
spectrum is zero and the integral of the absorption mode is positive: it is the phase of the integral
of the complex spectrum, which equals the phase of the first point of the FID. Points in the
chemical-shift intervals `exclude`, for example `[4.6 .. 4.9]`, are left out of the integral. No
first-order correction is applied.

The integral includes everything in the spectrum, so baseline offsets, broad background signals and
corrupted first FID points all bias the result. For multidimensional data the correction is
estimated from a single trace along `dim` (see [`estimatepoles`](@ref)) and applied to the whole
dataset. The correction is recorded in the dimension's `:phase` metadata.

See also [`autophase`](@ref), [`phase`](@ref).
"""
function autophase0(spec::NMRData, dim; exclude=())
    di = resolvedim(spec, dim)
    d = dims(spec, di)
    d isa FrequencyDimension ||
        throw(ArgumentError("autophase0: dimension $di ($(typeof(d))) is not a " *
                            "FrequencyDimension; Fourier transform it first"))
    isnothing(metadata(d, :mcindex)) &&
        throw(NMRToolsError("autophase0: dimension $di has no :mcindex metadata"))
    s = _referencetrace(spec, di)
    keep = [!any(r -> x in r, exclude) for x in data(d)]
    p0 = rad2deg(angle(sum(s[keep])))
    @debug "autophase0: p0 = $(round(p0; digits=1))°"
    return phase(spec, di, p0, 0.0)
end

function autophase0(spec::NMRData; kw...)
    fdims = _frequencydims(spec)
    isempty(fdims) &&
        throw(ArgumentError("autophase0: no frequency dimensions to phase; Fourier transform first"))
    return foldl((s, di) -> autophase0(s, di; kw...), fdims; init=spec)
end

autophase0(; kw...) = spec -> autophase0(spec; kw...)
autophase0(dim; kw...) = spec -> autophase0(spec, dim; kw...)
function autophase0(d::AbstractDict{Symbol,<:NMRData}, args...; kw...)
    return Dict(k => autophase0(v, args...; kw...) for (k, v) in d)
end
