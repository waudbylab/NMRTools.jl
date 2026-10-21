# Automatic processing #############################################################################

# dimension number of a dimension type, for per-dimension settings given by type
const _DIMNUMBER = Dict{Type,Int}(T1Dim => 1, F1Dim => 1, T2Dim => 2, F2Dim => 2,
                                  T3Dim => 3, F3Dim => 3, T4Dim => 4, F4Dim => 4)

"""
    _perdim(setting, di)

The value of a processing setting for array dimension `di`: `setting` itself, or, if it is a
dictionary, its entry for `di` keyed by dimension number or dimension type (`T1Dim` and `F1Dim` both
mean dimension 1). Returns `default` when a dictionary has no entry.
"""
function _perdim(setting, di; default=nothing)
    setting isa AbstractDict || return setting
    for (k, v) in setting
        (k isa Integer ? k : get(_DIMNUMBER, k, 0)) == di && return v
    end
    return default
end

"whether to apply solvent suppression: `solvent` is `true`, a shift, or `:auto` for ¹H in H₂O"
function _usesolventfilter(spec, solvent)
    solvent === :auto || return solvent !== false && !isnothing(solvent)
    sv = metadata(spec, :solvent)
    isnothing(sv) && return false
    nuc = metadata(dims(spec, 1), :nucleus)
    return occursin("h2o", lowercase(string(sv))) && nuc == H1
end

"""
    autoprocess(fid; window=Cos2Window(), zerofill=2, solvent=:auto, phase=:auto, rebuild=false,
                real=true)
    autoprocess(; kw...)

Process time-domain data from [`loadfid`](@ref NMRTools.NMRIO.loadfid) into a spectrum in one step.
The second form returns a function for use with `|>`. The steps, in order:

1. With `rebuild=true`, rebuild points missing at the start of a FID that begins after t = 0 (as
   `loadfid` reports for data modified by TopSpin's `cryoproc1d`), by backward
   [`linearpredict`](@ref).
2. Suppress the water signal in the direct dimension with [`solventfilter`](@ref).
3. In each Fourier dimension: [`apodize`](@ref), [`zerofill`](@ref) and [`ft`](@ref).
4. Phase each frequency dimension, direct dimension first, with [`autophase`](@ref) or given phases.
   A dimension in which `autophase` finds no lines is left unphased, with a warning. When the
   solvent was filtered from a 1D spectrum, its phase is estimated from the unfiltered spectrum,
   leaving out ±0.2 ppm around the solvent: the solvent's dispersion tails at the edges of the
   spectrum help the estimate, and the filter removes them.

The result is the real spectrum (`Float64`, real in every dimension, as from `realest`). Pass
`real=false` to keep the `Multicomplex` spectrum, for example to phase it further. Non-Fourier
dimensions (relaxation delays, offsets, …) are left unchanged.

# Keywords
- `window`: window function, for example `ExponentialWindow(2.0)` or
  `GeneralSineWindow(0.33, 1, 2)`.
- `zerofill`: zero-filling factor (`1` for none).
- `solvent`: `:auto` suppresses water in ¹H spectra whose solvent contains H₂O; `true` suppresses a
  signal at the carrier, a number suppresses one at that chemical shift (ppm), and `false` skips
  suppression.
- `phase`: `:auto` runs `autophase`; a tuple `(p0, p1)`, or a single `p0`, applies those phases
  (degrees); `nothing` skips phasing.
- `real`: return the real spectrum (default) rather than the `Multicomplex` one.
- `rebuild`: rebuild points missing after t = 0 by backward linear prediction. Missing points leave
  a ripple around strong lines; rebuilding them removes it, but can add a broad baseline distortion
  when the FID contains fast-decaying signals, so it is off by default.

`window`, `zerofill` and `phase` can also be given per dimension as a dictionary keyed by dimension
number or type (`T1Dim` and `F1Dim` both mean the first dimension). Dimensions missing from a
dictionary get the default.

# Examples
```julia
spec = loadfid("expt/1") |> autoprocess
spec = loadfid("expt/1") |> autoprocess(window=ExponentialWindow(1.0), solvent=false)
spec = autoprocess(loadfid("hsqc/1"); window=Dict(T1Dim => Cos2Window(), T2Dim => Cos2Window()),
                   phase=Dict(F2Dim => (0.0, 0.0)))
```

See also [`autophase`](@ref), [`apodize`](@ref), [`zerofill`](@ref), [`ft`](@ref).
"""
function autoprocess(fid::NMRData; window=Cos2Window(), zerofill=2, solvent=:auto, phase=:auto,
                     rebuild::Bool=false, real::Bool=true)
    spec = fid
    # 1. rebuild points missing after t = 0
    gap = -_groupdelaypoints(dims(spec, 1))
    if rebuild && gap > 0
        n = size(spec, 1)
        spec = linearpredict(spec, 1; direction=:backward, npoints=n + round(Int, gap))
    end

    # 2. solvent suppression in the direct dimension
    unfiltered = nothing
    if _usesolventfilter(spec, solvent)
        unfiltered = spec
        if solvent isa Real && !(solvent isa Bool)
            target = Float64(solvent)
            spec = solventfilter(spec, 1; shift=solvent)
        else
            target = metadata(dims(spec, 1), :offsetppm)
            spec = solventfilter(spec, 1)
        end
    end

    # 3. window, zero-fill and transform each Fourier dimension
    spec = _transformfourier(spec, window, zerofill)

    # 4. phase, direct dimension first
    for di in _fourierdims(spec)
        p = _perdim(phase, di; default=:auto)
        if p === :auto
            try
                if di == 1 && !isnothing(unfiltered) && length(_fourierdims(spec)) == 1
                    # The solvent's dispersion tails at the edges of a 1D spectrum are a precise
                    # phase reference, which the filter removes. Phase the direct dimension from
                    # the unfiltered spectrum instead, leaving out the solvent line itself.
                    # (Multidimensional spectra are phased from traces through their cross
                    # peaks, which do not use the edges.)
                    excl = [(target - SOLVENT_EXCLUDE_PPM) .. (target + SOLVENT_EXCLUDE_PPM)]
                    ref = autophase(_transformfourier(unfiltered, window, zerofill), 1; exclude=excl)
                    ph = metadata(dims(ref, 1), :phase)
                    spec = NMRBase.phase(spec, 1, ph.p0, ph.p1)
                else
                    spec = autophase(spec, di)
                end
            catch err
                err isa NMRToolsError || rethrow()
                @warn "autoprocess: dimension $di left unphased: $(err.message)"
            end
        elseif p isa Real
            spec = NMRBase.phase(spec, di, p)
        elseif p isa Tuple
            spec = NMRBase.phase(spec, di, p...)
        elseif !isnothing(p)
            throw(ArgumentError("autoprocess: phase must be :auto, a number, a tuple (p0, p1) " *
                                "or nothing, got $p"))
        end
    end
    return real ? realest(spec) : spec
end

autoprocess(; kw...) = fid -> autoprocess(fid; kw...)

"half-width (ppm) of the region around a filtered solvent left out of automatic phasing"
const SOLVENT_EXCLUDE_PPM = 0.2

"apodize, zero-fill and Fourier transform each Fourier time dimension, as in `autoprocess`"
function _transformfourier(spec, window, zerofill)
    for di in _fourierdims(spec)
        dims(spec, di) isa TimeDimension || continue
        w = _perdim(window, di; default=Cos2Window())
        isnothing(w) || (spec = apodize(spec, di, w))
        zf = _perdim(zerofill, di; default=2)
        zf > 1 && (spec = NMRBase.zerofill(spec, di; factor=zf))
        spec = ft(spec, di)
    end
    return spec
end
