# Raw time-domain loading ##########################################################################
#
# `loadfid` reads raw Bruker fid/ser files into time-domain `NMRData` whose element type is
# `Multicomplex`. Quadrature recombination, digital-filter removal and first-point scaling are all
# performed here so that the returned object is a clean time-domain signal ready for the processing
# pipeline (`apodize |> zerofill |> ft |> phase |> realest`).
#
# Handles 1D, 2D and 3D data with any mix of phase-sensitive (States, States-TPPI, echo/anti-echo)
# and pseudo (non-phase-sensitive) indirect dimensions. The multicomplex order equals the number of
# phase-sensitive dimensions; pseudo dimensions become X-dimensions. 4D+ raises a clear error.
#
# NOTE ON CONVENTIONS: several sign conventions in raw-data handling (imaginary sign, digital-filter
# phase-ramp direction, echo/anti-echo recombination) can only be pinned down against a reference
# processing (`loadpdata`). Where a choice has been made it is flagged with a comment so it can be
# validated by the round-trip tests in `test/`.

# FnMODE (acqu2s) → quadrature scheme symbol
const _FNMODE_SCHEMES = Dict(0 => :undefined,
                             1 => :qf,
                             2 => :qseq,
                             3 => :tppi,
                             4 => :states,
                             5 => :states_tppi,
                             6 => :echoantiecho)

# schemes that produce a phase-sensitive (hypercomplex) indirect dimension
const _PHASESENSITIVE = (:states, :states_tppi, :echoantiecho)

"""
    loadfid(filename; removefilter=true, firstpoint=nothing) -> NMRData

Load raw Bruker time-domain data (`fid` for 1D, `ser` for nD) as an [`NMRData`](@ref) with
`Multicomplex` elements and `TimeDimension` axes, ready for processing. `filename` is the `fid` or
`ser` file or the experiment folder.

Phase-sensitive indirect dimensions (States, States-TPPI, echo/anti-echo) are combined into extra
imaginary units, so a 2D experiment has `Multicomplex{Float64,2,4}` elements and even 1D data use
`Multicomplex{Float64,1,2}`. Each time dimension records its imaginary unit in `:mcindex`
(`T1Dim` → 1, `T2Dim` → 2, …). Non-Fourier dimensions such as a relaxation series become
`X2Dim`/`X3Dim` with `:mcindex == nothing`. Intensities match TopSpin's: the raw data are scaled by
`2^NC` as TopSpin scales them, and [`ft`](@ref NMRTools.NMRBase.ft) applies no normalisation.

# Keywords
- `removefilter`: account for the digital filter's group delay. The filter's build-up points stay
  in the FID, on a time axis that starts `GRPDLY` dwell times before zero, and the delay is recorded
  in the direct dimension's `:grpdly` metadata. [`ft`](@ref NMRTools.NMRBase.ft) removes it with a
  linear phase ramp, as TopSpin and NMRPipe do. With `false`, the time axis starts at zero and no
  delay is recorded.
- `firstpoint`: scaling factor for the t = 0 point of each time dimension, as a number for all
  dimensions or a dictionary keyed by dimension type, for example `Dict(T2Dim => 1.0)`. The default
  is 0.5, except in the direct dimension of baseopt data (`DIGMOD = 3`), where the spectrometer
  already gives the t = 0 point its half weight, and 1.0 is used.

FIDs modified by TopSpin's `cryoproc1d` (`DSPFVS` outside 10–23) have had their first points
removed. `loadfid` compares `TD` in `acqu` and `acqus` to find how many points after t = 0 are
missing, starts the time axis there, and warns. The missing points leave a ripple around strong
lines; [`linearpredict`](@ref NMRTools.NMRBase.linearpredict) with `direction=:backward` can rebuild
them.

Annotations ([`annotate`](@ref)) and sample metadata ([`addsampleinfo`](@ref)) are applied as for
[`loadnmr`](@ref).

# Examples
```julia
fid = loadfid("myexperiment/1")
spec = fid |> apodize(Cos2Window()) |> zerofill |> ft |> autophase |> realest
```

See also [`autoprocess`](@ref NMRTools.NMRBase.autoprocess), [`savefid`](@ref).
"""
function loadfid(filename; removefilter=true, firstpoint=nothing)
    filename = abspath(string(filename))

    # 1. locate the binary file and enclosing experiment folder
    if isdir(filename)
        expfolder = filename
        if isfile(joinpath(filename, "ser"))
            datafile = joinpath(filename, "ser")
        elseif isfile(joinpath(filename, "fid"))
            datafile = joinpath(filename, "fid")
        else
            throw(NMRToolsError("loadfid: no fid or ser file found in $filename"))
        end
    else
        datafile = filename
        expfolder = dirname(filename)
        isfile(datafile) || throw(SystemError("loadfid: cannot open $datafile"))
    end

    # 2. load text metadata (acqus, acqu2s, …, title, channels, pulse program)
    md = loadmetadata(expfolder)
    md[:format] = :bruker_fid
    md[:filename] = datafile
    acqusdic = get(md, :acqus, nothing)
    isnothing(acqusdic) &&
        throw(NMRToolsError("loadfid: no acqus metadata found in $expfolder"))

    # 3. read the binary data, recombine quadrature and build the NMRData
    spec = _buildfidspec(datafile, md, acqusdic; removefilter=removefilter,
                         firstpoint=firstpoint)

    # 4. initial noise estimate (refreshed by `ft` once in the frequency domain)
    spec = estimatenoise(spec)

    # 5. resolve pulse-programme annotations and sample information (as loadnmr does)
    spec = annotate(spec)
    spec = addsampleinfo(spec)

    return spec
end

# Binary reading ###################################################################################

"""
    _readfid(datafile, acqusdic, nrows) -> Matrix

Read a raw Bruker `fid`/`ser` file into a `(paddedtd, nrows)` matrix of the on-disk numeric type,
handling endianness (`BYTORDA`) and data type (`DTYPA`). The directly-detected dimension is
zero-padded on disk to a block boundary; the padded row length is recovered from the file size and
the number of stored FIDs (`nrows`), so no block-size assumption is needed. Padding is trimmed by
the caller.
"""
function _readfid(datafile, acqusdic, nrows)
    dtype = get(acqusdic, :dtypa, 0) == 2 ? Float64 : Int32
    bigendian = get(acqusdic, :bytorda, 0) == 1

    raw = Vector{dtype}(undef, filesize(datafile) ÷ sizeof(dtype))
    open(datafile) do f
        return read!(f, raw)
    end
    raw = bigendian ? ntoh.(raw) : ltoh.(raw)

    paddedtd = length(raw) ÷ nrows
    return reshape(raw, paddedtd, nrows)
end

"""
    _directcomplex(raw, td) -> Array{ComplexF64}

Trim the direct-dimension padding (to `td` real points) and deinterleave the real/imaginary points
into a complex array — the im1 quadrature of the directly detected dimension — with `td÷2` points
along the direct axis.
"""
function _directcomplex(raw, td)
    trimmed = selectdim(raw, 1, 1:td)
    re = Float64.(selectdim(trimmed, 1, 1:2:td))
    im_ = Float64.(selectdim(trimmed, 1, 2:2:td))
    return complex.(re, im_)
end

# Assembly #########################################################################################

const _TDIMS = (T1Dim, T2Dim, T3Dim, T4Dim)
const _XDIMS = (X1Dim, X2Dim, X3Dim, X4Dim)

# quadrature scheme for an indirect dimension; unknown/absent FnMODE → :undefined (treated as pseudo)
_scheme(aq) = get(_FNMODE_SCHEMES, Int(get(aq, :fnmode, 0)), :undefined)

"""
    _indirectstorageorder(acqusdic, ndim) -> Vector{Int}

Return the acqus-file indices (`2:ndim`) of the indirect dimensions, ordered from the fastest-varying
(innermost) loop to the slowest (outermost) as stored in the `ser` file. For 3D this depends on the
Bruker `AQSEQ` parameter: `AQSEQ=0` (3-2-1) nests `acqu2s` inside `acqu3s`; `AQSEQ=1` (3-1-2) swaps
them. 1D/2D have a single (or no) indirect dimension.
"""
function _indirectstorageorder(acqusdic, ndim)
    ndim < 3 && return collect(2:ndim)
    return Int(get(acqusdic, :aqseq, 0)) == 1 ? [3, 2] : [2, 3]
end

function _buildfidspec(datafile, md, acqusdic; removefilter, firstpoint)
    # per-dimension acquisition dictionaries (directly detected dimension first)
    acqdicts = Any[acqusdic]
    for k in 2:4
        aq = get(md, Symbol("acqu$(k)s"), nothing)
        isnothing(aq) && break
        push!(acqdicts, aq)
    end
    ndim = length(acqdicts)
    ndim <= 3 || throw(NMRToolsError("loadfid: $(ndim)D raw data is not yet supported"))

    # The ser file nests the indirect loops fastest→slowest; which indirect dimension is the inner
    # loop is set by the acquisition (Bruker AQSEQ for 3D), not fixed. Order the acquisition
    # dictionaries by array axis (axis 1 = direct, axis 2 = innermost indirect, …) so the reshape
    # puts each indirect dimension's size on the correct axis.
    axisacq = Any[acqusdic]
    for k in _indirectstorageorder(acqusdic, ndim)
        push!(axisacq, acqdicts[k])
    end

    # read the binary data and reshape to (paddedtd, td2, td3, …)
    tddirect = Int(acqusdic[:td])
    indirecttds = [Int(axisacq[k][:td]) for k in 2:ndim]
    nrows = isempty(indirecttds) ? 1 : prod(indirecttds)
    raw = _readfid(datafile, acqusdic, nrows)
    # (paddedtd, td2, td3, …); for 1D this collapses to a plain vector
    raw = reshape(raw, size(raw, 1), indirecttds...)
    C = _directcomplex(raw, tddirect)          # ComplexF64 (nc1, td2, td3, …), im1 = direct
    # raw data are stored scaled by 2^-NC; restoring it gives the same intensities as TopSpin
    nc = get(acqusdic, :nc, 0)
    nc isa Real && nc != 0 && (C .*= 2.0^nc)

    # recombine each phase-sensitive indirect dimension into an extra multicomplex order.
    # mcindex counts imaginary units in acquisition order (1 = direct); pseudo-dims get `nothing`.
    y = Multicomplex.(C)                        # order-1 array
    mcindex = Union{Int,Nothing}[1]
    nextunit = 2
    for k in 2:ndim
        if _scheme(axisacq[k]) in _PHASESENSITIVE
            y = _recombine_axis(y, k, _scheme(axisacq[k]))
            push!(mcindex, nextunit)
            nextunit += 1
        else
            push!(mcindex, nothing)
        end
    end

    # Digital-filter group delay (direct dimension only). The lead-in stays in the data at negative
    # times; `ft` removes the delay. First-point scaling acts on the t = 0 point.
    grpdly = removefilter ? _groupdelay(acqusdic, get(md, :experimentfolder, nothing)) : 0.0
    # In baseopt data the t = 0 point already carries its half weight, so by default it is not
    # scaled again; other dimensions and data use 0.5. If the FID starts after t = 0 (negative
    # delay: a modified FID), there is no t = 0 point to scale.
    if grpdly < 0
        firstpoint = merge(Dict{Type,Float64}(T2Dim => _firstpointfactor(T2Dim, firstpoint),
                                              T3Dim => _firstpointfactor(T3Dim, firstpoint)),
                           Dict{Type,Float64}(T1Dim => 1.0))
    elseif isnothing(firstpoint) && grpdly > 0 && _isbaseopt(acqusdic)
        firstpoint = Dict{Type,Float64}(T1Dim => 1.0)
    end
    fourieraxes = [k for k in 1:ndim if !isnothing(mcindex[k])]
    fpfactors = _scalefirstpoint!(y, [k => _TDIMS[mcindex[k]] for k in fourieraxes], firstpoint;
                                  firstindex=Dict(1 => max(round(Int, grpdly), 0) + 1))

    # reorder axes so Fourier dimensions come first (in acquisition order) and pseudo-dimensions
    # last, matching the dimension ordering used by loadnmr/loadpdata (for example the [1,3,2] permute when
    # a pseudo-dimension was acquired before a frequency one)
    pseudoaxes = [k for k in 1:ndim if isnothing(mcindex[k])]
    neworder = vcat(fourieraxes, pseudoaxes)
    neworder == collect(1:ndim) || (y = permutedims(y, neworder))

    # build dimensions in the reordered layout. For a Fourier dimension the new position equals its
    # :mcindex, so the dimension number (Tn) and imaginary unit stay aligned; pseudo-dims trail.
    dimlist = NMRDimension[]
    for (p, k) in enumerate(neworder)
        if isnothing(mcindex[k])
            push!(dimlist,
                  _XDIMS[p](1:size(y, p);
                            metadata=Dict{Symbol,Any}(:label => "", :pseudodim => true)))
        else
            d = _timedim(axisacq[k], p, mcindex[k]; delay=(k == 1 ? grpdly : 0.0))
            k > 1 && (d[:quadrature] = _scheme(axisacq[k]))
            d[:firstpoint] = fpfactors[k]
            k == 1 && (d[:grpdly] = grpdly)
            push!(dimlist, d)
        end
    end

    outmd = _mainmetadata(md, ndim, dimlist)
    outmd[:acqorder] = neworder          # storage-order axis of each array axis, for `savefid`
    return NMRData(y, Tuple(dimlist); metadata=outmd)
end

"""
    _recombine_axis(y, k, scheme) -> Array{<:Multicomplex}

Recombine the interleaved real/imaginary records of a phase-sensitive indirect dimension (array axis
`k`) into an extra multicomplex order, halving that axis. `scheme` selects States, States-TPPI or
echo/anti-echo (Rance–Kay) recombination. The new imaginary unit encodes that dimension's quadrature.
"""
function _recombine_axis(y, k, scheme)
    a = copy(selectdim(y, k, 1:2:size(y, k)))   # "cosine" / echo records
    b = copy(selectdim(y, k, 2:2:size(y, k)))   # "sine" / anti-echo records

    if scheme === :states_tppi
        # π phase shift on alternate increments
        for i in 2:2:size(a, k)
            selectdim(a, k, i:i) .*= -1
            selectdim(b, k, i:i) .*= -1
        end
    elseif scheme === :echoantiecho
        # Rance–Kay recombination (P = echo, N = anti-echo). The 90° shift is in the direct (im1)
        # plane; +im1 gives the correct indirect frequency sign.
        P, N = a, b
        a = P .+ N
        b = (P .- N) .* imN(1)
    end

    return Multicomplex.(a, b)   # order N → N+1
end

# Digital filter removal ###########################################################################

"""
    _groupdelay(acqusdic) -> Float64

The Bruker digital-filter group delay (points) of the direct dimension: the number of points of
filter build-up ("lead-in") that precede time zero in the stored FID. The data are not modified; the
lead-in stays in the FID at negative times and [`ft`](@ref NMRTools.NMRBase.ft) removes the delay with
a linear phase ramp after the Fourier transform. Which value applies is set by the acquisition
parameters:

- `DIGMOD = 0` (analogue filtering): no group delay.
- `DSPFVS` 20–23: `GRPDLY`. For baseopt acquisitions (`DIGMOD = 3`) it is an integer and the
  spectrometer has already corrected the initial points of the FID.
- `DSPFVS` 10–19 (older data): the delay should be looked up from `DECIM` and `DSPFVS`. This is not
  implemented; `GRPDLY` is used if valid, otherwise no delay is applied, with a warning.
- `DSPFVS` outside 10–23: not a firmware version Bruker defines. TopSpin's `cryoproc1d`, for example,
  sets `DSPFVS = 0` after rewriting the FID with its start removed: the filter build-up and some
  points after time zero. The number of points removed is found by comparing `TD` in `acqu` (the
  original) and `acqus`; the points beyond `GRPDLY` were after time zero, so the stored FID starts
  that many dwells after t = 0. This is returned as a negative delay, so that `ft` restores the time
  origin. A warning reports the missing points. If `acqu` does not show how many were removed, no
  delay is applied.
"""
function _groupdelay(acqusdic, expfolder=nothing)
    Int(get(acqusdic, :digmod, 1)) == 0 && return 0.0
    dspfvs = get(acqusdic, :dspfvs, nothing)
    if dspfvs isa Real && !(10 <= dspfvs <= 23)
        offset = _modifiedfidoffset(acqusdic, expfolder)
        if offset > 0
            @warn "loadfid: DSPFVS = $dspfvs indicates the FID has been modified (for example by TopSpin's " *
                  "cryoproc1d), and comparison with acqu shows its first $offset points after " *
                  "t = 0 were removed. The time axis starts at $offset dwell times; the missing " *
                  "points may leave a baseline distortion."
            return -Float64(offset)
        end
        @warn "loadfid: DSPFVS = $dspfvs is outside the range 10–23, which indicates the FID has " *
              "already been modified (for example by TopSpin's cryoproc1d, which removes the start of the " *
              "FID). No group delay will be removed. The first points of the FID may be missing; " *
              "they can be rebuilt by backward linear prediction."
        return 0.0
    end
    grpdly = get(acqusdic, :grpdly, -1.0)
    if isnothing(grpdly) || grpdly < 0
        @warn "loadfid: no valid GRPDLY found; no group delay will be removed (legacy " *
              "DSPFVS/DECIM handling is not yet implemented)"
        return 0.0
    end
    return Float64(grpdly)
end

"""
    _modifiedfidoffset(acqusdic, expfolder) -> Int

For a FID rewritten with its start removed (for example by `cryoproc1d`), the number of complex points
after time zero that were removed: the drop in `TD` from `acqu` (the original) to `acqus`, less the
group delay `GRPDLY`. Zero if this cannot be determined.
"""
function _modifiedfidoffset(acqusdic, expfolder)
    isnothing(expfolder) && return 0
    acqufile = joinpath(expfolder, "acqu")
    isfile(acqufile) || return 0
    original = try
        loadjdx(acqufile)
    catch err
        @debug "loadfid: could not read $acqufile" exception = err
        return 0
    end
    tdorig = get(original, :td, nothing)
    td = get(acqusdic, :td, nothing)
    grpdly = get(acqusdic, :grpdly, nothing)
    (tdorig isa Real && td isa Real && grpdly isa Real && grpdly >= 0) || return 0
    removed = (tdorig - td) ÷ 2
    return max(round(Int, removed - grpdly), 0)
end

"baseopt acquisition (`DIGMOD = 3`): the spectrometer has corrected the initial points of the FID"
_isbaseopt(acqusdic) = Int(get(acqusdic, :digmod, 1)) == 3

# First-point scaling ##############################################################################

"""
scale the t = 0 point of each axis in `axistype_pairs` (point `firstindex[axis]`, default 1); returns
the factor used for each axis
"""
function _scalefirstpoint!(y, axistype_pairs, firstpoint; firstindex=Dict{Int,Int}())
    nd = ndims(y)
    factors = Dict{Int,Float64}()
    for (di, T) in axistype_pairs
        f = _firstpointfactor(T, firstpoint)
        i0 = get(firstindex, di, 1)
        idx = ntuple(i -> i == di ? (i0:i0) : Colon(), nd)
        @views y[idx...] .*= f
        factors[di] = f
    end
    return factors
end

_firstpointfactor(::Type, ::Nothing) = 0.5
_firstpointfactor(T::Type, fp::AbstractDict) = get(fp, T, 0.5)
_firstpointfactor(::Type, fp::Real) = fp

# Dimension / metadata construction ################################################################

"""
    _timedim(aq, di, mcindex) -> TimeDimension

Construct the `TimeDimension` at array position `di` (naming it `T1Dim`…`T4Dim`) from a Bruker
acquisition-parameter dictionary `aq` (`acqus`, `acqu2s`, …), linking it to imaginary unit `mcindex`.
Frequencies are already in Hz (`replacefrequencies!` has run). The time axis is
`t = ((0:td-1) .- delay) / sw`, so points of digital-filter lead-in (`delay` points) have negative
times.
"""
function _timedim(aq, di, mcindex; delay=0.0)
    swhz = Float64(aq[:sw_h])
    bf = Float64(aq[:bf1])
    sf = Float64(aq[:sfo1])
    o1 = Float64(get(aq, :o1, sf - bf))
    offsetppm = o1 / bf * 1e6
    swppm = 1e6 * swhz / bf
    ncomplex = Int(aq[:td]) ÷ 2

    nucstr = string(get(aq, :nuc1, ""))
    nuc = try
        nucleus(nucstr)
    catch
        @debug "loadfid: could not parse nucleus '$nucstr'"
        nothing
    end

    dt = 1 / swhz
    t = range(-delay * dt; step=dt, length=ncomplex)

    dd = Dict{Symbol,Any}(:mcindex => mcindex,
                          :label => nucstr,
                          :nucleus => nuc,
                          :bf => bf,
                          :sf => sf,
                          :swhz => swhz,
                          :swppm => swppm,
                          :offsethz => o1,
                          :offsetppm => offsetppm,
                          :td => ncomplex,
                          :tdzf => ncomplex,
                          :npoints => ncomplex,
                          :region => missing,
                          :window => nothing,
                          :pseudodim => false,
                          :units => "s")
    return _TDIMS[di](t; metadata=dd)
end

"assemble the top-level NMRData metadata from the loaded text metadata and the constructed dims"
function _mainmetadata(md, ndim, axes)
    outmd = deepcopy(md)
    outmd[:ndim] = ndim
    outmd[:nuclei] = Set{Nucleus}()
    for ax in axes
        nuc = metadata(ax, :nucleus)
        isnothing(nuc) || push!(outmd[:nuclei], nuc)
    end
    return outmd
end
