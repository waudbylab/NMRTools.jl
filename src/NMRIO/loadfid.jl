# Raw time-domain loading ##########################################################################
#
# `loadfid` reads raw Bruker fid/ser files into time-domain `NMRData` whose element type is
# `Multicomplex`. Quadrature recombination, digital-filter removal and first-point scaling are all
# performed here so that the returned object is a clean time-domain signal ready for the processing
# pipeline (`apodize |> zerofill |> ft |> phase |> real`).
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

Load raw Bruker time-domain data (`fid` for 1D, `ser` for nD) into an [`NMRData`](@ref) with
`Multicomplex` element type and `TimeDimension` axes, ready for processing.

`filename` may point directly at a `fid`/`ser` file or at the enclosing experiment directory.

# Keywords
- `removefilter=true`: remove the Bruker digital-filter group delay (read from `GRPDLY`).
- `firstpoint=nothing`: first-point scaling factor for each time dimension (default `0.5`). Override
  per dimension with a dictionary keyed by dimension type, e.g. `firstpoint=Dict(T2Dim => 1.0)`.

The element type is `Multicomplex{Float64,N,2^N}` where `N` is the number of phase-sensitive
(Fourier) dimensions — even 1D data uses `Multicomplex{Float64,1,2}`. Each `TimeDimension` carries an
`:mcindex` linking it to its imaginary unit (`T1Dim → 1`, `T2Dim → 2`, …); non-Fourier axes such as a
relaxation series use `X2Dim`/`X3Dim` with `:mcindex == nothing`.

Annotations ([`annotate`](@ref)) and sample metadata ([`addsampleinfo`](@ref)) are applied as for
[`loadnmr`](@ref).

See also [`apodize`](@ref), [`zerofill`](@ref), [`ft`](@ref), [`phase`](@ref).
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

    # digital-filter removal (direct dimension only) and first-point scaling (acquisition order)
    removefilter && _removedigitalfilter!(y, acqusdic)
    fourieraxes = [k for k in 1:ndim if !isnothing(mcindex[k])]
    _scalefirstpoint!(y, [k => _TDIMS[mcindex[k]] for k in fourieraxes], firstpoint)

    # reorder axes so Fourier dimensions come first (in acquisition order) and pseudo-dimensions
    # last, matching the dimension ordering used by loadnmr/loadpdata (e.g. the [1,3,2] permute when
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
            d = _timedim(axisacq[k], p, mcindex[k])
            k > 1 && (d[:quadrature] = _scheme(axisacq[k]))
            push!(dimlist, d)
        end
    end

    outmd = _mainmetadata(md, ndim, dimlist)
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
    _removedigitalfilter!(y, acqusdic)

Remove the Bruker digital-filter group delay in place along the direct dimension of `y`. Uses
`GRPDLY` when available (modern data): Fourier transform the direct dimension, apply a first-order
phase ramp of `2π·GRPDLY` across the spectrum, then inverse transform. Legacy `(DSPFVS, DECIM)`
data without `GRPDLY` is left unshifted with a warning.
"""
function _removedigitalfilter!(y, acqusdic)
    grpdly = get(acqusdic, :grpdly, -1.0)
    if isnothing(grpdly) || grpdly < 0
        @warn "loadfid: no valid GRPDLY found; digital filter not removed (legacy DSPFVS/DECIM " *
              "handling is not yet implemented)"
        return y
    end
    grpdly == 0 && return y

    N = size(y, 1)
    M = eltype(y)
    imk = M(imN(1))

    NMRBase._mcfft!(y, 1, 1)                      # FT along direct dimension (im1)
    θ = 2π * grpdly .* (0:(N - 1)) ./ N           # first-order phase ramp (+ sign: to be validated)
    pf = @. cos(θ) + sin(θ) * imk
    y .= y .* reshape(pf, N, ntuple(_ -> 1, ndims(y) - 1)...)
    NMRBase._mcfft!(y, 1, 1; inverse=true)        # inverse FT
    return y
end

# First-point scaling ##############################################################################

function _scalefirstpoint!(y, axistype_pairs, firstpoint)
    nd = ndims(y)
    for (di, T) in axistype_pairs
        f = _firstpointfactor(T, firstpoint)
        idx = ntuple(i -> i == di ? (1:1) : Colon(), nd)
        @views y[idx...] .*= f
    end
    return y
end

_firstpointfactor(::Type, ::Nothing) = 0.5
_firstpointfactor(T::Type, fp::AbstractDict) = get(fp, T, 0.5)
_firstpointfactor(::Type, fp::Real) = fp

# Dimension / metadata construction ################################################################

"""
    _timedim(aq, di, mcindex) -> TimeDimension

Construct the `TimeDimension` at array position `di` (naming it `T1Dim`…`T4Dim`) from a Bruker
acquisition-parameter dictionary `aq` (`acqus`, `acqu2s`, …), linking it to imaginary unit `mcindex`.
Frequencies are already in Hz (`replacefrequencies!` has run). The time axis is `t = (0:td-1)/sw`.
"""
function _timedim(aq, di, mcindex)
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
    t = range(0.0; step=dt, length=ncomplex)

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
