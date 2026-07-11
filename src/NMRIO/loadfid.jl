# Raw time-domain loading ##########################################################################
#
# `loadfid` reads raw Bruker fid/ser files into time-domain `NMRData` whose element type is
# `Multicomplex`. Quadrature recombination, digital-filter removal and first-point scaling are all
# performed here so that the returned object is a clean time-domain signal ready for the processing
# pipeline (`apodize |> zerofill |> ft |> phase |> real`).
#
# Currently implemented: 1D, pseudo-2D (non-phase-sensitive indirect dimension) and phase-sensitive
# 2D (States, States-TPPI, echo/anti-echo). 3D+ raise a clear error until validated.
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
    _readfid(datafile, acqusdic) -> Matrix

Read a raw Bruker `fid`/`ser` file into a `(paddedtd, nrows)` matrix of the on-disk numeric type,
handling endianness (`BYTORDA`) and data type (`DTYPA`). The direct dimension is padded on disk to a
multiple of 1024 bytes; padding is retained here and trimmed by the caller.
"""
function _readfid(datafile, acqusdic)
    dtype = get(acqusdic, :dtypa, 0) == 2 ? Float64 : Int32
    bigendian = get(acqusdic, :bytorda, 0) == 1

    td = acqusdic[:td]                      # real points in the direct dimension
    bytes = sizeof(dtype)
    blockpts = 1024 ÷ bytes                 # 256 (Int32) or 128 (Float64)
    paddedtd = cld(td, blockpts) * blockpts

    raw = Vector{dtype}(undef, filesize(datafile) ÷ bytes)
    open(datafile) do f
        return read!(f, raw)
    end
    raw = bigendian ? ntoh.(raw) : ltoh.(raw)

    nrows = length(raw) ÷ paddedtd
    return reshape(raw, paddedtd, nrows)
end

"""
    _directcomplex(rawmatrix, td) -> Matrix{ComplexF64}

Trim the direct-dimension padding and deinterleave the real/imaginary points of every FID into a
`(td÷2, nrows)` complex matrix (im1 quadrature of the directly detected dimension).
"""
function _directcomplex(rawmatrix, td)
    raw = @view rawmatrix[1:td, :]
    re = Float64.(@view raw[1:2:end, :])
    im_ = Float64.(@view raw[2:2:end, :])
    return complex.(re, im_)
end

# Assembly #########################################################################################

function _buildfidspec(datafile, md, acqusdic; removefilter, firstpoint)
    rawmatrix = _readfid(datafile, acqusdic)
    zc = _directcomplex(rawmatrix, acqusdic[:td])   # (ncomplex1, nrows)

    # dimensionality from the presence of indirect acquXs parameter sets
    acqu2s = get(md, :acqu2s, nothing)
    acqu3s = get(md, :acqu3s, nothing)

    if isnothing(acqu2s)
        return _build1d(zc, md, acqusdic; removefilter=removefilter, firstpoint=firstpoint)
    elseif isnothing(acqu3s)
        return _build2d(zc, md, acqusdic, acqu2s; removefilter=removefilter,
                        firstpoint=firstpoint)
    else
        throw(NMRToolsError("loadfid: 3D and higher raw data are not yet supported " *
                            "(experiment has $(md[:ndim]) dimensions)"))
    end
end

function _build1d(zc, md, acqusdic; removefilter, firstpoint)
    size(zc, 2) == 1 ||
        @warn "loadfid: 1D fid file contains $(size(zc,2)) records; using the first only"
    y = Multicomplex.(vec(zc[:, 1]))      # order-1 multicomplex vector

    removefilter && _removedigitalfilter!(y, acqusdic)
    _scalefirstpoint!(y, [1 => T1Dim], firstpoint)

    xdim = _timedim(acqusdic, 1, T1Dim)
    outmd = _mainmetadata(md, 1, [xdim])
    return NMRData(y, (xdim,); metadata=outmd)
end

function _build2d(zc, md, acqusdic, acqu2s; removefilter, firstpoint)
    scheme = _FNMODE_SCHEMES[get(acqu2s, :fnmode, 0)]

    if scheme in _PHASESENSITIVE
        y = _recombine_indirect(zc, scheme)         # order-2 multicomplex (nc1, nc2)
        removefilter && _removedigitalfilter!(y, acqusdic)
        _scalefirstpoint!(y, [1 => T1Dim, 2 => T2Dim], firstpoint)

        xdim = _timedim(acqusdic, 1, T1Dim)
        ydim = _timedim(acqu2s, 2, T2Dim)
        ydim[:quadrature] = scheme
        outmd = _mainmetadata(md, 2, [xdim, ydim])
        return NMRData(y, (xdim, ydim); metadata=outmd)
    else
        # non-phase-sensitive indirect dimension → pseudo-2D (relaxation/diffusion/CEST series)
        y = Multicomplex.(zc)                        # order-1 multicomplex (nc1, nrows)
        removefilter && _removedigitalfilter!(y, acqusdic)
        _scalefirstpoint!(y, [1 => T1Dim], firstpoint)

        xdim = _timedim(acqusdic, 1, T1Dim)
        n2 = size(y, 2)
        ydim = X2Dim(1:n2; metadata=Dict{Symbol,Any}(:label => "", :pseudodim => true))
        outmd = _mainmetadata(md, 2, [xdim, ydim])
        return NMRData(y, (xdim, ydim); metadata=outmd)
    end
end

"""
    _recombine_indirect(zc, scheme) -> Matrix{<:Multicomplex}

Recombine the interleaved real/imaginary records of a phase-sensitive indirect dimension into a
bicomplex `(ncomplex1, ncomplex2)` matrix, where im1 encodes the direct-dimension quadrature and im2
the indirect-dimension quadrature.
"""
function _recombine_indirect(zc, scheme)
    a = zc[:, 1:2:end]          # "cosine" / echo records
    b = zc[:, 2:2:end]          # "sine" / anti-echo records

    if scheme === :states
        # nothing further to do
    elseif scheme === :states_tppi
        # π phase shift on alternate increments: negate every other complex point
        for k in axes(a, 2)
            if iseven(k)
                @views a[:, k] .*= -1
                @views b[:, k] .*= -1
            end
        end
    elseif scheme === :echoantiecho
        # Rance–Kay recombination (P = echo, N = anti-echo). The +90° combination conjugates the
        # indirect (im2) quadrature relative to −90°, giving the correct indirect frequency sign.
        P = copy(a)
        N = copy(b)
        a = P .+ N
        b = (P .- N) .* im
    end

    return Multicomplex.(a, b)   # Multicomplex(::Complex, ::Complex) → order-2 (rr,ir,ri,ii)
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
    imk = M(im1)

    NMRBase._mcfft!(y, 1, 1)                      # FT along direct dimension (im1)
    θ = 2π * grpdly .* (0:(N - 1)) ./ N           # first-order phase ramp (+ sign: to be validated)
    pf = @. cos(θ) + sin(θ) * imk
    y .= y .* reshape(pf, N, ntuple(_ -> 1, ndims(y) - 1)...)
    NMRBase._mcfft!(y, 1, 1; inverse=true)        # inverse FT
    return y
end

# First-point scaling ##############################################################################

function _scalefirstpoint!(y, dimtypes, firstpoint)
    nd = ndims(y)
    for (di, T) in dimtypes
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
    _timedim(acqusdic, mcindex, DimType) -> TimeDimension

Construct a `TimeDimension` from a Bruker acquisition-parameter dictionary (`acqus`, `acqu2s`, …).
Frequencies are already in Hz (`replacefrequencies!` has run). The time axis is `t = (0:td-1)/sw`.
"""
function _timedim(acqusdic, mcindex, DimType)
    swhz = Float64(acqusdic[:sw_h])
    bf = Float64(acqusdic[:bf1])
    sf = Float64(acqusdic[:sfo1])
    o1 = Float64(get(acqusdic, :o1, sf - bf))
    offsetppm = o1 / bf * 1e6
    swppm = 1e6 * swhz / bf
    ncomplex = acqusdic[:td] ÷ 2

    nucstr = string(get(acqusdic, :nuc1, ""))
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
    return DimType(t; metadata=dd)
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
