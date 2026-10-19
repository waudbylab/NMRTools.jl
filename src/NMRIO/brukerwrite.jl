# Writing Bruker data ##############################################################################
#
# `savefid` writes time-domain data as a raw Bruker `fid`/`ser`, undoing everything `loadfid` did
# (first-point scaling, digital-filter removal, quadrature recombination, axis reordering), so that
# `loadfid(savefid(path, loadfid(expt)))` reproduces the loaded data. `savenmr` writes
# frequency-domain data as Bruker processed data (`pdata/<procno>/1r`, `2rr`, …).
#
# Parameter files are copied from the experiment the data were loaded from when it is available
# (patching only the parameters that describe the written data), so that TopSpin sees the full
# acquisition parameter set. Otherwise minimal parameter files are generated from the dimension
# metadata.
#
# Sign convention: NMRTools spectra are the complex conjugate of TopSpin's in every imaginary unit
# (the real parts agree), so imaginary components are conjugated on writing processed data, as they
# are on reading.

# JCAMP-DX parameter files ---------------------------------------------------------------------------

_jdxvalue(v::AbstractString) = "<$v>"
_jdxvalue(v::Symbol) = "<$v>"
_jdxvalue(v::Integer) = string(v)
_jdxvalue(v::Real) = string(Float64(v))

"write a minimal JCAMP-DX parameter file from `key => value` pairs (keys as Bruker names)"
function _writejdx(filename, pairs)
    open(filename, "w") do io
        # the version at the end of the title line is read by `topspinversion`
        println(io, "##TITLE= Parameter file, written by NMRTools.jl, TopSpin 4.0.0")
        println(io, "##JCAMPDX= 5.0")
        println(io, "##DATATYPE= Parameter Values")
        println(io, "##ORIGIN= NMRTools.jl")
        println(io, "##OWNER= NMRTools.jl")
        for (k, v) in pairs
            println(io, "##\$", k, "= ", _jdxvalue(v))
        end
        return println(io, "##END=")
    end
    return filename
end

"""
    _patchjdx(text, pairs) -> String

Set single-valued parameters in the text of a JCAMP-DX parameter file, replacing the existing
`##\$KEY= value` line or inserting one before `##END=`.
"""
function _patchjdx(text::AbstractString, pairs)
    for (k, v) in pairs
        line = "##\$$k= $(_jdxvalue(v))"
        pattern = Regex("^##\\\$$k= .*\$", "m")
        if occursin(pattern, text)
            text = replace(text, pattern => line)
        else
            text = replace(text, "##END=" => line * "\n##END=")
        end
    end
    return text
end

function _patchjdxfile(filename, pairs)
    text = replace(read(filename, String), "\r\n" => "\n")
    write(filename, _patchjdx(text, pairs))
    return filename
end

# Experiment folders ---------------------------------------------------------------------------------

"source Bruker experiment folder of `spec`, or `nothing` if unknown or missing"
function _sourceexperiment(spec)
    src = metadata(spec, :experimentfolder)
    (src isa AbstractString && isfile(joinpath(src, "acqus"))) || return nothing
    return src
end

"""
    _copyparameters(src, dest)

Copy the parameter files of experiment folder `src` into `dest`: every regular file except the
binary data (`fid`, `ser`), plus the `lists` folder. Processed data (`pdata`) are not copied.
"""
function _copyparameters(src, dest)
    abspath(src) == abspath(dest) && return dest
    for f in readdir(src)
        path = joinpath(src, f)
        if isfile(path) && f ∉ ("fid", "ser")
            cp(path, joinpath(dest, f); force=true)
        elseif isdir(path) && f == "lists"
            cp(path, joinpath(dest, f); force=true)
        end
    end
    return dest
end

_mhz(hz) = isnothing(hz) ? 0.0 : hz * 1e-6

"Bruker FnMODE for a quadrature scheme"
_fnmode(scheme) = Dict(:qf => 1, :qseq => 2, :tppi => 3, :states => 4, :states_tppi => 5,
                       :echoantiecho => 6)[something(scheme, :states)]

"""
    _minimalacqus(dest, spec, axisdims, tds)

Write minimal `acqus`/`acqu2s`/… files into `dest` describing `spec`. `axisdims` are the dimensions
in storage order (direct first) and `tds` the corresponding TD values (points on disk).
"""
function _minimalacqus(dest, spec, axisdims, tds; dtypa=2, grpdly=0.0)
    for (k, (d, td)) in enumerate(zip(axisdims, tds))
        fourier = !isnothing(metadata(d, :mcindex))
        bf = something(metadata(d, :bf), 0.0)
        offsethz = something(metadata(d, :offsethz), 0.0)
        swhz = something(metadata(d, :swhz), 1.0)
        pairs = Pair{String,Any}["TD" => td,
                                 "SW_h" => swhz,
                                 "SW" => bf > 0 ? 1e6 * swhz / bf : 0.0,
                                 "BF1" => _mhz(bf),
                                 "SFO1" => _mhz(bf + offsethz),
                                 "O1" => offsethz,
                                 "NUC1" => string(something(metadata(d, :label), "off"))]
        if k == 1
            append!(pairs,
                    Pair{String,Any}["DTYPA" => dtypa, "BYTORDA" => 0, "NC" => 0,
                                     "GRPDLY" => grpdly,
                                     "DIGMOD" => grpdly > 0 ? 1 : 0, "AQ_mod" => 3,
                                     "AQSEQ" => 0, "NS" => something(metadata(spec, :ns), 1),
                                     "RG" => something(metadata(spec, :rg), 1),
                                     "PULPROG" => string(something(metadata(spec, :pulseprogram),
                                                                   "unknown")),
                                     "DATE" => round(Int, time())])
        else
            push!(pairs, "FnMODE" => fourier ? _fnmode(metadata(d, :quadrature)) : 0)
        end
        _writejdx(joinpath(dest, k == 1 ? "acqus" : "acqu$(k)s"), pairs)
    end
    return dest
end

"write a pdata title file from the spectrum's title or label"
function _writetitle(pdir, spec)
    title = something(metadata(spec, :title), metadata(spec, :label), "")
    write(joinpath(pdir, "title"), string(title) * "\n")
    return pdir
end

function _checkoverwrite(files, overwrite)
    overwrite && return nothing
    for f in files
        isfile(f) &&
            throw(ArgumentError("$f already exists; pass `overwrite=true` to replace it"))
    end
end

# savefid ----------------------------------------------------------------------------------------------

"""
    _splitunit(y, k, scheme) -> Array

Inverse of `_recombine_axis`: split the highest imaginary unit of `y` back into interleaved
real/imaginary records along array axis `k`, doubling that axis.
"""
function _splitunit(y, k, scheme)
    a = real.(y)      # strips the highest-order imaginary unit
    b = imag.(y)
    if scheme === :states_tppi
        for i in 2:2:size(a, k)
            selectdim(a, k, i:i) .*= -1
            selectdim(b, k, i:i) .*= -1
        end
    elseif scheme === :echoantiecho
        # a = P + N, b = (P − N)·im1  ⇒  P = (a − im1·b)/2, N = (a + im1·b)/2
        i1 = eltype(a)(imN(1))
        P = (a .- i1 .* b) ./ 2
        N = (a .+ i1 .* b) ./ 2
        a, b = P, N
    end
    sz = collect(size(a))
    sz[k] *= 2
    out = similar(a, sz...)
    selectdim(out, k, 1:2:sz[k]) .= a
    selectdim(out, k, 2:2:sz[k]) .= b
    return out
end

"""
    savefid(path, spec; dtype=Float64, overwrite=false) -> String
    savefid(path; kwargs...)

Write time-domain data as a raw Bruker experiment in the folder `path`, as `fid` (1D) or `ser` (nD)
with acquisition parameter files, and return the path of the data file.

The forms without the data argument return a function for use with `|>`.

`spec` is time-domain data such as that from [`loadfid`](@ref), possibly after time-domain
processing such as [`linearpredict`](@ref NMRTools.NMRBase.linearpredict). The changes `loadfid`
makes are undone, so the file has the layout TopSpin expects: t = 0 points are unscaled, the `2^NC`
intensity scaling is removed, phase-sensitive indirect dimensions are re-encoded with their original
quadrature scheme (States, States-TPPI or echo/anti-echo), and the dimensions are put back in
acquisition order. Digital-filter build-up points stay in the data, so they are written back
unchanged. `loadfid` on the written experiment returns the original data.

Parameter files (`acqus`, `acqu2s`, the pulse programme, lists, …) are copied from the experiment
the data were loaded from, with `TD`, `DTYPA`, `BYTORDA` and `GRPDLY` updated. Without a source
experiment, minimal parameter files are generated from the dimension metadata.

# Keywords
- `dtype`: `Float64` (`DTYPA = 2`, as written by Avance NEO spectrometers) or `Int32` (`DTYPA = 0`,
  readable by older TopSpin versions; data are rounded and must fit in 32 bits).
- `overwrite`: replace existing data in `path`.

# Examples
```julia
fid = loadfid("myexperiment/1")
savefid("extended/1", linearpredict(fid, T2Dim))
```

See also [`savenmr`](@ref).
"""
function savefid(path, spec::NMRData; dtype::Type=Float64, overwrite::Bool=false)
    dtype in (Float64, Int32) ||
        throw(ArgumentError("savefid: dtype must be Float64 or Int32, got $dtype"))
    ndim = ndims(spec)
    ndim <= 3 || throw(ArgumentError("savefid: $(ndim)D data are not supported"))
    for (i, d) in enumerate(dims(spec))
        if d isa FrequencyDimension
            throw(ArgumentError("savefid: dimension $i is a frequency dimension; inverse " *
                                "Fourier transform it first, or use savenmr"))
        end
    end
    isnothing(metadata(dims(spec, 1), :mcindex)) &&
        throw(ArgumentError("savefid: the first dimension must be the directly detected dimension"))

    path = abspath(string(path))
    datafile = joinpath(path, ndim == 1 ? "fid" : "ser")
    _checkoverwrite([datafile], overwrite)
    mkpath(path)

    # 1. undo first-point scaling (at t = 0) in array order. Digital-filter build-up points are
    # still in the data before t = 0, so the FID is already in the stored layout.
    y = copy(parent(spec))
    for (i, d) in enumerate(dims(spec))
        f = metadata(d, :firstpoint)
        (isnothing(f) || f == 0 || f == 1) && continue
        i0 = d isa TimeDimension ? argmin(abs.(collect(NMRBase.data(d)))) : 1
        selectdim(y, i, i0:i0) ./= f
    end
    grpdly = something(metadata(dims(spec, 1), :grpdly), 0.0)

    # 2. return to storage (acquisition) order
    acqorder = something(metadata(spec, :acqorder), collect(1:ndim))
    storagedims = Vector{Any}(undef, ndim)
    for (p, k) in enumerate(acqorder)
        storagedims[k] = dims(spec, p)
    end
    acqorder == collect(1:ndim) || (y = permutedims(y, invperm(acqorder)))

    # 3. re-encode phase-sensitive indirect dimensions, highest imaginary unit first
    mcindex = [metadata(d, :mcindex) for d in storagedims]
    for unit in maximum(something.(mcindex, 0)):-1:2
        k = findfirst(==(unit), mcindex)
        y = _splitunit(y, k, metadata(storagedims[k], :quadrature))
    end

    # 4. interleave the direct-dimension real/imaginary points and pad rows to 1024 bytes. Undo the
    # 2^NC scaling applied by loadfid when the source parameters (with their NC) are copied.
    src = _sourceexperiment(spec)
    nc = isnothing(src) ? 0 : get(something(metadata(spec, :acqus), Dict()), :nc, 0)
    c = [ComplexF64(realest(m), imag(m)) * 2.0^(-nc) for m in y]
    raw = zeros(Float64, 2 * size(c, 1), size(c)[2:end]...)
    selectdim(raw, 1, 1:2:size(raw, 1)) .= real.(c)
    selectdim(raw, 1, 2:2:size(raw, 1)) .= imag.(c)
    td = size(raw, 1)
    if ndim > 1
        blk = 1024 ÷ sizeof(dtype)
        padded = cld(td, blk) * blk
        if padded > td
            rawp = zeros(Float64, padded, size(raw)[2:end]...)
            selectdim(rawp, 1, 1:td) .= raw
            raw = rawp
        end
    end
    if dtype === Int32
        maximum(abs, raw) < typemax(Int32) ||
            throw(ArgumentError("savefid: data exceed the Int32 range; use dtype=Float64"))
        raw = round.(Int32, raw)
    end
    write(datafile, htol.(vec(raw)))

    # 5. parameter files
    tds = [td; [size(y, k) for k in 2:ndim]]
    dtypa = dtype === Float64 ? 2 : 0
    if isnothing(src)
        _minimalacqus(path, spec, storagedims, tds; dtypa, grpdly=max(grpdly, 0.0))
    else
        _copyparameters(src, path)
        acqusdic = metadata(spec, :acqus)
        # a negative delay marks a FID already modified in TopSpin (for example by cryoproc1d): its
        # parameters are copied unchanged so that loadfid recognises it again
        patch = Pair{String,Any}["TD" => td, "DTYPA" => dtypa, "BYTORDA" => 0]
        grpdly >= 0 && push!(patch, "GRPDLY" => grpdly)
        _patchjdxfile(joinpath(path, "acqus"), patch)
        # indirect storage axis k is described by acqu file `order[k-1]` (depends on AQSEQ)
        order = _indirectstorageorder(something(acqusdic, Dict()), ndim)
        for (j, fileindex) in enumerate(order)
            f = joinpath(path, "acqu$(fileindex)s")
            isfile(f) && _patchjdxfile(f, Pair{String,Any}["TD" => tds[j + 1]])
        end
    end
    pdir = joinpath(path, "pdata", "1")
    mkpath(pdir)
    isfile(joinpath(pdir, "title")) || _writetitle(pdir, spec)
    return datafile
end

# savenmr ----------------------------------------------------------------------------------------------

const _PDATAFILES = Dict(1 => ["1r", "1i"],
                         2 => ["2rr", "2ir", "2ri", "2ii"],
                         3 => ["3rrr"])

"window function → (WDW, LB, SSB) procs parameters"
_wdwpars(w::ExponentialWindow) = (1, w.lb, 0.0)
_wdwpars(::CosWindow) = (3, 0.0, 2.0)
_wdwpars(::Cos2Window) = (4, 0.0, 2.0)
_wdwpars(w::GeneralSineWindow) = (w.power == 2 ? 4 : 3, 0.0,
                                  w.offset >= 1 ? 0.0 : 1 / (1 - w.offset))
_wdwpars(::Any) = (0, 0.0, 0.0)

"procs parameters for dimension `d` of size `n`"
function _procspars(spec, d, n, nc, isfirst)
    pairs = Pair{String,Any}["SI" => n, "XDIM" => n, "STSI" => n, "STSR" => 0, "FTSIZE" => n]
    if d isa FrequencyDimension
        δ = NMRBase.data(d)
        bf = metadata(d, :bf)
        td = something(metadata(d, :td), n)
        wdw, lb, ssb = _wdwpars(metadata(d, :window))
        ph = something(metadata(d, :phase), (p0=0.0, p1=0.0))
        append!(pairs,
                Pair{String,Any}["FT_mod" => 6,
                                 "TDeff" => 2 * td,
                                 "SF" => _mhz(bf),
                                 "SW_p" => n * abs(δ[2] - δ[1]) * bf * 1e-6,
                                 "OFFSET" => δ[1],
                                 "AXNUC" => string(something(metadata(d, :label), "off")),
                                 "WDW" => wdw, "LB" => lb, "SSB" => ssb,
                                 "PHC0" => ph.p0, "PHC1" => ph.p1])
    else
        append!(pairs, Pair{String,Any}["FT_mod" => 0, "TDeff" => n, "SF" => 0.0, "SW_p" => 0.0,
                                        "OFFSET" => 0.0, "AXNUC" => "off"])
    end
    isfirst && append!(pairs, Pair{String,Any}["NC_proc" => nc, "DTYPP" => 0, "BYTORDP" => 0])
    return pairs
end

"""
    savenmr(path, spec; procno=1, overwrite=false) -> String
    savenmr(path; kwargs...)

Write a frequency-domain spectrum as Bruker processed data, and return the processed-data folder.
`path` is an experiment folder (data go in `path/pdata/<procno>`) or a processed-data folder
`…/pdata/<n>`. 1D, 2D and 3D spectra are supported, including pseudo-dimensions.

The forms without the data argument return a function for use with `|>`.

The real part is written to `1r`, `2rr` or `3rrr`. For multicomplex 1D and 2D spectra the imaginary
parts are written too (`1i`; `2ir`, `2ri`, `2ii`), conjugated to TopSpin's sign convention. Data are
stored as 32-bit integers, with the `NC_proc` exponent chosen to keep full precision. The processing
parameter files, both current (`proc`, `proc2`, …) and status (`procs`, `proc2s`, …), record the axes,
the window function and the phase correction applied (`PHC0`, `PHC1`).

If the experiment folder has no acquisition parameters, they are copied from the experiment the
spectrum was loaded from or, failing that, generated from the dimension metadata, so that
[`loadnmr`](@ref) can read the result.

# Keywords
- `procno`: processed-data number, used when `path` is an experiment folder.
- `overwrite`: replace existing processed data.

# Examples
```julia
spec = loadfid("myexperiment/1") |> autoprocess
savenmr("myexperiment/1", spec; procno=2)
```

See also [`savefid`](@ref).
"""
function savenmr(path, spec::NMRData; procno::Integer=1, overwrite::Bool=false)
    ndim = ndims(spec)
    ndim <= 3 || throw(ArgumentError("savenmr: $(ndim)D data are not supported"))
    for (i, d) in enumerate(dims(spec))
        if d isa TimeDimension && !isnothing(metadata(d, :mcindex))
            throw(ArgumentError("savenmr: dimension $i is a time dimension; Fourier transform it " *
                                "first, or use savefid"))
        end
    end
    dims(spec, 1) isa FrequencyDimension ||
        throw(ArgumentError("savenmr: the first dimension must be a frequency dimension"))

    path = rstrip(abspath(string(path)), '/')
    if basename(dirname(path)) == "pdata"
        pdir = path
        expdir = dirname(dirname(path))
    else
        expdir = path
        pdir = joinpath(path, "pdata", string(procno))
    end

    # components to write: real part, plus imaginary parts when the multicomplex order matches
    A = parent(spec)
    files = _PDATAFILES[ndim]
    comps = if eltype(A) <: Multicomplex && MulticomplexNumbers.order(first(A)) == ndim &&
               length(files) > 1
        # coefficients of 1, im1, (im2, im1·im2); conjugate every imaginary unit for TopSpin
        [[(isodd(count_ones(j - 1)) ? -1 : 1) * m.value[j] for m in A] for j in 1:(2^ndim)]
    else
        [_realpart(A)]
    end
    files = files[1:length(comps)]
    _checkoverwrite(joinpath.(pdir, files), overwrite)
    mkpath(pdir)

    # 32-bit integers with a power-of-two scaling, as TopSpin stores processed data
    maxabs = maximum(c -> maximum(abs, c), comps)
    nc = maxabs > 0 ? ceil(Int, log2(maxabs)) - 29 : 0
    for (f, c) in zip(files, comps)
        write(joinpath(pdir, f), htol.(round.(Int32, vec(c) .* 2.0^(-nc))))
    end

    # TopSpin reads both the current (proc, proc2, …) and status (procs, proc2s, …) parameters
    procfiles = ["proc", "proc2", "proc3"]
    for (i, d) in enumerate(dims(spec))
        pars = _procspars(spec, d, size(spec, i), nc, i == 1)
        for f in (procfiles[i], procfiles[i] * "s")
            _writejdx(joinpath(pdir, f), pars)
        end
    end
    _writetitle(pdir, spec)

    # acquisition parameters, needed by loadnmr
    if !isfile(joinpath(expdir, "acqus"))
        src = _sourceexperiment(spec)
        if isnothing(src)
            tds = [2 * something(metadata(d, :td), size(spec, i)) for (i, d) in enumerate(dims(spec))]
            _minimalacqus(expdir, spec, collect(dims(spec)), tds)
        else
            _copyparameters(src, expdir)
        end
    end
    return pdir
end

_realpart(A::AbstractArray{<:Multicomplex}) = realest.(A)
_realpart(A::AbstractArray{<:Complex}) = real.(A)
_realpart(A::AbstractArray{<:Real}) = A

# curried forms for ending a pipeline
savefid(path; kwargs...) = spec -> savefid(path, spec; kwargs...)
savenmr(path; kwargs...) = spec -> savenmr(path, spec; kwargs...)
