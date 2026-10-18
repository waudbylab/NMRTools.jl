"""
    loadnmr(filename, experimentfolder=nothing)

Main function for loading NMR data. `experimentfolder` contains the path to an experiment directory,
for identification of metadata, if the filename is not directly within an experiment.

Returns an `NMRData` structure, or throws an `NMRToolsError` is there is a problem.

# Examples

nmrPipe import:

```julia
loadnmr("exampledata/1D_1H/1/test.ft1");
loadnmr("exampledata/1D_19F/1/test.ft1");
loadnmr("exampledata/2D_HN/1/test.ft2");
loadnmr("exampledata/pseudo2D_XSTE/1/test.ft1");
loadnmr("exampledata/pseudo3D_HN_R2/1/ft/test%03d.ft2");
```

Bruker pdata import:

```julia
loadnmr("exampledata/1D_19F/1");
loadnmr("exampledata/1D_19F/1/");
loadnmr("exampledata/1D_19F/1/pdata/1");
loadnmr("exampledata/1D_19F/1/pdata/1/");
```

ucsf (Sparky) import:

```julia
loadnmr("exampledata/2D_HN/1/hmqc.ucsf");
loadnmr("exampledata/1D_19F/1/");
loadnmr("exampledata/1D_19F/1/pdata/1");
loadnmr("exampledata/1D_19F/1/pdata/1/");
```
"""
function loadnmr(filename; experimentfolder=nothing, allcomponents=false)
    # 0. replace filename with absolute path
    #    - necessary for checking samples that may be in a parent directory
    if filename isa Number
        # allow passing numbered experiment folders directly
        filename = string(filename)
    end
    filename = abspath(filename)
    if experimentfolder !== nothing
        experimentfolder = abspath(experimentfolder)
    end

    # 1. get format
    format, filename = getformat(filename)

    # raw time-domain data is handled by the dedicated loadfid entry point
    format == :bruker_fid && return loadfid(filename)

    # 2. get acqus metadata
    # TODO - this assumes we're dealing with Bruker data…
    aqdic = getacqusmetadata(format, filename, experimentfolder)

    # 3. load data
    if format == :nmrpipe
        spectrum = loadnmrpipe(filename)
    elseif format == :ucsf
        spectrum = loaducsf(filename)
    elseif format == :pdata
        spectrum = loadpdata(filename, allcomponents)
    else
        # unknown format
        throw(NMRToolsError("unknown file format for loadnmr\nfilename = " * filename))
    end
    _setmcindex!(spectrum)

    # 4. merge in acqus metadata
    merge!(metadata(spectrum), aqdic)
    try
        spectrum[:field] = 2π * metadata(spectrum, 1, :bf) /
                           (gyromagneticratio(H1) *
                            xi_ratio(metadata(spectrum, 1, :nucleus)))
    catch e
        @debug "Could not calculate field strength" exception=e
    end

    # 5. estimate the spectrum noise level
    spectrum = estimatenoise(spectrum)

    # 6. parse and resolve pulse programme annotations
    spectrum = annotate(spectrum)

    # 7. add sample information if available
    spectrum = addsampleinfo(spectrum)

    return spectrum
end

"""
    getformat(filename)

Take an input filename and return either :ucsf, :nmrpipe, :pdata (bruker processed),
:bruker_fid (raw fid/ser, handled by [`loadfid`](@ref)), or :unknown after checking
whether the filename matches any known format.
"""
function getformat(filename)
    # ucsf: match XXX.ucsf
    isucsf = occursin(r"\.ucsf$", filename)
    isucsf && return :ucsf, filename

    # NMRPipe: match XXX.ft, test.ft1, test.ft2, test.ft3, test.ft4, test001.ft1, test0001.ft3, ...
    ispipe = occursin(r"[a-zA-Z0-9]+(\%0[34]d)?\.ft[1234]?$", filename)
    ispipe && return :nmrpipe, filename

    # Bruker raw time-domain data: a path ending in /fid or /ser
    israwfid = occursin(r"/(fid|ser)$", filename)
    israwfid && return :bruker_fid, filename

    # Bruker: match pdata/1, 123/pdata/101, 1/pdata/23/, etc...
    # NB. This relies on identifying the pdata folder, so would fail if the working directory is already the pdata
    isbruker = occursin(r"pdata/[0-9]+/?$", filename)
    isbruker && return :pdata, filename

    # No pdata directory specified? Test if filename/pdata/1 exists
    # If it does, return the updated path filename/pdata/1
    if isdir(filename)
        testpath = joinpath(filename, "pdata", "1")
        isdir(testpath) && return :pdata, testpath
    end

    return :unknown, filename
end

"""
    loadnmr(expt::NMRExperiment; pdatadir=1, allcomponents=false) -> NMRData

Load the full NMR spectrum for an experiment previously discovered by [`scanexperiments`](@ref NMRTools.NMRIO.scanexperiments).
"""
function loadnmr(expt::NMRExperiment; pdatadir::Int=1, allcomponents::Bool=false)
    return loadnmr(joinpath(expt.path, "pdata", string(pdatadir)); allcomponents)
end

"""
    _setmcindex!(spec)

Link each frequency dimension of a loaded spectrum to its imaginary unit (`:mcindex`), in order:
`im1` for the direct dimension, `im2` for the next frequency dimension, and so on, matching the
components of data loaded with `allcomponents=true` and spectra produced by `ft`.
"""
function _setmcindex!(spec)
    k = 0
    for d in dims(spec)
        d isa FrequencyDimension || continue
        metadata(d)[:mcindex] = (k += 1)
    end
    return spec
end
