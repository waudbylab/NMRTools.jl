"""
    _pdatafiles(ndim)

Names of the pdata files of an `ndim`-dimensional spectrum, in the order of the multicomplex
components: file `k` (counting from 0) has an `i` for each set bit of `k`, the direct dimension
first, so `2ir` is the coefficient of `im1`.
"""
function _pdatafiles(ndim::Integer)
    return [string(ndim, join(isodd(k >> (d - 1)) ? 'i' : 'r' for d in 1:ndim))
            for k in 0:(2^ndim - 1)]
end

"""
    loadpdata(filename, allcomponents=false)

Load processed Bruker NMR data from a pdata folder.

# Arguments
- `filename`: Path to a pdata directory or a file within it
- `allcomponents`: If `true`, returns `Multicomplex` data with all components, one imaginary unit
  per frequency dimension (`im1` for the direct dimension), for example
  `2rr + 2ir*im1 + 2ri*im2 + 2ii*im1*im2` for 2D data (TopSpin's `2ir` is imaginary in the direct
  dimension, `2ri` in the indirect one). Default is `false` (returns real data only).

Imaginary components are conjugated relative to TopSpin's files, so that loaded data follow the
same convention as spectra produced by [`ft`](@ref NMRTools.NMRBase.ft): TopSpin's `PHC0`/`PHC1`
applied with [`phase`](@ref NMRTools.NMRBase.phase) give the same spectrum.
"""
function loadpdata(filename, allcomponents=false)
    # 1. get reference to pdata/X directory
    if isdir(filename)
        pdir = filename
    else
        pdir = dirname(filename)
    end
    isdir(pdir) ||
        throw(NMRToolsError("can't load bruker data, pdata directory $pdir does not exist"))

    # 2. get a list of pdata input files (which depends on the dimension of the spectrum)
    ndim = 0
    datafiles = []
    if isfile(joinpath(pdir, "1r"))
        ndim = 1
    elseif isfile(joinpath(pdir, "2rr"))
        ndim = 2
    elseif isfile(joinpath(pdir, "3rrr"))
        ndim = 3
    else
        throw(NMRToolsError("can't load bruker data, pdata directory $pdir does not contain binary data files of supported dimensions (1D, 2D, 3D)"))
    end
    ndim > 0 ||
        throw(NMRToolsError("can't load bruker data, pdata directory $pdir does not contain binary data files (1r/2rr/3rrr etc.)"))
    datafiles = _pdatafiles(ndim)

    # 2b. only include the realest component unless requested
    if !allcomponents
        datafiles = [first(datafiles)]
    end

    # 2c. check these files actually exist
    datafiles = map(x -> joinpath(pdir, x), datafiles)
    filter!(isfile, datafiles)

    # 3. read procs files (containing axis metadata)
    procsfiles = ["procs", "proc2s", "proc3s", "proc4s"]
    procsdics = [loadjdx(joinpath(pdir, procsfiles[i])) for i in 1:ndim]

    # 4. TODO parse procs into main axis metadata
    md = Dict{Symbol,Any}()  # main dictionary for metadata
    md[:ndim] = ndim
    # populate metadata for each dimension
    axesmd = []
    for i in 1:ndim
        dic = Dict{Symbol,Any}()
        procdic = procsdics[i]

        # add some data in a nice format
        if procdic[:ft_mod] == 0
            dic[:pseudodim] = true
            dic[:label] = ""
            dic[:npoints] = procdic[:tdeff] # NB use this and not SI to avoid regions of zeros - trim them out later
            dic[:val] = 1:dic[:npoints] # coordinates for this dimension are just 1 to N
        else # frequency domain
            dic[:pseudodim] = false
            dic[:label] = procdic[:axnuc]
            dic[:nucleus] = try
                nucleus(procdic[:axnuc]) # convert to Nucleus enum
            catch e
                @warn "unable to parse nucleus from $(procdic[:axnuc])"
                nothing
            end
            dic[:bf] = procdic[:sf]  # NB this includes the addition of SR
            # see bruker processing reference, p. 85
            # (SR = SF - BF1, where SF is proc par and BF1 is acqu par)

            # extleftppm = procdic[:offset] # edge of extracted region, in ppm
            extswhz = procdic[:sw_p] # sw of extracted region, in Hz
            extswppm = 1e6 * extswhz / dic[:bf]

            # full sw = ftsize / stsi * sw of extracted region
            dic[:swhz] = procdic[:ftsize] / procdic[:stsi] * procdic[:sw_p]
            dic[:swppm] = 1e6 * dic[:swhz] / dic[:bf]

            # chemical shift at point STSR = offset ppm
            # ppm per point = full sw / ftsize
            # midpoint = ftsize ÷ 2
            # true offset = offset + (STSR-midpoint)*ppmperpoint
            extoffset = procdic[:offset]
            ppmperpoint = dic[:swppm] / procdic[:ftsize]
            midpoint = procdic[:ftsize] ÷ 2
            offsetppm = extoffset + (procdic[:stsr] - midpoint) * ppmperpoint

            dic[:offsetppm] = offsetppm #dic[:offsethz] / dic[:bf]
            dic[:offsethz] = 1e-6 * offsetppm * dic[:bf] #procdic[:offset]*dic[:bf] - dic[:swhz]/2
            dic[:sf] = dic[:offsethz] + dic[:bf]

            # if procdic[:lpbin] == 0
            dic[:td] = procdic[:tdeff] ÷ 2  # number of COMPLEX points, after LP, before ZF
            # else
            #     dic[:td] = procdic[:lpbin] ÷ 2
            # end
            dic[:tdzf] = procdic[:ftsize] #÷ 2  # number of COMPLEX points
            dic[:aq] = dic[:td] / dic[:swhz]

            # get any extracted regions and number of points
            dic[:npoints] = procdic[:stsi]
            if procdic[:stsr] == 0
                dic[:region] = missing
            else
                dic[:region] = procdic[:stsr]:(procdic[:stsr] + procdic[:stsi] - 1)
            end

            # chemical shift values
            δ = LinRange(procdic[:offset], procdic[:offset] - extswppm, dic[:npoints] + 1)
            dic[:val] = δ[1:(end - 1)]

            # create a representation of the window function
            # calculate acquisition time = td / 2*sw
            w = procdic[:wdw]
            if w == 0
                window = NullWindow(dic[:aq])
            elseif w == 1
                window = ExponentialWindow(procdic[:lb], dic[:aq])
            elseif w == 2
                @warn("Gaussian window not yet implemented")
                window = UnknownWindow(dic[:aq])
            elseif w == 3
                ssb = procdic[:ssb]
                if ssb < 1
                    ssb = 1
                end
                window = SineWindow(1 - 1.0 / ssb, 1, 1, dic[:aq]) # offset, endpoint, power
            elseif w == 4
                ssb = procdic[:ssb]
                if ssb < 1
                    ssb = 1
                end
                window = SineWindow(1 - 1.0 / ssb, 1, 2, dic[:aq]) # offset, endpoint, power
            else
                window = UnknownWindow(dic[:aq])
            end
            dic[:window] = window
        end

        dic[:procs] = procsdics[i] # store procs file in main axis dictionary
        push!(axesmd, dic)
    end

    # create set of nuclei across all axes
    md[:nuclei] = Set{Nucleus}()
    for ax in axesmd
        if get(ax, :nucleus, nothing) !== nothing
            push!(md[:nuclei], ax[:nucleus])
        end
    end

    # 5. determine shape and submatrix size
    shape = [procs[:si] for procs in procsdics]
    submatrix = [procs[:xdim] for procs in procsdics]

    # 6. check data format (default to Int32)
    dtype = get(procsdics[1], :dtypp, 1) == 2 ? Float64 : Int32

    # 7. check endianness (default to little)
    endian = get(procsdics[1], :bytordp, 0) == 1 ? "b" : "l"

    # 8. read the data files
    dat = map(datafile -> readpdatabinary(datafile, shape, submatrix, dtype, endian),
              datafiles)

    # 9. combine into real or multicomplex output. Pseudo-dimensions have no imaginary files, so the
    #    files present are the components over the frequency dimensions, in the same order.
    if length(dat) == 1
        allcomponents &&
            @warn "unable to find imaginary data for import of $filename - returning real values only"
        y = first(dat)
    else
        ispow2(length(dat)) ||
            throw(NMRToolsError("incomplete set of imaginary data files for $filename"))
        # conjugate in every imaginary unit (NMRTools convention): the sign alternates with the
        # number of imaginary units in each component
        y = Multicomplex.((iseven(count_ones(k - 1)) ? d : -d for (k, d) in enumerate(dat))...)
    end

    # 10. scale
    y = one(Float64) * y
    scalepdata!(y, procsdics[1])

    # 11. trim out any regions of zero data
    y = y[[1:axesmd[i][:npoints] for i in 1:ndim]...]

    # 12. form NMRData and return
    if ndim == 1
        valx = axesmd[1][:val]
        delete!(axesmd[1], :val) # remove values from metadata to prevent confusion when slicing up

        xaxis = F1Dim(valx; metadata=axesmd[1])

        NMRData(y, (xaxis,); metadata=md)
    elseif ndim == 2
        val1 = axesmd[1][:val]
        val2 = axesmd[2][:val]
        delete!(axesmd[1], :val) # remove values from metadata to prevent confusion when slicing up
        delete!(axesmd[2], :val)

        xaxis = F1Dim(val1; metadata=axesmd[1])

        ax2 = axesmd[2][:pseudodim] ? X2Dim : F2Dim
        yaxis = ax2(val2; metadata=axesmd[2])

        NMRData(y, (xaxis, yaxis); metadata=md)
    elseif ndim == 3
        # rearrange data into a useful order - always place pseudo-dimension last - and generate axes
        # 1 is always direct, and should be placed first (i.e. 1 x x)
        # if is there is a pseudodimension, put that last (i.e. 1 y p)
        pdim = [axesmd[i][:pseudodim] for i in 1:3]

        val1 = axesmd[1][:val]
        val2 = axesmd[2][:val]
        val3 = axesmd[3][:val]
        delete!(axesmd[1], :val) # remove values from metadata to prevent confusion when slicing up
        delete!(axesmd[2], :val)
        delete!(axesmd[3], :val)

        if pdim[2]
            if pdim[3]
                # two pseudodimensions - keep default ordering
                xaxis = F1Dim(val1; metadata=axesmd[1])
                yaxis = X2Dim(val2; metadata=axesmd[2])
                zaxis = X3Dim(val3; metadata=axesmd[3])
            else
                # dimensions are x p y => we want ordering 1 3 2
                xaxis = F1Dim(val1; metadata=axesmd[1])
                yaxis = F2Dim(val3; metadata=axesmd[3])
                zaxis = X3Dim(val2; metadata=axesmd[2])
                y = permutedims(y, [1, 3, 2])
            end
        elseif pdim[3]
            # dimensions are x y p => we want ordering 1 2 3
            xaxis = F1Dim(val1; metadata=axesmd[1])
            yaxis = F2Dim(val2; metadata=axesmd[2])
            zaxis = X3Dim(val3; metadata=axesmd[3])
        else
            # no pseudodimension, use Z axis not Ti
            # dimensions are x y z => we want ordering 1 2 3
            xaxis = F1Dim(val1; metadata=axesmd[1])
            yaxis = F2Dim(val2; metadata=axesmd[2])
            zaxis = F3Dim(val3; metadata=axesmd[3])
        end

        NMRData(y, (xaxis, yaxis, zaxis); metadata=md)
    else
        throw(NMRToolsError("loadnmr cannot currently handle 4D+ experiments"))
    end
end

function scalepdata!(y, procs, reverse=false)
    if reverse
        y .*= 2.0^(-get(procs, :nc_proc, 1))
    else
        y .*= 2.0^(get(procs, :nc_proc, 1))
    end
end

function readpdatabinary(filename, shape, submatrix, dtype, endian)
    # preallocate header
    np = foldl(*, shape)
    y = zeros(dtype, np)

    # read the file
    open(filename) do f
        return read!(f, y)
    end
    convertor = endian == "b" ? ntoh : ltoh # get convertor function for big/little-to-host
    y = convertor(y)

    return reordersubmatrix(y, shape, submatrix)
end

function reordersubmatrix(yin, shape, submatrix, reverse=false)
    ndim = length(shape)
    if ndim == 1
        # do nothing to 1D data
        return yin
    end

    nsub = @. Int(ceil(shape / submatrix))
    nsubs = foldl(*, nsub)
    N = length(yin)
    sublength = foldl(*, submatrix)

    subshape = vcat(nsubs, submatrix)
    if reverse
        yout = zeros(eltype(yin), N)
    else
        yout = zeros(eltype(yin), shape...)
    end

    ci = CartesianIndices(ones(nsub...))
    for n in 1:nsubs
        idx = ci[n]
        slices = [(1 + (idx[i] - 1) * submatrix[i]):(idx[i] * submatrix[i]) for i in 1:ndim]
        vecslice = (1 + (n - 1) * sublength):(n * sublength)
        if reverse
            yout[vecslice] = yin[slices...]
        else
            yout[slices...] = yin[vecslice]
        end
    end

    return yout
end
