"""
    loadmetadata(experimentfolder) -> Dict{Symbol,Any}

Load all text-based metadata from a Bruker experiment folder without reading
any binary data (ser, fid, 2rr, etc.). Returns a metadata dictionary compatible
with NMRData metadata.
"""
function loadmetadata(experimentfolder::String)
    md = Dict{Symbol,Any}()
    md[:experimentfolder] = experimentfolder

    # parse the acqus file
    acqusfilename = joinpath(experimentfolder, "acqus")
    if isfile(acqusfilename)
        acqusmetadata = parseacqus(acqusfilename)
        md[:acqus] = acqusmetadata
        md[:topspin] = acqusmetadata[:topspin]
        md[:ns] = acqusmetadata[:ns]
        md[:rg] = acqusmetadata[:rg]
        md[:pulseprogram] = acqusmetadata[:pulprog]
        md[:acqusfilename] = acqusfilename
        md[:temperature] = get(acqusmetadata, :te, nothing)
        md[:solvent] = get(acqusmetadata, :solvent, "unknown")
        md[:channels] = _buildchannels(acqusmetadata)
        if haskey(acqusmetadata, :date)
            md[:date] = unix2datetime(acqusmetadata[:date])
        end
    else
        @warn "cannot locate acqus file in $experimentfolder - some metadata will be missing"
    end

    # parse the acquXs files, if they exist, and store in :acquXs dict entries
    for i in 2:4
        acquXsfilename = joinpath(experimentfolder, "acqu$(i)s")
        if isfile(acquXsfilename)
            acquXsmetadata = parseacqus(acquXsfilename, false) # don't parse aux files
            md[Symbol("acqu$(i)s")] = acquXsmetadata
        else
            break
        end
    end

    # load the title file
    titlefilename = joinpath(experimentfolder, "pdata", "1", "title")
    if isfile(titlefilename)
        title = read(titlefilename, String)
        title = replace(title, "\r\n" => "\n") # replace Windows line endings with Unix
        md[:title] = strip(title)
        # use first line of title for experiment label
        titleline1 = split(title, "\n")[1]
        md[:label] = strip(titleline1)
    else
        @debug "cannot locate title file in $experimentfolder"
    end

    return md
end

"""
    _buildchannels(acqusdic) -> Dict{Symbol,Any}

Build a dictionary of spectrometer channels from a parsed acqus dictionary, keyed
by channel label (`:f1`…`:f8`). Each entry is itself a dictionary describing the
nucleus and frequencies on that channel, independent of which dimensions were
detected. Only channels with an assigned nucleus (`nucN` present and not `"off"`)
are included.

Frequencies (`:bf`, `:sfo`) are in Hz — `replacefrequencies!` has already run by
the time the acqus dictionary is built.
"""
function _buildchannels(acqusdic)
    channels = Dict{Symbol,Any}()
    for n in 1:8
        nuc = get(acqusdic, Symbol("nuc$n"), nothing)
        (nuc isa AbstractString && lowercase(strip(nuc)) != "off") || continue

        bf = get(acqusdic, Symbol("bf$n"), nothing)
        sfo = get(acqusdic, Symbol("sfo$n"), nothing)

        ch = Dict{Symbol,Any}()
        ch[:label] = Symbol("f$n")
        ch[:nucleus] = try
            nucleus(nuc)
        catch
            @debug "unable to parse nucleus '$nuc' for channel f$n"
            nothing
        end
        ch[:bf] = bf
        ch[:sfo] = sfo
        if !isnothing(bf) && !isnothing(sfo)
            ch[:offsethz] = sfo - bf
            ch[:offsetppm] = (sfo - bf) * 1e6 / bf
        else
            ch[:offsethz] = nothing
            ch[:offsetppm] = nothing
        end

        channels[Symbol("f$n")] = ch
    end
    return channels
end

function getacqusmetadata(format, filename, experimentfolder=nothing)
    md = Dict{Symbol,Any}()
    md[:format] = format

    # try to locate the folder containing the acqus file
    if isnothing(experimentfolder)
        if format == :nmrpipe
            # assume that acqus is in same directory as test.ft2 files, and one up from ft/test%03d.ft2 files
            experimentfolder = dirname(filename)
            if occursin("%", basename(filename))
                # move up a directory
                experimentfolder = dirname(experimentfolder)
            end
        elseif format == :ucsf
            # assume that acqus is in same directory as ucsf file
            experimentfolder = dirname(filename)
        elseif format == :pdata
            # move up two directories from X/pdata/Y to X
            if isdir(filename)
                if isdirpath(filename)
                    # filename = X/pdata/Y/
                    # dirname = X/pdata/Y
                    experimentfolder = dirname(dirname(dirname(filename)))
                else
                    # filename = X/pdata/Y
                    # dirname = X/pdata
                    experimentfolder = dirname(dirname(filename))
                end
            else
                experimentfolder = dirname(dirname(dirname(filename)))
            end
        end
    else
        # if the passed experiment folder isn't a directory, use the enclosing directory
        if !isdir(experimentfolder)
            @warn "passed experiment folder $experimentfolder is not a directory"
            if isfile(experimentfolder)
                experimentfolder = dirname(experimentfolder)
            end
        end
    end

    md[:filename] = filename
    merge!(md, loadmetadata(experimentfolder))
    md[:format] = format  # restore format (loadmetadata doesn't set it)

    return md
end

function parseacqus(acqusfilename::String, auxfiles=true)
    # auxfiles flag is used for reading acqu2s, acqu3s etc.
    # where we don't want to parse aux files again

    dic = loadjdx(acqusfilename)

    if auxfiles
        # add the topspin version to the dictionary
        dic[:topspin] = topspinversion(acqusfilename)

        # parse the pulse program
        parsepulseprogram!(dic, dirname(acqusfilename))

        # parse pulseprogram lists
        parsepulseprogramlists!(dic, dirname(acqusfilename))

        # check for referenced files like vclist, fq1list, and load these in place of filename
        # (if not already loaded above)
        parseacqusauxfiles!(dic, dirname(acqusfilename))
    end

    return dic
end

"""
    replacepowers!(dic)

Ensure all powers in the acqus dictionary are stored as Power types.
"""
function replacepowers!(dic)
    if haskey(dic, :plw)
        # convert Dict{Int64, Float64} to Dict{Int64, Power}
        dic[:plw] = Dict(k => Power(v, :W) for (k, v) in dic[:plw])
        dic[:pl] = dic[:plw]  # pl is an alias for plw
    end
    if haskey(dic, :spw)
        # convert Dict{Int64, Float64} to Dict{Int64, Power}
        dic[:spw] = Dict(k => Power(v, :W) for (k, v) in dic[:spw])
        dic[:sp] = dic[:spw]  # sp is an alias for spw
    end
end

"""
    replacedurations!(dic)

Ensure all durations in the acqus dictionary are in seconds.
"""
function replacedurations!(dic)
    # ensure all durations are in seconds
    # for dictionaries in us: p, in, inf, inp, pcpd
    dics = [:in, :inf, :inp, :p, :pcpd]
    for d in dics
        if haskey(dic, d) # stored as microseconds
            dic[d] = Dict(k => v * 1e-6 for (k, v) in dic[d])
        end
    end
    # convert single parameters from us to s: DE
    pars = [:de]
    for p in pars
        if haskey(dic, p)
            dic[p] = dic[p] * 1e-6
        end
    end
end

function replacefrequencies!(dic)
    # ensure all frequencies are in Hz

    # convert single parameters from MHz to Hz
    pars = [:bf1, :bf2, :bf3, :bf4, :bf5, :bf6, :bf7, :bf8, :sf,
            :sfo1, :sfo2, :sfo3, :sfo4, :sfo5, :sfo6, :sfo7, :sfo8]
    for p in pars
        if haskey(dic, p)
            dic[p] = dic[p] * 1e6
        end
    end
end

"""
    parsetopspinversion(acqusfilename)

Return the TopSpin version as a VersionNumber (e.g. v"4.2.0").

This is obtained from the end of the first line of the acqus file, e.g.
```
##TITLE= Parameter file, TopSpin 4.2.0
##TITLE= Parameter file, TOPSPIN		Version 2.1
##TITLE= Parameter file, TOPSPIN		Version 3.2
```
"""
function topspinversion(acqusfilename)
    isfile(acqusfilename) ||
        throw(NMRToolsError("getting TopSpin version: $(acqusfilename) is not a valid acqus file"))

    firstline = readlines(acqusfilename)[1]
    version = split(firstline)[end]

    return VersionNumber(version)
end

function parsepulseprogram!(dic, basedir)
    if dic[:topspin] < v"4.0.8"
        pulprogfilename = joinpath(basedir, "pulseprogram")
        if isfile(pulprogfilename)
            prog = read(pulprogfilename, String)
            prog = replace(prog, "\r\n" => "\n") # replace Windows line endings with Unix
            dic[:pulseprogram_precomp] = prog
        end
    elseif dic[:topspin] < v"4.1.4"
        # TS 4.0.8 - pulseprogram / pulseprogram.precomp introduced
        pulprogfilename = joinpath(basedir, "pulseprogram")
        if isfile(pulprogfilename)
            prog = read(pulprogfilename, String)
            prog = replace(prog, "\r\n" => "\n") # replace Windows line endings with Unix
            dic[:pulseprogram_code] = prog
        end
        pulprogfilename = joinpath(basedir, "pulseprogram.precomp")
        if isfile(pulprogfilename)
            prog = read(pulprogfilename, String)
            prog = replace(prog, "\r\n" => "\n") # replace Windows line endings with Unix
            dic[:pulseprogram_precomp] = prog
        end
    else
        # TS 4.1.4 - lists directory introduced
        pulprogfilename = joinpath(basedir, "lists", "pp", dic[:pulprog])
        if isfile(pulprogfilename)
            prog = read(pulprogfilename, String)
            prog = replace(prog, "\r\n" => "\n") # replace Windows line endings with Unix
            dic[:pulseprogram_code] = prog
        end
        pulprogfilename = joinpath(basedir, "pulseprogram.precomp")
        if isfile(pulprogfilename)
            prog = read(pulprogfilename, String)
            prog = replace(prog, "\r\n" => "\n") # replace Windows line endings with Unix
            dic[:pulseprogram_precomp] = prog
        end
    end
    # convenience - link :pulsepgraom to :pulseprogram_precomp
    if haskey(dic, :pulseprogram_precomp)
        dic[:pulseprogram] = dic[:pulseprogram_precomp]
    end
end

function pulseprogram(spec::NMRData; precomp=true)
    if precomp
        acqus(spec, :pulseprogram_precomp)
    else
        acqus(spec, :pulseprogram_code)
    end
end

function pulseprogram(e::NMRExperiment; precomp=true)
    if precomp
        acqus(e, :pulseprogram_precomp)
    else
        acqus(e, :pulseprogram_code)
    end
end

"""
Example lines to parse:
define list<gradient> EA=<EA>
define list<gradient> EA3 = { 1.0000 0.8750 }
define list<pulse> taulist = <\\\$VPLIST>
define list<power> powerlist = <\\\$VALIST>
define list<gradient> diff=<Difframp>
define list<frequency> F19sat = <\\\$FQ1LIST>

NB gradient files have weird syntax!!
"""
function parsepulseprogramlists!(dic, basedir)
    haskey(dic, :pulseprogram_precomp) || return nothing
    pulprog = dic[:pulseprogram_precomp]
    lines = split(pulprog, '\n')

    for line in lines
        if occursin(r"^\s*define list<", line)
            # parse the line
            m = match(r"define list<(\w+)> (\w+)\s*=\s*(.+)", line)
            if m !== nothing
                listtype = m.captures[1]
                listname = m.captures[2]
                listvalue = strip(m.captures[3])
                parselist!(dic, basedir, listtype, listname, listvalue)
            end
        end
    end
end

function parselist!(dic, basedir, listtype, listname, listvalue)
    @debug "Parsing list definition: type='$listtype', name='$listname', value='$listvalue'"

    # Use regex to extract the actual value, handling comments and whitespace
    # Match <$PARAMNAME>, <FILENAME>, or {...} patterns
    if (m = match(r"^<\$(\w+)>", listvalue)) !== nothing
        # e.g. <$VCLIST> - dereference filename using vclist parameter
        paramname = lowercase(m.captures[1])
        # if topspin 3, use this as filename directly, otherwise lookup and load from lists directory
        if dic[:topspin] < v"4.1.4"
            filename = joinpath(basedir, paramname)
        else
            listfile = get(dic, Symbol(paramname), "")
            filename = joinpath(basedir, "lists", listdirectory(listtype), listfile)
        end
        isfile(filename) || return nothing
        lines = readlines(filename)
    elseif (m = match(r"^<([^>$]+)>", listvalue)) !== nothing
        # e.g. <EA> - direct filename reference
        listfile = m.captures[1]
        if dic[:topspin] < v"4.1.4"
            filename = joinpath(basedir, listfile)
        else
            filename = joinpath(basedir, "lists", listdirectory(listtype), listfile)
        end
        isfile(filename) || return nothing
        lines = readlines(filename)
    elseif (m = match(r"^\{([^}]*)\}", listvalue)) !== nothing
        # e.g. { 1.0000 0.8750 } or {0, 40, 20} - inline list
        inner = strip(m.captures[1])
        lines = split(inner, r"[,\s]+")
    else
        @debug "unsupported list value format '$listvalue' in pulseprogram"
        return nothing
    end
    return parselistdata!(dic, listtype, listname, lines)
end

function listdirectory(listtype)
    if listtype == "loopcounter"
        return "vc"
    elseif listtype == "delay"
        return "vd"
    elseif listtype == "pulse"
        return "vp"
    elseif listtype == "power"
        return "va"
    elseif listtype == "frequency"
        return "f1"
    elseif listtype == "gradient"
        return "gp"
    else
        @debug "unsupported list type '$listtype' in pulseprogram"
        return ""
    end
end

function parseacqusauxfiles!(dic, basedir)
    # TS 4.0.8 - pulseprogram / pulseprogram.precomp introduced
    # TS 4.1.4 - lists directory introduced
    if dic[:topspin] < v"4.1.4"
        _parseacqusauxfiles_TS3!(dic, basedir)
    else
        _parseacqusauxfiles_TS4!(dic, basedir)
    end
end

function _parseacqusauxfiles_TS3!(dic, basedir)
    # filenames aren't actually used before TS4 - e.g. vclist is always stored as vclist
    if get(dic, :vclist, "") != "" && get(dic, :vclist, "") isa AbstractString
        filename = joinpath(basedir, "vclist")
        if isfile(filename)
            dic[:vclist] = parsevclist(readlines(filename))
        end
    end

    if get(dic, :vdlist, "") != "" && get(dic, :vdlist, "") isa AbstractString
        filename = joinpath(basedir, "vdlist")
        if isfile(filename)
            dic[:vdlist] = parsevdlist(readlines(filename))
        end
    end

    if get(dic, :vplist, "") != "" && get(dic, :vplist, "") isa AbstractString
        filename = joinpath(basedir, "vplist")
        if isfile(filename)
            dic[:vplist] = parsevplist(readlines(filename))
        end
    end

    if get(dic, :valist, "") != "" && get(dic, :valist, "") isa AbstractString
        filename = joinpath(basedir, "valist")
        if isfile(filename)
            dic[:valist] = parsevalist(readlines(filename))
        end
    end

    for k in
        (:fq1list, :fq2list, :fq3list, :fq4list, :fq5list, :fq6list, :fq7list, :fq8list)
        get(dic, k, "") == "" && continue
        get(dic, k, "") isa AbstractString || continue

        filename = joinpath(basedir, string(k))
        isfile(filename) || continue

        dic[k] = parsefqlist(readlines(filename))
    end
end

function _parseacqusauxfiles_TS4!(dic, basedir)
    if get(dic, :vclist, "") isa AbstractString
        vclistfile = joinpath(basedir, "lists", "vc", get(dic, :vclist, ""))
        if isfile(vclistfile)
            dic[:vclist] = parsevclist(readlines(vclistfile))
        end
    end

    if get(dic, :vdlist, "") isa AbstractString
        vdlistfile = joinpath(basedir, "lists", "vd", get(dic, :vdlist, ""))
        if isfile(vdlistfile)
            dic[:vdlist] = parsevdlist(readlines(vdlistfile))
        end
    end

    if get(dic, :vplist, "") isa AbstractString
        vplistfile = joinpath(basedir, "lists", "vp", get(dic, :vplist, ""))
        if isfile(vplistfile)
            dic[:vplist] = parsevplist(readlines(vplistfile))
        end
    end

    if get(dic, :valist, "") isa AbstractString
        valistfile = joinpath(basedir, "lists", "va", get(dic, :valist, ""))
        if isfile(valistfile)
            dic[:valist] = parsevalist(readlines(valistfile))
        end
    end

    for k in
        (:fq1list, :fq2list, :fq3list, :fq4list, :fq5list, :fq6list, :fq7list, :fq8list)
        if get(dic, k, "") isa AbstractString
            fqlistfile = joinpath(basedir, "lists", "f1", get(dic, k, ""))
            if isfile(fqlistfile)
                dic[k] = parsefqlist(readlines(fqlistfile))
            end
        end
    end

    # lists/pp (pulseprogram)
    # lists/gp
    # lists/vc
    # lists/vd
    # lists/va
    # lists/vp
    # lists/f1 (frequencies)
    # referred to by name in the acqus file
end

function parselistdata!(dic, listtype, listname, lines)
    @debug "Parsing list data for listtype='$listtype', listname='$listname'"
    if listtype == "loopcounter"
        dic[Symbol(listname)] = parsevclist(lines)
    elseif listtype == "delay"
        dic[Symbol(listname)] = parsevdlist(lines)
    elseif listtype == "pulse"
        dic[Symbol(listname)] = parsevplist(lines)
    elseif listtype == "power"
        dic[Symbol(listname)] = parsevalist(lines)
    elseif listtype == "frequency"
        dic[Symbol(listname)] = parsefqlist(lines)
    elseif listtype == "gradient"
        dic[Symbol(listname)] = parsegplist(lines)
    else
        @debug "unsupported list type '$listtype' in pulseprogram"
        return nothing
    end
end

function parsevclist(lines)
    xint = tryparse.(Int, lines)
    if any(xint .== nothing)
        @debug "Unable to parse format of vclist"
        return lines
    end

    return xint
end

"return vdlist contents in seconds"
function parsevdlist(lines)
    # default unit for vdlist is seconds
    lines = replace.(lines, "u" => "e-6")
    lines = replace.(lines, "m" => "e-3")
    lines = replace.(lines, "s" => "")
    xf = tryparse.(Float64, lines)
    if any(xf .== nothing)
        @debug "Unable to parse format of vdlist"
        return lines
    end

    return xf
end

"return vplist contents in seconds"
function parsevplist(lines)
    # default unit for vplist in topspin is microseconds
    # so convert everything to microseconds first
    lines = map(lines) do line
        line = replace(line, "u" => "")
        line = replace(line, "m" => "e3")
        return line = replace(line, "s" => "e6")
    end

    xf = tryparse.(Float64, lines)

    if any(xf .== nothing)
        @debug "Unable to parse format of vplist"
        return lines
    end

    # finally convert microseconds to seconds
    return xf * 1e-6  # return vplist in seconds
end

"return valist contents as Powers"
function parsevalist(lines)
    # power unit must be specified on first lien
    powertoken = popfirst!(lines)

    # parse the rest of the list
    xf = tryparse.(Float64, lines)

    if any(xf .== nothing)
        @debug "Unable to parse format of valist"
        return lines
    end

    # convert to dB if needed
    if powertoken == "Watt"
        return Power.(xf, :W)
    else
        return Power.(xf, :dB)
    end
end

"""
    parsefqlist(lines)

Return contents of the specified fqlist.

fqlists can have several different formats:

first line   | reference   | unit
----------------------------------
             | sfo         | Hz
sfo hz       | sfo         | Hz
sfo ppm      | sfo         | ppm
bf hz        | bf          | Hz
bf ppm       | bf          | ppm
p            | sfo         | ppm
P            | bf          | ppm

"""
function parsefqlist(lines)
    if tryparse(Float64, lines[1]) !== nothing
        # first line is a number => no header line => sfo hz
        unit = :Hz
        relative = true
    else
        firstline = strip(popfirst!(lines))
        if firstline == "p"
            unit = :ppm
            relative = true
        elseif firstline == "P"
            unit = :ppm
            relative = false
        elseif lowercase(firstline) == "bf hz"
            unit = :Hz
            relative = false
        elseif lowercase(firstline) == "bf ppm"
            unit = :ppm
            relative = false
        elseif lowercase(firstline) == "sfo hz"
            unit = :Hz
            relative = true
        elseif lowercase(firstline) == "sfo ppm"
            unit = :ppm
            relative = true
        else
            @warn "Unable to parse format of fqlist: unrecognised header line" firstline
            return lines
        end
    end

    # parse the rest of the list
    xf = tryparse.(Float64, lines)

    if any(xf .== nothing)
        badlines = lines[xf .== nothing]
        @warn "Unable to parse format of fqlist: non-numeric entries" badlines
        return lines
    end

    fqlist = FQList(xf, unit, relative)

    return fqlist
end

function parsegplist(lines)
    # filter out comments (starting with #)
    lines = filter(line -> !startswith(strip(line), "#"), lines)

    # gradient list - unit is percent of max gradient
    xf = tryparse.(Float64, lines)

    if any(xf .== nothing)
        @debug "Unable to parse format of gradient list"
        return lines
    end

    return xf
end