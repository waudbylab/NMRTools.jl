# Peak lists, in the `peaklist.csv` format written and read by NMRAnalysis.jl (see its
# `docs/src/advanced/conventions.md`), and in Sparky's `.list` format.
#
#   # comment lines carry provenance
#   label,plane,x (ppm),y (ppm),lwx (Hz),lwy (Hz),intensity
#   G10,,8.2113,120.402,21.3,14.8,1.2e6
#
# `x` is the direct dimension and `y` the indirect one, matching the dimension order of an
# NMRData. Column names carry their units in parentheses, a blank `plane` means the row
# applies to every plane of a series, and a value that does not exist is written `NA`.
# NMRAnalysis reads the `label`, `x` and `y` columns and ignores the rest.

"""
    loadpeaks(filename) -> Vector{Peak{N}}
    loadpeaks(filename, spec) -> Vector{Peak{N}}

Read a peak list: either a `peaklist.csv` as written by [`savepeaks`](@ref) and by
NMRAnalysis.jl, a hand-made list of `label x y` lines, or a Sparky `.list` file.

CSV files are matched by column name, so column order and extra columns do not matter.
Linewidths and intensities are read where present and are otherwise `NaN`. A list in
which peaks were tracked plane by plane (with a `plane` column filled in) is read as one
position per peak, taking each peak's first row.

Sparky lists name their dimensions `w1`, `w2`, … with the indirect dimension first,
which is the reverse of NMRTools' order. Passing the spectrum `spec` lets the axis order
be checked against where the shifts actually fall, so a list written the other way round
is still read correctly.

# Examples
```julia
pks = loadpeaks("out/peaklist.csv")
pks = loadpeaks("hsqc.list", spec)
```
"""
loadpeaks(filename::AbstractString) = loadpeaks(filename, nothing)

function loadpeaks(filename::AbstractString, spec)
    isfile(filename) || throw(ArgumentError("no such peak list: $filename"))
    pks = _issparky(filename) ? _readsparky(filename, spec) : _readpeakcsv(filename)
    if !isnothing(spec) && !isempty(pks) && ndims(first(pks)) != ndims(spec)
        throw(DimensionMismatch("$filename holds $(ndims(first(pks)))D peaks but the " *
                                "spectrum is $(ndims(spec))D"))
    end
    @debug "Read $(length(pks)) peaks from $filename"
    return pks
end

"""
    savepeaks(filename, peaks)

Write a peak list. A filename ending in `.list` is written in Sparky format; anything else
is written as a `peaklist.csv` that NMRAnalysis.jl reads directly, with positions in ppm,
linewidths in Hz and intensities at full precision.

Peak labels identify peaks in NMRAnalysis, so they must be unique and may not contain
commas. Unlabelled peaks are written as `X1`, `X2`, … by their position in the list.

# Examples
```julia
pks = detectpeaks(spec)
savepeaks("peaklist.csv", pks)
savepeaks("hsqc.list", pks)
```
"""
function savepeaks(filename::AbstractString, pks::AbstractVector{<:Peak})
    labels = _savelabels(pks)
    mkpath(dirname(abspath(filename)))
    if lowercase(splitext(filename)[2]) == ".list"
        _writesparky(filename, pks, labels)
    else
        _writepeakcsv(filename, pks, labels)
    end
    @debug "Wrote $(length(pks)) peaks to $filename"
    return filename
end

function _savelabels(pks)
    labels = [isempty(label(p)) ? "X$i" : label(p) for (i, p) in enumerate(pks)]
    i = findfirst(l -> occursin(',', l) || any(isspace, l), labels)
    isnothing(i) ||
        throw(ArgumentError("peak label $(repr(labels[i])) contains a comma or space"))
    length(unique(labels)) == length(labels) ||
        throw(ArgumentError("peak labels must be unique: " *
                            join(unique(filter(l -> count(==(l), labels) > 1, labels)), ", ")))
    return labels
end

# ---- peaklist.csv -------------------------------------------------------------------

_csvvalue(x::Real) = isnan(x) ? "NA" : string(Float64(x))

function _writepeakcsv(filename, pks::AbstractVector{<:Peak{N}}, labels) where {N}
    names = NMRBase._peakcolumnnames(N)
    header = ["label"; "plane"; ["$n (ppm)" for n in names]; ["lw$n (Hz)" for n in names];
              "intensity"]
    open(filename, "w") do io
        println(io, "# NMRTools.jl $(pkgversion(parentmodule(@__MODULE__))) peak list")
        println(io, join(header, ","))
        for (p, l) in zip(pks, labels)
            row = [l; ""; _csvvalue.(collect(ppm(p))); _csvvalue.(collect(linewidth(p)));
                   _csvvalue(intensity(p))]
            println(io, join(row, ","))
        end
    end
end

_splitfields(line) = strip.(occursin(',', line) ? split(line, ',') : split(line))
_stripunit(name) = lowercase(strip(replace(name, r"\s*\(.*\)\s*$" => "")))
_parsevalue(s) = something(tryparse(Float64, s), NaN)

function _readpeakcsv(filename)
    lines = [strip(l) for l in eachline(filename)]
    lines = filter(l -> !isempty(l) && !startswith(l, '#'), lines)
    isempty(lines) && return Peak{0}[]
    names = _stripunit.(_splitfields(first(lines)))
    if "label" in names
        col = Dict(n => i for (i, n) in enumerate(names))
        rows = lines[2:end]
        N = something(findlast(n -> haskey(col, n), collect(NMRBase._peakcolumnnames(4))), 0)
        N == 0 && throw(ArgumentError("$filename has no x column"))
    else
        # a hand-made `label x y …` list without a header
        N = length(_splitfields(first(lines))) - 1
        N ≥ 1 || throw(ArgumentError("$filename has no peak positions"))
        dimnames = NMRBase._peakcolumnnames(N)
        col = Dict("label" => 1, (dimnames[n] => n + 1 for n in 1:N)...)
        rows = lines
    end
    dimnames = NMRBase._peakcolumnnames(N)
    value(fields, name) = haskey(col, name) && col[name] ≤ length(fields) ?
                          _parsevalue(fields[col[name]]) : NaN

    pks = Peak{N}[]
    seen = Set{String}()
    tracked = false
    for (k, line) in enumerate(rows)
        fields = _splitfields(line)
        l = String(fields[col["label"]])
        if l in seen && !isempty(l)
            tracked = true
            continue
        end
        position = ntuple(n -> value(fields, dimnames[n]), N)
        if any(isnan, position)
            @warn "Skipping a line of $filename without a valid position: $line"
            continue
        end
        push!(seen, l)
        push!(pks,
              Peak(position; label=l, intensity=value(fields, "intensity"),
                   linewidth=ntuple(n -> value(fields, "lw" * dimnames[n]), N)))
    end
    tracked && @warn "$filename tracks peaks plane by plane; using the first position of each"
    return pks
end

# ---- Sparky -------------------------------------------------------------------------
# Sparky writes a whitespace-delimited list with an `Assignment w1 w2 …` header, often
# followed by further columns such as `Data Height`. Unassigned peaks are labelled `?-?`.

function _issparky(filename)
    for line in eachline(filename)
        l = strip(line)
        (isempty(l) || startswith(l, '#')) && continue
        return !occursin(',', l) && lowercase(first(split(l))) == "assignment"
    end
    return false
end

function _readsparky(filename, spec)
    header = nothing
    entries = Tuple{String,Vector{Float64},Float64}[]
    for line in eachline(filename)
        l = strip(line)
        (isempty(l) || startswith(l, '#')) && continue
        fields = split(l)
        if isnothing(header)
            header = lowercase.(fields)
            continue
        end
        N = count(f -> occursin(r"^w\d+$", f), header)
        length(fields) ≥ N + 1 || continue
        w = tryparse.(Float64, fields[2:(N + 1)])
        any(isnothing, w) && continue
        height = "height" in header && length(fields) ≥ N + 2 ?
                 _parsevalue(fields[N + 2]) : NaN
        lab = String(fields[1])
        push!(entries, (occursin(r"^[?-]+$", lab) ? "" : lab, Float64.(w), height))
    end
    isempty(entries) && throw(ArgumentError("$filename contains no Sparky peaks"))

    reversed = _sparkyreversed(entries, spec)
    return [Peak(Tuple(reversed ? reverse(w) : w); label=l, intensity=h)
            for (l, w, h) in entries]
end

# Sparky's w1 is conventionally the indirect dimension, so its order is reversed relative to
# NMRTools. Given a spectrum, check that against how many peaks fall within the spectral
# ranges each way round; a tie keeps the convention.
function _sparkyreversed(entries, spec)
    isnothing(spec) && return true
    ranges = [extrema(data(d)) for d in dims(spec)]
    inside(w) = length(w) == length(ranges) &&
                all(r[1] ≤ v ≤ r[2] for (v, r) in zip(w, ranges))
    nrev = count(e -> inside(reverse(e[2])), entries)
    nfwd = count(e -> inside(e[2]), entries)
    nfwd > nrev && @info "Sparky list appears to have w1 as the direct dimension - reading it that way"
    return nfwd ≤ nrev
end

function _writesparky(filename, pks::AbstractVector{<:Peak{N}}, labels) where {N}
    hasheight = any(p -> !isnan(intensity(p)), pks)
    open(filename, "w") do io
        print(io, lpad("Assignment", 17))
        for n in 1:N
            print(io, lpad("w$n", 11))
        end
        hasheight && print(io, lpad("Data Height", 14))
        println(io, "\n")
        for (p, l) in zip(pks, labels)
            print(io, lpad(isempty(label(p)) ? "?-?" : l, 17))
            for δ in reverse(ppm(p))
                print(io, lpad(string(round(δ; digits=3)), 11))
            end
            hasheight && print(io, lpad(string(round(intensity(p); sigdigits=6)), 14))
            println(io)
        end
    end
end
