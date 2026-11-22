# Terminal display of NMRData ######################################################################
#
# DimensionalData prints an array as a box: a summary line, a `dims` block, a `metadata` block, and
# then the data. For NMRData the metadata dictionary holds the whole of acqus, the channel tables
# and so on, so printing it verbatim pushes the data off the screen. Here we replace only the
# metadata block, by overriding `DimensionalData.show_main` (the documented hook for this), with
# two lines of key fields. DimensionalData's dims block and data matrix are kept.

# DimensionalData lays out the metadata block by calling `isempty`, and `show` with the MIME type,
# on whatever it is given, so this small type stands in for the dictionary.
struct MetadataSummary
    lines::Vector{String}
end

Base.isempty(s::MetadataSummary) = isempty(s.lines)

# DimensionalData indents the first line itself, and ends the last
Base.show(io::IO, ::MIME"text/plain", s::MetadataSummary) = join(io, s.lines, "\n  ")

# Entries left at their defaults (`nothing` or empty text) say nothing, so are not shown
_isset(x) = true
_isset(::Nothing) = false
_isset(::Missing) = false
_isset(x::AbstractString) = !isempty(x)

# Cut text to fit a width, marking the cut
function _shorten(s::AbstractString, width::Integer)
    textwidth(s) <= width && return String(s)
    return string(first(s, max(width - 1, 0)), "…")
end

_firstline(x::AbstractString) = first(eachsplit(x, '\n')) # title files have many lines
_firstline(x) = x

# `name=value` for a field that is set, and nothing otherwise
_field(name, value) = _isset(value) ? string(name, '=', value) : nothing
_field(name, value, unit) = _isset(value) ? string(name, '=', value, ' ', unit) : nothing
function _field(name, value::Real, unit)
    rounded = round(value; sigdigits=3)
    return string(name, '=', isinteger(rounded) ? Int(rounded) : rounded, ' ', unit)
end

# Key fields on at most two lines, each cut to `width`
function _summarylines(A::AbstractNMRData, width::Integer)
    label = metadata(A, :label)
    pulseprogram = metadata(A, :pulseprogram)
    date = metadata(A, :date)
    # Line 1 describes the sample, line 2 the acquisition. This also keeps the two lines of
    # similar length (about 75 and 85 characters at most). Sample and annotation status are
    # always shown.
    line1 = (_field("label", _isset(label) ? _shorten(_firstline(label), 20) : nothing),
             _field("sample", hassample(A)),
             _field("temperature", metadata(A, :temperature), "K"),
             _field("solvent", metadata(A, :solvent)))
    line2 = (_field("date",
                    date isa Dates.DateTime ? Dates.format(date, "yyyy-mm-dd") : date),
             _field("field", metadata(A, :field), "T"),
             _field("ns", metadata(A, :ns)),
             _field("pulseprogram",
                    _isset(pulseprogram) ? _shorten(pulseprogram, 20) : nothing),
             _field("annotations", hasannotations(A)))

    lines = [join(filter(!isnothing, line), "; ") for line in (line1, line2)]
    return [_shorten(line, width) for line in lines if !isempty(line)]
end

function DD.show_main(io::IO, mime::MIME"text/plain", A::AbstractNMRData;
                      blockwidth=0, displaywidth=displaysize(io)[2])
    # leave room for the box edge and indentation either side
    summary = MetadataSummary(_summarylines(A, displaywidth - 6))

    # as DimensionalData.show_main: measure the metadata block first, so that the box drawn
    # around the dimensions is wide enough for it
    _, blockwidth, _ = DD.print_metadata_block(IOBuffer(), mime, summary; blockwidth,
                                               displaywidth)
    lines_t, blockwidth, displaywidth, separatorwidth, istop = DD.print_top(io, mime, A;
                                                                            blockwidth,
                                                                            displaywidth)
    lines_m, blockwidth, stilltop = DD.print_metadata_block(io, mime, summary; blockwidth,
                                                            displaywidth, separatorwidth,
                                                            istop)
    return lines_t + lines_m, blockwidth, stilltop
end
