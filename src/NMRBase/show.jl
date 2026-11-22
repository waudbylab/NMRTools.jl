# Terminal display of NMRData ######################################################################
#
# DimensionalData prints an array as a box: a summary line, a `dims` block, a `metadata` block, and
# then the data. For NMRData the metadata dictionary holds the whole of acqus, the channel tables
# and so on, so printing it verbatim pushes the data off the screen. Here we replace only the
# metadata block, by overriding `DimensionalData.show_main` (the documented hook for this), and keep
# DimensionalData's dims block and data matrix.

const METADATADISPLAYMODES = (:none, :summary, :keys, :full)

const METADATADISPLAY = Ref(:summary)

"""
    metadatadisplay()

Return the mode that sets how the metadata of an [`NMRData`](@ref) are shown at the REPL:
`:summary` (the default), `:keys`, `:none` or `:full`. See [`metadatadisplay!`](@ref).
"""
metadatadisplay() = METADATADISPLAY[]

"""
    metadatadisplay!(mode)

Set how the metadata of an [`NMRData`](@ref) are shown at the REPL, and return the mode. The
`mode` is one of

- `:summary` (default): a few key entries (label, pulse programme, solvent, temperature, field,
  scans, noise, date and file), and a count of the remainder.
- `:keys`: the names of all entries, without their values.
- `:none`: no metadata block, only the dimensions and the data.
- `:full`: every entry and value, as DimensionalData.jl prints them. Entries such as `:acqus` are
  long, so this can fill the screen.

The setting applies to the whole session. To change it for one display only, set the
`:nmrmetadata` key of an `IOContext`:

```julia
show(IOContext(stdout, :nmrmetadata => :full), MIME"text/plain"(), spec)
```

All the metadata remain available through [`metadata`](@ref) whatever the display mode.

See also [`metadatadisplay`](@ref).
"""
function metadatadisplay!(mode::Symbol)
    mode in METADATADISPLAYMODES ||
        throw(ArgumentError("unknown metadata display mode $(repr(mode)): use one of $(join(repr.(METADATADISPLAYMODES), ", "))"))
    METADATADISPLAY[] = mode
    return mode
end

# The metadata block, formatted as rows of `key => value` text. DimensionalData lays out this block
# by calling `isempty`, and `show` with the MIME type, on whatever it is given, so this small type
# stands in for the dictionary. Rows are formatted to the display width when it is constructed,
# because DimensionalData measures the width without the IOContext that the final `show` receives.
struct MetadataSummary
    rows::Vector{Pair{String,String}}
end

Base.isempty(s::MetadataSummary) = isempty(s.rows)

function Base.show(io::IO, ::MIME"text/plain", s::MetadataSummary)
    keywidth = maximum(textwidth ∘ first, s.rows)
    for (i, (key, value)) in enumerate(s.rows)
        # DimensionalData indents the first line itself, and ends the last
        i > 1 && print(io, "\n  ")
        printstyled(io, rpad(key, keywidth + 2); color=:light_black)
        print(io, value)
    end
end

# Keys shown in the summary, in order, with the units appended to their values
const SUMMARYKEYS = (:label, :pulseprogram, :solvent, :temperature, :field, :ns, :noise,
                     :date,
                     :filename)
const SUMMARYUNITS = Dict(:temperature => " K", :field => " T")

_summaryvalue(x) = sprint(show, x; context=:compact => true)
_summaryvalue(x::AbstractString) = first(eachsplit(x, '\n')) # title files have many lines
_summaryvalue(x::AbstractFloat) = string(round(x; sigdigits=5))
_summaryvalue(x::Dates.DateTime) = Dates.format(x, "yyyy-mm-dd HH:MM")

# Cut text to fit a width. Paths are most informative at the end, other text at the start.
function _truncate(s::AbstractString, width::Integer, keep::Symbol)
    textwidth(s) <= width && return s
    width <= 1 && return "…"
    chars = collect(s)
    return keep === :end ? string("…", String(chars[(end - width + 2):end])) :
           string(String(chars[1:(width - 1)]), "…")
end

# Entries of the form `key => value`, one per row. `width` is the room available for the value.
function _summaryrow(md, key::Symbol, width::Integer)
    value = string(_summaryvalue(md[key]), get(SUMMARYUNITS, key, ""))
    return string(key) => _truncate(value, width, key === :filename ? :end : :start)
end

# Entries left at their defaults (`nothing` or empty text) say nothing, so are not shown
_isset(x) = true
_isset(::Nothing) = false
_isset(::Missing) = false
_isset(x::AbstractString) = !isempty(x)

# `width` is the text width available inside the metadata block, including the key column
function _metadatarows(::Val{:summary}, md, width::Integer)
    present = filter(key -> haskey(md, key) && _isset(md[key]), collect(SUMMARYKEYS))
    isempty(present) && return Pair{String,String}[]
    valuewidth = width - maximum(length ∘ string, present) - 2
    rows = [_summaryrow(md, key, valuewidth) for key in present]
    nother = length(md) - length(present)
    nother > 0 && push!(rows, "…" => string(nother, " more entries: see metadata(spec)"))
    return rows
end

function _metadatarows(::Val{:keys}, md, width::Integer)
    isempty(md) && return Pair{String,String}[]
    names = sort!(map(string, collect(keys(md))))
    label = string(length(names), " entries")
    available = width - textwidth(label) - 2
    lines = String[]
    line = ""
    for name in names
        candidate = isempty(line) ? name : string(line, ", ", name)
        if !isempty(line) && textwidth(candidate) > available
            push!(lines, string(line, ","))
            line = name
        else
            line = candidate
        end
    end
    push!(lines, line)
    return [(i == 1 ? label : "") => l for (i, l) in enumerate(lines)]
end

_metadatarows(::Val{:none}, md, width::Integer) = Pair{String,String}[]

function _metadatarows(::Val{M}, md, width::Integer) where {M}
    return throw(ArgumentError("unknown metadata display mode $(repr(M)): use one of $(join(repr.(METADATADISPLAYMODES), ", "))"))
end

function DD.show_main(io::IO, mime::MIME"text/plain", A::AbstractNMRData;
                      blockwidth=0, displaywidth=displaysize(io)[2])
    mode = get(io, :nmrmetadata, metadatadisplay())
    mode === :full &&
        return invoke(DD.show_main, Tuple{Any,Any,DD.AbstractBasicDimArray}, io, mime, A;
                      blockwidth, displaywidth)

    # leave room for the box edge and indentation either side
    summary = MetadataSummary(_metadatarows(Val(mode), metadata(A), displaywidth - 6))

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
