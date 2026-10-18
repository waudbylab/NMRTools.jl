# Assigned chemical shifts from NMR-STAR files and the BMRB.
#
# Only the `_Atom_chem_shift` loops are read. An NMR-STAR file is a sequence of tokens:
# bare words, words quoted with ' or " (a quote closes only when followed by whitespace),
# and multi-line text fields between lines beginning with a semicolon.

const BMRB_URL = "https://bmrb.io/ftp/pub/bmrb/entry_directories/bmr%ID%/bmr%ID%_3.str"

"""
    loadbmrb(id::Integer; list=nothing) -> Vector{ChemicalShift}
    loadbmrb(filename; list=nothing) -> Vector{ChemicalShift}
    loadbmrb(id_or_filename, dims...; list=nothing) -> Vector{Peak{N}}

Download the assigned chemical shifts of BMRB entry `id`, or read them from a local
NMR-STAR (v3) file.

Given dimensions as well, the shifts are assembled straight into the peak list of an
experiment, exactly as [`peaklist`](@ref NMRTools.NMRBase.peaklist) does: `H1, N15` for a
¹⁵N HSQC, `H1, C13` for a ¹³C HSQC, or atom names such as `"H", "N", "CO-1"` for an HNCO.

An entry may hold several shift lists, e.g. under different sample conditions. The first is
used unless `list` selects another by its ID. Shifts of nuclei NMRTools does not know
about (such as ¹⁷O) are skipped.

# Examples
```julia
loadbmrb(15000, H1, N15)                 # ¹⁵N HSQC peak list, side chains included
loadbmrb(15000, H1, C13)                 # ¹³C HSQC peak list
loadbmrb(15000, "H", "N")                # backbone amides only
loadbmrb(15000, "CO-1", "N")             # CON

shifts = loadbmrb(15000)                 # the shifts themselves, to filter or reuse
shifts = loadbmrb("bmr15000_3.str")      # a local NMR-STAR file
```
"""
loadbmrb(filename::AbstractString; kwargs...) = _readshifts(read(filename, String); kwargs...)

function loadbmrb(id::Integer; kwargs...)
    url = replace(BMRB_URL, "%ID%" => string(id))
    @info "Downloading BMRB entry $id"
    text = String(take!(Downloads.download(url, IOBuffer())))
    return _readshifts(text; kwargs...)
end

function loadbmrb(source::Union{Integer,AbstractString}, dim::Union{AbstractString,Nucleus},
                  dims::Union{AbstractString,Nucleus}...; kwargs...)
    return peaklist(loadbmrb(source; kwargs...), dim, dims...)
end

function _readshifts(text; list=nothing)
    tags, rows = _starloop(text, "_Atom_chem_shift.")
    isempty(tags) && throw(NMRToolsError("no assigned chemical shifts found"))
    function col(name)
        i = findfirst(==("_Atom_chem_shift.$name"), tags)
        isnothing(i) && throw(NMRToolsError("shift list has no $name column"))
        return i
    end
    ilist = findfirst(==("_Atom_chem_shift.Assigned_chem_shift_list_ID"), tags)
    ientity = findfirst(==("_Atom_chem_shift.Entity_assembly_ID"), tags)
    iseq, icomp, iatom = col("Seq_ID"), col("Comp_ID"), col("Atom_ID")
    itype, iiso, ival = col("Atom_type"), col("Atom_isotope_number"), col("Val")
    ierr = findfirst(==("_Atom_chem_shift.Val_err"), tags)

    lists = isnothing(ilist) ? ["1"] : unique(r[ilist] for r in rows)
    chosen = isnothing(list) ? first(lists) : string(list)
    chosen in lists ||
        throw(ArgumentError("no shift list $chosen; this entry has lists $(join(lists, ", "))"))
    if isnothing(list) && length(lists) > 1
        @warn "Entry has $(length(lists)) shift lists ($(join(lists, ", "))); using list " *
              "$chosen. Select another with the `list` keyword."
    end

    number(s) = something(tryparse(Float64, s), NaN)
    shifts = ChemicalShift[]
    for r in rows
        isnothing(ilist) || r[ilist] == chosen || continue
        δ = tryparse(Float64, r[ival])
        residue = tryparse(Int, r[iseq])
        nuc = get(NMRBase.NUCLEUS_LOOKUP, r[iiso] * uppercase(r[itype]), nothing)
        if any(isnothing, (δ, residue, nuc))
            @debug "Skipping shift" row = r
            continue
        end
        entity = isnothing(ientity) ? 1 : something(tryparse(Int, r[ientity]), 1)
        push!(shifts,
              ChemicalShift(entity, residue, r[icomp], r[iatom], nuc, δ,
                            isnothing(ierr) ? NaN : number(r[ierr])))
    end
    @debug "Read $(length(shifts)) chemical shifts"
    return shifts
end

# The tags and rows of every loop whose tags start with `prefix`, concatenated.
function _starloop(text, prefix)
    tokens = _startokens(text)
    tags = String[]
    rows = Vector{String}[]
    i = 1
    while i ≤ length(tokens)
        if tokens[i] == "loop_" && i < length(tokens) && startswith(tokens[i + 1], prefix)
            looptags = String[]
            i += 1
            while i ≤ length(tokens) && startswith(tokens[i], "_")
                push!(looptags, tokens[i])
                i += 1
            end
            values = String[]
            while i ≤ length(tokens) && tokens[i] != "stop_"
                push!(values, tokens[i])
                i += 1
            end
            n = length(looptags)
            length(values) % n == 0 ||
                throw(NMRToolsError("malformed NMR-STAR loop: $(length(values)) values " *
                                    "for $n tags"))
            isempty(tags) && append!(tags, looptags)
            if looptags == tags
                append!(rows, [values[(k + 1):(k + n)] for k in 0:n:(length(values) - 1)])
            else
                @debug "Skipping an NMR-STAR loop with different tags"
            end
        end
        i += 1
    end
    return tags, rows
end

function _startokens(text)
    tokens = String[]
    lines = split(replace(text, "\r\n" => "\n"), '\n')
    k = 1
    while k ≤ length(lines)
        line = lines[k]
        if startswith(line, ';')
            # multi-line text field, closed by a line that is just ';'
            field = [line[2:end]]
            k += 1
            while k ≤ length(lines) && !startswith(lines[k], ';')
                push!(field, lines[k])
                k += 1
            end
            push!(tokens, join(field, "\n"))
        else
            for m in eachmatch(r"'(?:[^']|'(?=\S))*'(?=\s|$)|\"(?:[^\"]|\"(?=\S))*\"(?=\s|$)|#.*$|\S+", line)
                t = m.match
                startswith(t, '#') && break
                if length(t) ≥ 2 && first(t) in ('\'', '"') && last(t) == first(t)
                    t = t[2:(end - 1)]
                end
                push!(tokens, String(t))
            end
        end
        k += 1
    end
    return tokens
end
