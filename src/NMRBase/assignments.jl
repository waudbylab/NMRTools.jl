"""
    ChemicalShift

One assigned chemical shift: an atom of a residue, its nucleus, and its shift.

| field         | meaning                                                    |
|:--------------|:-----------------------------------------------------------|
| `entity`      | which molecule in the assembly (BMRB `Entity_assembly_ID`) |
| `residue`     | residue number (BMRB `Seq_ID`)                             |
| `residuename` | residue type, e.g. `"LEU"`                                 |
| `atom`        | IUPAC atom name, e.g. `"HD11"`                             |
| `nucleus`     | e.g. `H1`, `C13`                                           |
| `shift`       | chemical shift (ppm)                                       |
| `error`       | uncertainty (ppm), `NaN` if not reported                   |

A shift list is a `Vector{ChemicalShift}`, read with
[`loadshifts`](@ref NMRTools.NMRIO.loadshifts) and turned into the peaks expected in a given
experiment with [`peaklist`](@ref). Use `filter` to select residues before building peaks,
e.g. `filter(s -> s.residuename in ("ILE", "LEU", "VAL"), shifts)`.
"""
struct ChemicalShift
    entity::Int
    residue::Int
    residuename::String
    atom::String
    nucleus::Nucleus
    shift::Float64
    error::Float64
end

ppm(s::ChemicalShift) = s.shift
nucleus(s::ChemicalShift) = s.nucleus
label(s::ChemicalShift) = _residuecode(s.residuename) * string(s.residue) * s.atom

Base.show(io::IO, s::ChemicalShift) = print(io, "ChemicalShift(", label(s), ", ", s.shift, ")")

"""
    peaklist(shifts, atoms...) -> Vector{Peak{N}}
    peaklist(shifts, H1, X) -> Vector{Peak{2}}

The peaks expected in an `N`-dimensional experiment, assembled from a list of assigned
chemical shifts (see [`loadshifts`](@ref NMRTools.NMRIO.loadshifts)). Dimensions are
given in the order of the spectrum, direct dimension first.

**By atom name.** Each argument names an atom, optionally with a residue offset relative
to the residue that the peak belongs to: `"CA-1"` is the Cα of the preceding residue. A
peak is produced for each residue where every atom has an assigned shift. `"CO"` and `"HN"`
are accepted for the BMRB names `"C"` and `"H"`.

**By nucleus.** Passing `H1` and another nucleus gives one peak for every directly bonded
pair of atoms with assigned shifts: every N–H for an ¹⁵N HSQC, including side-chain NH and
NH₂ groups, or every C–H for a ¹³C HSQC. Methyl and methylene protons with identical
shifts give a single peak.

Peaks are labelled by residue and atom, e.g. `L23H-N`, `K22C-L23N` or `L23HD1-CD1`.
Positions come from the shift list, while intensities and linewidths are `NaN`.

# Examples
```julia
shifts = loadshifts(15000)              # download BMRB entry 15000

peaklist(shifts, "H", "N")              # backbone amide ¹H,¹⁵N HSQC
peaklist(shifts, H1, N15)               # all NH groups, including side chains
peaklist(shifts, H1, C13)               # ¹H,¹³C HSQC
peaklist(shifts, "C-1", "N")            # CON: C′ of residue i-1 with N of residue i
peaklist(shifts, "H", "N", "C-1")       # HNCO
peaklist(shifts, "H", "N", "CA")        # the i peaks of an HNCA
peaklist(shifts, "H", "N", "CA-1")      # and the i-1 peaks

ilv = filter(s -> s.residuename in ("ILE", "LEU", "VAL"), shifts)
peaklist(ilv, H1, C13)                  # ILV methyl (and other CH) peaks
```
"""
function peaklist(shifts::AbstractVector{ChemicalShift}, atoms::AbstractString...)
    isempty(atoms) && throw(ArgumentError("peaklist needs at least one atom"))
    specs = map(_atomspec, atoms)
    table = _shifttable(shifts)
    multientity = length(unique(first.(keys(table)))) > 1
    pks = Peak{length(specs)}[]
    for key in sort(collect(keys(table)))
        entity, residue = key
        found = map(specs) do (atom, offset)
            residueatoms = get(table, (entity, residue + offset), nothing)
            isnothing(residueatoms) ? nothing : get(residueatoms, atom, nothing)
        end
        any(isnothing, found) && continue
        push!(pks, Peak(Tuple(ppm.(found)); label=_peaklabel(collect(found), multientity)))
    end
    return pks
end

function peaklist(shifts::AbstractVector{ChemicalShift}, a::Nucleus, b::Nucleus)
    (a == H1) ⊻ (b == H1) ||
        throw(ArgumentError("peaklist by nucleus builds X–H correlations, so exactly one " *
                            "nucleus must be H1; use atom names for other experiments"))
    heavy = a == H1 ? b : a
    table = _shifttable(shifts)
    multientity = length(unique(first.(keys(table)))) > 1
    pks = Peak{2}[]
    for key in sort(collect(keys(table)))
        atoms = table[key]
        # group the protons by the heavy atom they are bonded to
        groups = Dict{String,Vector{ChemicalShift}}()
        for h in values(atoms)
            nucleus(h) == H1 || continue
            x = _bondedatom(h.atom, atoms)
            (isnothing(x) || nucleus(atoms[x]) != heavy) && continue
            push!(get!(groups, x, ChemicalShift[]), h)
        end
        for x in sort(collect(keys(groups)))
            X = atoms[x]
            protons = sort(groups[x]; by=h -> h.atom)
            for hs in _degenerate(protons)
                h = length(hs) == 1 ? only(hs) : _pseudoatom(hs)
                pair = a == H1 ? [h, X] : [X, h]
                push!(pks, Peak(ppm(pair[1]), ppm(pair[2]); label=_peaklabel(pair, multientity)))
            end
        end
    end
    return pks
end

# "CA-1" → ("CA", -1); BMRB names carbonyl carbons C and amide protons H
function _atomspec(s::AbstractString)
    m = match(r"^([A-Za-z0-9'\"]+?)([+-]\d+)?$", strip(s))
    isnothing(m) && throw(ArgumentError("cannot read atom $(repr(s)); use e.g. \"CA\" or \"CA-1\""))
    atom = uppercase(m.captures[1])
    atom = get(Dict("CO" => "C", "HN" => "H"), atom, atom)
    offset = isnothing(m.captures[2]) ? 0 : parse(Int, m.captures[2])
    return (atom, offset)
end

# (entity, residue) => atom name => shift
function _shifttable(shifts)
    table = Dict{Tuple{Int,Int},Dict{String,ChemicalShift}}()
    for s in shifts
        get!(table, (s.entity, s.residue), Dict{String,ChemicalShift}())[s.atom] = s
    end
    return table
end

# The heavy atom a proton is bonded to, following IUPAC naming: the proton H<x> is bonded to
# the heavy atom named <element><x> (HE1 → NE1 in Trp, HD1 → CD1 in Phe, H1′ → C1′), or,
# for methyl and methylene protons, <element><x> less its final digit or prime (HB2 → CB,
# HD11 → CD1, H5″ → C5′). The backbone amide proton H is bonded to N.
function _bondedatom(h::AbstractString, atoms)
    h == "H" && return haskey(atoms, "N") ? "N" : nothing
    startswith(h, "H") || return nothing
    suffix = h[nextind(h, 1):end]
    candidates = [suffix]
    length(suffix) > 1 && (isdigit(last(suffix)) || last(suffix) == '\'') &&
        push!(candidates, chop(suffix))
    for c in candidates, element in ("C", "N", "O", "S", "P")
        name = element * c
        haskey(atoms, name) && return name
    end
    return nothing
end

# Split protons on one heavy atom into groups with identical shifts.
function _degenerate(protons)
    groups = Vector{ChemicalShift}[]
    for h in protons
        i = findfirst(g -> ppm(first(g)) == ppm(h), groups)
        isnothing(i) ? push!(groups, [h]) : push!(groups[i], h)
    end
    return groups
end

# A methyl's three protons, HD11/HD12/HD13, are named together as HD1.
function _pseudoatom(hs)
    h = first(hs)
    return ChemicalShift(h.entity, h.residue, h.residuename, chop(h.atom), h.nucleus,
                         h.shift, h.error)
end

# Residue and atom names, with the residue repeated only where it changes: L23H-N,
# K22C-L23N, L23H-N-K22C. Multi-entity assemblies are prefixed with the entity, 2:L23H-N.
function _peaklabel(found, multientity)
    parts = String[]
    previous = nothing
    for s in found
        r = (s.entity, s.residue)
        push!(parts, r == previous ? s.atom : label(s))
        previous = r
    end
    l = join(parts, "-")
    return multientity ? "$(first(found).entity):$l" : l
end

const RESIDUE_CODES = Dict("ALA" => "A", "ARG" => "R", "ASN" => "N", "ASP" => "D",
                           "CYS" => "C", "GLN" => "Q", "GLU" => "E", "GLY" => "G",
                           "HIS" => "H", "ILE" => "I", "LEU" => "L", "LYS" => "K",
                           "MET" => "M", "PHE" => "F", "PRO" => "P", "SER" => "S",
                           "THR" => "T", "TRP" => "W", "TYR" => "Y", "VAL" => "V")

_residuecode(name) = get(RESIDUE_CODES, uppercase(name), name)
