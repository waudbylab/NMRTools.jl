# Example data come from NMRExampleData.jl, which downloads each dataset the first time it is
# used. The BMRB shifts of ubiquitin are small enough to ship with NMRTools itself.

const BMRB_UBIQUITIN = joinpath(@__DIR__, "..", "..", "data", "bmr17769_3.str")
const BMRB_UBIQUITIN_TITLE = "Ubiquitin assigned chemical shifts (BMRB 17769)"

isbmrbubiquitin(name) = lowercase(name) == "bmrb_ubiquitin"

"""
    exampledata(name; fid=false, path=false, kwargs...)

Load an example dataset from [NMRExampleData.jl](https://waudbylab.org/NMRExampleData.jl/stable/),
downloading it the first time it is used.

Names ignore case, and can be the full name, an alias (such as the names used by earlier
versions of NMRTools, like `"2D_HN"`), or the start of the name or of any part of it. Where
a shortened name matches several datasets, the one added first is used. Call `exampledata()`
to choose from the available datasets, or see the
[catalogue](https://waudbylab.org/NMRExampleData.jl/stable/catalogue/).

The first time a dataset is used in a session, who to credit, what to cite and the licence
are shown. [`exampleinfo`](@ref) gives them at any time.

`BMRB_ubiquitin` returns the assigned chemical shifts of ubiquitin (BMRB entry 17769) as a
`Vector{ChemicalShift}`, which can be compared with the ubiquitin example spectra.

# Arguments
- `name`: name of the example dataset
- `fid`: if `true`, load the raw time-domain data (`fid` or `ser`) with
  [`loadfid`](@ref NMRTools.NMRIO.loadfid) instead of the processed spectrum
- `path`: if `true`, return the path of the dataset instead of loading it
- other keywords are passed to `loadfid` or [`loadnmr`](@ref NMRTools.NMRIO.loadnmr), for
  example `firstpoint` or `allcomponents`

A dataset with several experiments, such as a titration or a relaxation series, returns a
vector, in the order given by the dataset.

# Examples
```julia
spec = exampledata("2D_HN")                 # processed spectrum
fid = exampledata("2D_HN"; fid=true)        # raw data
fid = exampledata("2D_HN"; fid=true, firstpoint=Dict(T2Dim => 1.0))
dir = exampledata("2D_HN"; path=true)       # experiment directory
shifts = exampledata("BMRB_ubiquitin")      # assigned chemical shifts
spec = exampledata()                        # choose from a menu
```
"""
function exampledata(name::AbstractString; fid::Bool=false, path::Bool=false, kwargs...)
    if isbmrbubiquitin(name)
        fid && _nofid(name)
        return path ? BMRB_UBIQUITIN : loadbmrb(BMRB_UBIQUITIN; kwargs...)
    end
    info = exampleinfo(name)
    p = examplepath(info.name; fid)
    _showcredit(info)
    path && return p
    loadone(q) = fid ? loadfid(_examplefidpath(info.name, q); kwargs...) :
                 _loadexample(q; kwargs...)
    return p isa AbstractVector ? map(loadone, p) : loadone(p)
end

"""
    exampledata()

Choose an example dataset from a menu and load it. Outside an interactive session, list the
available datasets instead.
"""
function exampledata()
    names = [examples(); "BMRB_ubiquitin"]
    titles = [[exampleinfo(n).title for n in examples()]; BMRB_UBIQUITIN_TITLE]
    width = maximum(length, names)
    entries = [rpad(n, width) * "  " * t for (n, t) in zip(names, titles)]
    if isinteractive() && stdin isa Base.TTY
        menu = RadioMenu(entries; pagesize=min(20, length(entries)))
        i = request("Choose an example dataset:", menu)
        return i == -1 ? nothing : exampledata(names[i])
    end
    println("Example datasets, to load with exampledata(name):\n")
    foreach(e -> println("  ", e), entries)
    return nothing
end

# processed data in NMR-STAR format are BMRB shifts; anything else is a spectrum
_loadexample(path; kwargs...) = endswith(path, ".str") ? loadbmrb(path; kwargs...) :
                                loadnmr(path; kwargs...)

# experiment directory holding the raw data of an example data set
function _examplefidpath(name, p)
    dir = isdir(p) ? p : dirname(p)
    while !isfile(joinpath(dir, "acqus")) && dirname(dir) != dir
        dir = dirname(dir)
    end
    isfile(joinpath(dir, "fid")) || isfile(joinpath(dir, "ser")) || _nofid(name)
    return dir
end

_nofid(name) = throw(ArgumentError("example data $name has no raw time-domain data"))

# datasets whose credit has been shown in this session
const CREDITED = Set{String}()

function _showcredit(info)
    info.name in CREDITED && return
    push!(CREDITED, info.name)
    lines = ["Example data $(info.name): $(info.title)"]
    isempty(info.credit) || push!(lines, "Credit: " * join(info.credit, ", "))
    isempty(info.citation) || push!(lines, "Cite: " * info.citation)
    if isempty(info.licence) || info.licence == "None"
        push!(lines, "Licence: none granted. Ask those credited before using these data " *
                     "beyond trying out the software.")
    else
        push!(lines, "Licence: " * info.licence)
    end
    @info join(lines, "\n")
    return
end
