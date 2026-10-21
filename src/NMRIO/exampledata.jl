const EXAMPLEDATA_NAMES = ["1D_1H", "1D_19F", "1D_19F_titration", "2D_HN", "2D_HN_titration",
                           "pseudo2D_T2", "pseudo2D_XSTE", "pseudo3D_HN_R2", "pseudo3D_kinetics",
                           "3D_HNCA", "3D_HNCA_pipe", "3D_HNCA_ucsf", "BMRB_ubiquitin"]

function exampledata()
    return println("""
       Usage: exampledata("1D_19F_titration")

       Available example data:
       $(join("- " .* EXAMPLEDATA_NAMES, "\n"))
       """)
end

"""
    exampledata(name; fid=false, path=false, kwargs...)

Load example NMR data provided with NMRTools.jl.

Call `exampledata()` to list the available data sets. The ubiquitin data sets (`1D_1H`,
`2D_HN`, `3D_HNCA`, `3D_HNCA_pipe` and `3D_HNCA_ucsf`) can be compared with the assigned
chemical shifts of ubiquitin in `BMRB_ubiquitin` (BMRB entry 17769), which loads as a
`Vector{ChemicalShift}`.

# Arguments
- `name`: name of the example data set
- `fid`: if `true`, load the raw time-domain data (`fid` or `ser`) with
  [`loadfid`](@ref NMRTools.NMRIO.loadfid) instead of the processed spectrum
- `path`: if `true`, return the path of the data set instead of loading it
- other keywords are passed to `loadfid` or [`loadnmr`](@ref NMRTools.NMRIO.loadnmr), for
  example `firstpoint` or `allcomponents`

Series (`1D_19F_titration`, `2D_HN_titration`) return a vector.

# Examples
```julia
spec = exampledata("2D_HN")                 # processed spectrum
fid = exampledata("2D_HN"; fid=true)        # raw data
fid = exampledata("2D_HN"; fid=true, firstpoint=Dict(T2Dim => 1.0))
dir = exampledata("2D_HN"; path=true)       # experiment directory
shifts = exampledata("BMRB_ubiquitin")      # assigned chemical shifts
```
"""
function exampledata(name; fid=false, path=false, kwargs...)
    p = _examplepath(name)
    path && return p
    name == "BMRB_ubiquitin" && return fid ? _nofid(name) : loadbmrb(p; kwargs...)
    load = fid ? (q -> loadfid(_examplefidpath(name, q); kwargs...)) : (q -> loadnmr(q; kwargs...))
    return p isa AbstractVector ? map(load, p) : load(p)
end

function _examplepath(name)
    if name == "1D_1H"
        artifact"1D_1H"
    elseif name == "1D_19F"
        artifact"1D_19F"
    elseif name == "1D_19F_titration"
        [joinpath(artifact"1D_19F_titration", "$i") for i in 1:11]
    elseif name == "2D_HN"
        artifact"2D_HN"
    elseif name == "2D_HN_titration"
        [joinpath(artifact"2D_HN_titration", "$i", "test.ft2") for i in 1:11]
    elseif name == "pseudo2D_T2"
        artifact"pseudo2D_T2"
    elseif name == "pseudo2D_XSTE"
        artifact"pseudo2D_XSTE"
    elseif name == "pseudo3D_HN_R2"
        artifact"pseudo3D_HN_R2"
    elseif name == "pseudo3D_kinetics"
        artifact"pseudo3D_kinetics"
    elseif name == "3D_HNCA"
        artifact"3D_HNCA"
    elseif name == "3D_HNCA_pipe"
        joinpath(artifact"3D_HNCA_pipe", "smile", "test%03d.ft3")
    elseif name == "3D_HNCA_ucsf"
        joinpath(artifact"3D_HNCA_ucsf", "hnca.ucsf")
    elseif name == "BMRB_ubiquitin"
        joinpath(@__DIR__, "..", "..", "data", "bmr17769_3.str")
    else
        throw(NMRToolsError("example data $name not found. Call exampledata() to list available data"))
    end
end

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
