"""
    NMRSample(path::String) -> NMRSample

Load an NMR sample JSON file and return an `NMRSample` descriptor.
"""
function NMRSample(path::String)
    isfile(path) || throw(NMRToolsError("sample file not found: $path"))
    return NMRSample(abspath(path), loadsample(path))
end

"""
    addsampleinfo(expt::NMRData) -> NMRData

Search for and load sample metadata from JSON files in the experiment folder.
Returns the NMRData object with sample information added to metadata under `:sample`
as an [`NMRSample`](@ref).

Sample files are JSON files following the schema at https://github.com/waudbylab/nmr-sample-schema.
The function matches samples based on creation/ejection timestamps against the experiment date.

# Examples
```julia
spec = loadnmr("experiment/1/pdata/1")
spec = addsampleinfo(spec)  # Load matching sample metadata
sample(spec, "notes")  # Access sample information
```

See also [`sample`](@ref), [`hassample`](@ref).
"""
function addsampleinfo(expt::NMRData)
    experimentfolder = expt[:experimentfolder]
    isnothing(experimentfolder) && return expt
    sampledate = expt[:date]
    isnothing(sampledate) && return expt

    candidates = _scan_sample_files(dirname(experimentfolder))
    for (path, tcreation, tejection) in candidates
        if tcreation <= sampledate <= tejection
            expt[:sample] = NMRSample(path)
            return expt
        end
    end
    return expt
end

"""
    _addsampleinfo_metadata!(md::Dict{Symbol,Any})

Search for a matching sample JSON in the experiment parent directory and store
the result as an `NMRSample` under `md[:sample]`. Operates on a raw metadata dict.
"""
function _addsampleinfo_metadata!(md::Dict{Symbol,Any})
    experimentfolder = get(md, :experimentfolder, nothing)
    isnothing(experimentfolder) && return md

    sampledate = get(md, :date, nothing)
    isnothing(sampledate) && return md

    candidates = _scan_sample_files(dirname(experimentfolder))
    for (path, tcreation, tejection) in candidates
        if tcreation <= sampledate <= tejection
            md[:sample] = NMRSample(path)
            return md
        end
    end
    return md
end

"""
    scansamples(dir; recursive=false) -> Vector{NMRSample}

Scan `dir` for NMR sample JSON files and return a vector of [`NMRSample`](@ref)
descriptors sorted by creation timestamp. Files that cannot be parsed or have an
unsupported schema version are skipped silently.

By default only the immediate directory is searched (not subdirectories), since
sample files are typically stored flat in a single directory.
"""
function scansamples(dir::AbstractString; recursive::Bool=false)::Vector{NMRSample}
    dir = abspath(dir)
    isdir(dir) || throw(NMRToolsError("not a directory: $dir"))

    candidates = _scan_sample_files(dir; recursive)
    samples = NMRSample[]
    for (path, _, _) in candidates
        try
            push!(samples, NMRSample(path))
        catch e
            @debug "Skipping sample file" path exception = e
        end
    end
    return samples
end

"""
    findsample(expt::NMRExperiment) -> Union{NMRSample, Nothing}
    findsample(expt::NMRExperiment, dir::String) -> Union{NMRSample, Nothing}
    findsample(expt::NMRExperiment, candidates::Vector{NMRSample}) -> Union{NMRSample, Nothing}

Find the sample that matches `expt` based on experiment date falling within the
sample's creation–ejection timestamp window.

- With no second argument, searches the parent directory of the experiment folder.
- With a directory path, searches that directory for JSON sample files.
- With a pre-scanned `Vector{NMRSample}`, filters without any disk I/O.
"""
function findsample(expt::NMRExperiment)::Union{NMRSample,Nothing}
    experimentfolder = expt[:experimentfolder]
    isnothing(experimentfolder) && return nothing
    return findsample(expt, dirname(experimentfolder))
end

function findsample(expt::NMRExperiment, dir::AbstractString)::Union{NMRSample,Nothing}
    return findsample(expt, scansamples(dir))
end

function findsample(expt::NMRExperiment,
                    candidates::Vector{NMRSample})::Union{NMRSample,Nothing}
    sampledate = expt[:date]
    isnothing(sampledate) && return nothing
    for s in candidates
        tcreation, tejection = _sample_timestamps(s)
        isnothing(tcreation) && continue
        tcreation <= sampledate <= tejection && return s
    end
    return nothing
end

# Internal helpers #################################################################################

function _scan_sample_files(dir::AbstractString; recursive::Bool=false)
    # Returns list of (path, tcreation, tejection) for valid sample JSON files, sorted by tcreation
    results = Tuple{String,DateTime,DateTime}[]
    _collect_sample_files!(results, dir, recursive)
    sort!(results; by=x -> x[2])
    return results
end

function _collect_sample_files!(results, dir, recursive)
    for entry in readdir(dir; join=true)
        if isdir(entry) && recursive
            _collect_sample_files!(results, entry, recursive)
        elseif isfile(entry) && endswith(entry, ".json")
            try
                j = loadsample(entry)
                m = get(j, "metadata", nothing)
                isnothing(m) && continue
                v = VersionNumber(get(m, "schema_version", "0"))
                v >= v"0.0.3" || continue

                creation = get(m, "created_timestamp", nothing)
                isnothing(creation) && continue
                tcreation = DateTime(creation, dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")

                ejection = get(m, "ejected_timestamp", nothing)
                tejection = isnothing(ejection) ?
                            DateTime("3000-01-01T00:00:00.000Z",
                                     dateformat"yyyy-mm-ddTHH:MM:SS.sssZ") :
                            DateTime(ejection, dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")

                push!(results, (entry, tcreation, tejection))
            catch e
                @debug "Skipping sample entry due to parse error" entry exception = e
            end
        end
    end
end

function _sample_timestamps(s::NMRSample)
    m = get(s.metadata, "metadata", nothing)
    isnothing(m) && return (nothing, nothing)
    creation = get(m, "created_timestamp", nothing)
    isnothing(creation) && return (nothing, nothing)
    tcreation = DateTime(creation, dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")
    ejection = get(m, "ejected_timestamp", nothing)
    tejection = isnothing(ejection) ?
                DateTime("3000-01-01T00:00:00.000Z", dateformat"yyyy-mm-ddTHH:MM:SS.sssZ") :
                DateTime(ejection, dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")
    return (tcreation, tejection)
end
