# default constructor: dictionary of Symbol => Any #################################################

# struct NMRMetadata <: DimensionalData.Dimensions.LookupArrays.AbstractMetadata
# Metadata() = Metadata(Dict{Symbol,Any}())

# getters ##########################################################################################

# metadata accessor functions
"""
    metadata(nmrdata, key)
    metadata(nmrdata, dim, key)
    metadata(nmrdimension, key)

Return the metadata for specified key, or `nothing` if not found. Keys are passed as symbols.

# Examples (spectrum metadata)
- `:ns`: number of scans
- `:ds`: number of dummy scans
- `:rg`: receiver gain
- `:ndim`: number of dimensions
- `:title`: spectrum title (contents of title pdata file)
- `:filename`: spectrum filename
- `:pulseprogram`: title of pulse program used for acquisition
- `:experimentfolder`: path to experiment
- `:noise`: RMS noise level

# Examples (dimension metadata)
- `:pseudodim`: flag indicating non-frequency domain data
- `:npoints`: final number of (real) data points in dimension (after extraction)
- `:td`: number of complex points acquired
- `:tdzf`: number of complex points when FT executed, including LP and ZF
- `:bf`: base frequency, in Hz
- `:sf`: carrier frequency, in Hz
- `:offsethz`: carrier offset from bf, in Hz
- `:offsetppm`: carrier offset from bf, in ppm
- `:swhz`: spectrum width, in Hz
- `:swppm`: spectrum width, in ppm
- `:region`: extracted region, expressed as a range in points, otherwise missing
- `:window`: `WindowFunction` indicating applied apodization
- `:referenceoffset`: referencing (in ppm) applied to the dimension

See also [`estimatenoise`](@ref).
"""
function metadata end

metadata(A::AbstractNMRData, key::Symbol) = get(metadata(A), key, nothing)
metadata(A::AbstractNMRData, dim) = metadata(dims(A, dim))
metadata(A::AbstractNMRData, dim, key::Symbol) = get(metadata(A, dim), key, nothing)
metadata(d::NMRDimension, key::Symbol) = get(metadata(d), key, nothing)

# NMRExperiment metadata accessors
metadata(e::NMRExperiment) = e.metadata
metadata(e::NMRExperiment, key::Symbol) = get(e.metadata, key, nothing)

Base.getindex(A::AbstractNMRData, key::Symbol) = metadata(A, key)
Base.getindex(A::AbstractNMRData, dim, key::Symbol) = metadata(A, dim, key)
Base.getindex(d::NMRDimension, key::Symbol) = metadata(d, key)
Base.setindex!(A::AbstractNMRData, v, key::Symbol) = setindex!(metadata(A), v, key)  #(A[key] = v  =>  metadata(A)[key] = v)
function Base.setindex!(A::AbstractNMRData, v, dim, key::Symbol)
    return setindex!(metadata(A, dim), v, key)
end  #(A[dim, key] = v  =>  metadata(A, dim)[key] = v)
Base.setindex!(d::NMRDimension, v, key::Symbol) = setindex!(metadata(d), v, key)  #(d[key] = v  =>  metadata(d)[key] = v)
Base.getindex(e::NMRExperiment, key::Symbol) = metadata(e, key)
Base.setindex!(e::NMRExperiment, v, key::Symbol) = (e.metadata[key] = v)

"""
    units(nmrdata)
    units(nmrdata, dim)
    units(nmrdimension)

Return the physical units associated with an `NMRData` structure or an `NMRDimension`.
"""
function units end
units(A::AbstractNMRData) = metadata(A, :units)
units(A::AbstractNMRData, dim) = metadata(A, dim, :units)
units(d::NMRDimension) = get(metadata(d), :units, nothing)

"""
    label(nmrdata)
    label(nmrdata, dim)
    label(nmrdimension)

Return a short label associated with an `NMRData` structure or an `NMRDimension`.
By default, for a spectrum this is obtained from the first line of the title file.
For a frequency dimension, this is normally something of the form `1H chemical shift (ppm)`.

See also [`label!`](@ref).
"""
function label end
label(A::AbstractNMRData) = metadata(A, :label)
label(e::NMRExperiment) = metadata(e, :label)
label(A::AbstractNMRData, dim) = metadata(A, dim, :label)
label(d::NMRDimension) = get(metadata(d), :label, nothing)

"""
    label!(nmrdata, labeltext)
    label!(nmrdata, dim, labeltext)
    label!(nmrdimension, labeltext)

Set the label associated with an `NMRData` structure or an `NMRDimension`.

See also [`label`](@ref).
"""
function label! end
label!(A::AbstractNMRData, labeltext::AbstractString) = (A[:label] = labeltext)
label!(A::AbstractNMRData, dim, labeltext::AbstractString) = (A[dim, :label] = labeltext)
label!(d::NMRDimension, labeltext::AbstractString) = (d[:label] = labeltext)

"""
    acqus(nmrdata)
    acqus(nmrdata, key)
    acqus(nmrdata, key, index)

Return data from a Bruker acqus file, or `nothing` if it does not exist.
Keys can be passed as symbols or strings. If no key is specified, a dictionary
is returned representing the entire acqus file.

If present, the contents of auxilliary files such as `vclist` and `vdlist` can
be accessed using this function.

# Examples
```julia-repl
julia> acqus(expt, :pulprog)
"zgesgp"
julia> acqus(expt, "TE")
276.9988
julia> acqus(expt, :p, 1)
9.2
julia> acqus(expt, "D", 1)
0.1
julia> acqus(expt, :vclist)
11-element Vector{Int64}:
[...]
```

See also [`metadata`](@ref).
"""
acqus(A::AbstractNMRData) = metadata(A, :acqus)
acqus(e::NMRExperiment) = metadata(e, :acqus)
function acqus(A::Union{AbstractNMRData,NMRExperiment}, key::Symbol)
    return ismissing(acqus(A)) ? missing : get(acqus(A), key, missing)
end
function acqus(A::Union{AbstractNMRData,NMRExperiment}, key::String)
    return acqus(A, Symbol(lowercase(key)))
end
acqus(A::Union{AbstractNMRData,NMRExperiment}, key, index) = acqus(A, key)[index]

"""
    annotations(nmrdata)
    annotations(nmrdata, key)
    annotations(nmrdata, key, index)
    annotations(nmrdata, key1, key2)
    annotations(nmrdata, keys...)

Return annotation data from pulse programme, or `nothing` if it does not exist.
Keys can be passed as strings or symbols. If no key is specified, the entire
annotations dictionary is returned.

Dotted names are automatically split and dereferenced, e.g. `annotations(spec, "cest.duration")`
is equivalent to `annotations(spec, "cest", "duration")`.

This function provides nested access to annotations stored in `metadata[:annotations]`.
Multiple keys can be chained to access nested dictionaries or specific array elements.

# Examples
```julia-repl
julia> annotations(spec, "title")
"19F CEST"

julia> annotations(spec, "dimensions")
2-element Vector{String}:
 "cest.offset"
 "f1"

julia> annotations(spec, "dimensions", 1)
"cest.offset"

julia> annotations(spec, "cest", "duration")
0.5

julia> annotations(spec, "cest.duration")
0.5

julia> annotations(spec, :reference_pulse)
1-element Vector{Dict{String, Any}}:
 Dict("channel" => "f1", "pulse" => 9.2e-6, "power" => -3.0)

julia> annotations(spec, "calibration.duration.start")
0.001
```

See also [`metadata`](@ref), [`acqus`](@ref).
"""
annotations(x) = metadata(x, :annotations)

# Variadic version: handles any number of keys
function annotations(A, keys...)
    annot = annotations(A)
    isnothing(annot) && return nothing

    # Expand dotted notation and navigate
    current = annot
    for key in keys
        if key isa Integer
            # Array indexing
            if current isa AbstractVector && 1 <= key <= length(current)
                current = current[key]
            else
                return nothing
            end
        else
            # String/symbol key - may contain dots
            key_parts = (k -> split(k, '.'))(string(key))

            for part in key_parts
                if current isa AbstractDict
                    current = get(current, part, nothing)
                    isnothing(current) && return nothing
                else
                    return nothing
                end
            end
        end
    end

    return current
end

"""
    hasannotations(nmrdata) -> Bool

Check if the NMRData object has any annotations.
"""
function hasannotations(spec)
    annot = annotations(spec)
    return !isnothing(annot) && !isempty(annot)
end

"""
    sample(nmrdata)
    sample(nmrdata, keys...)

Return the sample metadata associated with the NMRData object, or `nothing` if not available.
Sample metadata follows the schema defined at https://github.com/waudbylab/nmr-sample-schema.

When keys are provided, navigates nested dictionaries using the given keys. Keys can be strings
or symbols and are automatically converted to lowercase strings for consistent access.

# Examples
```julia-repl
julia> sample(spec)
Dict{String, Any} with 6 entries:
  "nmr_tube" => Dict{String, Any}("diameter"=>"5 mm", "type"=>"regular", "sample_volume_uL"=>600)
  "buffer"   => Dict{String, Any}("solvent"=>"10% D2O", "reference_unit"=>"%w/v", "chemical_shift_reference"=>"none")
  [...]

julia> sample(spec, "notes")
"here's a note!"

julia> sample(spec, :people, :users)
1-element Vector{Any}:
 "Chris"

julia> sample(spec, "sample", "components")
2-element Vector{Any}:
 Dict{String, Any}("name" => "HEWL", "isotopic_labelling" => "unlabelled", "unit" => "mM", "concentration" => 10)
 Dict{String, Any}("name" => "gadodiamide", "isotopic_labelling" => "unlabelled", "unit" => "mM", "concentration" => 1)

julia> sample(spec, :sample, :components, :name)
2-element Vector{Any}:
 "HEWL"
 "gadodiamide"
```

When the final key is applied to an array (such as `components`), the key is mapped over each element and a filtered array is returned.

For `isotopic_labelling` and `solvent`, if the stored value is `"custom"`, the corresponding `custom_labelling` or `custom_solvent` field is returned automatically.

See also [`metadata`](@ref), [`hassample`](@ref).
"""
function sample(spec)
    return metadata(spec, :sample)
end

sample(s::NMRSample) = s

# Variadic version: handles any number of keys, works for NMRData/NMRExperiment and NMRSample directly
function sample(spec, keys::Union{String,Symbol}...)
    sample_data = sample(spec)
    isnothing(sample_data) && return nothing
    return sample(sample_data, keys...)
end

function sample(s::NMRSample, keys::Union{String,Symbol}...)
    return _sample_traverse(s.metadata, keys...)
end

# Fields whose value "custom" should resolve to a sibling custom_* field
const _CUSTOM_FIELD_MAP = Dict{String,String}(
    "isotopic_labelling" => "custom_labelling",
    "solvent" => "custom_solvent",
)

function _sample_traverse(current, keys::Union{String,Symbol}...)
    isempty(keys) && return current
    key_str = lowercase(string(first(keys)))
    remaining = Base.tail(keys)
    if current isa AbstractArray
        results = [_sample_traverse(el, keys...) for el in current]
        return filter(!isnothing, results)
    end
    current isa AbstractDict || return nothing
    value = get(current, key_str, nothing)
    isnothing(value) && return nothing
    if value == "custom"
        custom_key = get(_CUSTOM_FIELD_MAP, key_str, nothing)
        if !isnothing(custom_key)
            value = get(current, custom_key, value)
        end
    end
    return _sample_traverse(value, remaining...)
end

"""
    hassample(nmrdata) -> Bool

Check if the NMRData object has any associated sample metadata.
"""
function hassample(spec)
    return !isnothing(sample(spec))
end

"""
    samplefile(x) -> Union{String, Nothing}

Return the path to the JSON sample file associated with `x`, or `nothing` if no
sample has been matched. Works with [`NMRData`](@ref) and [`NMRExperiment`](@ref).

See also [`sample`](@ref), [`hassample`](@ref).
"""
samplefile(x) = hassample(x) ? sample(x).path : nothing
samplefile(s::NMRSample) = s.path

"""
    referencepulse(spec, nucleus) -> (pulse_length, power)

Get the reference pulse calibration for a given nucleus from pulse programme annotations.
Returns a tuple of (pulse_length, power) or `nothing` if not found.

The reference pulse is the calibrated 90° hard pulse used as the basis for all power 
calculations on that channel.

# Arguments
- `spec`: NMRData object with annotations
- `nucleus`: Nucleus symbol or string (e.g., :19F, "19F", :1H, "1H")

# Returns
- `(pulse_length, power)`: Tuple where pulse_length is in μs and power is a Power object
- `nothing`: If no reference pulse found for the specified nucleus

# Examples
```julia-repl
julia> referencepulse(spec, "19F")
(13.29, Power(-9.03 dB, 8.0 W))

julia> referencepulse(spec, :1H)
(9.2, Power(-3.0 dB, 32.0 W))
```

See also [`annotations`](@ref).
"""
function referencepulse(spec::AbstractNMRData, nucleus)
    ref_pulses = annotations(spec, "reference_pulse")
    isnothing(ref_pulses) && return nothing

    nucleus_str = string(nucleus)

    for pulse in ref_pulses
        if get(pulse, "channel", nothing) == nucleus_str
            pulse_length = get(pulse, "duration", nothing)
            power = get(pulse, "power", nothing)

            if !isnothing(pulse_length) && !isnothing(power)
                return (pulse_length, power)
            end
        end
    end

    return nothing
end

# # Metadata for NMRData #############################################################################

function defaultmetadata(::Type{<:AbstractNMRData})
    defaults = Dict{Symbol,Any}(:title => "",
                                :label => "",
                                :filename => nothing,
                                :format => nothing,
                                :lastserialised => nothing,
                                :pulseprogram => nothing,
                                :ns => nothing,
                                :rg => nothing,
                                :noise => nothing)

    return Metadata{NMRData}(defaults)
end

# Metadata for Dimensions ##########################################################################

function defaultmetadata(T::Type{<:NMRDimension})
    defaults = Dict{Symbol,Any}(:label => "",
                                :units => nothing)
    return Metadata{T}(defaults)
end

function defaultmetadata(::Type{<:FrequencyDimension})
    defaults = Dict{Symbol,Any}(:label => "",
                                :coherence => nothing,
                                :bf => nothing,
                                :offsethz => nothing,
                                :offsetppm => nothing,
                                :swhz => nothing,
                                :swppm => nothing,
                                :td => nothing,
                                :tdzf => nothing,
                                :npoints => nothing,
                                :region => nothing,
                                :units => nothing,
                                :window => nothing,
                                :referenceoffset => 0.0,
                                :mcindex => nothing)
    return Metadata{FrequencyDimension}(defaults)
end

function defaultmetadata(::Type{<:TimeDimension})
    defaults = Dict{Symbol,Any}(:label => "",
                                :units => nothing,
                                :window => nothing,
                                :mcindex => nothing)
    return Metadata{TimeDimension}(defaults)
end

function defaultmetadata(::Type{<:TrelaxDim})
    defaults = Dict{Symbol,Any}(:label => "Relaxation time",
                                :units => nothing,
                                :delay_type => nothing)
    return Metadata{TrelaxDim}(defaults)
end

function defaultmetadata(::Type{<:TkinDim})
    defaults = Dict{Symbol,Any}(:label => "Time elapsed",
                                :units => nothing)
    return Metadata{TkinDim}(defaults)
end

function defaultmetadata(::Type{<:UnknownDimension})
    defaults = Dict{Symbol,Any}(:label => "",
                                :units => nothing)
    return Metadata{UnknownDimension}(defaults)
end

function defaultmetadata(::Type{<:OffsetDimension})
    defaults = Dict{Symbol,Any}(:label => "Offset",
                                :units => "ppm",
                                :nucleus => nothing)
    return Metadata{OffsetDimension}(defaults)
end

function defaultmetadata(::Type{<:FieldDimension})
    defaults = Dict{Symbol,Any}(:label => "Field strength",
                                :units => "Hz")
    return Metadata{FieldDimension}(defaults)
end

# Metadata entry definitions and help ##############################################################
function metadatahelp(A::Type{T}) where {T<:Union{NMRData,NMRDimension}}
    m = Metadata{A}()
    return Dict(map(kv -> kv[1] => metadatahelp(kv[1]), m))
end

function metadatahelp(key::Symbol)
    return get(metadatahelp(), key, "No information is available on this entry.")
end

function metadatahelp()
    d = Dict(:title => "Experiment title (may be multi-line)",
             :label => "Short experiment or axis description",
             :filename => "Original filename or template",
             :format => ":NMRPipe or :pdata",
             :pulseprogram => "Pulse program (PULPROG) from acqus file",
             :ns => "Number of scans",
             :rg => "Receiver gain",
             :noise => "RMS noise level",
             :units => "Axis units",
             :window => "Window function",
             :mcindex => "Index of imaginary component associated with axis (for multicomplex data)")
    return d
end