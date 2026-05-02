"""
    NMR Pulse Programme Annotation Parser for NMRTools.jl
    
    Parses semantic annotations from NMR pulse programme strings into structured dictionaries.
    Supports schema ≥ v0.0.2.
"""

"""
    annotate(spec::NMRData) -> NMRData

Parse pulse programme annotations and store them in the spectrum metadata.
Returns the NMRData object with annotations added to metadata.

Annotations are embedded within pulse programme comments marked with `;@ ` at the start
of lines, following YAML format with schema version 0.0.2.

# Examples
```julia
spec = loadnmr("experiment/1/pdata/1")
spec = annotate(spec)  # Parse and store annotations
annotations(spec, "title")  # Access parsed annotations
```

See also [`annotations`](@ref), [`hasannotations`](@ref).
"""
function annotate(spec::NMRData)
    pp = pulseprogram(spec)
    if isnothing(pp) || ismissing(pp) || length(pp) == 0
        return spec
    end
    parsed_annotations = parse_annotations(pp)

    # Check schema version
    schema_version = get(parsed_annotations, "schema_version", nothing)
    if isnothing(schema_version)
        return spec
    end

    # Parse and validate schema version (support ≥ 0.0.2)
    try
        version = VersionNumber(schema_version)
        min_version = v"0.0.2"
        if version < min_version
            @warn "Pulse programme uses unsupported schema version $schema_version. Minimum supported version is v$min_version. Annotations may not parse correctly."
            return spec
        end
    catch e
        @warn "Failed to parse schema version '$schema_version': $e. Annotations may not parse correctly."
        return spec
    end

    # Resolve parameter references to actual values
    resolve_parameter_references!(parsed_annotations, spec)

    # Resolve programmatic list patterns
    resolve_programmatic_lists!(parsed_annotations, spec)

    # Reverse dimensions to match data order
    if haskey(parsed_annotations, "dimensions")
        parsed_annotations["dimensions"] = reverse(parsed_annotations["dimensions"])
    end

    spec[:annotations] = parsed_annotations

    # Apply dimension transformations based on annotations
    spec = _apply_dimension_annotations(spec)

    return spec
end

"""
    _annotate_metadata!(md::Dict{Symbol,Any})

Parse pulse programme annotations from a raw metadata dict and store them under `:annotations`.
Only parses the raw YAML — skips parameter resolution and programmatic list expansion, which
require loaded data. Suitable for use with `NMRExperiment` descriptors.
"""
function _annotate_metadata!(md::Dict{Symbol,Any})
    acqus_dict = get(md, :acqus, nothing)
    isnothing(acqus_dict) && return md

    pp = get(acqus_dict, :pulseprogram_precomp, nothing)
    if isnothing(pp) || ismissing(pp) || isempty(pp)
        return md
    end

    parsed_annotations = parse_annotations(pp)

    schema_version = get(parsed_annotations, "schema_version", nothing)
    isnothing(schema_version) && return md

    # Parse and validate schema version (support ≥ 0.0.2)
    try
        version = VersionNumber(schema_version)
        min_version = v"0.0.2"
        if version < min_version
            @warn "Pulse programme uses unsupported schema version $schema_version. Minimum supported version is v$min_version."
            return md
        end
    catch e
        @warn "Failed to parse schema version '$schema_version': $e."
        return md
    end

    md[:annotations] = parsed_annotations
    return md
end

"""
    parse_annotations(content::String) -> Dict{String, Any}

Parse pulse programme annotations from content string. Annotations are embedded 
within pulse programme comments marked with `;@ ` at the start of lines.

The function:
1. Filters out lines starting with `;@ `
2. Strips the `;@ ` prefix from each line
3. Parses the resulting YAML document into a dictionary
"""
function parse_annotations(content::String)::Dict{String,Any}
    # Filter lines starting with ";@ " and strip the prefix
    annotation_lines = String[]

    for line in split(content, '\n')
        if startswith(line, ";@ ")
            # Strip the ";@ " prefix (4 characters)
            stripped_line = line[4:end]
            push!(annotation_lines, stripped_line)
        end
    end

    # If no annotation lines found, return empty dictionary
    if isempty(annotation_lines)
        return Dict{String,Any}()
    end

    # Join lines back into YAML content
    yaml_content = join(annotation_lines, '\n')

    # Parse YAML content
    try
        return YAML.load(yaml_content)
    catch e
        @debug "Failed to parse pulse programme annotations as YAML: $e"
        return Dict{String,Any}()
    end
end

"""
    resolve_parameter_references!(annotations::Dict{String, Any}, spec::NMRData)

Recursively traverse the annotations dictionary and resolve parameter references.
References are bare parameter names like p1, pl1, d20, VALIST, etc.

Parameters are resolved using the acqus() function to get values from the acqus file.

# Arguments
- `annotations::Dict{String, Any}`: Annotations dictionary to modify in-place
- `spec::NMRData`: NMRData object containing acqus metadata

# Examples
```julia
# Before resolution: {"pulse": "p1", "power": "pl1", "duration": "VDLIST"}
# After resolution: {"pulse": 9.2, "power": Power(...), "duration": [0.01, 0.02, ...]}
```
"""
function resolve_parameter_references!(annotations::Dict, spec::NMRData)
    _resolve_recursive!(annotations, spec)
    return annotations
end

function _resolve_recursive!(obj::Dict, spec::NMRData)
    for (key, value) in obj
        obj[key] = _resolve_recursive!(value, spec)
    end
    return obj
end

function _resolve_recursive!(obj::Vector, spec::NMRData)
    for i in eachindex(obj)
        obj[i] = _resolve_recursive!(obj[i], spec)
    end
    return obj
end

function _resolve_recursive!(obj::String, spec::NMRData)
    # Check if this string matches a parameter pattern
    resolved_value = _resolve_parameter(obj, spec)
    return isnothing(resolved_value) ? obj : resolved_value
end

function _resolve_recursive!(obj, spec::NMRData)
    # For any other type (numbers, booleans, etc.), return as-is
    return obj
end

"""
    _resolve_parameter(param_str::String, spec::NMRData) -> Union{Nothing, Any}

Attempt to resolve a parameter string using the acqus data.
Returns the resolved value or nothing if not a recognized parameter.

Handles parameter patterns (case-insensitive):
- Indexed parameters: p1, pl1, d20, cnst8, gpz6, etc.
- Nucleus assignments: f1, f2, f3, f4
- Any other parameter names that exist in acqus (including lists)
"""
function _resolve_parameter(param_str::String, spec::NMRData)
    param = strip(param_str)

    # Nucleus assignments: f1 => nuc1, f2 => nuc2, etc.
    m = match(r"^f(\d+)$", param)
    if !isnothing(m)
        index = parse(Int, m.captures[1])
        nuc = acqus(spec, Symbol("nuc$index"))
        # return nucleus(nuc) # causes bugs when trying to replace a string with a Nucleus type
        return nuc
    end

    # Define indexed parameter types (strings map to symbols in acqus)
    indexed_params = ["p", "pl", "sp", "spnam", "d", "cnst", "gpnam", "gpx", "gpy", "gpz"]

    # Try to match indexed parameters
    for prefix in indexed_params
        m = match(Regex("^$(prefix)(\\d+)\$"), param)
        isnothing(m) && continue

        index = parse(Int, m.captures[1])
        symbol = Symbol(prefix)
        return acqus(spec, symbol, index)
    end

    # Check if this parameter exists directly in acqus (handles lists and other parameters)
    # Convert to symbol for acqus lookup
    param_symbol = Symbol(param)
    direct_value = acqus(spec, param_symbol)

    if !isnothing(direct_value) && !ismissing(direct_value)
        return direct_value
    end

    # No match found - return nothing to keep original string
    return nothing
end

"""
    resolve_programmatic_lists!(annotations::Dict{String, Any}, spec::NMRData)

Resolve programmatic list patterns in annotations to actual vectors.

Programmatic lists are specified as dictionaries with one of two patterns:

**Pattern 1: Type-based (linear/log)**
- `type`: "linear" or "log" (logarithmic)
- `start`: starting value
- `step`: step size (for linear with start/step)
- `end`: ending value (for linear or log with start/end)

**Pattern 2: Counter-based**
- `counter`: array of values (typically integers from a loop counter)
- `scale`: scale factor to multiply each counter element by

For type-based patterns, the function:
1. Identifies fields containing programmatic list patterns ({type, start, step/end})
2. Finds the dimension number for that field from the dimensions annotation
3. Gets the number of points in that dimension
4. Replaces the pattern with a correctly sized vector

For counter-based patterns, the function:
1. Multiplies each element of the counter array by the scale value
2. Returns the resulting vector

# Examples
```julia
# Linear spacing with start/step
# {type: linear, start: p9, step: p9} with 10 points → [p9, 2*p9, 3*p9, ...]

# Linear spacing with start/end
# {type: linear, start: 0.001, end: 0.1} with 10 points → [0.001, 0.012, ..., 0.1]

# Logarithmic spacing with start/end
# {type: log, start: 0.001, end: 0.1} with 10 points → [0.001, 0.0016, ..., 0.1]

# Counter-based pattern
# {counter: [0, 1, 2, 3], scale: 0.01} → [0.0, 0.01, 0.02, 0.03]
```
"""
function resolve_programmatic_lists!(annotations::Dict, spec::NMRData)
    # Get dimensions array to map field names to dimension numbers
    dimensions = get(annotations, "dimensions", nothing)
    if isnothing(dimensions)
        return annotations
    end

    # Find all programmatic list patterns and resolve them
    _resolve_programmatic_recursive!(annotations, annotations, spec, dimensions)
    return annotations
end

"""
    _resolve_programmatic_recursive!(obj, root_annotations, spec, dimensions)

Recursively traverse annotations to find and resolve programmatic list patterns.
"""
function _resolve_programmatic_recursive!(obj::Dict, root_annotations::Dict, spec::NMRData,
                                          dimensions::Vector)
    for (key, value) in obj
        if _is_programmatic_pattern(value)
            # This is a programmatic list pattern - resolve it
            resolved = _resolve_programmatic_pattern(key, value, root_annotations, spec,
                                                     dimensions)
            if !isnothing(resolved)
                obj[key] = resolved
            end
        else
            # Recursively process nested structures
            _resolve_programmatic_recursive!(value, root_annotations, spec, dimensions)
        end
    end
    return obj
end

function _resolve_programmatic_recursive!(obj::Vector, root_annotations::Dict,
                                          spec::NMRData, dimensions::Vector)
    # Don't process vectors themselves, only their elements if they're dicts
    for i in eachindex(obj)
        if obj[i] isa Dict
            _resolve_programmatic_recursive!(obj[i], root_annotations, spec, dimensions)
        end
    end
    return obj
end

function _resolve_programmatic_recursive!(obj, root_annotations::Dict, spec::NMRData,
                                          dimensions::Vector)
    # For any other type, return as-is
    return obj
end

"""
    _is_programmatic_pattern(value) -> Bool

Check if a value is a programmatic list pattern.
Must be a Dict with either:
- "type" key and either "step" or "end" key (linear/log patterns)
- "counter" and "scale" keys (counter-based pattern)
"""
function _is_programmatic_pattern(value)
    if !(value isa Dict)
        return false
    end

    # Pattern 1: type-based (linear/log) - must have "type" and either "step" or "end"
    has_type = haskey(value, "type")
    has_step_or_end = haskey(value, "step") || haskey(value, "end")
    is_type_pattern = has_type && has_step_or_end

    # Pattern 2: counter-based - must have both "counter" and "scale"
    has_counter = haskey(value, "counter")
    has_scale = haskey(value, "scale")
    is_counter_pattern = has_counter && has_scale

    return is_type_pattern || is_counter_pattern
end

"""
    _resolve_programmatic_pattern(key, pattern, root_annotations, spec, dimensions)

Resolve a single programmatic list pattern to a vector.
"""
function _resolve_programmatic_pattern(key::String, pattern::Dict, root_annotations::Dict,
                                       spec::NMRData, dimensions::Vector)
    # Check if this is a counter-based pattern
    if haskey(pattern, "counter") && haskey(pattern, "scale")
        counter = get(pattern, "counter", nothing)
        scale = get(pattern, "scale", nothing)

        if isnothing(counter)
            @debug "Counter-based pattern missing 'counter' value for: $key"
            return nothing
        end

        if isnothing(scale)
            @debug "Counter-based pattern missing 'scale' value for: $key"
            return nothing
        end

        # Counter must be a vector or array
        if !(counter isa AbstractVector)
            @debug "Counter-based pattern 'counter' must be an array for: $key"
            return nothing
        end

        # Generate list by multiplying each counter element by scale
        return [c * scale for c in counter]
    end

    # Otherwise, handle type-based patterns (linear/log)
    # Find which dimension this key corresponds to
    dim_path = _find_dimension_path(key, root_annotations)
    if isnothing(dim_path)
        @debug "Could not find dimension path for programmatic list key: $key"
        return nothing
    end

    # Find dimension number in the dimensions array
    dim_index = findfirst(==(dim_path), dimensions)
    if isnothing(dim_index)
        @debug "Could not find dimension for programmatic list: $dim_path not in dimensions array"
        return nothing
    end

    # Get number of points in this dimension (before reversal)
    npoints = size(spec, length(dimensions) - dim_index + 1)

    # Get pattern parameters
    pattern_type = get(pattern, "type", "linear")
    start_val = get(pattern, "start", nothing)
    step_val = get(pattern, "step", nothing)
    end_val = get(pattern, "end", nothing)

    if isnothing(start_val)
        @debug "Programmatic list pattern missing 'start' value for: $key"
        return nothing
    end

    # Generate the list based on pattern type
    if pattern_type == "linear"
        if !isnothing(step_val)
            # Linear with start/step
            return [start_val + (i - 1) * step_val for i in 1:npoints]
        elseif !isnothing(end_val)
            # Linear with start/end
            return collect(range(start_val, end_val; length=npoints))
        else
            @debug "Programmatic list pattern must have either 'step' or 'end': $key"
            return nothing
        end
    elseif pattern_type == "log"
        if !isnothing(end_val)
            # Logarithmic spacing
            return exp.(range(log(start_val), log(end_val); length=npoints))
        else
            @debug "Logarithmic programmatic list pattern requires 'end' value: $key"
            return nothing
        end
    else
        @debug "Unknown programmatic list type: $pattern_type (expected 'linear' or 'log')"
        return nothing
    end
end

"""
    _find_dimension_path(key, annotations, prefix="")

Find the dotted path to a key in the annotations dictionary.
Returns the path that should match entries in the dimensions array.
"""
function _find_dimension_path(key::String, annotations::Dict, prefix::String="")
    # Check if this key exists at the current level
    if haskey(annotations, key)
        return isempty(prefix) ? key : "$prefix.$key"
    end

    # Recursively search in nested dictionaries
    for (k, v) in annotations
        if v isa Dict && k != "dimensions"  # Don't search in dimensions array itself
            new_prefix = isempty(prefix) ? k : "$prefix.$k"
            result = _find_dimension_path(key, v, new_prefix)
            if !isnothing(result)
                return result
            end
        end
    end

    return nothing
end

# Dimension application logic ################################################################

"""
    _apply_dimension_annotations(spec::NMRData) -> NMRData

Apply semantic dimension types based on the dimensions array in annotations.
The dimensions array provides explicit ordering - index i in the array
corresponds to dimension i in the data. Returns new NMRData with updated dimensions.
"""
function _apply_dimension_annotations(spec::NMRData)
    dimensions = NMRBase.annotations(spec, :dimensions)
    isnothing(dimensions) && return spec

    # dimensions[i] specifies what dimension i represents
    for (dim_index, dim_spec) in enumerate(dimensions)
        spec = _apply_single_dimension(spec, dim_index, dim_spec)
    end

    return spec
end

function _apply_single_dimension(spec, dim_index, dim_spec::String)
    parts = split(dim_spec, ".")

    # Skip simple dimensions like "f1", "f2" - already correct
    length(parts) == 1 && return spec

    block_name = parts[1]  # e.g., "relaxation", "diffusion", "cest", "r1rho"
    param_name = parts[2]  # e.g., "duration", "offset", "g", "power"

    # Get the parameter block from annotations
    list = NMRBase.annotations(spec, dim_spec)
    isnothing(list) && return spec

    # Adjust list size to match dimension if needed
    if list isa AbstractVector
        ndata = size(spec, dim_index)
        nlist = length(list)
        if nlist != ndata
            if nlist > ndata
                @warn "Annotation list for $dim_spec has $nlist entries but dimension $dim_index has $ndata points; truncating"
                list = list[1:ndata]
            else
                @warn "Annotation list for $dim_spec has $nlist entries but dimension $dim_index has $ndata points; extending cyclically"
                list = list[mod.(0:(ndata - 1), nlist) .+ 1]
            end
            spec[:annotations][block_name][param_name] = list
        end
    end

    # Dispatch to specific handler based on parameter type
    if param_name == "duration"
        return _apply_duration(spec, dim_index, list, block_name)
    elseif param_name == "g"
        return _apply_gradient(spec, dim_index, list, block_name)
    elseif param_name == "offset"
        return _apply_offset(spec, dim_index, list, block_name)
    elseif param_name == "power"
        return _apply_power(spec, dim_index, list, block_name)
    end

    return spec
end

function _apply_duration(spec, dim_index, values, block_name)
    # All duration arrays use TrelaxDim
    label_text = if block_name == "relaxation"
        "Relaxation delay"
    elseif block_name == "r1rho"
        "Spinlock duration"
    elseif block_name == "calibration"
        "Pulse duration"
    else
        "Delay"
    end

    newdim = TrelaxDim(values)
    spec = replacedimension(spec, dim_index, newdim)
    label!(spec, dim_index, label_text)
    spec[dim_index, :units] = "s"
    spec[dim_index, :delay_type] = Symbol(block_name)
    return spec
end

function _apply_gradient(spec, dim_index, values, block_name)
    block = NMRBase.annotations(spec, block_name)
    gmax = get(block, "gmax", nothing)
    if isnothing(gmax)
        @warn "No gmax specified for diffusion gradient axis"
        gmax = 0.55  # Default assumption
    end
    gvals = gmax .* values

    # Select appropriate gradient dimension type based on index
    DimConstructor = if dim_index == 1
        G1Dim
    elseif dim_index == 2
        G2Dim
    elseif dim_index == 3
        G3Dim
    else
        G4Dim
    end

    newdim = DimConstructor(gvals)
    spec = replacedimension(spec, dim_index, newdim)
    label!(spec, dim_index, "Gradient strength")
    spec[dim_index, :units] = "T m⁻¹"
    return spec
end

function _apply_offset(spec, dim_index, values, block_name)
    block = NMRBase.annotations(spec, block_name)
    channel = get(block, "channel", nothing)

    # Convert FQList to ppm if needed
    if values isa FQList
        freq_dim_index = _find_frequency_dim_for_channel(spec, channel)
        if !isnothing(freq_dim_index)
            freq_dim = dims(spec, freq_dim_index)
            offset_values = ppm(values, freq_dim)
        else
            offset_values = data(values)  # Fall back to raw values
        end
    else
        offset_values = values
    end

    newdim = OffsetDim(offset_values)
    spec = replacedimension(spec, dim_index, newdim)
    label!(spec, dim_index, "Saturation offset")
    spec[dim_index, :units] = "ppm"
    if !isnothing(channel)
        spec[dim_index, :nucleus] = channel
    end
    return spec
end

function _apply_power(spec, dim_index, values, block_name)
    block = NMRBase.annotations(spec, block_name)
    channel = get(block, "channel", nothing)

    # Convert Power values to Hz using reference pulse
    if !isnothing(channel)
        ref_pulse = referencepulse(spec, channel)
        if isnothing(ref_pulse)
            @debug "No reference pulse found for channel $channel, cannot convert power to Hz"
            return spec
        end
        ref_duration, ref_power = ref_pulse

        # Convert Power array to Hz
        if !isempty(values) && first(values) isa Power
            field_hz = [hz(p, ref_power, ref_duration, 90.0) for p in values]
        else
            field_hz = values  # Assume already in Hz
        end
    else
        field_hz = values
    end

    newdim = SpinlockDim(field_hz)
    spec = replacedimension(spec, dim_index, newdim)
    label!(spec, dim_index, "Spinlock field")
    spec[dim_index, :units] = "Hz"
    return spec
end

"""
    _find_frequency_dim_for_channel(spec, channel) -> Union{Int, Nothing}

Find the frequency dimension index that corresponds to the given channel.
Channel can be "f1", "f2", etc. or a nucleus string like "1H", "19F".
"""
function _find_frequency_dim_for_channel(spec, channel::AbstractString)
    channel = String(channel)  # Convert SubString to String
    # If channel is "f1", "f2", look up the nucleus
    m = match(r"^f(\d+)$", channel)
    if !isnothing(m)
        nuc_index = parse(Int, m.captures[1])
        nuc_str = acqus(spec, Symbol("nuc$nuc_index"))
        if !isnothing(nuc_str) && !ismissing(nuc_str)
            channel = nuc_str
        end
    end

    # Find dimension with matching nucleus
    for (i, d) in enumerate(dims(spec))
        if d isa FrequencyDimension
            dim_nuc = metadata(d, :nucleus)
            if !isnothing(dim_nuc) && string(dim_nuc) == channel
                return i
            end
        end
    end
    return nothing
end

_find_frequency_dim_for_channel(spec, channel::Nothing) = nothing