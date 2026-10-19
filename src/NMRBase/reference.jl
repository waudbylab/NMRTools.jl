# Chemical shift referencing functions for NMR spectra

# Common aqueous solvents for auto-detection
const AQUEOUS_SOLVENTS = ["d2o", "h2o", "h2o+d2o", "water"]

"""
    isaqueous(spec) -> Bool

Determine whether a spectrum was acquired in an aqueous solvent.

# Arguments
- `spec`: NMRData object

Checks against the metadata field `:solvent`.
"""
function isaqueous(spec)
    solvent = metadata(spec, :solvent)
    if isnothing(solvent)
        @warn "Could not determine solvent from metadata, assuming organic (TMS reference)"
        return false
    end
    # Check if solvent matches any aqueous pattern
    solvent_lower = lowercase(string(solvent))
    if solvent_lower in AQUEOUS_SOLVENTS
        return true
    end
    return false
end

"""
    watershift(temperature_K)

Calculate the chemical shift of water at a given temperature using the empirical formula:

    δ(H₂O) = 7.83 − T / 96.9

where T is the temperature in Kelvin.

# Arguments
- `temperature_K`: Temperature in Kelvin

# Returns
- Chemical shift of water in ppm

# References
Bruker uses 4.7 ppm as the default water reference when locking (with SR=0).
"""
watershift(temperature_K) = 7.83 - temperature_K / 96.9

"""
    reference(spec; indirect=true, aqueous=:auto, temperature=nothing)
    reference(; kwargs...)

Reference a spectrum to water based on the sample temperature.

The forms without the data argument return a function for use with `|>`.

This method is used when no explicit reference shift is provided. It calculates the
expected water chemical shift from the temperature and applies a correction relative
to the Bruker default of 4.7 ppm.

# Arguments
- `spec`: NMRData object to reference
- `indirect=true`: If true, also reference heteronuclear dimensions using Xi ratios
- `temperature=nothing`: Temperature in Kelvin. If `nothing`, reads from metadata.

# Returns
- New NMRData object with referenced chemical shifts

Only works for aqueous solvents (D2O, H2O, H2O+D2O).

See also [`reference(spec, dim, pair)`](@ref), [`xi_ratio`](@ref).
"""
function reference(spec; indirect=true, temperature=nothing)
    # Check if aqueous
    if !isaqueous(spec)
        throw(NMRToolsError("reference: Water referencing only available for aqueous solvents. " *
                            "Use reference(spec, dim, old => new) for explicit referencing."))
    end

    # Get temperature
    temp = temperature
    if isnothing(temp)
        temp = metadata(spec, :temperature)
        if isnothing(temp)
            throw(NMRToolsError("reference: Temperature not found in metadata. " *
                                "Provide temperature explicitly using temperature=T keyword."))
        end
    end

    # Find 1H dimension
    h1_dim = finddim(spec, H1)
    if isnothing(h1_dim)
        throw(NMRToolsError("reference: No 1H dimension found. " *
                            "Water referencing requires a 1H dimension."))
    end

    # Calculate water shift at this temperature
    expected_water = watershift(temp)
    bruker_default = 4.7  # Bruker locks to 4.7 ppm by default

    # The offset to apply
    offset = expected_water - bruker_default

    @info "Water referencing" temperature = temp expected_δ_water = round(expected_water;
                                                                          digits=3) correction_ppm = round(offset;
                                                                                                           digits=4)

    # Apply referencing using the internal function
    return _reference_with_offset(spec, h1_dim, offset; indirect=indirect, aqueous=true)
end

"""
    reference(spec, dim, pair::Pair; indirect=true, aqueous=:auto)
    reference(dim, pair::Pair; kwargs...)

Reference a spectrum dimension by specifying an old and new chemical shift.

The forms without the data argument return a function for use with `|>`.

# Arguments
- `spec`: NMRData object to reference
- `dim`: Dimension to reference (Int, dimension type like `F1Dim`, or `Nucleus` like `H1`)
- `pair`: A pair `old_shift => new_shift` specifying the reference point
- `indirect=true`: If true, reference other frequency dimensions using Xi ratios
- `aqueous=:auto`: Solvent type for Xi ratio selection (`:auto`, `true`, or `false`)

# Examples
```julia
# Reference 1H dimension
spec2 = reference(spec, 1, -0.12 => 0.0)
spec2 = reference(spec, F1Dim, -0.12 => 0.0)
spec2 = reference(spec, H1, -0.12 => 0.0)

# Reference without indirect referencing of other dimensions
spec2 = reference(spec, H1, -0.12 => 0.0; indirect=false)
```

See also [`shiftdim`](@ref), [`xi_ratio`](@ref).
"""
function reference(spec, dim, pair::Pair{<:Number,<:Number}; indirect=true, aqueous=:auto)
    old_shift, new_shift = pair
    offset = new_shift - old_shift

    dim_no = resolvedim(spec, dim)

    # Get nucleus info for this dimension
    dim_obj = dims(spec, dim_no)
    nuc = metadata(dim_obj, :nucleus)
    nuc_str = isnothing(nuc) ? "dim $dim_no" : string(nuc)

    @info "Referencing $nuc_str" old_shift = old_shift new_shift = new_shift offset = round(offset;
                                                                                            digits=4)

    return _reference_with_offset(spec, dim_no, offset; indirect=indirect, aqueous=aqueous)
end

"""
    reference(spec, dims, pairs; indirect=true, aqueous=:auto)
    reference(dims, pairs; kwargs...)

Reference multiple spectrum dimensions simultaneously.

The forms without the data argument return a function for use with `|>`.

# Arguments
- `spec`: NMRData object to reference
- `dims`: Collection of dimensions to reference
- `pairs`: Collection of `old_shift => new_shift` pairs (same length as dims)
- `indirect=true`: If true, reference remaining frequency dimensions using Xi ratios
- `aqueous=:auto`: Solvent type for Xi ratio selection

# Examples
```julia
# Reference both 1H and 15N dimensions
spec2 = reference(spec, [F1Dim, F2Dim], [4.7 => 4.8, 120.0 => 119.5])
spec2 = reference(spec, (H1, N15), (4.7 => 4.8, 120.0 => 119.5))
```

See also [`xi_ratio`](@ref).
"""
function reference(spec, dims_refs, pairs; indirect=true, aqueous=:auto)
    if length(dims_refs) != length(pairs)
        throw(ArgumentError("Number of dimensions ($(length(dims_refs))) must match number of pairs ($(length(pairs)))"))
    end

    result = spec
    referenced_dims = Int[]

    # First pass: apply direct references (with indirect=false to avoid duplicate work)
    for (dim_ref, pair) in zip(dims_refs, pairs)
        dim_no = resolvedim(result, dim_ref)
        push!(referenced_dims, dim_no)

        old_shift, new_shift = pair
        offset = new_shift - old_shift

        dim_obj = dims(result, dim_no)
        nuc = metadata(dim_obj, :nucleus)
        nuc_str = isnothing(nuc) ? "dim $dim_no" : string(nuc)

        @info "Referencing $nuc_str" old_shift = old_shift new_shift = new_shift offset = round(offset;
                                                                                                digits=4)

        result = shiftdim(result, dim_no, offset)
    end

    # Second pass: apply indirect referencing if enabled
    if indirect
        result = _apply_indirect_referencing(result, referenced_dims; aqueous=aqueous)
    end

    return result
end

"""
    _reference_with_offset(spec, dim_no, offset; indirect=true, aqueous=:auto)

Internal function to apply a reference offset and optionally indirect referencing.
"""
function _reference_with_offset(spec, dim_no, offset; indirect=true, aqueous=:auto)
    # Apply the direct offset
    result = shiftdim(spec, dim_no, offset)

    # Apply indirect referencing if requested
    if indirect
        result = _apply_indirect_referencing(result, [dim_no]; aqueous=aqueous)
    end

    return result
end

"""
    _apply_indirect_referencing(spec, referenced_dims; aqueous=:auto)

Apply indirect referencing using Xi ratios to all frequency dimensions not in referenced_dims.

The effective 1H frequency is back-calculated from the referenced dimension(s), and Xi ratios
are used to calculate the expected frequencies for other nuclei.
"""
function _apply_indirect_referencing(spec, referenced_dims; aqueous=:auto)
    if aqueous == :auto
        is_aqueous = isaqueous(spec)
    else
        is_aqueous = aqueous
    end
    ref_compound = is_aqueous ? "DSS" : "TMS"

    # Find all frequency dimensions
    freq_dims = Int[]
    for (i, d) in enumerate(dims(spec))
        if d isa FrequencyDimension
            push!(freq_dims, i)
        end
    end

    # Determine which dimensions still need referencing
    unreferenced = setdiff(freq_dims, referenced_dims)
    if isempty(unreferenced)
        return spec
    end

    # We need to back-calculate the effective 1H reference frequency
    # from one of the referenced dimensions
    ref_dim = referenced_dims[1]
    ref_axis = dims(spec, ref_dim)
    ref_nuc = metadata(ref_axis, :nucleus)

    if isnothing(ref_nuc)
        @warn "Cannot apply indirect referencing: nucleus not known for dimension $ref_dim"
        return spec
    end

    ref_xi = xi_ratio(ref_nuc; aqueous=is_aqueous)
    if isnothing(ref_xi)
        @warn "Cannot apply indirect referencing: Xi ratio not defined for $(string(ref_nuc))"
        return spec
    end

    # Get the base frequency of the referenced dimension
    # After referencing, the carrier frequency (sf) gives us the spectrometer frequency
    # at the carrier position. We use bf (base frequency) for the calculation.
    ref_bf = metadata(ref_axis, :bf)

    # The effective 1H base frequency (in Hz)
    # Xi = ν(X) / ν(1H), so ν(1H) = ν(X) / Xi
    effective_h1_bf = ref_bf / ref_xi

    @info "Indirect referencing using $ref_compound Xi ratios" reference_nucleus = string(ref_nuc) effective_1H_bf_MHz = round(effective_h1_bf /
                                                                                                                               1e6;
                                                                                                                               digits=3)

    result = spec

    # Apply indirect referencing to remaining dimensions
    for dim_no in unreferenced
        dim_axis = dims(result, dim_no)
        target_nuc = metadata(dim_axis, :nucleus)

        if isnothing(target_nuc)
            @warn "Skipping dimension $dim_no: nucleus not known"
            continue
        end

        target_xi = xi_ratio(target_nuc; aqueous=is_aqueous)
        if isnothing(target_xi)
            @warn "Skipping $(string(target_nuc)): Xi ratio not defined"
            continue
        end

        # Calculate the expected base frequency for this nucleus
        expected_bf = effective_h1_bf * target_xi

        # Current base frequency
        current_bf = metadata(dim_axis, :bf)

        # The frequency difference translates to a ppm offset
        # Δν = expected_bf - current_bf (in Hz)
        # Δppm = Δν / current_bf * 1e6
        freq_diff = expected_bf - current_bf
        offset_ppm = freq_diff / current_bf * 1e6

        if abs(offset_ppm) > 0.0001  # Only apply if offset is significant
            @info "  Indirect referencing $(string(target_nuc))" expected_bf_MHz = round(expected_bf /
                                                                                         1e6;
                                                                                         digits=6) current_bf_MHz = round(current_bf /
                                                                                                                          1e6;
                                                                                                                          digits=6) offset_ppm = round(offset_ppm;
                                                                                                                                                       digits=4)
            result = shiftdim(result, dim_no, offset_ppm)
        else
            @info "  $(string(target_nuc)) already correctly referenced (offset < 0.0001 ppm)"
        end
    end

    return result
end

# curried forms for piping
reference(; kwargs...) = spec -> reference(spec; kwargs...)
reference(dim, pair::Pair{<:Number,<:Number}; kwargs...) = spec -> reference(spec, dim, pair; kwargs...)
function reference(dims::Union{Tuple,AbstractVector}, pairs::Union{Tuple,AbstractVector}; kwargs...)
    return spec -> reference(spec, dims, pairs; kwargs...)
end
