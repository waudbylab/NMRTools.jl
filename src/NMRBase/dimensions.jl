"""
    NMRDimension

Abstract supertype for all axes used in NMRData objects.

See also [`FrequencyDimension`](@ref) and [`NonFrequencyDimension`](@ref).
"""
abstract type NMRDimension{T} <: DimensionalData.Dimension{T} end

"""
    FrequencyDimension <: NMRDimension

Abstract supertype for frequency dimensions used in NMRData objects.
Concrete types `F1Dim`, `F2Dim`, `F3Dim` and `F4Dim` are generated for
use in creating objects.

See also [`NonFrequencyDimension`](@ref).
"""
abstract type FrequencyDimension{T} <: NMRDimension{T} end

"""
    NonFrequencyDimension <: NMRDimension

Abstract supertype for non-frequency dimensions used in NMRData objects.
Sub-types include [`TimeDimension`](@ref), [`GradientDimension`](@ref), 
and [`UnknownDimension`](@ref).

See also [`FrequencyDimension`](@ref).
"""
abstract type NonFrequencyDimension{T} <: NMRDimension{T} end

"""
    TimeDimension <: NonFrequencyDimension <: NMRDimension

Abstract supertype for time dimensions used in NMRData objects.
Concrete types `T1Dim`, `T2Dim`, `T3Dim` and `T4Dim` are generated for
time-domains representing frequency evolution, and `TrelaxDim` and
`TkinDim` are generated for representing relaxation and real-time
kinetics.
"""
abstract type TimeDimension{T} <: NonFrequencyDimension{T} end
# abstract type QuadratureDimension{T} <: NMRDimension{T} end

"""
    UnknownDimension <: NonFrequencyDimension <: NMRDimension

Abstract supertype for unknown, non-frequency dimensions used in NMRData objects.
Concrete types `X1Dim`, `X2Dim`, `X3Dim` and `X4Dim` are generated for
use in creating objects.
"""
abstract type UnknownDimension{T} <: NonFrequencyDimension{T} end

"""
    GradientDimension <: NonFrequencyDimension <: NMRDimension

Abstract supertype for gradient-encoded dimensions used in NMRData objects.
Concrete types `G1Dim`, `G2Dim`, `G3Dim` and `G4Dim` are generated for
use in creating objects.
"""
abstract type GradientDimension{T} <: NonFrequencyDimension{T} end

"""
    OffsetDimension <: NonFrequencyDimension <: NMRDimension

Abstract supertype for frequency offset dimensions used in NMRData objects.
Used for CEST saturation offsets, R1rho off-resonance offsets, etc.
Concrete type `OffsetDim` is generated for use in creating objects.
"""
abstract type OffsetDimension{T} <: NonFrequencyDimension{T} end

"""
    FieldDimension <: NonFrequencyDimension <: NMRDimension

Abstract supertype for RF field strength dimensions used in NMRData objects.
Used for spinlock power in R1rho dispersion experiments.
Concrete type `SpinlockDim` is generated for use in creating objects.
"""
abstract type FieldDimension{T} <: NonFrequencyDimension{T} end

# override DimensionalData.Dimensions macro to generate default metadata
macro NMRdim(typ::Symbol, supertyp::Symbol, args...)
    return NMRdimmacro(typ, supertyp, args...)
end
function NMRdimmacro(typ, supertype, name::String=string(typ))
    esc(quote
            Base.@__doc__ struct $typ{T} <: $supertype{T}
                val::T
            end
            function $typ(val::AbstractArray; kw...)
                v = values(kw)
                # e.g. values(kw) = (metadata = NoMetadata(),)
                if :metadata ∉ keys(kw)
                    # if no metadata defined, define it
                    # alternatively, if there is valid metadata, merge in the defaults
                    @debug "Creating default dimension metadata"
                    v = merge((metadata=defaultmetadata($typ),), v)
                elseif v[:metadata] isa Metadata{$typ}
                    @debug "Merging dimension metadata with defaults"
                    v2 = merge(v, (metadata=defaultmetadata($typ),))
                    merge!(v2[:metadata].val, v[:metadata].val)
                    v = v2
                elseif v[:metadata] isa Dict
                    @debug "Merging metadata dictionary with defaults"
                    md = v[:metadata]
                    v = merge(v, (metadata=defaultmetadata($typ),))
                    merge!(v[:metadata].val, md)
                else
                    # if NoMetadata (or an invalid type), define the correct default metadata
                    @debug "Dimension metadata is NoMetadata - replace with defaults"
                    v = merge(v, (metadata=defaultmetadata($typ),))
                end
                val = AutoLookup(val, v)
                return $typ{typeof(val)}(val)
                # @show tmpdim = $typ{typeof(val)}(val)
                # @show newlookup = DimensionalData.Dimensions._format(tmpdim, axes(tmpdim,1))
                # return $typ{typeof(newlookup)}(newlookup)
            end
            function $typ(val::T) where {T<:DimensionalData.Dimensions.Lookups.Lookup}
                # HACK - this would better be replaced with a call to DD.format in the function above
                # e.g.
                # DimensionalData.Dimensions.format(DimensionalData.Lookups.val(axH), DimensionalData.Lookups.basetypeof(axH), Base.OneTo(11))
                return $typ{T}(val)
            end
            $typ() = $typ(:)
            Dimensions.name(::Type{<:$typ}) = $(QuoteNode(Symbol(name)))
            Dimensions.key2dim(::Val{$(QuoteNode(typ))}) = $typ()
        end)
end

@NMRdim F1Dim FrequencyDimension
@NMRdim F2Dim FrequencyDimension
@NMRdim F3Dim FrequencyDimension
@NMRdim F4Dim FrequencyDimension
@NMRdim T1Dim TimeDimension
@NMRdim T2Dim TimeDimension
@NMRdim T3Dim TimeDimension
@NMRdim T4Dim TimeDimension
@NMRdim TrelaxDim TimeDimension
@NMRdim TkinDim TimeDimension
# @NMRdim Q1Dim QuadratureDimension
# @NMRdim Q2Dim QuadratureDimension
# @NMRdim Q3Dim QuadratureDimension
# @NMRdim Q4Dim QuadratureDimension
@NMRdim X1Dim UnknownDimension
@NMRdim X2Dim UnknownDimension
@NMRdim X3Dim UnknownDimension
@NMRdim X4Dim UnknownDimension
@NMRdim G1Dim GradientDimension
@NMRdim G2Dim GradientDimension
@NMRdim G3Dim GradientDimension
@NMRdim G4Dim GradientDimension
@NMRdim OffsetDim OffsetDimension
@NMRdim SpinlockDim FieldDimension
# @NMRdim SpatialDim NMRDimension

# Getters ########
"""
    data(nmrdimension)

Return the numerical data associated with an NMR dimension.
"""
data(d::NMRDimension) = d.val.data

"""
    ppm(axis)

Return the ppm values for points along a frequency axis.
"""
ppm(ax::FrequencyDimension) = data(ax)

"""
    hz(axis)

Return the offsets (in Hz) for points along a frequency axis.
"""
hz(ax::FrequencyDimension) = ax[:bf] * (data(ax) .- ax[:offsetppm]) * 1e-6

"""
    ppm(offset, axis)

Return the chemical shifts for a given offset or list of offsets along a frequency axis.
"""
function ppm(offset::Union{Number,AbstractArray{<:Number}}, ax::FrequencyDimension)
    return (offset .+ ax[:offsethz]) * 1e6 / ax[:bf]
end

"""
    hz(δ, axis)

Return the offset (in Hz) for a chemical shift (or list of shifts) on a frequency axis.
"""
hz(δ::Union{Number,AbstractArray{<:Number}}, ax::FrequencyDimension) = ax[:bf] *
                                                                       (δ .- ax[:offsetppm]) *
                                                                       1e-6

"""
    finddim(spec, nucleus::Nucleus)

Find the dimension index in an NMRData object that corresponds to the given nucleus.
Returns the dimension index (1-based) or `nothing` if not found.

# Examples
```julia
spec = loadnmr("hsqc.ft2")
finddim(spec, H1)   # Returns 1 for a 1H-15N HSQC
finddim(spec, N15)  # Returns 2 for a 1H-15N HSQC
```

See also [`shiftdim`](@ref), [`reference`](@ref).
"""
function finddim(spec, nuc::Nucleus)
    for (i, d) in enumerate(dims(spec))
        if d isa FrequencyDimension
            dim_nuc = metadata(d, :nucleus)
            if dim_nuc == nuc
                return i
            end
        end
    end
    return nothing
end

"""
    resolvedim(spec, dim_ref) -> Int

Resolve a dimension reference to a dimension index. Accepts:
- `Int`: dimension index (returned as-is)
- Dimension type (e.g., `F1Dim`): finds first matching dimension
- `Nucleus` (e.g., `H1`, `C13`): finds dimension with matching nucleus

Returns the dimension index or throws an error if not found.
"""
function resolvedim(spec, dim_ref)
    if dim_ref isa Int
        return dim_ref
    elseif dim_ref isa Nucleus
        dim_no = finddim(spec, dim_ref)
        if isnothing(dim_no)
            throw(ArgumentError("resolvedim: No dimension found for nucleus $dim_ref"))
        end
        return dim_no
    else
        # Assume it's a dimension type like F1Dim
        dim_no = findfirst(d -> d isa dim_ref, dims(spec))
        if isnothing(dim_no)
            throw(ArgumentError("resolvedim: Dimension $dim_ref not found in NMRData object"))
        end
        return dim_no
    end
end

# Mapping from dimension index to FrequencyDimension constructor
const FDIM_CONSTRUCTORS = Dict(1 => F1Dim, 2 => F2Dim, 3 => F3Dim, 4 => F4Dim)

"""
    shiftdim(data::NMRData, dim_ref, offset)

Add an offset to a frequency dimension in an NMRData object. The dimension can be
specified as a numerical index, a dimension type like `F1Dim`, or a `Nucleus` enum
(e.g., `H1`, `C13`).

The metadata is copied using `replacedimension`, and an entry is added or updated
in the dimension metadata to record the offset change.

# Examples
```julia
spec2 = shiftdim(spec, 1, 0.5)       # Using dimension index
spec2 = shiftdim(spec, F1Dim, 0.5)   # Using dimension type
spec2 = shiftdim(spec, H1, 0.5)      # Using nucleus
```

See also [`reference`](@ref), [`finddim`](@ref).
"""
function shiftdim(spec, dim_ref, offsetppm)
    dim_no = resolvedim(spec, dim_ref)

    # check that dim_no is valid (within the range of dims)
    if dim_no > length(spec.dims) || dim_no < 1
        throw(ArgumentError("shiftdim: Dimension index $dim_no is out of range"))
    end
    # check that the dimension is a FrequencyDimension
    if !(dims(spec, dim_no) isa FrequencyDimension)
        throw(ArgumentError("shiftdim: Dimension $dim_ref is not a FrequencyDimension"))
    end

    dim = dims(spec, dim_no)

    new_data = data(dim) .+ offsetppm
    md = deepcopy(metadata(dim))

    # Get the appropriate FrequencyDimension constructor
    DimConstructor = get(FDIM_CONSTRUCTORS, dim_no, nothing)
    if isnothing(DimConstructor)
        throw(ArgumentError("shiftdim: Unsupported dimension index $dim_no (max 4)"))
    end
    newdim = DimConstructor(new_data; metadata=md)

    new_spec = replacedimension(spec, dim_no, newdim)

    if :referenceoffset ∉ keys(metadata(new_spec, dim_no))
        metadata(new_spec, dim_no)[:referenceoffset] = offsetppm
    else
        metadata(new_spec, dim_no)[:referenceoffset] += offsetppm
    end

    # adjust other metadata
    metadata(new_spec, dim_no)[:offsetppm] += offsetppm
    metadata(new_spec, dim_no)[:offsethz] += 1e-6 * offsetppm *
                                             metadata(new_spec, dim_no, :bf)
    metadata(new_spec, dim_no)[:sf] += 1e-6 * offsetppm * metadata(new_spec, dim_no, :bf)

    return new_spec
end