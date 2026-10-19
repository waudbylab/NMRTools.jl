"""
    AbstractNMRData <: DimensionalData.AbstractDimArray

Abstract supertype for objects that wrap an array of NMR data, and metadata about
its contents.

`AbstractNMRData`s inherit from [`AbstractDimArray`]($DDarraydocs)
from DimensionalData.jl. They can be indexed as regular Julia arrays or with
DimensionalData.jl [`Dimension`]($DDdimdocs)s.
"""
abstract type AbstractNMRData{T,N,D,A} <: AbstractDimArray{T,N,D,A} end

"""
    NMRSample

A lightweight wrapper around an NMR sample JSON file. Holds the file path and
parsed metadata. Access fields via [`sample`](@ref).

See also: [`scansamples`](@ref NMRTools.NMRIO.scansamples), [`findsample`](@ref NMRTools.NMRIO.findsample).
"""
struct NMRSample
    path::String
    metadata::Dict{String,Any}
end

"""
    NMRExperiment

A lightweight, metadata-only descriptor for a Bruker NMR experiment folder.
Does not load any binary data. Upgrade to a full [`NMRData`](@ref) via `loadnmr`.

See also: [`scanexperiments`](@ref NMRTools.NMRIO.scanexperiments), [`findsample`](@ref NMRTools.NMRIO.findsample), [`loadnmr`](@ref NMRTools.NMRIO.loadnmr).
"""
struct NMRExperiment
    path::String
    metadata::Dict{Symbol,Any}
end

function Base.show(io::IO, e::NMRExperiment)
    lbl = get(e.metadata, :label, nothing)
    pp = get(e.metadata, :pulseprogram, nothing)
    d = get(e.metadata, :date, nothing)
    print(io, "NMRExperiment(", repr(e.path))
    isnothing(lbl) || isempty(lbl) || print(io, " \"$lbl\"")
    isnothing(pp) || print(io, " [$pp]")
    isnothing(d) || print(io, " $(Dates.format(d, "yyyy-mm-dd"))")
    return print(io, ")")
end

function Base.show(io::IO, s::NMRSample)
    lbl = get(get(s.metadata, "sample", Dict()), "label", nothing)
    print(io, "NMRSample(", repr(s.path))
    isnothing(lbl) || print(io, " \"$lbl\"")
    return print(io, ")")
end

# Base methods #####################################################################################

function Base.:(==)(A::AbstractNMRData{T,N}, B::AbstractNMRData{T,N}) where {T,N}
    return size(A) == size(B) && all(A .== B)
end

Base.parent(A::AbstractNMRData) = A.data
# NOTE In the future, specialise for file-based data (see Rasters.jl/array.jl)
Base.Array(A::AbstractNMRData) = Array(parent(A))
Base.collect(A::AbstractNMRData) = collect(parent(A))
# TODO define similar() to inherit metadata appropriately

# Interface methods ################################################################################

"""
    missingval(x)

Returns the value representing missing data in the dataset.
"""
function missingval end
missingval(x) = missing
missingval(A::AbstractNMRData) = A.missingval

# DimensionalData methods ##########################################################################

# Rebuild types of AbstractNMRData
# Note: metadata is deep-copied by default to ensure sliced/copied data have independent metadata.
# See GitHub issue #25: https://github.com/waudbylab/NMRTools.jl/issues/25
function DD.rebuild(X::A, data, dims::Tuple, refdims, name,
                    metadata=deepcopy(metadata(X)),
                    missingval=missingval(X)) where {A<:AbstractNMRData}
    # HACK use A.name.wrapper to return the type (i.e. constructor) stripped of parameters.
    # This may not be a stable feature. See discussions:
    # - https://discourse.julialang.org/t/get-generic-constructor-of-parametric-type/57189/5
    # - https://github.com/JuliaObjects/ConstructionBase.jl/blob/master/src/constructorof.md
    return (A.name.wrapper)(data, dims, refdims, name, metadata, missingval)
end

function DD.rebuild(A::AbstractNMRData;
                    data=parent(A), dims=dims(A), refdims=refdims(A), name=name(A),
                    metadata=deepcopy(metadata(A)), missingval=missingval(A))
    return rebuild(A, data, dims, refdims, name, metadata, missingval)
end

# NOTE In the future, specialise for file-based data (see Rasters.jl/array.jl)
function DD.modify(f, A::AbstractNMRData)
    newdata = f(parent(A))
    size(newdata) == size(A) ||
        throw(DimensionMismatch("$f returns an array with size $(size(newdata)) when the original size was $(size(A))"))
    return rebuild(A, newdata)
end

# function DD.DimTable(As::Tuple{<:AbstractNMRData,Vararg{<:AbstractNMRData}}...)
#     return DD.DimTable(DimStack(map(read, As...)))
# end

# Concrete implementation ##########################################################################

"""
    NMRData <: AbstractNMRData
    NMRData(A::AbstractArray{T,N}, dims; kw...)
    NMRData(A::AbstractNMRData; kw...)

A generic [`AbstractNMRData`](@ref) for NMR array data. It holds memory-backed arrays.

# Keywords

- `dims`: `Tuple` of `NMRDimension`s for the array.
- `name`: `Symbol` name for the array, which will also retreive named layers if `NMRData`
    is used on a multi-layered file like a NetCDF.
- `missingval`: value representing missing data, normally detected from the file. Set manually
    when you know the value is not specified or is incorrect. This will *not* change any
    values in the NMRData, it simply assigns which value is treated as missing.

# Internal Keywords

In some cases it is possible to set these keywords as well.

- `data`: can replace the data in an `AbstractNMRData`

- `refdims`: `Tuple of` position `Dimension`s the array was sliced from, defaulting to `()`.
"""
struct NMRData{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me,Mi} <:
       AbstractNMRData{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
end

function NMRData(A::AbstractArray, dims::Tuple;
                 refdims=(), name=Symbol(""), metadata=defaultmetadata(NMRData),
                 missingval=missing)
    return NMRData(A, Dimensions.format(dims, A), refdims, name, metadata, missingval)
end

function NMRData(A::AbstractArray{<:Any,1}, dims::Tuple{<:Dimension,<:Dimension,Vararg};
                 kw...)
    return NMRData(reshape(A, map(length, dims)), dims; kw...)
end

# function NMRData(table, dims::Tuple; name=first(_not_a_dimcol(table, dims)), kw...)
#     Tables.istable(table) ||
#         throw(ArgumentError("First argument to `NMRData` is not a table or other known object: $table"))
#     isnothing(name) && throw(UndefKeywordError(:name))
#     cols = Tables.columns(table)
#     A = reshape(cols[name], map(length, dims))
#     return NMRData(A, dims; name, kw...)
# end

NMRData(A::AbstractArray; dims, kw...) = NMRData(A, dims; kw...)

# Note: metadata is deep-copied by default to ensure copied data have independent metadata.
# See GitHub issue #25: https://github.com/waudbylab/NMRTools.jl/issues/25
function NMRData(A::AbstractDimArray;
                 data=parent(A), dims=dims(A), refdims=refdims(A),
                 name=name(A), metadata=deepcopy(metadata(A)), missingval=missingval(A),
                 kw...)
    return NMRData(data, dims; refdims, name, metadata, missingval, kw...)
end

# function NMRData(filename::AbstractString, dims::Tuple{<:Dimension,<:Dimension,Vararg};
#                  kw...)
#     return NMRData(filename; dims, kw...)
# end

# TODO ??
# DD.dimconstructor(::Tuple{<:Dimension{<:AbstractProjected},Vararg{<:Dimension}}) = NMRData

# Getters ##########################################################################################
"""
    data(nmrdata)

Return the numerical data associated with the specified NMRData.
"""
data(A::NMRData) = A.data

"""
    data(nmrdata, dim)

Return the numerical data associated with the specified dimension.
"""
data(A::NMRData, dim) = data(dims(A, dim))

# Methods ##########################################################################################
"""
    scale(d::AbstractNMRData)

Return a scaling factor for the data combining the number of scans, receiver gain, and, if
specified, the sample concentration.
```math
\\mathrm{scale} = \\mathrm{ns} \\cdot \\mathrm{rg} \\cdot \\mathrm{conc}
```
"""
function scale(d::AbstractNMRData)
    # get ns, rg, and concentration, with safe defaults in case missing
    ns = get(metadata(d), :ns, 1)
    rg = get(metadata(d), :rg, 1)
    conc = get(metadata(d), :conc, 1)
    return ns * rg * conc
end

# Set pseudo-dimension data

"""
    setkinetictimes(A::NMRData, [dimnumber], tvals, units="")
    setkinetictimes([dimnumber], tvals, units="")

Return a new NMRData with a kinetic time axis containing the passed values (and optionally, units).
If a dimension number is specified, that dimension will be replaced with a `TkinDim`. If not, the
function will search for a unique non-frequency dimension, and replace that. If there are multiple
non-frequency dimensions, the dimension number must be specified explicitly, and the function will throw
an error.

The forms without the data argument return a function for use with `|>`.
"""
function setkinetictimes() end

function setkinetictimes(A::NMRData, tvals::AbstractVector, units="")
    hasnonfrequencydimension(A) ||
        throw(NMRToolsError("cannot set time values: data does not have a non-frequency dimension"))

    nonfreqdims = isa.(dims(A), NonFrequencyDimension)
    sum(nonfreqdims) == 1 ||
        throw(NMRToolsError("multiple non-frequency dimensions are present - ambiguous command"))

    olddimnumber = findfirst(nonfreqdims)

    return setkinetictimes(A::NMRData, olddimnumber, tvals, units)
end

function setkinetictimes(A::NMRData, dimnumber::Integer, tvals::AbstractVector, units="")
    newdim = TkinDim(tvals)
    newA = replacedimension(A::NMRData, dimnumber, newdim)

    label!(newA, newdim, "Time elapsed")
    newA[newdim, :units] = units

    return newA
end

"""
    setrelaxtimes(A::NMRData, [dimnumber], tvals, units="")
    setrelaxtimes([dimnumber], tvals, units="")

Return a new NMRData with a relaxation time axis containing the passed values (and optionally, units).
If a dimension number is specified, that dimension will be replaced with a `TrelaxDim`. If not, the
function will search for a unique non-frequency dimension, and replace that. If there are multiple
non-frequency dimensions, the dimension number must be specified explicitly, and the function will throw
an error.

The forms without the data argument return a function for use with `|>`.
"""
function setrelaxtimes() end

function setrelaxtimes(A::NMRData, tvals::AbstractVector, units="")
    hasnonfrequencydimension(A) ||
        throw(NMRToolsError("cannot set time values: data does not have a non-frequency dimension"))

    nonfreqdims = isa.(dims(A), NonFrequencyDimension)
    sum(nonfreqdims) == 1 ||
        throw(NMRToolsError("multiple non-frequency dimensions are present - ambiguous command"))

    olddimnumber = findfirst(nonfreqdims)

    return setrelaxtimes(A::NMRData, olddimnumber, tvals, units)
end

function setrelaxtimes(A::NMRData, dimnumber::Integer, tvals::AbstractVector, units="")
    newdim = TrelaxDim(tvals)
    newA = replacedimension(A::NMRData, dimnumber, newdim)

    label!(newA, newdim, "Relaxation time")
    newA[newdim, :units] = units

    return newA
end

"""
    setgradientlist(A::NMRData, [dimnumber], relativegradientlist, Gmax=nothing)
    setgradientlist([dimnumber], relativegradientlist, Gmax=nothing)

Return a new NMRData with a gradient axis containing the passed values. If no maximum strength is specified,
a default gradient strength of 0.55 T m⁻¹ will be set, but a warning raised for the user.

The forms without the data argument return a function for use with `|>`.

If a dimension number is specified, that dimension will be replaced. If not, the function will search for
a unique non-frequency dimension, and replace that. If there are multiple non-frequency dimensions, the
dimension number must be specified explicitly, and the function will throw an error.
"""
function setgradientlist() end

function setgradientlist(A::NMRData, dimnumber::Integer,
                         relativegradientlist::AbstractVector, Gmax=nothing)
    if isnothing(Gmax)
        @warn("a maximum gradient strength of 0.55 T m⁻¹ is being assumed - this is roughly correct for modern Bruker systems but calibration is recommended")
        gvals = 0.55 * relativegradientlist
    else
        gvals = Gmax * relativegradientlist
    end

    if dimnumber == 1
        newdim = G1Dim(gvals)
    elseif dimnumber == 2
        newdim = G2Dim(gvals)
    elseif dimnumber == 3
        newdim = G3Dim(gvals)
    end
    newA = replacedimension(A::NMRData, dimnumber, newdim)

    label!(newA, newdim, "Gradient strength")
    newA[newdim, :units] = "T m⁻¹"

    return newA
end

function setgradientlist(A::NMRData, relativegradientlist::AbstractVector, Gmax=nothing)
    hasnonfrequencydimension(A) ||
        throw(NMRToolsError("cannot set gradient values: data does not have a non-frequency dimension"))

    unknowndims = isa.(dims(A), UnknownDimension)
    sum(unknowndims) == 1 ||
        throw(NMRToolsError("multiple unknown dimensions are present - ambiguous command"))

    olddimnumber = findfirst(unknowndims)

    return setgradientlist(A, olddimnumber, relativegradientlist, Gmax)
end

"""
    setoffsets(A::NMRData, dimnumber, offsets, units="ppm")

Return a new NMRData with an offset axis (for CEST, R1rho off-resonance, etc.).
Dimension number must be specified explicitly.
"""
function setoffsets(A::NMRData, dimnumber::Integer, offsets::AbstractVector, units="ppm")
    newdim = OffsetDim(offsets)
    newA = replacedimension(A, dimnumber, newdim)
    label!(newA, newdim, "Offset")
    newA[newdim, :units] = units
    return newA
end

"""
    setspinlockfield(A::NMRData, dimnumber, fields, units="Hz")

Return a new NMRData with a spinlock field strength axis.
Dimension number must be specified explicitly.
"""
function setspinlockfield(A::NMRData, dimnumber::Integer, fields::AbstractVector,
                          units="Hz")
    newdim = SpinlockDim(fields)
    newA = replacedimension(A, dimnumber, newdim)
    label!(newA, newdim, "Spinlock field")
    newA[newdim, :units] = units
    return newA
end

"""
    replacedimension(nmrdata, olddimnumber, newdim)
    replacedimension(olddimnumber, newdim)

Return a new NMRData, in which the numbered axis is replaced by a new `Dimension`.

The forms without the data argument return a function for use with `|>`.
"""
function replacedimension(A::NMRData, olddimnumber, newdim)
    olddim = dims(A, olddimnumber)
    length(olddim) == length(newdim) ||
        throw(DimensionMismatch("replacedimension: size of old ($(length(olddim))) and new ($(length(newdim))) dimensions are not compatible"))
    merge!(metadata(newdim).val, metadata(olddim).val)

    olddims = Vector{NMRDimension}([dims(A)...])
    olddims[olddimnumber] = newdim
    newdims = tuple(olddims...)

    return NMRData(A.data, newdims; name=A.name, refdims=A.refdims,
                   metadata=deepcopy(A.metadata))
end

"""
    decimate(data, n dims=1)
    decimate(n)

Decimate NMR data into n-point averages along the specified dimension.
Note that data are *averaged* and not *summed*. Noise metadata is not updated.

The forms without the data argument return a function for use with `|>`.
"""
function decimate(expt::NMRData, n, dims=1)
    np = size(expt, dims) ÷ n
    out = selectdim(expt, dims, 1:n:(np * n)) # preallocate decimated NMRData

    sz = [size(out)...]
    insert!(sz, dims, n)

    y = selectdim(data(expt), dims, 1:(np * n))
    y = reshape(y, sz...)
    y = reduce(+, y; dims=dims) / n

    outsz = [size(out)...]
    outsz[dims] = np
    data(out) .= reshape(y, outsz...)

    return out
end

function decimate(signal, n, dims=1)
    np = size(signal, dims) ÷ n
    y = selectdim(signal, dims, 1:(np * n))

    sz = [size(signal)...]
    sz[dims] = np
    insert!(sz, dims, n)

    y = reshape(y, sz...)
    y = sum(y; dims=dims) / n

    outsz = [size(signal)...]
    outsz[dims] = np
    return reshape(y, outsz...)
end

"""
    stack(expts::Vector{NMRData})

Combine a collection of equally-sized NMRData into one larger array, by arranging them
along a new dimension, of type `UnknownDimension`.

Throws a `DimensionMismatch` if data are not of compatible shapes.
"""
function Base.stack(expts::Vector{D}) where {D<:AbstractNMRData{T,N}} where {T,N}
    # construct new data
    newdata = stack(data, expts)

    # construct new dimensions
    n = length(expts)
    if N == 1
        newdim = X2Dim(1:n)
    elseif N == 2
        newdim = X3Dim(1:n)
    elseif N == 3
        newdim = X4Dim(1:n)
    end
    newdims = [dims(expts[1])..., newdim]
    # push!(newdims, newdim)
    newdims = tuple(newdims...)

    # check and warn if dimensions don't match
    allequal(map(dims, expts)) || warn("stack: dimensions of expts 1 and $i are not equal")

    # construct new metadata, based on metadata for first spectrum
    newmd = deepcopy(metadata(first(expts)))
    newmd[:ndim] = N + 1

    # check and warn if ns don't match between experiments
    ns = map(e -> e[:ns], expts)
    allequal(ns) || @warn "stack: experiments do not have same ns" ns
    rg = map(e -> e[:rg], expts)
    allequal(rg) || @warn "stack: experiments do not have same rg" rg

    return NMRData(newdata, newdims; metadata=newmd)
end

# Arithmetic #######################################################################################

# BUG - it would be nice to propagate the noise, but at the moment metadata is shallow-copied during
# arithmetic operations so the operation changes the original input values as well.

# function Base.:+(a::NMRData, b::NMRData)
#     c = invoke(+, Tuple{DimensionalData.AbstractDimArray, DimensionalData.AbstractDimArray}, a, b)
#     c = NMRData(c.data, dims(c), name=c.name, refdims=c.refdims, metadata=copy(metadata(c)))
#     delete!(metadata(c), :noise) # delete old entry
#     c[:noise] = sqrt(get(metadata(a), :noise, 0)^2 + get(metadata(b), :noise, 0)^2)
#     return c
# end

# function Base.:-(a::NMRData, b::NMRData)
#     c = invoke(-, Tuple{DimensionalData.AbstractDimArray, DimensionalData.AbstractDimArray}, a, b)
#     c = NMRData(c.data, dims(c), name=c.name, refdims=c.refdims, metadata=copy(metadata(c)))
#     delete!(metadata(c), :noise) # delete old entry
#     c[:noise] = sqrt(get(metadata(a), :noise, 0)^2 + get(metadata(b), :noise, 0)^2)
#     return c
# end

# function Base.:*(a::Number, b::NMRData)
#     c = invoke(*, Tuple{Number, DimensionalData.AbstractDimArray}, a, b)
#     c = NMRData(c.data, dims(c), name=c.name, refdims=c.refdims, metadata=copy(metadata(c)))
#     delete!(metadata(c), :noise) # delete old entry
#     c[:noise] = a * get(metadata(b), :noise, 0)
#     return c
# end

# function Base.:/(a::NMRData, b::Number)
#     c = invoke(/, Tuple{DimensionalData.AbstractDimArray, Number}, a, b)
#     c = NMRData(c.data, dims(c), name=c.name, refdims=c.refdims, metadata=copy(metadata(c)))
#     delete!(metadata(c), :noise) # delete old entry
#     c[:noise] = get(metadata(a), :noise, 0) / b
#     return c
# end

# Traits ###########################################################################################
"""
    @traitdef HasNonFrequencyDimension{D}

A trait indicating whether the data object has a non-frequency domain dimension.

# Example

```julia
@traitfn f(x::X) where {X; HasNonFrequencyDimension{X}} = "This spectrum has a non-frequency domain dimension!"
@traitfn f(x::X) where {X; SimpleTraits.Not{HasNonFrequencyDimension{X}}} = "This is a pure frequency-domain spectrum!"
```
"""
@traitdef HasNonFrequencyDimension{D}

# adapted from expansion of @traitimpl
# this is a pure function because it only depends on the type definition of NMRData
Base.@pure function SimpleTraits.trait(t::Type{HasNonFrequencyDimension{D}}) where {D<:NMRData{T,
                                                                                               N,
                                                                                               A}} where {T,
                                                                                                          N,
                                                                                                          A}
    if any(map(dim -> (typeintersect(dim, FrequencyDimension) == Union{}), A.parameters))
        HasNonFrequencyDimension{D}
    else
        SimpleTraits.Not{HasNonFrequencyDimension{D}}
    end
end

"""
    hasnonfrequencydimension(spectrum)

Return true if the spectrum contains a non-frequency domain dimension.

# Example

```julia
julia> y2=loadnmr("exampledata/2D_HN/test.ft2");
julia> hasnonfrequencydimension(y2)
false
julia> y3=loadnmr("exampledata/pseudo3D_HN_R2/ft/test%03d.ft2");
julia> hasnonfrequencydimension(y3)
true
```
"""
function hasnonfrequencydimension(spectrum::NMRData{T,N,A}) where {T,N,A}
    return any(map(dim -> (typeintersect(dim, FrequencyDimension) == Union{}),
                   A.parameters))
end

# curried forms for piping
setkinetictimes(tvals::AbstractVector, units::AbstractString="") = A -> setkinetictimes(A, tvals, units)
function setkinetictimes(dimnumber::Integer, tvals::AbstractVector, units::AbstractString="")
    return A -> setkinetictimes(A, dimnumber, tvals, units)
end
setrelaxtimes(tvals::AbstractVector, units::AbstractString="") = A -> setrelaxtimes(A, tvals, units)
function setrelaxtimes(dimnumber::Integer, tvals::AbstractVector, units::AbstractString="")
    return A -> setrelaxtimes(A, dimnumber, tvals, units)
end
setgradientlist(list::AbstractVector, Gmax=nothing) = A -> setgradientlist(A, list, Gmax)
function setgradientlist(dimnumber::Integer, list::AbstractVector, Gmax=nothing)
    return A -> setgradientlist(A, dimnumber, list, Gmax)
end
replacedimension(olddimnumber, newdim::NMRDimension) = A -> replacedimension(A, olddimnumber, newdim)
decimate(n::Integer) = A -> decimate(A, n)
