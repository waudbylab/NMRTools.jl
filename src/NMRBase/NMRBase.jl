module NMRBase

const EXPERIMENTAL = """
    WARNING: This feature is experimental. It may change in future versions, and may
    not be 100% reliable in all cases. Please file github issues if problems occur.
    """

# DimensionalData documentation urls
const DDdocs = "https://rafaqz.github.io/DimensionalData.jl/stable/api"
const DDdimdocs = joinpath(DDdocs, "#DimensionalData.Dimension")
const DDarraydocs = joinpath(DDdocs, "#DimensionalData.AbstractDimensionalArray")
# const DDabssampleddocs = joinpath(DDdocs, "#DimensionalData.AbstractSampled")
# const DDsampleddocs = joinpath(DDdocs, "#DimensionalData.Sampled")
# const DDlocusdocs = joinpath(DDdocs, "#DimensionalData.Locus")
# const DDselectordocs = joinpath(DDdocs, "#DimensionalData.Selector")
# const DDtidocs = joinpath(DDdocs, "#DimensionalData.Ti")

using Dates: Dates, DateTime, Date, unix2datetime
using LinearAlgebra
using MulticomplexNumbers
using Optim
using Peaks
using SimpleTraits
using SpecialFunctions
using Statistics

# using Reexport: Reexport
# Reexport.@reexport using DimensionalData
using DimensionalData
using DimensionalData.Lookups,
      DimensionalData.Dimensions
const DD = DimensionalData

# using DimensionalData: Name, NoName
# using .Dimensions: StandardIndices, DimTuple, Dimension
# using .Lookups: LookupTuple

import DimensionalData.refdims_title
import .Lookups: metadata, set, _set, rebuild, basetypeof, Metadata,
                 order, span, sampling, locus, val, bounds, hasselection, units,
                 SelectorOrInterval, ForwardOrdered
import .Dimensions: dims, refdims, lookup, hasdim, label, _astuple

include("exceptions.jl")
include("nuclei.jl")
include("coherences.jl")
include("isotopedata.jl")
include("power.jl")
include("dimensions.jl")
include("nmrdata.jl")
include("windows.jl")
include("metadata.jl")
include("noise.jl")
include("peaks.jl")
include("pcls.jl")
include("reference.jl")

macro exportinstances(enum)
    eval = GlobalRef(Core, :eval)
    return :($eval($__module__, Expr(:export, map(Symbol, instances($enum))...)))
end

# Exceptions
export NMRToolsError

# Nuclei
export Nucleus
@exportinstances Nucleus
export nucleus
export spin
export gyromagneticratio
export xi_ratio

# coherences
export Coherence
export SQ, MQ
export coherenceorder

# power
export Power
export db, watts
export hz

# NMRData
export AbstractNMRData
export NMRData
export NMRExperiment
export NMRSample
# Selectors
export Selector, IntSelector, ArraySelector
export At, Between, Touches, Contains, Near, Where, All, ..
# getter methods
export data, parent, dims, refdims, lookup, bounds, missingval
# traits
export hasnonfrequencydimension, HasNonFrequencyDimension
# Dimension/Lookup primitives
export dimnum, hasdim, hasselection, otherdims
# utils
export refdims_title
export set, rebuild, reorder, modify, broadcast_dims, broadcast_dims!, ForwardOrdered
# NMR properties
export scale
export estimatenoise
export decimate

# Dimensions
export NMRDimension
export FrequencyDimension
export NonFrequencyDimension
export TimeDimension
export GradientDimension
export UnknownDimension
export OffsetDimension
export FieldDimension
export F1Dim, F2Dim, F3Dim, F4Dim
export T1Dim, T2Dim, T3Dim, T4Dim
export TrelaxDim, TkinDim
export G1Dim, G2Dim, G3Dim, G4Dim
export X1Dim, X2Dim, X3Dim, X4Dim
export OffsetDim, SpinlockDim
# export GradientDim
# export SpatialDim
export ppm, hz
export replacedimension
export setkinetictimes
export setrelaxtimes
export setgradientlist
export setoffsets
export setspinlockfield
export shiftdim
export finddim, resolvedim

# Chemical shift referencing
export reference
export isaqueous, watershift

# Metadata
export AbstractMetadata, Metadata, NoMetadata
export metadata
export defaultmetadata
export metadatahelp
export label, label!, units
export acqus
export annotations, hasannotations
export sample, hassample, samplefile
export referencepulse

# Window functions
export WindowFunction
export NullWindow
export UnknownWindow
export ExponentialWindow
export SineWindow
export GeneralSineWindow, CosWindow, Cos²Window
export GaussWindow
export GeneralGaussWindow, LorentzToGaussWindow
export lineshape
export apod

export LineshapeComplexity
export RealLineshape
export ComplexLineshape

# peaks and peak detection
export AbstractPeak, SimplePeak
export detectpeaks

# phase-constrained least squares
export pcls

end
