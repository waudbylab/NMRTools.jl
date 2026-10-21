# Dimensions

Selectors (`At`, `Near`, `..`, `Where`, …) and functions such as `dims`, `rebuild` and `set` are re-exported from DimensionalData.jl; see the [DimensionalData.jl documentation](https://rafaqz.github.io/DimensionalData.jl/stable/).

```@index
Pages = ["dimensions.md"]
```

## Data

```@docs
NMRData
AbstractNMRData
data
missingval
```

## Dimension types

```@docs
NMRDimension
FrequencyDimension
NonFrequencyDimension
TimeDimension
GradientDimension
UnknownDimension
OffsetDimension
FieldDimension
HasNonFrequencyDimension
hasnonfrequencydimension
```

## Working with dimensions

```@docs
replacedimension
setrelaxtimes
setkinetictimes
setgradientlist
setoffsets
setspinlockfield
shiftdim
finddim
resolvedim
ppm
hz
```
