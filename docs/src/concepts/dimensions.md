# Dimensions

Each axis of an [`NMRData`](@ref) has a dimension type, which says what the axis measures. The type
determines how processing functions treat the axis, how plots label it, and which selectors make
sense on it.

## The dimension types

```
NMRDimension
├── FrequencyDimension          F1Dim  F2Dim  F3Dim  F4Dim
└── NonFrequencyDimension
    ├── TimeDimension           T1Dim  T2Dim  T3Dim  T4Dim   (frequency evolution)
    │                           TrelaxDim                    (relaxation delays)
    │                           TkinDim                      (real time, kinetics)
    ├── GradientDimension       G1Dim  G2Dim  G3Dim  G4Dim   (diffusion)
    ├── OffsetDimension         OffsetDim                    (CEST and other offsets)
    ├── FieldDimension          SpinlockDim                  (spin-lock fields)
    └── UnknownDimension        X1Dim  X2Dim  X3Dim  X4Dim   (anything else)
```

- **Frequency dimensions** (`F1Dim`–`F4Dim`) have chemical shift coordinates in ppm. They are
  produced by Fourier transformation, or loaded directly from processed data.
- **Time dimensions for frequency evolution** (`T1Dim`–`T4Dim`) have coordinates in seconds. They
  are what [`loadfid`](@ref NMRTools.NMRIO.loadfid) returns, and what [`ft`](@ref) turns into
  frequency dimensions.
- **Other non-frequency dimensions** describe the extra axis of a pseudo-2D or pseudo-3D series:
  relaxation delays, kinetic time points, gradient strengths, offsets or spin-lock fields. They are
  never Fourier transformed.
- **Unknown dimensions** (`X1Dim`–`X4Dim`) are non-frequency axes whose meaning NMRTools could not
  work out from the acquisition parameters. A pseudo-2D relaxation series is loaded with an
  `X2Dim`, and you replace it with a `TrelaxDim` once you know the delays.

The number in a name gives the order of the dimension in the data: `F1Dim` (or `T1Dim`) is the
direct dimension, `F2Dim` the first indirect dimension, and so on.

## From time to frequency

`ft` replaces each `TnDim` with the matching `FnDim`, keeping its metadata. Non-Fourier
dimensions are left alone. For a 2D experiment:

```
loadfid   →  (T1Dim, T2Dim)       time domain, Multicomplex{Float64,2,4}
ft        →  (F1Dim, F2Dim)       frequency domain, Multicomplex{Float64,2,4}
realest   →  (F1Dim, F2Dim)       frequency domain, Float64
```

and for a pseudo-2D relaxation series:

```
loadfid         →  (T1Dim, X2Dim)
setrelaxtimes   →  (T1Dim, TrelaxDim)
ft              →  (F1Dim, TrelaxDim)
```

Processing functions use the dimension types to decide what to act on: called without a dimension,
[`apodize`](@ref), [`zerofill`](@ref) and [`ft`](@ref) act on every `TimeDimension` that has an
imaginary unit, and skip the rest.

## Referring to a dimension

Functions that take a dimension accept its type (`F1Dim`), its number (`1`), or, for frequency
dimensions, its nucleus (`H1`, `N15`). [`resolvedim`](@ref) converts any of these to a number and
[`finddim`](@ref) finds the dimension of a given nucleus:

```@example dims
using NMRTools
spec = exampledata("2D_HN")
finddim(spec, N15)
```

## Changing a dimension

[`replacedimension`](@ref) swaps one dimension for another of a different type with new
coordinates. The `set*` functions do this for the common cases:

| Function | New dimension |
|---|---|
| [`setrelaxtimes`](@ref) | `TrelaxDim` |
| [`setkinetictimes`](@ref) | `TkinDim` |
| [`setgradientlist`](@ref) | `G1Dim`… |
| [`setoffsets`](@ref) | `OffsetDim` |
| [`setspinlockfield`](@ref) | `SpinlockDim` |

[`shiftdim`](@ref) moves the coordinates of a frequency dimension by a fixed number of ppm.

## Related pages

- [Multicomplex numbers and quadrature](multicomplex.md): how each time and frequency dimension is
  linked to its imaginary unit
- [Process step by step](../howto/process.md): processing pseudo-2D data
