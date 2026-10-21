# Metadata

[`NMRData`](@ref) objects carry metadata on acquisition and processing, filled in automatically
when data are loaded and updated by processing functions.

## Two levels of metadata

Metadata are attached at two levels:

1. **Spectrum metadata** describe the experiment as a whole: number of scans, receiver gain, pulse
   programme, title, temperature, the noise level (estimated when the spectrum is loaded),
   acquisition parameters from the `acqus` files, and the contents of auxiliary files such as
   vclists and vdlists.
2. **Axis metadata** describe one dimension: number of points, the original time-domain size
   (before zero-filling, linear prediction or extraction of a region), carrier and base
   frequencies, spectral width, nucleus, window function, phase correction and the dimension's
   imaginary unit (`:mcindex`).

Both are dictionaries keyed by symbols such as `:ns` or `:window`. The standard keys are listed in
[Metadata keys](../reference/metadatakeys.md).

## Looking up metadata

Spectrum metadata are looked up with a key; axis metadata with a dimension and a key. The
dimension can be given by type, number or (for frequency dimensions) nucleus:

```julia
metadata(spec, :ns)            # spectrum metadata
spec[:ns]                      # the same

metadata(spec, F1Dim, :bf)     # axis metadata
spec[F1Dim, :bf]               # the same
spec[1, :bf]                   # by dimension number

metadata(spec)                 # the whole spectrum dictionary
metadata(spec, F1Dim)          # the whole dictionary for one axis
```

```@example meta
using NMRTools
spec = exampledata("2D_HN")
spec[:pulseprogram], spec[F2Dim, :nucleus], spec[F1Dim, :window]
```

Acquisition parameters are kept as parsed from the `acqus` files and read with
[`acqus`](@ref), which accepts Bruker parameter names in either case; see
[Find acquisition parameters](../howto/acqus.md).

## Metadata types

Where a parameter has structure beyond a single number, NMRTools stores it as a typed object
instead of a raw value:

| Parameter | Type | Reference |
|---|---|---|
| apodization applied in processing | [`WindowFunction`](@ref) | [Window functions](../reference/windows.md) |
| power levels (`PLW`, `PL`) | [`Power`](@ref) | [Power levels](../reference/power.md) |
| frequency lists (`fq1list`, …) | [`FQList`](@ref) | [Frequency lists](../reference/fqlists.md) |
| sample description | [`NMRSample`](@ref) | [Sample metadata schema](../reference/samples.md) |
| spectrometer channels | `Dict` per channel | [Find acquisition parameters](../howto/acqus.md) |
| pulse programme annotations | `Dict` | [Pulse programme annotations](../experimental/annotations.md) |

Times are stored in seconds and frequencies in Hz throughout; see
[Time and frequency units](units.md).

## How processing updates metadata

Processing functions record what they did on the dimension they acted on:

- [`apodize`](@ref) sets `:window`;
- [`zerofill`](@ref) updates `:tdzf`, and [`linearpredict`](@ref) updates `:td`;
- [`ft`](@ref) carries the time-dimension metadata over to the new frequency dimension and sets
  `:grpdly` to 0;
- [`phase`](@ref) adds to `:phase`, the total correction applied so far;
- [`phase`](@ref), [`extract`](@ref) and [`ft`](@ref) re-estimate the spectrum's `:noise`.

!!! warning "Shared metadata"
    Slices and copies of an `NMRData` share the same metadata dictionary as the original, so
    changing metadata on one changes it on all of them. Use `deepcopy` to make an independent copy.
