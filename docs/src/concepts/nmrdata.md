# The NMRData model

Every spectrum, FID or series in NMRTools is an [`NMRData`](@ref): an array of data points,
a set of dimensions describing its axes, and a dictionary of metadata. Loading, processing,
plotting and peak picking all take and return `NMRData`.

```@example nmrdata
using NMRTools
spec = exampledata("2D_HN")
size(spec), eltype(spec)
```

## An array with axes

`NMRData` is a subtype of `AbstractDimArray` from
[DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl/stable/). It behaves like an
ordinary Julia array: it can be indexed, sliced, broadcast and passed to functions that expect an
`AbstractArray`. The difference is that each axis carries its coordinates, so the array also
answers questions about chemical shift, time or gradient strength:

```@example nmrdata
dims(spec)
```

Indexing by position and indexing by value use the same square brackets. Value-based selectors
come from DimensionalData: `Near(x)` for the nearest point, `At(x)` for an exact match, and
`a .. b` for a closed interval:

```@example nmrdata
spec[Near(8.5), 120 .. 121]
```

The result of slicing is another `NMRData`, with its axes cut to match. A slice through a 2D
spectrum is a 1D spectrum that can be plotted, phased or peak-picked like any other.

`data(spec)` returns the plain array underneath, and `data(spec, dim)` the plain vector of axis
coordinates.

## Operations preserve structure

Arithmetic and broadcasting return `NMRData` with the same axes and metadata:

```@example nmrdata
spec2 = 2 * spec .+ 1
typeof(spec2) == typeof(spec)
```

Reductions over a dimension, such as `sum(spec; dims=F1Dim)`, keep the remaining dimensions. To
build a new `NMRData` with different data but the same axes and metadata, use `rebuild`:

```julia
rebuild(spec, newdata)
```

!!! warning "Shared metadata"
    Slices and copies of an `NMRData` share the same metadata dictionary as the original, so
    changing metadata on one changes it on all of them. Use `deepcopy` to make an independent
    copy. Acquisition metadata (`:acqus`, `:td`, …) describe the experiment as acquired and are not
    updated when data are sliced.

## Element types

The element type records what kind of data an `NMRData` holds:

| Element type | Meaning | Produced by |
|---|---|---|
| `Float64` | real spectrum | `loadnmr` (default), `autoprocess`, `real`, `realest` |
| `Multicomplex{Float64,1,2}` | complex in one dimension | `loadfid` on a 1D or pseudo-2D, `loadnmr(...; allcomponents=true)` on 1D data |
| `Multicomplex{Float64,2,4}` | complex in two dimensions | `loadfid` on a 2D |
| `Multicomplex{Float64,3,8}` | complex in three dimensions | `loadfid` on a 3D |

Processing keeps the element type unchanged from `loadfid` to the final phase correction, and
`realest` at the end gives a real spectrum. See
[Multicomplex numbers and quadrature](multicomplex.md).

## Metadata

Metadata are attached at two levels: to the spectrum as a whole (`spec[:ns]`) and to each
dimension (`spec[F1Dim, :bf]`). See [Metadata](metadata.md).

## Related pages

- [Dimensions](dimensions.md): the dimension types and what they mean
- [Working with spectra](../tutorials/spectra.md): a tutorial on indexing, axes and metadata
- [DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl/stable/): selectors,
  `rebuild`, `set` and other functions re-exported by NMRTools
